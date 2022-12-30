import datetime
import json
import logging
from pathlib import Path

import pandas as pd
import numpy as np
from typing import List, Sequence, Set, Union
import plotly

import multiprocessing
from dataclasses import dataclass, field

from lib.utils import NB_Utils

class StepFailedError(Exception):
    ...
    pass

@dataclass
class MAGpy_Cache:
    """
    Cache dataclass for PyMAG class.
    """
    iterations:   int       = field(default=0)
    iter_identifiers: List[Set[str]] = field(default_factory=list)
    iter_detection_mse: List[float] = field(default_factory=list)

    identifiers_dropped:      Set[str] = field(default_factory=set)
    identifiers_outliers:   Set[str] =  field(default_factory=set)

    def to_dict(cls):
        outdict = {
            'identifiers_dropped' : list(cls.identifiers_dropped),
            'identifiers_outliers' : list(cls.identifiers_dropped)
        }
        outdict.update({
            'iterations ':{
                i : {
                'mse': cls.iter_detection_mse[i],
                'identifiers': list(cls.iter_identifiers[i])
                }
                for i in range(cls.iterations+1)
            }
        })
        return outdict


    def __repr__(self):
        dict_repr = ', '.join(
            f'{k}:{type(v)}'
            for k, v in filter(
                lambda item: not item[0].startswith('_'),
                self.__dict__.items()
            )
        )
        return f'{self.__class__.__name__}({dict_repr})'


class MAGpy:
    """
    Implementation of NegBinom optimisation of kmer catalogue.
    The implementation is based on the original work MAGinator by Trine Zachariasen.

    The class implements an catalogue refiner which takes a count matrix with rows
    of identifiers and columns of samples and attempts to determine a subset of identifiers which
    have the best fit to a detection saturation curve.

    In essense the refinement attempts to improve the fit to the detection curve by choosing those
    identifiers which count best fit a NegativeBinomial distribution.

    ...

    Attributes
    ----------
    df_counts : pd.DataFrame
        Dataframe containing the counts from fp_counts subset to identifier in fp_meta

    Static Methods:
    -----
    refine_all:
        Takes a list of counts and metafiles and instantiates a PyMAG class for each and runs refinement.
        Main interface unless doing manual tuning.
    Class Methods:
    -----
    preflight:
        Runs input check, initiates cache with a "hot start"
    run_refinement:
        Runs the iterative refinement (multiple class to refine_step() )
    generate_minimal_output:
        Writes count matrix for refined identifier set and mse iterations to output dir.
    generate_verbose_output:
        Writes detection saturation plots and dumps cache to output.
    
    
    """
    def __init__(self, fp_counts: Path, fp_meta: Path, max_threads=2, catalogue_size:int=500, logfile=None, verbosity="DEBUG", rng_seed=2112):
        """
        Constructs all the necessary attributes for the person object.
        Checks the inputs are compliant.
        
        Parameters
        -------
        fp_counts : Path
            Path to .tsv file.
            index = identifiers, cols = samples
        fp_meta : Path
            Path to tsv file describing the identifiers.
            expected cols: identifier, prevalence, lenght
            NOTE: no headers.
        catalogue_size : int
            Size of generated catalogue of identifiers
        max_threads : int
            Maximum number of pool workers to run
        """

        
        self.fp_counts = Path(fp_counts)
        self.fp_meta   = Path(fp_meta)
        self.catalogue_size = catalogue_size

        self.rng_state = np.random.RandomState(rng_seed)

        if (safe_threads := multiprocessing.cpu_count()-1) < max_threads:
            self.log.warning(f"Reducing max-threads to (cpu_count-1 = {safe_threads})")
            max_threads = safe_threads
        self.max_threads = max_threads



        for input_file in [self.fp_counts, self.fp_meta]:
            if not input_file.is_file():
                raise FileExistsError(input)

        self.logfile    = Path(logfile) if logfile else None
        self.loglvl     = verbosity

        try:
            self._check_input_compliance()
        except Exception as err:
            self.log.exception("Failed input compliance check [TERMINATING]")
            raise err
        
    ## Front matter
    @property
    def log_setup(self):
        self.runname = self.fp_meta.stem
        return dict(
                name = self.__class__.__name__+"."+self.runname,
                level = logging.getLevelName(self.loglvl),
                log_file = self.logfile
            )    
    @property
    def log(self):
        if not hasattr(self, "_log"):
            setup = self.log_setup
            logger = logging.getLogger(setup["name"])
                      
           
            F = "[%(asctime)s %(name)s:%(funcName)s]%(levelname)s: %(message)s"
            formatter = logging.Formatter(F, datefmt='%d-%b-%y %H:%M:%S')
            console = logging.StreamHandler()
            console.setFormatter(formatter)
            console.setLevel(setup["level"])
            handlers = [console]
            
            if setup["log_file"]:
                Path(setup["log_file"]).parent.mkdir(parents=True, exist_ok=True)
                file_handler = logging.FileHandler(setup["log_file"])
                file_handler.setFormatter(formatter)
                file_handler.setLevel(setup["level"])
                handlers.append(file_handler)

            logger.handlers = handlers
            logger.setLevel(setup["level"])
            self._log = logger
        return self._log

        
    @property
    def cache(self):
        if not hasattr(self, "_cache"):
            self.log.debug("initialise cache")
            self._cache = MAGpy_Cache()
        return self._cache
    def reset_cache(self):
        try:
            delattr(self, "_cache")
        except AttributeError:
            self.log.debug("No cache to reset.")

    @property
    def df_counts_all(self):
        if not hasattr(self, "_df_counts_all"):
            self._df_counts_all = pd.read_csv(self.fp_counts, sep="\t", index_col=0)

            dims = self._df_counts_all.shape
            if dims[1] == 1:
                raise IOError(f"Bad format for fp_counts, expected .tsv file but loaded 1 column file. -> {self.fp_counts}")
            self.log.debug(f"Loaded counts of dimensions {dims}-> {self.fp_counts}")
        return self._df_counts_all

    @property
    def df_counts(self):
        if not hasattr(self, "_df_counts"):
            self._df_counts = self.df_counts_all.loc[self.df_meta["identifier"],:]
            #clear up memory
            delattr(self, "_df_counts_all")
        return self._df_counts

    @property
    def df_meta(self):
        if not hasattr(self, "_df_meta"):
            self.log.debug("loading meta")
            self._df_meta = pd.read_csv(self.fp_meta, sep="\t", names=["identifier","prevalence", "lenght"])
        return self._df_meta

    def _check_input_compliance(self) -> None:
        """
        Checks if identifiers given in fp_counts and fp_meta match.
        Raises Runtime error if they do not.
        Requires all identifiers in fp_meta to also be in fp_count but not reverse.
          

        """
        self.log.info("Checking compliance of input")
        meta_ids = self.df_meta["identifier"]
        set_meta_ids = set(meta_ids)
        if not len(meta_ids) == len(set_meta_ids):
            raise RuntimeError(f"count matrix contains duplicate identifiers")

        count_ids = self.df_counts_all.index
        set_count_ids = set(count_ids)
        if not len(count_ids) == len(set_count_ids):
            raise RuntimeError(f"count matrix contains duplicate identifier")

        if (missing_count:= set_meta_ids - set_count_ids):
            raise RuntimeError(f"Identifiers found in meta which is not found in count data -> {missing_count}")

        if (missing_meta := set_count_ids - set_meta_ids):
            self.log.warning(f"Following ids where not found in meta data (Allowed)) -> {len(missing_meta)}")

        self.log.info(f"Number of idenifiers ({len(set_meta_ids)})")
        self.log.info(f"Number of samples:   ({len(self.df_counts.columns)})")

    def detect_sample_wide_outliers(self, keep_percentile: int=80, hard_freq_limit=0.5) -> Set[str]:
        """
        Uses tukeys IQR method to call outliers per sample.
        Overall outliers are called based on the number of samples in which they where call as outlier.

        By default the 20% identifiers with highest number of within sample outlier calls are consideres as 
        sample wide outliers.
        
        Parameters
        -----
        keep_percentile: int
            The percentile of least outlying identifiers to keep (default: 80)
        hard_freq_limit: float
            Maximum freq of outlier occourences overwrites the percentile value if smaller (default: 0.5)

        Returns
        -----
        Outliers : Set[str]
            Set of identifiers called as overall outliers
        identifier_outlier_occcourence: List[int]
            list of outlier occourences per identifier
        """
        with multiprocessing.Pool(self.max_threads) as p:
            #Run outlier detection per sample and return to row=id, col=sample format.
            per_sample_outliers = np.array(p.map(NB_Utils.tukeys_outlier_detection, self.df_counts.values.T)).T
        # Summarise how often each id is an outlier
        # We could consider weighing the number by whether they were an outlier
        # In a sample with many outliers - but this will no be implemented at this point.
        identifiers = self.df_counts.index
        # Mark the ??% with most outlier occourences.
        identifier_outlier_occ = per_sample_outliers.sum(1)
        px_outlier_occ = np.percentile(identifier_outlier_occ, keep_percentile)
        
        hard_limit = int(hard_freq_limit * self.df_counts.shape[1])
        if hard_limit < px_outlier_occ:
            self.log.debug(f"Hard limit hit: hard_limit=({hard_limit:}) < percentile_val=({px_outlier_occ}).")
            max_occ = hard_limit
        else:
            self.log.debug(f"Using ({keep_percentile}) percentile cutoff = ({px_outlier_occ})")
            max_occ = px_outlier_occ

        outlier_identifier = set(identifiers[identifier_outlier_occ > max_occ])
        self.log.info(f"Outliers: ({len(outlier_identifier)})/({len(identifiers)}) have outlier occ > ({max_occ})")
        return outlier_identifier, identifier_outlier_occ


    def calculate_mean_rank_mse(self, identifiers: List[str]) -> List[float]:
        """
        Calculates the mean rank mse for each identifier.

        Each id is ranked within each sample and then the average rank is found
        for each identifier across all samples.

        Parameters:
        -------
        identifiers : List[str]
            list of identifiers in df_count

        Returns:
        -------
        mean_rank_mse (List[float]): list of floats of the average rank of mse across samples.
            Low values indicate the identifier generally performs well.
        """
        
        df_count = self.df_counts.loc[identifiers,:]
        
        sample_list = df_count.values.T
        with multiprocessing.Pool(self.max_threads) as p:
            rank_matrices, failed_fit =  zip(*p.map(NB_Utils.counts2ranked_pearson, sample_list))


        if any(failed_fit):
            self.log.debug(f"Fit warning, not used in ranking: {sum(failed_fit)}/{len(rank_matrices)} samples")
            #self.log.debug(f"Fit warning, not used in ranking: {df_count.columns.values[list(failed_fit)]}")
            rank_matrices = [matrix for matrix, has_failed in zip(rank_matrices, failed_fit) if not has_failed]
        identifier_to_sample_ranks = np.array(rank_matrices).T
        identifier_mean_rank = identifier_to_sample_ranks.mean(1)
        return identifier_mean_rank

    def calculate_detection_mse(self, identifiers: List[str]) -> float:
        """
        Calculate the mse between expected and observed detection vs. total count curves.

        Parameters
        ------
        identifiers : list[str]

        Returns
        mse : float
        """
        counts = self.df_counts.loc[identifiers,:].values
        mse = NB_Utils.calculate_detection_mse(counts)
        return mse

    def preflight(self, max_retries: int = 50) -> None:
        """
        Initial steps before the full refining produce.

        Prepares identifiers and finds outliers.
        Attempts a few initial seeds to get a good one.

        Paraters
        -----
        max_retries : int
            Number of retries is attempted if a random selection of identifiers fail

        Setting
        -----
        Updates cache
        """
        self.reset_cache()

        identifiers_outliers, outlier_occ    = self.detect_sample_wide_outliers(keep_percentile=80)
        self.df_meta["is_outlier"]   = self.df_meta.identifier.isin(identifiers_outliers)
        self.df_meta["outlier_occ"]  = outlier_occ
        max_accepted_outlier_occ = self.df_meta.query("not is_outlier")["outlier_occ"].max()
        self.df_meta["outlier_degree"] = self.df_meta.outlier_occ / max_accepted_outlier_occ

        self.df_meta["identifier_desireability"] = self.df_meta["prevalence"] + (1 - self.df_meta["outlier_degree"])

        self.cache.identifiers_outliers.update(identifiers_outliers)
        
        df_meta_eligible = self.df_meta.loc[~self.df_meta["identifier"].isin(identifiers_outliers)]

        #Generate 10 starting identifier compositions.
        starter_compositions = []
        starter_compositions_mse = []
        starter_compositions_ranks = []
        retries = 0
        while (len(starter_compositions) < 10) and (retries < max_retries):
            df_random_subset = df_meta_eligible.sample(n=self.catalogue_size, weights="identifier_desireability", replace=False, random_state=self.rng_state)
            identifiers_random = df_random_subset.identifier.values

            try:
                ranks_identifier_random = self.calculate_mean_rank_mse(identifiers=identifiers_random)
                identifier_random_detection_mse = self.calculate_detection_mse(identifiers=identifiers_random)
            except Exception as ee:
                retries += 1
                self.log.exception(f"Exception during fitting, re-sampling: {retries}/{max_retries}")
                pass
            else:
                starter_compositions.append(identifiers_random)
                starter_compositions_mse.append(identifier_random_detection_mse)
                starter_compositions_ranks.append(ranks_identifier_random)
                self.log.debug(f"Added a valid starter composition {len(starter_compositions)}/10")

        best_iteration = np.argmin(starter_compositions_mse)

        self.log.info(f"Initial identifiers found. Initial MSE: {starter_compositions_mse[best_iteration]:.2f}")

        self.cache.iterations = 0
        self.cache.iter_identifiers.append(set(starter_compositions[best_iteration]))
        self.cache.iter_detection_mse.append(starter_compositions_mse[best_iteration])

    def get_substitute_identifiers(self, n_substitutes:int, black_list:Set[str]= set(),  retry_identifiers:bool=True):
        """
        Find a set of possible substitute identifiers.
        The possible substitutes are drawn from a pool of candidates.

        Parameters
        -----
        n_substitutes : int
            Number of identifiers required.
        black_list : Sequence[str]
            blacklisted identifiers.
        allow_restries : bool
            identifiers_previously attempted (not outlier) identifiers to be retried.
        
        Returns
        -----
        identifiers : set
            Set of random identifiers
        """
        identifier_candidates  = set(self.df_counts.index) 
        identifier_candidates -= self.cache.iter_identifiers[self.cache.iterations]
        identifier_candidates -= self.cache.identifiers_outliers
        if not retry_identifiers:
            identifier_candidates -= self.cache.identifiers_dropped

        self.log.debug(f"Selecting from ({len(identifier_candidates)}) candidates")

        df_meta_eligible = self.df_meta.loc[self.df_meta["identifier"].isin(identifier_candidates)]
        df_random_subset = df_meta_eligible.sample(
            n=n_substitutes,
            weights="identifier_desireability", 
            replace=False, 
            random_state=self.rng_state)
        return set(df_random_subset.identifier)

        

    def refine_step(self, step_size:int = 10, max_attempts: int=5, min_improvement : float = 0.1) -> None:
        """
        Do one iteration of the refinement process.
        Calculates the per identifier pearson residual and drops the worst (step_size) part.
        
        It then introduces new identifiers and checks if the MSE is improved.
            This step is done (attempts) times until improvement.
        if improved the iteration is stored in cache.

        Parameters
        -----
        step_size : int (]0:100[)
            proportion of identifiers to exchange must be between >0 and <100
        attempts : int
            Number attempts to include new kmers to get and improvement.
        min_improvement : float
            Minimum fractional improvement required 0.1=10% smaller mse.
        
        Setting
        -----
        Updates cache
        """
        cur_identifiers = self.cache.iter_identifiers[self.cache.iterations]
        cur_mse         = self.cache.iter_detection_mse[self.cache.iterations]
        self.log.debug(f"Stepping with current MSE ({cur_mse:.2f}), stepsize: {step_size}%")

        #locks order
        cur_identifiers_list = np.array(list(cur_identifiers))

        cur_ranks     = self.calculate_mean_rank_mse(cur_identifiers_list)
        max_kept_rank = np.percentile(cur_ranks, 100-step_size)
        identifiers_keep = cur_identifiers_list[cur_ranks<=max_kept_rank]
        identifiers_drop = cur_identifiers_list[cur_ranks>max_kept_rank]
        self.log.debug(f"Keeping ({len(identifiers_keep)}) identifiers")

        n_new_identifiers = self.catalogue_size - len(identifiers_keep) #ensure we dont drop a id from numerical stuff in max_kept_rank

        attempts = 0
        new_mse  = None
        new_identifiers = set()
        while attempts < max_attempts:
            attempt_candidates  = self.get_substitute_identifiers(n_substitutes=n_new_identifiers, retry_identifiers=True)
            attempt_identifiers =  np.append(identifiers_keep, list(attempt_candidates))
            attempt_mse         = self.calculate_detection_mse(attempt_identifiers)

            perc_change = (attempt_mse-cur_mse)/cur_mse
            self.log.debug(f"Attempt: ({attempts+1})/({max_attempts}) - MSE: {attempt_mse:.2f} - Change: {perc_change*100:.1f}%")

            if perc_change < -min_improvement: #The change in signed fraction must be greater than 
                new_mse = attempt_mse
                new_identifiers = attempt_identifiers
                break
            attempts += 1
        if new_mse:
            self.log.info(f"Succesful step: New_MSE: {new_mse:.2f}: Change: {perc_change*100:.1f}%")
            self.cache.iterations +=1
            self.cache.iter_identifiers.append(set(new_identifiers))
            self.cache.iter_detection_mse.append(new_mse)
            self.cache.identifiers_dropped.update(identifiers_drop)
        else:
            #todo. custome error
            err_msg = f"Failed improving with stepsize = {step_size}%, minimum improvement = {min_improvement}%"
            #self.log.error(err_msg)
            raise StepFailedError(err_msg)

    def run_refinement(self, step_sizes:List[int]=[40,35,30,25,20,15,10,5,2], retries:int=10, min_improvement:float=0.05) -> None:
        """
        Main function running the full improvement cycle.
        Will run through the step_sizes and try (retries) times to iterratively improve the detection MSE.

        Parameters
        -----
        step_sizes : List[int] = [40,35,30,25,20,15,10,5,2]
            List of step sizes (% worst identifiers to drop).
            Recommended large -> small 
        retries : int = 10
            Number of tries for each step_size
        min_improvement : float = 0.05
            Minimum relative MSE improvement for keeping a step. 

        Setting
        -----
        Updates cache
        """

        for ss in step_sizes:
            step_size_has_succeded = False
            for step_iter in range(retries):
                try:
                    self.refine_step(step_size=ss, max_attempts=retries, min_improvement=min_improvement)
                except StepFailedError:
                    pass
                else:
                    step_size_has_succeded= True
            if step_size_has_succeded:
                pass
            else:
                self.log.warning(f"stepsize: ({ss}%) Failed.")

        total_tries = len(step_sizes)*retries
        final_iter = self.cache.iterations
        final_mse  = self.cache.iter_detection_mse[final_iter]
        mse_improvement = ((final_mse - self.cache.iter_detection_mse[0]) / self.cache.iter_detection_mse[0])*100
        self.log.info(f"Succesful iterations: ({final_iter}/{total_tries}): Final MSE: ({final_mse:.2f}) Improvement: {mse_improvement:.2f}%")

    def generate_minimal_output(self, output: Path):
        output = Path(output)
        if not output.is_dir():
            raise FileNotFoundError(output)
        final_iter = self.cache.iterations
        final_identifiers = self.cache.iter_identifiers[final_iter]

        final_counts = self.df_counts.loc[list(final_identifiers), :]
        mses = self.cache.iter_detection_mse

        final_counts.to_csv(output/"counts_refined.tsv", sep="\t", index=True)
        
        pd.DataFrame({
            'iteration': list(range(final_iter+1)),
            'mse': mses
        }).to_csv(output/"mse.tsv", sep="\t", index=False)
        self.log.info(f"extracted counts (for refined identifiers) and mse overview -> {output}")

    def generate_identifier_csv(self, output_file: Path):
        df_kmers = pd.DataFrame(
            {
                i: list(self.cache.iter_identifiers[i])
                for i in range(self.cache.iterations+1)
            }
        ).rename(columns={0:"init", self.cache.iterations:"best"})
        df_kmers["random"] = self.rng_state.choice(self.df_counts.index, size=self.catalogue_size, replace=False)
        df_kmers.to_csv(output_file, index=False)

    def generate_detection_curve(self, identifiers) -> plotly.graph_objects.Figure:
        """
        """
        df_plot = self.df_counts[identifiers, :]
        fig = NB_Utils.plot_expected_detection_curve(df_plot)
        return fig
    
    def generate_verbose_output(self, output: Path) -> None:
        """
        """
        #Generate figure
        output = Path(output)
        if not output.is_dir():
            raise FileNotFoundError(output)
        final_iteration = self.cache.iterations

        initial_identifiers = self.cache.iter_identifiers[0]
        final_identifiers   = self.cache.iter_identifiers[final_iteration]

        random_identifiers = self.rng_state.choice(self.df_counts.index, size=self.catalogue_size, replace=False)
        
        fig_random = NB_Utils.plot_expected_detection_curve(self.df_counts.loc[list(random_identifiers),:].values)
        fig_raw    = NB_Utils.plot_expected_detection_curve(self.df_counts.values)
        fig_ini    = NB_Utils.plot_expected_detection_curve(self.df_counts.loc[list(initial_identifiers),:].values)
        fig_final  = NB_Utils.plot_expected_detection_curve(self.df_counts.loc[list(final_identifiers),:].values)

        fig_random.write_image(output/"detection_curve_random.png", scale=2)
        fig_raw.write_image(output/"detection_curve_raw.png", scale=2)
        fig_ini.write_image(output/"detection_curve_initial.png", scale=2)
        fig_final.write_image(output/"detection_curve_final.png", scale=2)

        #Dump cache
        dump_file = output/"pymag_cache.json"
        dump_file.write_text(json.dumps(self.cache.to_dict(), indent=True))

        self.log.info(f"Plots and cache dumped -> {output}")

        
    @staticmethod
    def refine_all(
        fp_count_files: Union[Path, List[Path]], 
        fp_meta_files: List[Path], 
        output: Path, 
        verbose_output=True, 
        kw_pymag: dict={}, 
        kw_refinement: dict={}) -> None:
        """
        Static method to run PyMAG refinement of catalogue for multiple catalogues.
        Instanciates the PyMAG class for each catalogue and runs refinement.

        Usage:
            PyMAG.refine_all(fp_count_files=[...], fp_meta_files=[...], ...)

        Parameters
        -----
        fp_count_files : Path or List[Path]
            Either a concatted countfile containing counts for all identifiers descriped in the metafiles
            or individual countfiles matching the fp_meta_files parameter.
            index = identifiers, cols = samples
        fp_meta_fields : List[Path]
            Note filenames from metafiles are used for subdirs - ensure uniqueness
            metafiles describing the identifiers - one per refined catalogue
            expected cols: identifier, prevalence, lenght
            NOTE: no headers.
        verbose_output : bool = True
            If False only refined identifiers are printed.
            If True plots and reports are generated.
        output : Path
            Top output directory, for each catalogue a subdir is created.
        kw_pymag:
            other key word arguments passed to PyMAG initialiser.
        kw_refinement:
            other keyword arguments passed to PyMAG run_refinement method.

        Returns
        -----
        None
            Output are written to output directory subfolders.
"""
        #Input control
        if not isinstance(fp_count_files, list):
            fp_count_files = [Path(fp_count_files) for x in range(len(fp_meta_files))]
        elif len(fp_count_files)==1:
            fp_count_files = [Path(fp_count_files[0]) for x in range(len(fp_meta_files))]
        else:
            fp_count_files = [Path(f) for f in fp_count_files]
        fp_meta_files = [Path(f) for f in fp_meta_files]
        output = Path(output)

        for file in fp_count_files+fp_meta_files:
            if not file.is_file():
                raise FileNotFoundError(file)
        
        if output.is_dir():
            logging.warn("Output directory already exists - may overwrite existing files. Consider cleaning first.")
        else:
            output.mkdir(exist_ok=True)

        run_overview = []

        for fp_count, fp_meta in zip(fp_count_files, fp_meta_files):
            name_catalogue = fp_meta.stem
            output_catalogue = output / "catalogues" /name_catalogue
            
            output_catalogue.mkdir(parents=True, exist_ok=True)
            output_identifier_overview = output / "screened_flat" / (name_catalogue+"_kmers.csv")
            output_identifier_overview.parent.mkdir(exist_ok=True)

            logfile_catalogue = output_catalogue / "log.txt"

            #initialization
            try:
                pymag = MAGpy(
                    fp_counts=fp_count, 
                    fp_meta=fp_meta, 
                    logfile=logfile_catalogue,
                    **kw_pymag)
            except Exception as err:
                run_overview.append(
                    {'catalogue':name_catalogue, 'status': 'Failed: Initialisation', 'error': err}
                )
                continue
            #Preflight
            try:
                pymag.preflight()
            except Exception as err:
                run_overview.append(
                    {'catalogue':name_catalogue, 'status': 'Failed: Preflight', 'error': err}
                )
                continue

            #Refinement
            try:
                pymag.run_refinement(**kw_refinement)
            except Exception as err:
                run_overview.append(
                    {'catalogue':name_catalogue, 'status': 'Failed: Refinement', 'error': err}
                )
                continue
            else:
                #minimum output:
                pymag.generate_minimal_output(output_catalogue)
                pymag.generate_identifier_csv(output_identifier_overview)

                if verbose_output:
                    pymag.generate_verbose_output(output_catalogue)
                run_overview.append(
                    {'catalogue':name_catalogue, 'status': 'Success: Output Generated', 'error': None}
                )
        time_now = datetime.datetime.now().strftime("%H:%M:%S")
        time_now_str = time_now.replace(":","")
        pd.DataFrame(
            run_overview
        ).sort_values("catalogue").to_csv(output/f"run_overview_{time_now_str}.tsv", sep="\t", index=False)

            
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("--count-files", required=True, type=Path, nargs="+")
    parser.add_argument("--meta-files", required=True, type=Path, nargs="+")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--full-output", action="store_true")

    pymag_init = parser.add_argument_group('Arguments for PyMAG initializer')
    pymag_init.add_argument("--max-threads", type=int, default=10)
    pymag_init.add_argument("--catalogue-size", type=int, default=500)
    pymag_init.add_argument("--verbosity", type=str, default="INFO")
    pymag_init.add_argument("--rng-seed", type=int, default=2112)

    pymag_refine = parser.add_argument_group('Arguments for PyMAG refiner')
    pymag_refine.add_argument("--step-sizes", type=int, nargs="+", default = [40,35,30,25,20,15,10,5,2])
    pymag_refine.add_argument("--retries", type=int, default=15)
    pymag_refine.add_argument("--min-improvement", type=float, default=0.05)

    args = parser.parse_args()

    kw_init = ["max_threads","catalogue_size","verbosity","rng_seed"]
    kw_refine = ["step_sizes","retries","min_improvement"]

    kw_init_args = {
        kw:getattr(args, kw)
        for kw in kw_init
    }
    kw_refine_args = {
        kw:getattr(args, kw)
        for kw in kw_refine
    }

    print("Number of catalogues: ",len(args.meta_files))

    MAGpy.refine_all(
        fp_count_files=args.count_files, 
        fp_meta_files=args.meta_files, 
        output=args.output, 
        verbose_output=args.full_output,
        kw_pymag=kw_init_args,
        kw_refinement=kw_refine_args
        )



        









        