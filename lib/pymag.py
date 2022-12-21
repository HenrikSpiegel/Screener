import logging
from pathlib import Path
import warnings

import pandas as pd
import numpy as np
from typing import List, Sequence, Set, Union
import plotly
import plotly.express as px
import plotly.graph_objects as go

import scipy
import statsmodels.api as sm
import multiprocessing


class NB_Utils:
    """
    Namespace class for static methods related to fitting count data to Negetive Binomials.
    """
    @staticmethod
    def expected_detection(G, N) -> float:
        """
        Trines formula (1)
        Given the total counts N and G different identifier then
        the probability of not observing a specific identifier is given by
        P0 = ((G-1)/G)**N, and thus for a given N the expected number of 
        observed genes are given by (1-P0)*G.

        Parameters
        ------
        G : int
            Number of different identifiers
        N : Int
            Number of total observations of all identifiers

        Returns
        ------
        Expected number of detected identifiers : float
        """
        P0 = ((G-1)/G)**N
        return (1-P0)*G

    @staticmethod
    def calculate_detection_mse(count_matrix: np.ndarray) -> float:
        """
        Trines Model implementation step (2)
        Calculates the MSE between the expected and observed realionship
        between distinct observed identifiers and total counts.

        Parameters
        ------
        count_matrix : np.ndarray
            matrix containing counts (rows=identifiers, columns=samples)
        
        Returns
        ------
        mse : float

        """
        # determine kmer detection:
        total_counts = count_matrix.sum(axis=0)
        observed_identifier = (count_matrix>0).sum(axis=0)

        G = count_matrix.shape[0]
        expected = [NB_Utils.expected_detection(G, Nj) for Nj in total_counts]

        mse = np.mean((observed_identifier-expected)**2)
        return mse

    @staticmethod
    def plot_expected_detection_curve(count_matrix: np.ndarray) -> plotly.graph_objs.Figure:
        """
        Calculates and plots the expected vs observed relationship 
        between total assigned counts and number of observed members
        of a single catalogue.
        
        Parameters:
        -----
        count_matrix : np.ndarray
            Array with rows of identifers and columns of samples
        
        Returns:
        -----
        fig : plotly.graph_objs.Figure
        
        """    
        total_counts = count_matrix.sum(axis=0)
        observed_kmers = (count_matrix>0).sum(axis=0)

        G = count_matrix.shape[0]
        expected = [NB_Utils.expected_detection(G, Nj) for Nj in total_counts]

        mse = np.mean((observed_kmers-expected)**2)
        
        total_range = list(range(0, total_counts.max()))
        expected = [NB_Utils.expected_detection(G, Nj) for Nj in total_range]
        fig = px.scatter(
            color=["Observed" for x in total_counts],
            x=total_counts,
            y=observed_kmers,
            log_x=True,
            title = f'Relationship between number of observed members and total counts assigned to catalogue<br><sup> MSE (expected vs. observed): {mse:.2F}',
            labels = {
                'x': 'Total counts assigned to catalogue',
                'y': 'Total number of observed catalogue members',
                'color':''
            }
        )
        fig.add_trace(
            go.Scatter(
                x=total_range,
                y=expected,
                name="Expected"
            )
        )
        return fig

    @staticmethod
    def tukeys_outlier_detection(counts: np.ndarray, tukeys_const: float = 1.5):
        """
        Perform outlier detection based on tukeys IQR method.
        Outliers are defined as points outside of the q1 or q3 by iqr*tukeys_const or more.
        tukeys_range = (q3-q1)*tukeys_const
        outlier = (counts < (q1-tukeys_range)) | ((q3+tukeys_range) < counts)
        
        Parameters
        -----
        tukeys_const : float
            Weight for IQR (default is 1.5 which is also default for boxplot whiskers)
        
        Returns
        -----
        is_outlier : List[bool]
            list of boolean, true = outlier.
        """
        q1, q3 = np.percentile(counts, [25,75])
        tukeys_range = (q3-q1)*tukeys_const
        outlier = (counts < (q1-tukeys_range)) | ((q3+tukeys_range) < counts)
        return outlier

    @staticmethod
    def counts2ranked_pearson(endog: List[int], model="NB") -> List[int]:
        """
        Fit endog (counts) and get pearson residual per entry in endog.
        
        Parameters
        -----
        endog : List[int]
            count values
        
        Returns
        -----
        ranked_pearson : List[float]
            Ranked pearson correlation (matches order of input.)
        has_converged : bool
            Whether the fitter has converged. 
        """
        if model == "NB":
            m = sm.NegativeBinomial(endog, np.ones_like(endog))
        else:
            raise NotImplementedError("Currently only simply 'NB' is implemented")
        nobs = len(endog)
        m_fitted = m.fit(disp=0)

        statsmodels_major_version = int(sm.__version__.split(".")[1])
        if statsmodels_major_version >= 14:
            expec_bincount = m_fitted.predict(which="prob").mean(0)*nobs
            has_converged = m_fitted.converged
        else:
            #Distribution prediction is not included for all models in statsmodels 13.*
            mu = np.exp(m_fitted.params[0])
            over_disp = m_fitted.params[1]
            p = 1/(1+mu*over_disp)
            n = mu*p/(1-p)
            bin_range = list(range(0, endog.max()+1))
            expec_bin_freq = scipy.stats.nbinom.pmf(bin_range, n, p)
            expec_bincount = expec_bin_freq*nobs
            has_converged = m_fitted._results.__dict__["mle_retvals"]["converged"]
            
        
        obs_bincount = np.bincount(endog)

        pearson_res_per_bin = (obs_bincount-expec_bincount)/np.sqrt(expec_bincount)
        pearson_res_per_identifer = pearson_res_per_bin[endog] #(expand from bins to obs)

        rank_pearson_res_per_identifier = scipy.stats.rankdata(pearson_res_per_identifer)
        return rank_pearson_res_per_identifier, has_converged

from dataclasses import dataclass, field

@dataclass
class PyMAG_Cache:
    """
    Cache dataclass for PyMAG class.
    """
    iterations:   int       = field(default=0)
    iter_identifiers: List[Set[str]] = field(default_factory=list)
    iter_detection_mse: List[float] = field(default_factory=list)

    identifiers_dropped:      Set[str] = field(default_factory=set)
    identifiers_outliers:   Set[str] =  field(default_factory=set)

    
    def __repr__(self):
        dict_repr = ', '.join(
            f'{k}:{type(v)}'
            for k, v in filter(
                lambda item: not item[0].startswith('_'),
                self.__dict__.items()
            )
        )

        return f'{self.__class__.__name__}({dict_repr})'


class PyMAG:
    """
    Implementation of NegBinom optimisation of kmer catalogue.
    The implementation is based on the original work MAGinator by Trine Zachariasen.

    The class implements an catalogue refiner which attempts take a count matrix with rows
    of identifiers and columns of samples and determine a subset of identifiers which have
    the best fit to a NegativeBinomial distribution.

    ...

    Attributes
    ----------
    df_counts : pd.DataFrame
        Dataframe containing the counts from fp_counts subset to identifier in fp_meta
    
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
        return dict(
                name = self.__class__.__name__,
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
            self._cache = PyMAG_Cache()
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
            self.log.warning(f"Following ids where not found in count data (Allowed)) -> {missing_meta}")

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
            rank_matrices, has_converged =  zip(*p.map(NB_Utils.counts2ranked_pearson, sample_list))
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

    def preflight(self, max_retries: int = 10) -> None:
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

        self.df_meta["identifier_desireability"] = self.df_meta["prevalence"] - self.df_meta["outlier_degree"]

        self.cache.identifiers_outliers = identifiers_outliers
        
        df_meta_eligible = self.df_meta.loc[~self.df_meta["identifier"].isin(identifiers_outliers)]

        #Generate 10 starting identifier compositions.
        starter_compositions = []
        starter_compositions_mse = []
        starter_compositions_ranks = []
        retries = 0
        while (len(starter_compositions) < 10) and (retries < max_retries):
            df_random_subset = df_meta_eligible.sample(n=self.catalogue_size, weights="identifier_desireability", replace=False, random_state=self.rng_state)
            identifiers_random = df_random_subset.identifier.values

            #with warnings.catch_warnings():
            #warnings.filterwarnings('once')
            try:
                ranks_identifier_random = self.calculate_mean_rank_mse(identifiers=identifiers_random)
                identifier_random_detection_mse = self.calculate_detection_mse(identifiers=identifiers_random)
                #NB_Utils.calculate_detection_mse(self.df_counts.loc[identifiers_random,:])
            # except Warning as we:
            #     retries += 1
            #     self.log.exception(f"Warning from fitter, re-sampling: {retries}/{max_retries}")
            #     pass
            except Exception as ee:
                retries += 1
                self.log.exception(f"Exception during fitting, re-sampling: {retries}/{max_retries}")
                pass
            else:
                starter_compositions.append(identifiers_random)
                starter_compositions_mse.append(identifier_random_detection_mse)
                starter_compositions_ranks.append(ranks_identifier_random)
                self.log.debug(f"Added a valid starter composition {len(starter_compositions)}/10")
            #finally:
            #    warnings.filterwarnings()
        best_iteration = np.argmin(starter_compositions_mse)

        self.log.info(f"Initial identifiers found. Initial MSE: {starter_compositions_mse[best_iteration]:.2f}")

        self.cache.iterations = 0
        self.cache.iter_identifiers.append(set(starter_compositions[best_iteration]))
        self.cache.iter_detection_mse.append(starter_compositions_mse[best_iteration])

    def _get_substitute_identifiers(self, n_substitutes:int, black_list:Set[str]= set(),  retry_identifiers:bool=True):
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
        """
        cur_identifiers = self.cache.iter_identifiers[self.cache.iterations]
        cur_mse         = self.cache.iter_detection_mse[self.cache.iterations]
        self.log.debug(f"Stepping with current MSE ({cur_mse:.2f}), stepsize: {step_size}%")

        #locks order
        cur_identifiers_list = np.array(list(cur_identifiers))

        cur_ranks     = self.calculate_mean_rank_mse(cur_identifiers_list)
        max_kept_rank = np.percentile(cur_ranks, 100-step_size)
        identifiers_keep = cur_identifiers_list[cur_ranks<=max_kept_rank]
        self.log.debug(f"Keeping ({len(identifiers_keep)}) identifiers")

        n_new_identifiers = self.catalogue_size - len(identifiers_keep) #ensure we dont drop a id from numerical stuff in max_kept_rank

        attempts = 0
        new_mse  = None
        new_identifiers = set()
        while attempts < max_attempts:
            attempt_candidates  = self._get_substitute_identifiers(n_substitutes=n_new_identifiers, retry_identifiers=True)
            attempt_identifiers =  np.append(identifiers_keep, list(attempt_candidates))
            attempt_mse         = self.calculate_detection_mse(attempt_identifiers)
            self.log.debug(f"Attempt: ({attempts+1})/({max_attempts}) - MSE: {attempt_mse:.2f}")

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
        else:
            #todo. custome error
            raise RuntimeError("Failed stepping")






    
        



        









        