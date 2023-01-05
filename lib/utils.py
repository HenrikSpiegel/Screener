import numpy as np
from typing import List, Sequence, Set, Union
import plotly
import plotly.express as px
import plotly.graph_objects as go
from statsmodels.tools.sm_exceptions import HessianInversionWarning, ConvergenceWarning
import warnings

import scipy
import statsmodels.api as sm


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
            title = f'Relationship between unique and total member count.<br><sup> MSE (expected vs. observed): {mse:.2F}',
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
        fit_warning : bool
            Whether a warning was raised during fitting. 
        """
        if model == "NB":
            m = sm.NegativeBinomial(endog, np.ones_like(endog))
        else:
            raise NotImplementedError("Currently only simply 'NB' is implemented")
        nobs = len(endog)

        #Note we ignore the warnings to reduce clutter in stderr but handle them below by returning a flag.
        
        warnings.simplefilter('ignore',category=HessianInversionWarning)
        warnings.simplefilter('ignore',category=ConvergenceWarning)
        warnings.simplefilter('ignore',category=RuntimeWarning)
        try:
            m_fitted = m.fit(disp=0)
        except Exception:
            return np.zeros_like(endog), True

        warnings.simplefilter('default',category=HessianInversionWarning)
        warnings.simplefilter('default',category=ConvergenceWarning)
        warnings.simplefilter('default',category=RuntimeWarning)

        hessianfail = m_fitted._results.__dict__["normalized_cov_params"] is None
        is_converged = m_fitted._results.__dict__["mle_retvals"]["converged"]

        failed_fit =  hessianfail or (not is_converged)

        statsmodels_major_version = int(sm.__version__.split(".")[1])
        if statsmodels_major_version >= 14:
            expec_bincount = m_fitted.predict(which="prob").mean(0)*nobs

        else:
            #Distribution prediction is not included for all models in statsmodels 13.*
            mu = np.exp(m_fitted.params[0])
            over_disp = m_fitted.params[1]
            p = 1/(1+mu*over_disp)
            n = mu*p/(1-p)
            bin_range = list(range(0, endog.max()+1))
            expec_bin_freq = scipy.stats.nbinom.pmf(bin_range, n, p)
            expec_bincount = expec_bin_freq*nobs

        obs_bincount = np.bincount(endog)

        pearson_res_per_bin = (obs_bincount-expec_bincount)/np.sqrt(expec_bincount)
        pearson_res_per_identifer = pearson_res_per_bin[endog] #(expand from bins to obs)

        rank_pearson_res_per_identifier = scipy.stats.rankdata(np.abs(pearson_res_per_identifer))
        return rank_pearson_res_per_identifier, failed_fit