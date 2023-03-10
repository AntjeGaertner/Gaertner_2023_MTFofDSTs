import numpy as np
import pandas as pd
import statsmodels.stats.api as sms
import math
import scipy.optimize as opt
import sys

from classes import snag_MTF_func as sg


   
    
def ODR_linear(x, y, min_x = 0, max_x = 100):
    import scipy.odr as odr
    from scipy import stats

    def f(B, x):
        '''Linear function y = m*x + b'''
        # B is a vector of the parameters.
        # x is an array of the current x values.
        # x is in the same format as the x passed to Data or RealData.
        #
        # Return an array in the same format as y passed to Data or RealData.
        return B[0]*x + B[1]

    linear = odr.Model(f)

    mydata = odr.Data(x, y, wd=1./np.power(np.std(x),2), we=1./np.power(np.std(y),2))

    myodr = odr.ODR(mydata, linear, beta0=[0, 2.])

    myoutput = myodr.run()

    ss_res = myoutput.sum_square
    ss_tot = np.sum((y-np.mean(y))**2)
    r_squared = 1 - (ss_res / ss_tot)
    #print('R2=',np.round(r_squared,2))

    N = len(y)

    # Degrees of freedom are n_samples - n_parameters
    df = N - 2  # equivalently, df = odr_out.iwork[10]
    beta_0 = 0  # test if slope is significantly different from zero
    t_stat = (myoutput.beta[0] - beta_0) / myoutput.sd_beta[0]  # t statistic for the slope parameter
    p_val = stats.t.sf(np.abs(t_stat), df) * 2
    #print('Recovered equation: y={:3.2f}x + {:3.2f}, t={:3.2f}, p={:.2e}'.format(myoutput.beta[0],
    #                                                                             myoutput.beta[1],
    #                                                                             t_stat, p_val))
    beta = myoutput.beta
    X = np.linspace(min_x, max_x, max_x+np.abs(min_x))
    Y = f(beta, X)

    return (X, Y, beta, df, r_squared, t_stat, p_val)


def draw_bs_pairs_odr(x, y, size=1, min_x = 0, max_x = 100):
    """
    #--- Perform pairs bootstrap for orthogonal distance regression
    
    #- Input
    x        -> x data, original
    y        -> y data, original
    size     -> number of bootstrap replicates
    min_x    -> minimum X used for Y data generation
    max_x    -> maximum X used for Y data generation
    
    
    #- Output
    X        -> np.array of all X values in the range of min_x to max_x
    Y        -> np.array of all Y values in the range of min_x to max_x
    """

    # Set up array of indices to sample from: inds
    inds = np.arange(len(x))

    # Initialize replicates: bs_slope_reps, bs_intercept_reps
    X         = np.empty((size, max_x+np.abs(min_x)))
    Y         = np.empty((size, max_x+np.abs(min_x)))
    beta      = np.empty((size, 2))
    df        = np.empty(size)
    r_squared = np.empty(size)
    t_stat    = np.empty(size)
    p_val     = np.empty(size)

    # Generate replicates
    for i in range(size):
        bs_inds = np.random.choice(inds, len(inds))
        bs_x, bs_y = x[bs_inds], y[bs_inds]
        X[i], Y[i], beta[i], df[i], r_squared[i], t_stat[i], p_val[i] = ODR_linear(bs_x, bs_y, min_x, max_x)

    return X, Y
    
def draw_bs_pairs_linear(x, y, size=1, min_x = 0, max_x = 100):
    """
    #--- Perform pairs bootstrap for linear regression
    
    #- Input
    x        -> x data, original
    y        -> y data, original
    size     -> number of bootstrap replicates
    min_x    -> minimum X used for Y data generation
    max_x    -> maximum X used for Y data generation
    
    
    #- Output
    X        -> np.array of all X values in the range of min_x to max_x
    Y        -> np.array of all Y values in the range of min_x to max_x
    """

    # Set up array of indices to sample from: inds
    inds = np.arange(len(x))

    # Initialize replicates: bs_slope_reps, bs_intercept_reps
    X         = np.empty((size, max_x+np.abs(min_x)))
    Y         = np.empty((size, max_x+np.abs(min_x)))
    beta      = np.empty((size, 2))
    df        = np.empty(size)
    r_squared = np.empty(size)
    t_stat    = np.empty(size)
    p_val     = np.empty(size)

    # Generate replicates
    for i in range(size):
        bs_inds = np.random.choice(inds, len(inds))
        bs_x, bs_y = x[bs_inds], y[bs_inds]
        
        popt, pcov = opt.curve_fit(sg.mod_linear_origin, bs_x, bs_y,
                                   maxfev= 100000
                                   )
        X[i] = np.linspace(min_x, max_x, max_x+np.abs(min_x))
        Y[i] = sg.mod_linear_origin(X[i], *popt)
    return X, Y


# base code
import numpy as np
import seaborn as sns
from statsmodels.tools.tools import maybe_unwrap_results
from statsmodels.graphics.gofplots import ProbPlot
from statsmodels.stats.outliers_influence import variance_inflation_factor
import matplotlib.pyplot as plt
from typing import Type
import statsmodels

style_talk = 'seaborn-talk'    #refer to plt.style.available

class Linear_Reg_Diagnostic():
    """
    Diagnostic plots to identify potential problems in a linear regression fit.
    Mainly,
        a. non-linearity of data
        b. Correlation of error terms
        c. non-constant variance
        d. outliers
        e. high-leverage points
        f. collinearity

    Author:
        Prajwal Kafle (p33ajkafle@gmail.com, where 3 = r)
        Does not come with any sort of warranty.
        Please test the code one your end before using.
    """

    def __init__(self,
                 results: Type[statsmodels.regression.linear_model.RegressionResultsWrapper]) -> None:
        """
        For a linear regression model, generates following diagnostic plots:

        a. residual
        b. qq
        c. scale location and
        d. leverage

        and a table

        e. vif

        Args:
            results (Type[statsmodels.regression.linear_model.RegressionResultsWrapper]):
                must be instance of statsmodels.regression.linear_model object

        Raises:
            TypeError: if instance does not belong to above object

        Example:
        >>> import numpy as np
        >>> import pandas as pd
        >>> import statsmodels.formula.api as smf
        >>> x = np.linspace(-np.pi, np.pi, 100)
        >>> y = 3*x + 8 + np.random.normal(0,1, 100)
        >>> df = pd.DataFrame({'x':x, 'y':y})
        >>> res = smf.ols(formula= "y ~ x", data=df).fit()
        >>> cls = Linear_Reg_Diagnostic(res)
        >>> cls(plot_context="seaborn-paper")

        In case you do not need all plots you can also independently make an individual plot/table
        in following ways

        >>> cls = Linear_Reg_Diagnostic(res)
        >>> cls.residual_plot()
        >>> cls.qq_plot()
        >>> cls.scale_location_plot()
        >>> cls.leverage_plot()
        >>> cls.vif_table()
        """

        if isinstance(results, statsmodels.regression.linear_model.RegressionResultsWrapper) is False:
            raise TypeError("result must be instance of statsmodels.regression.linear_model.RegressionResultsWrapper object")

        self.results = maybe_unwrap_results(results)

        self.y_true = self.results.model.endog
        self.y_predict = self.results.fittedvalues
        self.xvar = self.results.model.exog
        self.xvar_names = self.results.model.exog_names

        self.residual = np.array(self.results.resid)
        influence = self.results.get_influence()
        self.residual_norm = influence.resid_studentized_internal
        self.leverage = influence.hat_matrix_diag
        self.cooks_distance = influence.cooks_distance[0]
        self.nparams = len(self.results.params)

    def __call__(self, plot_context='seaborn-paper'):
        # print(plt.style.available)
        with plt.style.context(plot_context):
            fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10,10))
            self.residual_plot(ax=ax[0,0])
            self.qq_plot(ax=ax[0,1])
            self.scale_location_plot(ax=ax[1,0])
            self.leverage_plot(ax=ax[1,1])
            plt.show()

        self.vif_table()
        return fig, ax


    def residual_plot(self, ax=None):
        """
        Residual vs Fitted Plot

        Graphical tool to identify non-linearity.
        (Roughly) Horizontal red line is an indicator that the residual has a linear pattern
        """
        if ax is None:
            fig, ax = plt.subplots()

        sns.residplot(
            x=self.y_predict,
            y=self.residual,
            lowess=True,
            scatter_kws={'alpha': 0.5},
            line_kws={'color': 'red', 'lw': 1, 'alpha': 0.8},
            ax=ax)

        # annotations
        residual_abs = np.abs(self.residual)
        abs_resid = np.flip(np.sort(residual_abs))
        abs_resid_top_3 = abs_resid[:3]
        for i, _ in enumerate(abs_resid_top_3):
            ax.annotate(
                i,
                xy=(self.y_predict[i], self.residual[i]),
                color='C3')

        ax.set_title('Residuals vs Fitted', fontweight="bold")
        ax.set_xlabel('Fitted values')
        ax.set_ylabel('Residuals')
        return ax

    def qq_plot(self, ax=None):
        """
        Standarized Residual vs Theoretical Quantile plot

        Used to visually check if residuals are normally distributed.
        Points spread along the diagonal line will suggest so.
        """
        if ax is None:
            fig, ax = plt.subplots()

        QQ = ProbPlot(self.residual_norm)
        QQ.qqplot(line='45', alpha=0.5, lw=1, ax=ax)

        # annotations
        abs_norm_resid = np.flip(np.argsort(np.abs(self.residual_norm)), 0)
        abs_norm_resid_top_3 = abs_norm_resid[:3]
        for r, i in enumerate(abs_norm_resid_top_3):
            ax.annotate(
                i,
                xy=(np.flip(QQ.theoretical_quantiles, 0)[r], self.residual_norm[i]),
                ha='right', color='C3')

        ax.set_title('Normal Q-Q', fontweight="bold")
        ax.set_xlabel('Theoretical Quantiles')
        ax.set_ylabel('Standardized Residuals')
        return ax

    def scale_location_plot(self, ax=None):
        """
        Sqrt(Standarized Residual) vs Fitted values plot

        Used to check homoscedasticity of the residuals.
        Horizontal line will suggest so.
        """
        if ax is None:
            fig, ax = plt.subplots()

        residual_norm_abs_sqrt = np.sqrt(np.abs(self.residual_norm))

        ax.scatter(self.y_predict, residual_norm_abs_sqrt, alpha=0.5);
        sns.regplot(
            x=self.y_predict,
            y=residual_norm_abs_sqrt,
            scatter=False, ci=False,
            lowess=True,
            line_kws={'color': 'red', 'lw': 1, 'alpha': 0.8},
            ax=ax)

        # annotations
        abs_sq_norm_resid = np.flip(np.argsort(residual_norm_abs_sqrt), 0)
        abs_sq_norm_resid_top_3 = abs_sq_norm_resid[:3]
        for i in abs_sq_norm_resid_top_3:
            ax.annotate(
                i,
                xy=(self.y_predict[i], residual_norm_abs_sqrt[i]),
                color='C3')
        ax.set_title('Scale-Location', fontweight="bold")
        ax.set_xlabel('Fitted values')
        ax.set_ylabel(r'$\sqrt{|\mathrm{Standardized\ Residuals}|}$');
        return ax

    def leverage_plot(self, ax=None):
        """
        Residual vs Leverage plot

        Points falling outside Cook's distance curves are considered observation that can sway the fit
        aka are influential.
        Good to have none outside the curves.
        """
        if ax is None:
            fig, ax = plt.subplots()

        ax.scatter(
            self.leverage,
            self.residual_norm,
            alpha=0.5);

        sns.regplot(
            x=self.leverage,
            y=self.residual_norm,
            scatter=False,
            ci=False,
            lowess=True,
            line_kws={'color': 'red', 'lw': 1, 'alpha': 0.8},
            ax=ax)

        # annotations
        leverage_top_3 = np.flip(np.argsort(self.cooks_distance), 0)[:3]
        for i in leverage_top_3:
            ax.annotate(
                i,
                xy=(self.leverage[i], self.residual_norm[i]),
                color = 'C3')

        xtemp, ytemp = self.__cooks_dist_line(0.5) # 0.5 line
        ax.plot(xtemp, ytemp, label="Cook's distance", lw=1, ls='--', color='red')
        xtemp, ytemp = self.__cooks_dist_line(1) # 1 line
        ax.plot(xtemp, ytemp, lw=1, ls='--', color='red')

        ax.set_xlim(0, max(self.leverage)+0.01)
        ax.set_title('Residuals vs Leverage', fontweight="bold")
        ax.set_xlabel('Leverage')
        ax.set_ylabel('Standardized Residuals')
        ax.legend(loc='upper right')
        return ax

    def vif_table(self):
        """
        VIF table

        VIF, the variance inflation factor, is a measure of multicollinearity.
        VIF > 5 for a variable indicates that it is highly collinear with the
        other input variables.
        """
        vif_df = pd.DataFrame()
        vif_df["Features"] = self.xvar_names
        vif_df["VIF Factor"] = [variance_inflation_factor(self.xvar, i) for i in range(self.xvar.shape[1])]

        print(vif_df
                .sort_values("VIF Factor")
                .round(2))
                
        
    def return_vif_table(self):
        """
        VIF table

        VIF, the variance inflation factor, is a measure of multicollinearity.
        VIF > 5 for a variable indicates that it is highly collinear with the
        other input variables.
        """
        vif_df = pd.DataFrame()
        vif_df["Features"] = self.xvar_names
        vif_df["VIF Factor"] = [variance_inflation_factor(self.xvar, i) for i in range(self.xvar.shape[1])]

        return vif_df


    def __cooks_dist_line(self, factor):
        """
        Helper function for plotting Cook's distance curves
        """
        p = self.nparams
        formula = lambda x: np.sqrt((factor * p * (1 - x)) / x)
        x = np.linspace(0.001, max(self.leverage), 50)
        y = formula(x)
        return x, y
        
def calc_pixels_cm_tiff_export_preview(dpi):
    """
    Calculate the pixels per cm to reach a specific dpi when exporting tiff images from preview.
    
    Also for purpose of exporting images / graphics from power point.
        1. Save slide / presentation as postscript (.ps) by pressing command + p
        2. open <file>.ps in preview and export as tiff (croppig is preserved in tiffs from preview)
        3. define pixels / cm in the export dialogue
    """
    
    dpi_factor_cm = 0.4 # this is approximate
    
    pixels_per_cm = dpi * dpi_factor_cm
    
    return pixels_per_cm

import itertools
from itertools import chain
from itertools import combinations
def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


def unique_models_fixed_interaction(var_group, interaction=False):
    
    var_group_mixed = var_group.copy()
    
    if interaction:
        for word, vals in var_group.items(): # Including interaction effects
            vals_new = [i+' + ' for i in vals] + [i+' * ' for i in vals]
            var_group_mixed[word] = vals_new
    else:
        for word, vals in var_group.items(): # Only fixed effects
            vals_new = [i+' + ' for i in vals]
            var_group_mixed[word] = vals_new

    keys, values = zip(*var_group_mixed.items())
    permutations_list = [list(v) for v in itertools.product(*values)]
    
    set_uni_model_list = []

    for n,l in enumerate(permutations_list):
        set1 = [list(i) for i in powerset(l)][1:]

        all_effects = [''.join(i) for i in set1]

        set_uni_model = list(np.unique(np.array(all_effects)))

        set_uni_model_list = set_uni_model_list + set_uni_model

    uni_model_list = list(np.unique(np.array(set_uni_model_list)))  
    
    for i,j in enumerate(uni_model_list):
        #print('\n',j)
        if type(j) != list:
            j = j[:-3]
            uni_model_list[i] = j
        else:
            #j = np.array(j)
            j[-1] = j[-1][:-3]
            uni_model_list[i] = j

    uni_model_list = list(np.unique(np.array(uni_model_list))) 
    
    return uni_model_list


