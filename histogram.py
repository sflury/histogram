'''
dependencies:
    numpy 1.17.3
    scipy 1.4.1

author:
    Sophia Flury 2022.06.07
'''
import numpy as np
from scipy.special import ndtr,ndtri,erf,gammaincc,gammaln,iv
from scipy.optimize import newton
'''
Name
    bernoulli_rvs

Purpose
    Sample Bernoulli random variates of a given shape for a given probability of
    success.
'''
def bernoulli_rvs(p,n_samp=int(1e4)):
    # sample Bernoulli variates using uniform random variates
    rvs_u = np.random.uniform(size=n_samp)
    # assume no successful trials at first
    rvs_b = np.zeros(rvs_u.shape)
    # if rvs_u < p, then bin trial is successful
    rvs_b[np.where(rvs_u < p)] = 1.
    return rvs_b
'''
Name
    poisson_gehrels

Purpose
    Compute the poisson confidence intervals for counts n. Solves Gehrels (1986)
    eqns 1 and 2 for lower and upper limits on n (p_l and p_u, respectively) by
    root-finding.

Arguments
    :n (*float*): total counts
    :p_conf (*float* or *np.ndarray*): two-sided probability intervals for which
                            to compute confidence intervals

Returns:
    :conf_l (*np.ndarray*): lower confidences on p1=n1/n corresponding to p_conf
    :conf_u (*np.ndarray*): upper confidences on p1=n1/n corresponding to p_conf
'''
def poisson_gehrels(n,p_conf):
    # check for array
    if not hasattr(p_conf,'__len__'):
        p_conf = np.array([p_conf])
    # upper quantiles from two-sided confidence intervals
    q_u = 0.5+0.5*p_conf
    # starting arrays for lower and upper confidences in p
    conf_l = np.zeros(p_conf.shape) # if n=0, then lower confidence is zero
    conf_u = 1.841*np.ones(p_conf.shape)  # if n=0, then upper confidence is 1.841
    # if n1 = 0, then pl = 0
    if n != 0 :
        # all counts from 0 to n1-1
        x = np.arange(0,n,1)//1
        # for each interval
        for i in range(len(p_conf)):
            # define function taking difference of confidence limit and
            # Gehrels (1986) Eqn 2
            func = lambda p_l: np.sum(np.exp(x*np.log(p_l)-(p_l+gammaln(x+1))))-q_u[i]
            # set equal to zero and solve for pl using root-finding
            p_l = newton(func,max([n-np.sqrt(n),0.5]),tol=2**-32)
            if p_l > 0:
                conf_l[i] = p_l
    # all counts from 0 to n1
    x = np.arange(0,n+1,1)//1
    # for each interval
    for i in range(len(p_conf)):
        # define function taking difference of confidence limit and
        # Gehrels (1986) Eqn 1
        func = lambda p_u: np.sum(np.exp(x*np.log(p_u)-(p_u+gammaln(x+1))))-1+q_u[i]
        # set equal to zero and solve for pl using root-finding
        conf_u[i] = newton(func,max([n+np.sqrt(n),2]),tol=2**-32)
    return conf_l,conf_u
'''
Name
    binom_gehrels

Purpose
    Compute the binomial confidence intervals for a binomial distribution for
    counts n1 that are some fraction p1 of a total number of counts n=n1+n2
    such that p1 = n1/n. Solves Gehrels (1986) eqns 16 and 17 for lower and
    upper limits on p1 (p_l and p_u, respectively) by root-finding.

Arguments
    :n1 (*float*): number of counts of subset of n
    :n (*float*): total counts
    :p_conf (*float* or *np.ndarray*): two-sided probability intervals for which
                            to compute confidence intervals

Returns:
    :conf_l (*np.ndarray*): lower confidences on p1=n1/n corresponding to p_conf
    :conf_u (*np.ndarray*): upper confidences on p1=n1/n corresponding to p_conf
'''
def binom_gehrels(n1,n,p_conf):
    # fraction p1 = n1/(n1+n2)
    p1 = n1/n
    # check for array
    if not hasattr(p_conf,'__len__'):
        p_conf = np.array([p_conf])
    # upper quantiles from confidence intervals
    q_u = 0.5+0.5*p_conf
    # starting arrays for lower and upper confidences in p
    conf_l = np.zeros(p_conf.shape) # if n1=0, then lower confidence is zero
    conf_u = np.ones(p_conf.shape)  # if n1=n, then upper confidence is one
    # if n1 = 0, then pl = 0
    if n1 != 0:
        # all counts from 0 to n1-1
        x = np.arange(0,n1,1)
        # binomial coefficient n choose x = n!/(x!(n-x)!)
        # and k! = Gamma(k+1), using ln(Gamma) to prevent overflow
        bincoef = np.exp(gammaln(n+1)-gammaln(x+1)-gammaln(n-x+1))
        # for each interval
        for i in range(len(p_conf)):
            # define function taking difference of confidence limit and
            # Gehrels (1986) Eqn 17
            func = lambda p_l: np.sum(bincoef*p_l**x*(1-p_l)**(n-x))-q_u[i]
            # set equal to zero and solve for pl using root-finding
            p_l = newton(func,p1,tol=2**-50)
            if p_l > 0:
                conf_l[i] = p_l
    # if n1 = n, then pu = 1
    if n1 != n:
        # all counts from 0 to n1
        x = np.arange(0,n1+1,1)
        # binomial coefficient n choose x = n!/(x!(n-x)!)
        # and k! = Gamma(k+1), using ln(Gamma) to prevent overflow
        bincoef = np.exp(gammaln(n+1)-gammaln(x+1)-gammaln(n-x+1))
        # for each interval
        for i in range(len(p_conf)):
            # define function taking difference of confidence limit and
            # Gehrels (1986) Eqn 16
            func = lambda p_u: np.sum(bincoef*p_u**x*(1-p_u)**(n-x))-1+q_u[i]
            # set equal to zero and solve for pl using root-finding
            p_u = newton(func,p1,tol=2**-50)
            if p_u < 1.:
                conf_u[i] = p_u
    return conf_l,conf_u

'''
Name
    hist_error

Purpose
    Compute the "true" histogram with uncertainties given a set of measurements,
    their uncertainties, and bins.

Arguments
    :x (*np.ndarray*): array of measurements

Returns
    :cts (*np.ndarray*): the "counts" of x in each bin
    :center (*np.ndarray*): the center of each bin
    :cts_lower (*np.ndarray*): the lower uncertainty in cts
    :cts_upper (*np.ndarray*): the upper uncertainty in cts
'''
# function to get error in histogram bin counts
#
def hist_error(x,x_err=None,bins=None):
    # define bins if not provided
    if np.any(bins == None):
        bins = np.linspace(min(x)//1,max(x)//1+1,int(np.sqrt(len(x))))
    # number of bins
    n_bins = int(len(bins)-1)
    # Poisson binomial method if errors provided
    if not np.any(x_err == None):
        # normal PMF for variable mean and standard deviation
        pmf = lambda x,mu,sig: np.exp(-0.5*((x-mu)/sig)**2)/np.sqrt(2*np.pi)
        pmf_ij = np.array([pmf(bin_i,x,x_err) for bin_i in bins])
        # normal CDF for variable mean and standard deviation
        cdf = lambda x,mu,sig: 0.5+0.5*erf((x-mu)/(abs(sig)*np.sqrt(2)))
        # array of probabilities p_ij that jth measurement x landed in ith bin
        # given the ith bin limits [b_i-1,b_i] and the uncertainty s_j in x_j
        # following P(b_i-1 < x_j < b_i) = cdf(x_j|b_i,s_j) - cdf(x_j|b_i-1,s_j)
        cdf_ij = np.array([cdf(bin_i,x,x_err) for bin_i in bins])
        prb_ij = cdf_ij[1:]-cdf_ij[:-1]
        # array of bin edge expected values
        # from integral u * pmf(u|x,x_err) du from -inf to bin_edge
        cent_ij = np.array([np.sqrt(2*np.pi)*x*cdf_ij[i]-x_err*pmf_ij[i]\
                            for i in range(n_bins+1)])/np.sqrt(2*np.pi)
        # bin center is integral u * pmf(u|x,x_err) du from bin_low to bin_up,
        # i.e., difference of bin edge expected values, normalized by bin counts
        cent = np.sum(cent_ij[1:]-cent_ij[:-1],axis=1)/np.sum(prb_ij,axis=1)
        # check for iffy values
        prb_ij[np.where(prb_ij<0)] = 0
        prb_ij[np.where(prb_ij>1)] = 1
        prb_ij[np.where(np.isnan(prb_ij))] = 0.5
        # mean counts is sum of probabilities
        counts = prb_ij.sum(axis=1)
        # variance is sum of p*(1-p)
        var = np.sum(prb_ij*(1-prb_ij),axis=1)
        p_lo = counts-np.sqrt(var)
        p_up = counts+np.sqrt(var)
    # Gehrels (1986) method if no errors provided
    else:
        # histogram
        counts,bins = np.nanhistogram(x,bins=bins)
        # bin centers
        cent = (bins[:-1]+bins[1:])/2.
        # confidences
        p_lo = np.zeros(n_bins)
        p_up = np.zeros(n_bins)
        for i in range(n_bins):
            p_lo[i],p_up[i] = poisson_gehrels(counts[i]//1,0.68268949)
    # return
    return counts,cent,counts-p_lo,p_up-counts
'''
Name
    calc_uniform_bins

Purpose
    Determine the edges for bins which are uniformly populated at confidence
    interval specified according to the normal cumulative distribution function.

Arguments
    :x (*np.ndarray*): array of measurements

Keyword Arguments
    :conf (*float*): confidence interval for probit function, default is 0.95.

Returns
    :bins (*np.ndarray*): the bin edges
'''
def calc_uniform_bins(x,bin_conf=0.95):
    # optimal number of uniformly populated bins
    n_bins = int((2*(2*float(len(x))**2/ndtri(bin_conf))**0.2)//1+1)
    # number of measurements per bin
    n_vals = int(len(x)/n_bins)
    # rank-ordered measurements
    x_sort = np.sort(x)
    # array of bin edges
    bins = np.zeros(n_bins)
    # intermediate bin edges are between n_vals*i and n_vals*i+1
    for i in range(1,n_bins):
        if n_vals*i < len(x) :
            bins[i] = (x_sort[n_vals*i]+x_sort[n_vals*i+1])/2.
    # upper limit is highest xvalue
    bins[-1] = x_sort[-1]
    # adjust lower and upper limits based on range of data set
    bins[0] = x_sort[0]-abs(x_sort[-1]-x_sort[0])
    bins[-1] = x_sort[-1]+abs(x_sort[-1]-x_sort[0])
    return bins
'''
Name
    frac_dist

Purpose
    Compute the fraction p1 of measurements x which satisfy some criterion and
    determine 1-sigma confidence limits following Gehrels (1986).

Arguments
    :x (*np.ndarray*): array of measurements
    :ind (*np.ndarray*): array of indices or booleans which select a subset of
                        x by some specified criterion

Returns
    :p1 (*np.ndarray*): the fraction of x subset in each bin
    :center (*np.ndarray*): the expected center of each bin from the
                        expected value of each x within the bin range
    :p1_lower (*np.ndarray*): the lower uncertainty in p1
    :p1_upper (*np.ndarray*): the upper uncertainty in p1
'''
def frac_dist(x,inds,x_err=None,bins=None,bin_conf=0.99999,n_samp=10000):
    # calculate number of bins
    if np.any(bins == None):
        bins = int(np.sqrt(len(x)))-1
    # number of bins
    if hasattr(bins,'__len__'):
        n_bins = len(bins)-1
    else:
        n_bins = np.copy(bins)
    # histograms
    thist,bins = np.nanhistogram(x,bins=bins)
    shist,bins = np.nanhistogram(x[inds],bins=bins)
    # bin centers
    cent = (bins[:-1]+bins[1:])/2.
    # if x uncertainties provided, do an MC simulation of Bernoulli trials
    if not np.any(x_err == None):
        # normal PMF for variable mean and standard deviation
        pmf = lambda x,mu,sig: np.exp(-0.5*((x-mu)/sig)**2)/np.sqrt(2*np.pi)
        pmf_ij = np.array([pmf(bin_i,x,x_err) for bin_i in bins])
        # normal CDF for variable mean and standard deviation
        cdf = lambda x,mu,sig: 0.5+0.5*erf((x-mu)/(abs(sig)*np.sqrt(2)))
        cdf_ij = np.array([cdf(bin_i,x,x_err) for bin_i in bins])
        # array of probabilities p_ij that jth measurement x landed in ith bin
        # given the ith bin limits [b_i-1,b_i] and the uncertainty s_j in x_j
        # following P(b_i-1 < x_j < b_i) = cdf(x_j|b_i,s_j) - cdf(x_j|b_i-1,s_j)
        prb_ij = cdf_ij[1:]-cdf_ij[:-1]
        # array of bin centers
        # from integral u * pmf(u|x,x_err) du from -inf to bin_edge
        cent_ij = np.array([np.sqrt(2*np.pi)*x*cdf_ij[i]-x_err*pmf_ij[i]\
                            for i in range(n_bins+1)])/np.sqrt(2*np.pi)
        cent = np.sum(cent_ij[1:]-cent_ij[:-1],axis=1)/np.sum(prb_ij,axis=1)
        # check for iffy values
        prb_ij[np.where(prb_ij<0)] = 0
        prb_ij[np.where(prb_ij>1)] = 1
        prb_ij[np.where(np.isnan(prb_ij))] = 0.5
        # mean counts is sum of probabilities, so use sums to get ratio
        p1 = prb_ij[:,inds].sum(axis=1)/prb_ij.sum(axis=1)
        # MC simulation to get uncertainties
        rvs = np.zeros((n_bins,len(x),10000))
        for i in range(n_bins):
            for j in range(len(x)):
                # sample Bernoulli random variates
                rvs[i,j] = bernoulli_rvs(prb_ij[i,j],n_samp=len(rvs[i,j]))
        # ratio of subset successes to total successes from variates
        p1_samp = rvs[:,inds].sum(axis=1)/rvs.sum(axis=1)
        # use quantiles to get confidence limits
        p_lo,p_up = np.quantile(p1_samp,[0.1587,0.8413],axis=1)
        # in case of nans, revert to Gehrels (1986) for upper limits
        for i in range(n_bins):
            if p_lo[i] <= 0.01 or not np.isfinite(p_lo[i]):
                p_lo[i] = binom_gehrels(shist[i],thist[i],0.68268949)[0]
            if p_up[i] >= 0.99 or not np.isfinite(p_up[i]):
                p_up[i] = binom_gehrels(shist[i],thist[i],0.68268949)[1]
    # if no uncertainties, use Gehrels (1986)
    else:
        # confidences
        p_lo = np.zeros(n_bins)
        p_up = np.zeros(n_bins)
        for i in range(n_bins):
            p_lo[i],p_up[i] = binom_gehrels(shist[i],thist[i],0.68268949)
        p1 = shist/thist
    # return to user
    return p1,cent,p1-p_lo,p_up-p1
