import numpy as np
import time
import scipy.stats as stats
from scipy.integrate import quad

def plwc_pval(x, alpha, lambd, gof):
    """
    Finds p-value for the tail-conditional power-law with exponential cut-off fit using a KS test.
    
    Input:
        x            ndarray, ndim = 1, dtype = integer
        alpha        float, exponent of the power-law
        lambd        float, rate of the exponential cut-off
        xmin         int, starting point for the power-law tail, must be >= 1
        gof          float, goodness of fit statistic (Kolmogorov-Smirnov)

    Output:
        p            p-value of the returned fit (reject power-law with exponential cut-off hypothesis for p < 0.1)
    """
    # Set desired precision level in p-value
    eps = 0.1
    xmin = np.min(x)
    num_resamps = 1000
    bootstraps = np.zeros(num_resamps)
    n = len(x)
    
    # Separate the tail and the head
    tailinds = x >= xmin
    xtail = x[tailinds]
    xhead = x[~tailinds]
    ntail = len(xtail)
    nhead = len(xhead)
    ptail = float(ntail) / n
    
    # Define the PDF for the power-law with exponential cut-off
    def pdf_powerlaw_exp(val, alpha, lambd):
        return val**(-alpha) * np.exp(-lambd * val)

    # Normalize the distribution for the tail (x >= xmin)
    def normalization_constant(alpha, lambd, xmin):
        return quad(lambda t: pdf_powerlaw_exp(t, alpha, lambd), xmin, np.inf)[0]

    Z = normalization_constant(alpha, lambd, xmin)

    # Define the CDF for the tail-conditional power-law with exponential cut-off
    def cdf_powerlaw_exp(val):
        # Handle scalar and array inputs
        val = np.atleast_1d(val)
        cdf_vals = np.zeros_like(val)

        for i, v in enumerate(val):
            if v < xmin:
                cdf_vals[i] = 0.0
            else:
                integral, _ = quad(lambda t: pdf_powerlaw_exp(t, alpha, lambd), xmin, v)
                cdf_vals[i] = integral / Z

        return cdf_vals if len(cdf_vals) > 1 else cdf_vals[0]

    # Semi-parametric bootstrap
    starttime = time.time()
    for resamp_ind in range(num_resamps):
        # Non-parametric bootstrap from the head of x
        nnewhead = n
        while nnewhead >= n:
            nnewhead = np.sum(np.random.rand(n) > ptail)
        headinds = np.array([np.floor(nhead * np.random.rand(nnewhead))], dtype=int)
        newhead = xhead[headinds][0]
        nnewtail = n - nnewhead
        
        # Parametric bootstrap for the power-law with exponential cut-off tail
        rtail = np.sort(np.random.rand(nnewtail))
        newtail = np.zeros(nnewtail, dtype=int)
        indrtail = 0
        indnewtail = 0
        xval = xmin
        while indrtail < len(rtail):
            while indrtail < len(rtail) and rtail[indrtail] <= cdf_powerlaw_exp(xval):
                newtail[indnewtail] = xval
                indnewtail += 1
                indrtail += 1
            xval += 1
            if indnewtail >= nnewtail:
                break
        
        # Combine the head and the tail
        newx = np.concatenate((newhead, newtail))
        if np.all(newx == np.zeros_like(newx)):
            import pdb; pdb.set_trace()

        # Perform the Kolmogorov-Smirnov test using the empirical CDF and the power-law with exponential cut-off CDF
        _, newgof = stats.kstest(newx, lambda t: cdf_powerlaw_exp(t))
        
        # Calculate current p-value (handling first iteration case)
        if resamp_ind == 0:
            current_p = 0.0  # No previous bootstrap to compare to in the first iteration
        else:
            current_p = np.sum(bootstraps[:resamp_ind] >= gof[0]) / float(resamp_ind)
    
        # Store goodness-of-fit statistic
        bootstraps[resamp_ind] = newgof
        
        # Early stopping if taking too long
        if time.time() - starttime > 500:
            if resamp_ind > num_resamps / 20.:
                if current_p < 0.05 or current_p > 0.5:
                    print(f"current p = {current_p}   elapsed time = {time.time() - starttime}")
                    return current_p
    
    # Final p-value
    p = np.sum(bootstraps >= gof[0]) / float(num_resamps)
    print(f"p = {p:.3f}   elapsed time = {time.time() - starttime}")
    return p

def exp_pval(x, lambd, gof):
    """ Finds p-value for the tail-conditional exponential fit using a KS test. This is based on
    Aaron's plpva.m Matlab code (http://tuvalu.santafe.edu/~aaronc/powerlaws/).

    Input:
        x            ndarray, ndim = 1, dtype = integer
        lam          float, exponent on x, must be > 1
        xmin         int, starting point for the exponential tail, must be >= 1
        gof          float, goodness of fit statistic (Kolmogorov-Smirnov)


    Output:
        p            p-value of the returned fit (reject exp hypothesis for p<0.1)
    """
    # set desired precision level in p-value
    eps = 0.1
    num_resamps = 1000
    bootstraps = np.zeros(num_resamps)
    n = len(x)
    xmin = np.min(x)

    # Separate the tail and the head
    tailinds = x>=xmin 
    xtail = x[tailinds]
    xhead = x[~tailinds]
    ntail = len(xtail)
    nhead = len(xhead)
    ptail = float(ntail)/n


    # precompute the exponential CDF
    cdf = lambda val: 1 - np.exp(-lambd * val)

    # Semi-parametric bootstrap
    starttime = time.time()
    for resamp_ind in range(num_resamps):
        # Non-parametric bootstrap from the head of x
        nnewhead = n
        while nnewhead >= n:
            nnewhead = np.sum(np.random.rand(n) > ptail)
        headinds = np.array([np.floor(nhead * np.random.rand(nnewhead))], dtype=int)
        newhead = xhead[headinds][0]
        nnewtail = n - nnewhead
        
        # Parametric bootstrap for the exponential tail
        rtail = np.sort(np.random.rand(nnewtail))
        newtail = np.zeros(nnewtail, dtype=int)
        indrtail = 0
        indnewtail = 0
        xval = xmin
        while indrtail < len(rtail):
            while indrtail < len(rtail) and rtail[indrtail] <= cdf(xval):
                newtail[indnewtail] = xval
                indnewtail += 1
                indrtail += 1
            xval += 1
            if indnewtail >= nnewtail:
                break
        
        # Combine the head and the tail
        newx = np.concatenate((newhead, newtail))
        if np.all(newx == np.zeros_like(newx)):
            import pdb; pdb.set_trace()

        # Fit this new sample
        _, newgof = stats.kstest(newx, 'expon', args=(xmin, 1.0 / lambd))
        
        # Calculate current p-value (handling first iteration case)
        if resamp_ind == 0:
            current_p = 0.0  # No previous bootstrap to compare to in the first iteration
        else:
            current_p = np.sum(bootstraps[:resamp_ind] >= gof[0]) / float(resamp_ind)
        
        # Store gof statistic
        bootstraps[resamp_ind] = newgof
        
        # Early stopping if taking too long
        if time.time() - starttime > 500:
            if resamp_ind > num_resamps / 20.:
                if current_p < 0.05 or current_p > 0.5:
                    print(f"current p = {current_p}   elapsed time = {time.time()-starttime}")
                    return current_p
    
    # Final p-value
    p = np.sum(bootstraps >= gof[0]) / float(num_resamps)
    print(f"p = {p:.3f}   elapsed time = {time.time()-starttime}")
    return p

def stretched_exp_pval(x, lambd, beta, gof):
    """
    Finds p-value for the tail-conditional stretched exponential fit using a KS test.
    
    Input:
        x            ndarray, ndim = 1, dtype = integer
        lambd        float, scale parameter (λ) of the stretched exponential
        beta         float, shape parameter (β) of the stretched exponential
        xmin         int, starting point for the stretched exponential tail, must be >= 1
        gof          float, goodness of fit statistic (Kolmogorov-Smirnov)

    Output:
        p            p-value of the returned fit (reject stretched exponential hypothesis for p < 0.1)
    """
    # Set desired precision level in p-value
    eps = 0.1
    num_resamps = 1000
    bootstraps = np.zeros(num_resamps)
    n = len(x)
    
    # Separate the tail and the head
    xmin = np.min(x)
    tailinds = x >= xmin
    xtail = x[tailinds]
    xhead = x[~tailinds]
    ntail = len(xtail)
    nhead = len(xhead)
    ptail = float(ntail) / n
    
    # Define the CDF for the tail-conditional stretched exponential distribution
    def cdf_stretched_exp(val):
        return (1 - np.exp(-lambd * val**beta)) / (1 - np.exp(-lambd * xmin**beta))
    
    # Semi-parametric bootstrap
    starttime = time.time()
    for resamp_ind in range(num_resamps):
        # Non-parametric bootstrap from the head of x
        nnewhead = n
        while nnewhead >= n:
            nnewhead = np.sum(np.random.rand(n) > ptail)
        headinds = np.array([np.floor(nhead * np.random.rand(nnewhead))], dtype=int)
        newhead = xhead[headinds][0]
        nnewtail = n - nnewhead
        
        # Parametric bootstrap for the stretched exponential tail
        rtail = np.sort(np.random.rand(nnewtail))
        newtail = np.zeros(nnewtail, dtype=int)
        indrtail = 0
        indnewtail = 0
        xval = xmin
        while indrtail < len(rtail):
            while indrtail < len(rtail) and rtail[indrtail] <= cdf_stretched_exp(xval):
                newtail[indnewtail] = xval
                indnewtail += 1
                indrtail += 1
            xval += 1
            if indnewtail >= nnewtail:
                break
        
        # Combine the head and the tail
        newx = np.concatenate((newhead, newtail))
        if np.all(newx == np.zeros_like(newx)):
            import pdb; pdb.set_trace()

        # Fit this new sample with the stretched exponential distribution
        def stretched_exp_cdf(val, lambd, beta):
            return 1 - np.exp(-lambd * val**beta)
        
        # Perform KS test using the stretched exponential CDF
        _, newgof = stats.kstest(newx, lambda t: stretched_exp_cdf(t, lambd, beta))
        
        # Calculate current p-value (handling first iteration case)
        if resamp_ind == 0:
            current_p = 0.0  # No previous bootstrap to compare to in the first iteration
        else:
            current_p = np.sum(bootstraps[:resamp_ind] >= gof[0]) / float(resamp_ind)
        
        # Store gof statistic
        bootstraps[resamp_ind] = newgof
        
        # Early stopping if taking too long
        if time.time() - starttime > 500:
            if resamp_ind > num_resamps / 20.:
                if current_p < 0.05 or current_p > 0.5:
                    print(f"current p = {current_p}   elapsed time = {time.time() - starttime}")
                    return current_p
    
    # Final p-value
    p = np.sum(bootstraps >= gof[0]) / float(num_resamps)
    print(f"p = {p:.3f}   elapsed time = {time.time() - starttime}")
    return p

def lognormal_pval(x, mu, sigma, gof):
    """
    Finds p-value for the tail-conditional log-normal fit using a KS test.
    
    Input:
        x            ndarray, ndim = 1, dtype = integer
        mu           float, mean of the logarithm of the variable for log-normal
        sigma        float, standard deviation of the logarithm of the variable for log-normal
        xmin         int, starting point for the log-normal tail, must be >= 1
        gof          float, goodness of fit statistic (Kolmogorov-Smirnov)

    Output:
        p            p-value of the returned fit (reject log-normal hypothesis for p < 0.1)
    """
    # set desired precision level in p-value
    eps = 0.1
    xmin = np.min(x)
    num_resamps = 1000
    bootstraps = np.zeros(num_resamps)
    n = len(x)
    
    # Separate the tail and the head
    tailinds = x >= xmin
    xtail = x[tailinds]
    xhead = x[~tailinds]
    ntail = len(xtail)
    nhead = len(xhead)
    ptail = float(ntail) / n
    
    # Define the CDF for the tail-conditional log-normal distribution
    def cdf_lognormal(val):
        return stats.lognorm.cdf(val, s=sigma, scale=np.exp(mu)) / (1 - stats.lognorm.cdf(xmin, s=sigma, scale=np.exp(mu)))
    
    # Semi-parametric bootstrap
    starttime = time.time()
    for resamp_ind in range(num_resamps):
        # Non-parametric bootstrap from the head of x
        nnewhead = n
        while nnewhead >= n:
            nnewhead = np.sum(np.random.rand(n) > ptail)
        headinds = np.array([np.floor(nhead * np.random.rand(nnewhead))], dtype=int)
        newhead = xhead[headinds][0]
        nnewtail = n - nnewhead
        
        # Parametric bootstrap for the log-normal tail
        rtail = np.sort(np.random.rand(nnewtail))
        newtail = np.zeros(nnewtail, dtype=int)
        indrtail = 0
        indnewtail = 0
        xval = xmin
        while indrtail < len(rtail):
            while indrtail < len(rtail) and rtail[indrtail] <= cdf_lognormal(xval):
                newtail[indnewtail] = xval
                indnewtail += 1
                indrtail += 1
            xval += 1
            if indnewtail >= nnewtail:
                break
        
        # Combine the head and the tail
        newx = np.concatenate((newhead, newtail))
        if np.all(newx == np.zeros_like(newx)):
            import pdb; pdb.set_trace()

        # Fit this new sample with the log-normal distribution
        _, newgof = stats.kstest(newx, 'lognorm', args=(sigma, 0, np.exp(mu)))
        
        # Calculate current p-value (handling first iteration case)
        if resamp_ind == 0:
            current_p = 0.0  # No previous bootstrap to compare to in the first iteration
        else:
            current_p = np.sum(bootstraps[:resamp_ind] >= gof[0]) / float(resamp_ind)
        
        # Store gof statistic
        bootstraps[resamp_ind] = newgof
        
        # Early stopping if taking too long
        if time.time() - starttime > 500:
            if resamp_ind > num_resamps / 20.:
                if current_p < 0.05 or current_p > 0.5:
                    print(f"current p = {current_p}   elapsed time = {time.time() - starttime}")
                    return current_p
    
    # Final p-value
    p = np.sum(bootstraps >= gof[0]) / float(num_resamps)
    print(f"p = {p:.3f}   elapsed time = {time.time() - starttime}")
    return p
