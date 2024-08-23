import numpy as np
from scipy.optimize import brentq
import time
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
import astropy.units as u
import pandas as pd

from scipy.stats.mstats import mquantiles
from scipy.stats import norm
import matplotlib.pyplot as plt

L=1350 #parsecs

# 1 - define the probability function
def prob(r,o,s):
    p = (r**2*np.exp(-r/L)/(s))*np.exp((-1/(2*(s)**2))*(o-(1/r))**2)
    if not np.isfinite(p):
        return -np.inf
    return p

# 2 - define the MCMC sampler:
def sampler(lnprob,mu_init=0,o=0,s=0,proposal_width=0.5,nsamples=50,prob=[0.05, 0.5, 0.95]):
    '''Samples the parameter space for a 1-parameter function.  It is specific to the problem at
       hand and not generalizable as written to other functions.
    Args:
        lnprob (func): the probability function
        mu_init (float): initial guess for critial point of function
        o (float): parallax for probability function
        s (float): parallax uncertainty
        proposal_width (float): how large should each jump be (try to get acceptance rate ~30-40%)
        nsamples (int): how many samples to test
        prob (array): quantiles to return.  Default is 5%, 50%, and 95%
    Returns
        posterior (array): The chain of samples
        mquantiles (array): the compiuted quantiles
        acceptance rate (float)
    '''
    # Begin with initial guess
    mu_current = mu_init
    # Initialize the posterior chain
    posterior = [mu_current]
    # Initialize the acceptance rate tracker
    yes_accept = 0
    for i in range(nsamples):
        # Propose a new value for mu from a normal distribution centered at the current
        # value of mu and spread by the proposal width
        mu_proposal = norm(mu_current,proposal_width).rvs()
        # Determine the likelihood of the new proposal
        prob_of_proposal = lnprob(mu_proposal,o,s)
        # Determine the likelihood of the old mu value
        prob_current = lnprob(mu_current,o,s)
        # Determing the likelihood ratio of new to old:
        p_accept = prob_of_proposal/prob_current
        # Random "dice roll"
        dice = np.random.rand()
        # Accept the new proposal if it is more likely than the dice roll
        accept = dice < p_accept
        if accept:
            # If accepted move to the new value.  Otherwise stay on the old value
            mu_current = mu_proposal
            yes_accept = yes_accept+1 #for tracking acceptance rate
        # Add the current value to the posterior chain and repeat.
        posterior.append(mu_current)
    print('posterior = ',posterior)
    print('prob = ',prob)
    print('mquantiles(posterior, prob=prob) = ',mquantiles(posterior, prob=prob))
    print('yes_accept = ',yes_accept)
    print('nsamples = ',nsamples)
    return posterior,mquantiles(posterior, prob=prob),float(yes_accept)/float(nsamples)

""" input parameter k: gaia table """
""" return distances in pc """
def get_gaia_distance(k, display=False):
    # Set scale length for prior:
    L=1350 #parsecs

    start = time.time()
    # convert to arcsec and add in zero-point shift
    parallaxes = k['parallax'].to_numpy()
    parallaxes = np.array([float(x) for x in parallaxes])
    parallax_errors = k['parallax_error'].to_numpy()
    parallax_errors = np.array([float(x) for x in parallax_errors])
    gdists = []
    fwhm_los = []
    fwhm_his=[]
    posts=[]
    quants=[]
    accept_rates=[]
    for i in range(parallaxes.shape[0]):
        parallax = parallaxes[i]
        parallax_error = parallax_errors[i]
        omega,sigma = (parallax+0.029)/1000,parallax_error/1000
        print('omega = ',omega,', sigma = ',sigma)
        gdist = np.array([])

        print('Computing distances')
        count=0
        f = parallax_error/parallax

        # Initialize the 95% CI arrays:
        dist_95ci_lo,dist_95ci_hi = np.array([]),np.array([])

        rmax = 1e6
        fwhm_lo,fwhm_hi = np.array([]),np.array([])
        # establish the coefficients of the mode-finding polynomial:
        coeff = np.array([(1./L),(-2),((omega)/((sigma)**2)),-(1./((sigma)**2))])
        # use numpy to find the roots:
        g = np.roots(coeff)
        # Find the number of real roots:
        reals = np.isreal(g)
        realsum = np.sum(reals)
        # If there is one real root, that root is the  mode:
        if realsum == 1:
            gd = np.real(g[np.where(reals)[0]])
        # If all roots are real:
        elif realsum == 3:
            if omega >= 0:
                # Take the smallest root:
                gd = np.min(g)
            elif omega < 0:
                # Take the positive root (there should be only one):
                gd = g[np.where(g>0)[0]]
        gdist = np.append(gdist,gd)
        rmode = gdist[0]
        M = (rmode**2*np.exp(-rmode/L)/sigma)*np.exp((-1./(2*(sigma)**2))*(omega-(1./rmode))**2)
        lo = brentq(lambda x: 2*np.log(x)-(x/L)-(((omega-(1./x))**2)/(2*sigma**2)) \
                +np.log(2)-np.log(M)-np.log(sigma), 0.001, rmode)
        hi = brentq(lambda x: 2*np.log(x)-(x/L)-(((omega-(1./x))**2)/(2*sigma**2)) \
                +np.log(2)-np.log(M)-np.log(sigma), rmode, rmax)
        fwhm_lo,fwhm_hi = np.append(fwhm_lo,lo),np.append(fwhm_hi,hi)
        if f<0.1:
            proposal_width=25
        elif f>0.1 and f<0.25:
            proposal_width=100
        else:
            proposal_width=200
#        print('prob = ',prob)
#        print('gdist = ',gdist)
#        print('omega = ',omega)
#        print('sigma = ',sigma)
#        print('proposal_width = ',proposal_width)
#        STOP
        post,quant,accept_rate = sampler(prob,mu_init=gdist[0],o=omega,s=sigma,nsamples=5000,proposal_width=proposal_width)
        dist_95ci_lo,dist_95ci_hi = np.append(dist_95ci_lo,quant[0]),np.append(dist_95ci_hi,quant[2])

        if display:
            plt.subplot(111)
            r = np.linspace(1,4400,999)
            print('r = ',r)
            post1 = (r**2*np.exp(-r/L)/sigma)*np.exp((-1/(2*(sigma)**2))*(omega-(1/r))**2)
            print('post1 = ',post1)
            plt.plot(r,post1/np.max(post1),label='f = {0}'.format(np.round(f,decimals=3)))
            plt.axvline(x=dist_95ci_lo)
            plt.axvline(x=dist_95ci_hi)
            plt.axvline(x=fwhm_lo,ls=':')
            plt.axvline(x=fwhm_hi,ls=':')
            plt.legend()
            plt.ylabel('$P*\;(r \mid \omega , \sigma_{\omega})$')
            plt.show()
        print('Finished distances for ',gdist.shape[0], 'sources')

        end = time.time()
        print('This took: ',(end - start),' s')
        gdists.append(gdist)
        fwhm_los.append(fwhm_lo)
        fwhm_his.append(fwhm_hi)
        posts.append(post)
        quants.append(quant)
        accept_rates.append(accept_rate)
    return gdists,fwhm_los,fwhm_his,posts,quants,accept_rates


def readGaiaMainTable(ra_deg, dec_deg, rad_deg, row_limit=10000000):
    Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"  # Reselect Data Release 3, default
    Gaia.ROW_LIMIT = row_limit  # Ensure the default row limit.
    coord = SkyCoord(ra=ra_deg, dec=dec_deg, unit=(u.degree, u.degree), frame='icrs')
    j = Gaia.cone_search_async(coord, radius=u.Quantity(rad_deg, u.deg))
    r = j.get_results()
    return r

raDeg = 244.9175# degrees
decDeg = -49.2331# degrees
gaia_stars = pd.read_csv('/Users/xxz/Desktop/LSR-labintern/all_GAIA_stars_in_area.csv')
print('gaia = ',gaia_stars)
print('dir(gaia) = ',dir(gaia_stars))
print('gaia_stars.colnames = ',gaia_stars.keys())
print('v_rad = ',gaia_stars['radial_velocity'])
print('Teff = ',gaia_stars['teff_gspphot'])
print('parallax = ',gaia_stars['parallax'])
dist,fwhm_lo,fwhm_hi,post,quant,accept_rate = get_gaia_distance(gaia_stars,display=False)
print('dist = ',dist)
print('fwhm_lo = ',fwhm_lo)
print('fwhm_hi = ',fwhm_hi)
print('post = ',post)
print('quant = ',quant)
print('accept_rate = ',accept_rate)

