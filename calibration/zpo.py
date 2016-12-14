import os
import sys
sys.path.insert(0, os.path.abspath('..'))
import subprocess
import time

import mysql.connector
from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.table import Table,operations

from sdf import photometry
from sdf import model
from sdf import spectrum
from sdf import plotting
from sdf import fitting
from sdf import filter
from sdf import config as cfg

import subprocess
import time
import mysql.connector
from scipy.optimize import minimize
from sdf import config as cfg
import numpy as np
def chisq_all(zpos):
    """Run SED fitting for everything and return photometry chisq."""

    # update the ZPOs
    zp_filt = np.array(['UJ','BJ',
                        'US','VS','BS',
                        'HP',
                        'RC','IC',
#                        'BT','VT',
#                        '2MJ','2MH',
#                        'WISE3P4'
                        ])

    f = open('/Users/grant/astro/projects/sdf/sdf/calibration/zpos.txt','w')
    f.write( ' '.join(zp_filt) )
    f.write('\n')
    f.write( ' '.join([str(s)+' ' for s in zpos]) )
    f.close()

    # clear everything
    subprocess.run(['python3','cleanup.py','--sample','zpo_cal_'])

    # run all the fitting, Popen doesn't wait for end, run does
    cmd = ['python3','fit.py','-w','-p',
           '--no-spectra',
           '--sample','zpo_cal_','--subset','public']
    subprocess.Popen(cmd,stdin=None,stdout=None,stderr=None)
    time.sleep(2)
    subprocess.Popen(cmd,stdin=None,stdout=None,stderr=None)
    time.sleep(2)
    subprocess.Popen(cmd,stdin=None,stdout=None,stderr=None)
    time.sleep(2)
    subprocess.Popen(cmd,stdin=None,stdout=None,stderr=None)
    time.sleep(2)
    subprocess.Popen(cmd,stdin=None,stdout=None,stderr=None)
    time.sleep(2)
    subprocess.Popen(cmd,stdin=None,stdout=None,stderr=None)
    time.sleep(2)
#    subprocess.Popen(cmd,stdin=None,stdout=None,stderr=None)
#    time.sleep(2)
    subprocess.run(cmd)
    
    # wait for everything to finish
    time.sleep(20)

    # get the db output
    cnx = mysql.connector.connect(user=cfg.mysql['user'],
                                  password=cfg.mysql['passwd'],
                                  host=cfg.mysql['host'],
                                  database=cfg.mysql['db_results'])
    cursor = cnx.cursor(buffered=True)
    cursor.execute("SELECT SUM(chi*chi) FROM sdb_samples.zpo_cal_ "
                   "LEFT JOIN sdb_results.model on sdbid=id "
                   "LEFT JOIN phot USING (id) WHERE filter "
                   "REGEXP('UJ_BJ|BJ_VJ|STROMC1|STROMM1|BS_YS|VJ|2MKS') "
                   "AND obs_upperlim=0")
    return float( cursor.fetchall()[0][0] )


zpo0 = [0.04050373,  0.0217984 ,  1.29187515,  0.23172751,  0.03709482,
        0.03903796,  0.04579623,  0.01117473 ]
res = minimize(chisq_all,zpo0,method='Nelder-Mead',options={'maxiter':10})


def mvb_zpo(fname,flam=None):
    """Get the Mann & von Braun ZPO for the CALSPEC Vega spectrum.
        
        Values are scaled so that the ZPO at Johnson V is 0.027.
        """
    
    f_vega = filter.Filter.get(fname)
    
    # this assumes the SVO has the correct zero points
    svn = filter.filter_info.filters[fname]['svo_name']
    f_svo = filter.Filter.svo_get(svn)
    zpo_flux = f_svo.zero_point / f_vega.zero_point
    zpo = 2.5 * np.log10( zpo_flux ) + 0.027 - 0.064605028
    print(fname)
    print('SVO: ',f_svo.zero_point)
    print('Vega:',f_vega.zero_point)
    print('ZPO: ',zpo)
    print('')
    
    # this uses the pivot wavelength to convert F_lam
    fnu = flam * u.Unit('erg/(cm**2*s*AA)').\
        to('Jy',equivalencies=u.spectral_density(f_vega.pivot_wavelength()*u.micron))
    
    zpo_flux = fnu / f_vega.zero_point
    zpo = 2.5 * np.log10( zpo_flux ) #+ 0.027 - 0.04777335
    
    
    print('MvB: ',fnu)
    print('Vega:',f_vega.zero_point)
    print('ZPO: ',zpo)


"""This was supposed to be a way to set zero point offsets using some
    set of reference stars, in particular Gaia benchmarks, but it didn't
    look like it was going to work very well so has been shelved
    """

def chisq_norm(norm,*args):
    
    
    par0,obs,mod = args
    res = minimize(fitting.chisq,
                   np.append(par0,par),
                   args=( (obs,), (mod,) ),
                   method='Nelder-Mead')
                   print(res['x'])
                   return res['fun']


def chisq_one(obs,mod,pars):
    """Get the chisq for obs given a model and best normalisation."""
    
    filt = np.invert( filter.iscolour(tuple(obs.filters)) )
    norm = 0.0
    par = np.append(pars,norm)
    mod_fnujy = mod.fnujy(par)[:len(obs.fnujy)]
    norm = np.log10( np.median(obs.fnujy[filt]/mod_fnujy[filt]) )
    chisq = fitting.chisq(np.append(pars,norm),(obs,),(mod,))
    return chisq,np.append(pars,norm)

# get the info for the targets
t1 = Table.read('calibration/gaia_benchmark/gb.csv')
t2 = Table.read('calibration/gaia_benchmark/ids.csv')
t = operations.join(t1,t2)
tok = [1,4,7,8,13,15,17,19,20,21,24,28,29,31,32,35,36,37]
keep = np.zeros(len(t),dtype=bool)
for i in tok:
    keep[i] = True
t = t[keep]

# set up filter multipliers, exclude 2MASS K and Johnson V
zp_filt = np.array(['UJ','BJ','RC','IC','US',
                    'VS','BS','2MJ','2MH',
                    'BT','VT','HP'])
filt_keep = np.append(zp_filt,['STROMM1','STROMC1','BS_YS','UJ_BJ','BJ_VJ',
                      'VJ_IC','VJ_RC'])
zp_fac = np.zeros(len(zp_filt))+0.03


# grab the photometry and the models
obs = []
mod = []
par0s = []
allfilt = []
chisq = 0.0

for i,id in enumerate(t['sdbid']):

    obs.append( photometry.Photometry.read_sdb_file('calibration/gaia_benchmark/'+id+'-rawphot.txt',keep_filters=filt_keep) )
    mod1,mod2 = model.get_models((obs[-1],),('kurucz_m',))
    mod.append(mod1[0][0])

    par0s.append( np.array( [t[i]['Teff'],t[i]['logg'],t[i]['[Fe/H]']] ) )

    chi,par = chisq_one(obs[i],mod[i],par0s[i])
    chisq += chi
    resid,_,_ = fitting.residual( par,(obs[i],),(mod[i],) )

    plt.plot(filter.mean_wavelength(obs[i].filters),resid+i*5)
#    plt.plot(filter.mean_wavelength(obs[i].filters),
#             obs[i].fnujy/mod[i].fnujy(par)[:len(obs[i].fnujy)]+i)
    plt.plot(filter.mean_wavelength(obs[i].filters),np.zeros(len(obs[i].filters))+i*5,'--')

for f in filt_keep:
    plt.text(filter.mean_wavelength([f]),-.5,f,rotation=90)

plt.show()
