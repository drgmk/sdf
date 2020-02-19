'''Analytic routines for debris disks.'''

import numpy as np

from . import photometry
from . import filter
from . import utils

class BB_Disk(object):
    '''A blackbody disk class.
        
    Takes multiple temperatures, the purpose being for use to show
    disk properties in parameter spaces such as fractional luminosity
    vs. temperature.
    
    Parameters
    ----------
    lstar : float
        Stellar luminosity in Solar units.
    tstar : float
        Stellar effective temperature in Kelvin.
    distance : float
        Stellar distance in parsec.
    wavelengths : 1-D array, optional
        Vector of wavelengths.
    temperatures : 1-D array, optional
        Vector of temperatures.

    .. todo:: distance not actually needed for calibration limited, fix.

    .. todo:: don't use a for loop over temperatures,
              fix utils.bnu_wav_micron instead.
    '''

    def __init__(self,wavelengths=None,temperatures=None,
                 lstar=None,tstar=None,distance=None):
        '''Initialise, default T=100K, Omega=1.0'''
        
        if wavelengths is None:
            self.wavelengths = 10**np.linspace(-1,4,1000)
        else:
            self.wavelengths = wavelengths
        
        if temperatures is None:
            self.temperatures = 10**np.linspace(1,3,1000)
        else:
            self.temperatures = temperatures

        self.lstar = lstar
        self.tstar = tstar
        self.distance = distance


    def blackbody_radii(self):
        '''Return the blackbody radii.'''
        
        return (278.3/self.temperatures)**2 * self.lstar**0.5


    def radiance(self):
        '''Return radiance, in W / m^2 / sr.'''
        
        return 5.67e-8 * self.temperatures**4 / np.pi


    def f_limits(self,lim_waves,flux_limits=None,r_limits=None,
                 stellar_flux=None,fwhm=None,lstar_1pc=None):
        '''Return fractional luminosity limits.
            
        This routine implements Wyatt (2008) equations 8 and 11.
        
        Parameters
        ----------
        lim_waves : numpy.ndarray
            Array of wavelengths at which limits apply.
        flux_limits : numpy.ndarray, optional
            Array of flux limits.
        r_limits : numpy.ndarray, optional
            Array of calibration limits (F_disk/F_star).
        stellar_flux : numpy.ndarray, optional
            Array of stellar fluxes at lim_waves.
        fwhm : numpy.ndarray, optional
            Array of spatial resolutions at lim_waves, affects flux
            limited observations if disk is resolved.
        lstar_1pc : float
            L_star at 1pc, used for flux limits when distance unknown.

        One of flux_limits or r_limits must be given. If both, they must
        have the same length, and correspond to the wavelengths given.
        Likewise for stellar_flux and fwhm.
        '''

        if flux_limits is not None and r_limits is not None:
            if len(flux_limits) != len(r_limits):
                raise RuntimeError(
                       'flux_limits must be same length as r_limits')

        # sensitivity limit
        if flux_limits is not None:
            
            slims = np.zeros((len(self.temperatures),len(flux_limits)))

            for i,temp in enumerate(self.temperatures):

                if self.distance is not None:
                    slims[i,:] = 3.4e9 * flux_limits * self.distance**2 / \
                                self.blackbody_radii()[i]**2 / \
                                utils.bnu_wav_micron(lim_waves,temp)
                else:
                    # distance independent calculation, 2487305. is
                    # pc^2/Lsun, haven't tracked down the 4 yet
                    ldisk_1pc = 4 * 5.6704e-8 * flux_limits * 2487305. * \
                        temp**4 / utils.bnu_wav_micron(lim_waves,temp)
                    slims[i,:] = ldisk_1pc / lstar_1pc

                # apply correction for resolved disks
                if self.distance is not None and fwhm is not None:
                    fwhm_fact = 2 * self.blackbody_radii()[i] / self.distance / fwhm
                    resolved = fwhm_fact > 1.0
                    slims[i,resolved] *= fwhm_fact[resolved]

        
        # calibration limit, use actual stellar flux if given
        if r_limits is not None:
            
            if stellar_flux is not None:
                if len(stellar_flux) != len(r_limits):
                    raise RuntimeError(
                           'Stellar flux ({}) must have same '
                           'length as r_limits ({})'.format(
                                                        len(stellar_flux),
                                                        len(r_limits)
                                                            )
                                       )
                fstar = stellar_flux
            else:
                fstar = 1.77 * utils.bnu_wav_micron(lim_waves,self.tstar) * \
                            self.lstar / self.tstar**4 / self.distance**2

            clims = np.zeros((len(self.temperatures),len(r_limits)))
            for i,temp in enumerate(self.temperatures):
                clims[i,:] = 6e9/1.77 * r_limits * fstar / \
                                utils.bnu_wav_micron(lim_waves,temp) * \
                                (self.distance/self.blackbody_radii()[i])**2

        if flux_limits is not None and r_limits is not None:
            return np.minimum(slims,clims)
        elif flux_limits is not None:
            return slims
        elif r_limits is not None:
            return clims
        else:
            raise RuntimeError('Need to pass flux_limits or r_limits')


    def f_limits_from_result(self,r,min_wavelength=8.0,
                             skip_filters=[],keep_filters=None):
        '''Derive fractional luminosity limits from an sdf result object.
            
        Rather than worry about flux vs. calibration limited, just do 
        the calculation assuming flux limited by calculating the flux
        limit for each observed filter (whether it was an upper limit
        or not).
        
        Parameters
        ----------
        r : sdf.result.Result
            Result object with photometry.
        min_wavelength : float, optional
            Exclude filters with a mean wavelength shorter than this.
        skip_filters : list, optional
            List of filters to skip.
        keep_filters : list, optional
            List of filters to keep, applied after skip_filters.
        '''

        waves = np.array([])
        filters = np.array([])
        f_lim = np.array([])
        f_star = np.array([])

        # get stellar luminosity at 1pc if no distance
        lstar = None
        if self.distance is None:
            lstar = 0.0
            if hasattr(r,'star'):
                for s in r.star:
                    lstar += s['lstar_1pc']
    
            if lstar == 0.0:
                raise utils.SdfError('dont have lstar_1pc or distance')
        
        for p in r.obs:
            if not isinstance(p,photometry.Photometry):
                continue
        
            ok = np.invert(p.ignore)
            # loop to grab correct stellar photometry
            for i,filt in enumerate(p.filters[ok]):
            
                new_wave = p.mean_wavelength()[ok][i]
                if (filter.iscolour(filt) or 
                    new_wave < min_wavelength or
                    filt in skip_filters):
                    continue
            
                if keep_filters is not None:
                    if filt not in keep_filters:
                        continue
                        
                waves = np.append(waves,new_wave)
                filters = np.append(filters,filt)
                filt_i = np.where(filt == np.array(r.all_filters))[0]
                f_star = np.append(f_star,r.all_star_phot[filt_i])
                
                if p.upperlim[ok][i]:
                    f_lim = np.append(f_lim,p.fnujy[ok][i])
                else:
                    # 1sigma uncertainty, observed and star in quadrature
                    unc = np.sqrt(
              p.e_fnujy[ok][i]**2 + \
              0.25*(r.all_star_phot_1sig_lo[filt_i] + r.all_star_phot_1sig_hi[filt_i])**2
                                  )
                    f_lim = np.append(f_lim,3*unc)

        lims = self.f_limits(waves,flux_limits=f_lim,
                             stellar_flux=f_star,lstar_1pc=lstar)
        return lims,filters

