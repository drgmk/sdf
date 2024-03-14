"""Analytic routines for debris disks."""

import numpy as np

from . import photometry
from . import filter
from . import utils


class BB_Disk(object):
    """A blackbody disk class.
        
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
    """

    def __init__(self, wavelengths=None, temperatures=None,
                 lstar=None, tstar=None, distance=None):
        """Initialise, default T=100K, Omega=1.0"""
        
        if wavelengths is None:
            self.wavelengths = 10**np.linspace(-1, 4, 1000)
        else:
            self.wavelengths = wavelengths
        
        if temperatures is None:
            self.temperatures = 10**np.linspace(1, 3, 1000)
        else:
            self.temperatures = temperatures

        self.lstar = lstar
        self.tstar = tstar
        self.distance = distance

    def blackbody_radii(self):
        """Return the blackbody radii."""
        return (278.3/self.temperatures)**2 * self.lstar**0.5

    def radiance(self):
        """Return radiance, in W / m^2 / sr."""
        return 5.67e-8 * self.temperatures**4 / np.pi

    def f_limits(self, lim_waves, flux_limits=None, r_limits=None,
                 stellar_flux=None, fwhm=None, lstar_1pc=None,
                 resolved='detection'):
        """Return fractional luminosity limits.
            
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
        resolved : str, optional
            Whether limit requires detection, radial, or fully resolved.
            
        Resolved keyword allows for disks to be resolved. Default is just
        detection, in which case s/n ~ 1/r, since n_beams increases as r^2,
        but noise is sqrt(n_beams) and signal is constant. For radial,
        (azimuthal averaging) signal decreases as 1/r and noise per
        radial bin increases as r, so s/n ~ 1/r^1.5. For fully noise is
        constant but signal decreases as r^2, so s/n ~ 1/r^2.

        One of flux_limits or r_limits must be given. If both, they must
        have the same length, and correspond to the wavelengths given.
        Likewise for stellar_flux and fwhm.
        """

        if flux_limits is not None and r_limits is not None:
            if len(flux_limits) != len(r_limits):
                raise RuntimeError(
                       'flux_limits must be same length as r_limits')

        # sensitivity limit
        if flux_limits is not None:
            
            slims = np.zeros((len(self.temperatures), len(flux_limits)))

            for i, temp in enumerate(self.temperatures):

                if self.distance is not None:
                    slims[i, :] = 3.4e9 * flux_limits * self.distance**2 / \
                                self.blackbody_radii()[i]**2 / \
                                utils.bnu_wav_micron(lim_waves, temp)
                else:
                    # distance independent calculation, 2487305. is
                    # pc^2/Lsun, haven't tracked down the 4 yet
                    ldisk_1pc = 4 * 5.6704e-8 * flux_limits * 2487305. * \
                        temp**4 / utils.bnu_wav_micron(lim_waves, temp)
                    slims[i, :] = ldisk_1pc / lstar_1pc

                # apply correction for resolved disks
                if self.distance is not None and fwhm is not None:
                    fwhm_fact = 2 * self.blackbody_radii()[i] / self.distance / fwhm
                    if resolved == 'radial':
                        fwhm_fact = fwhm_fact**1.5
                    elif resolved == 'fully':
                        fwhm_fact = fwhm_fact**2
                    elif resolved == 'detection':
                        pass
                    else:
                        raise utils.SdfError('Resolved needs to be "detection/radial/fully"')
                    resolvedi = fwhm_fact > 1.0
                    slims[i, resolvedi] *= fwhm_fact[resolvedi]

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
                fstar = 1.77 * utils.bnu_wav_micron(lim_waves, self.tstar) * \
                            self.lstar / self.tstar**4 / self.distance**2

            clims = np.zeros((len(self.temperatures), len(r_limits)))
            for i, temp in enumerate(self.temperatures):
                clims[i, :] = 6e9/1.77 * r_limits * fstar / \
                                utils.bnu_wav_micron(lim_waves, temp) * \
                                (self.distance/self.blackbody_radii()[i])**2

        if flux_limits is not None and r_limits is not None:
            return np.minimum(slims, clims)
        elif flux_limits is not None:
            return slims
        elif r_limits is not None:
            return clims
        else:
            raise RuntimeError('Need to pass flux_limits or r_limits')

    def f_limits_from_result(self, r, min_wavelength=8.0, sn=3,
                             x={}, x_det={},
                             skip_filters=[], keep_filters=None):
        """Derive fractional luminosity limits from an sdf result object.
        
        Also derive fractional luminosities and signal to noise of excess
        detections. Return low and high limits, expect to plot these
        with pyplot.fill_between and something like:
        
        ax.fill_between(temps, det_lo[:, i], det_hi[:, i],
                        where=(det_lo[:, i]<det_hi[:, i]), alpha=0.25)

        Account for long wavelength grain inefficiency with X factor,
        used per filter, e.g. {'WAV850':4}.

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
        sn : float, optional
            S/N at which detection significant, used only for detections.
        x : dict
            X factor to increase limits by: {filter, X}
        x_det : dict
            X factor to increase upper detection limit by: {filter, X}
        skip_filters : list, optional
            List of filters to skip.
        keep_filters : list, optional
            List of filters to keep, applied after skip_filters.
        """

        waves = np.array([])
        filters = np.array([])
        f_lim = np.array([])
        f_det = np.array([])
        e_det = np.array([])
        f_star = np.array([])

        # get stellar luminosity at 1pc if no distance
        lstar = None
        if self.distance is None:
            lstar = 0.0
            if hasattr(r, 'star'):
                for s in r.star:
                    lstar += s['lstar_1pc']
    
            if lstar == 0.0:
                raise utils.SdfError('dont have lstar_1pc or distance')
        
        for p in r.obs:
            if not isinstance(p, photometry.Photometry):
                continue
        
            ok = np.invert(p.ignore)
            # loop to grab correct stellar photometry
            for i, filt in enumerate(p.filters[ok]):
            
                new_wave = p.mean_wavelength()[ok][i]
                if (filter.iscolour(filt) or
                        new_wave < min_wavelength or
                        filt in skip_filters):
                    continue
            
                if keep_filters is not None:
                    if filt not in keep_filters:
                        continue
                        
                waves = np.append(waves, new_wave)
                filters = np.append(filters, filt)
                filt_i = np.where(filt == np.array(r.all_filters))[0]
                f_star = np.append(f_star, r.all_star_phot[filt_i])
                
                fac = 1
                if filt in x.keys():
                    fac = x[filt]

                if p.upperlim[ok][i]:
                    f_lim = np.append(f_lim, p.fnujy[ok][i]*fac)
                    f_det = np.append(f_det, 0)
                    e_det = np.append(e_det, 0)
                else:
                    # 1sigma uncertainty, observed and star in quadrature
                    unc = np.sqrt(
                        p.e_fnujy[ok][i]**2 +
                        0.25*(r.all_star_phot_1sig_lo[filt_i] + r.all_star_phot_1sig_hi[filt_i])**2
                                  )
                    f_lim = np.append(f_lim, 3*unc*fac)
                    f_det = np.append(f_det, p.fnujy[ok][i] - f_star[-1])
                    e_det = np.append(e_det, unc)

        lims = self.f_limits(waves, flux_limits=f_lim,
                             stellar_flux=f_star, lstar_1pc=lstar)
        dets = self.f_limits(waves, flux_limits=f_det,
                             stellar_flux=f_star, lstar_1pc=lstar)

        ok = e_det > 0
        sn_dets = np.zeros(lims.shape[1])
        sn_dets[ok] = f_det[ok] / e_det[ok]

        # now compute limit ranges for detections, first get ranges
        det_lo = np.zeros(lims.shape)
        det_hi = lims.copy()
        both_hi = lims.copy()
        for i in range(lims.shape[1]):
            if sn_dets[i] > sn:
                fac = 1
                if filters[i] in x_det.keys():
                    fac = x_det[filters[i]]
                det_lo[:, i] = dets[:, i]*(1-sn/sn_dets[i])
                det_hi[:, i] = dets[:, i]*(fac+sn/sn_dets[i])
                both_hi[:, i] = np.max([[det_hi[:, i]], [lims[:, i]]], axis=0)

        # now adjust high limit based on other limits
        for i in range(lims.shape[1]):
            other = np.arange(lims.shape[1]) != i
            det_hi[:, i] = np.min(np.hstack((both_hi[:, other], det_hi[:, i].reshape((-1, 1)))), axis=1)
        
        return lims, det_lo, det_hi, sn_dets, filters

    def f_limits_togrid(self, lims, f=None):
        """Return boolean grid in f - r/T space indicating detectability.
        
        Sum multiple of these to get the grid that shows how many of the
        systems it was possible to detect a disk for.
        
        Parameters
        ----------
        lims : array
            Array of f limits (i.e. n_temperatures x n_lim).
        f : array, optional
            Array of f to use in grid.
        """
        
        if f is None:
            f = 10**np.linspace(-7, -1, 100)
            
        fs, _ = np.meshgrid(f, self.temperatures)
        return fs > np.min(lims, axis=1), f
