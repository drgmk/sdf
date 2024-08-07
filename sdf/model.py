from functools import lru_cache
import os
import glob
import copy

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.interpolate import interp1d, RegularGridInterpolator
import astropy.units as u
from astropy.io import fits
from astropy.table import Table

from . import convolve
from . import photometry
from . import spectrum
from . import filter
from . import utils
from . import config as cfg


class Model(object):
    """Basic model class.
        
    PhotModel and SpecModel are derived from this, the only real
    difference is that the former has convolved flux data for an array
    of filters, and the latter has flux data for an array of
    wavelengths.
    
    Many of the methods and functions required are similar for both, so
    are contained in the base class with checking to ensure the right
    thing is done. """

    @lru_cache(maxsize=8)
    def read_file(file):
        """Read a model file."""

        # parameter names
        fh = fits.open(file)
        keywords = fh[0].header
        nparam = keywords['NPARAM']

        # see what type of model we have
        type_ = keywords['SDFTYPE']
        if type_ == 'PhotModel':
        
            self = PhotModel()
            # get the filter names
            dat = fh[2].data
            fs = np.array(dat, dtype=str)
            for i in range(len(fs)):
                fs[i] = fs[i].strip()
            self.filters = np.array(fs)
            
        elif type_ == 'SpecModel':
        
            self = SpecModel()
            # get the wavelengths
            dat = fh[2].data
            self.wavelength = np.array(dat, dtype=dat.dtype[0])
        
        else:
            raise utils.SdfError("file {} not a PhotModel or SpecModel, \
                           is a {}".format(file, type_))

        self.name = keywords['NAME']
        self.parameters = [keywords['PARAM'+str(i)] for i in range(nparam)]

        # parameter ranges, assume that all values have the same dtype
        # so it's OK if the resulting ndarray does
        d = {}
        for i, par in enumerate(self.parameters):
            j = i+3
            if str.upper(par) != fh[j].name:
                raise utils.SdfError("{}th parameter {} not equal to HDU\
                               with name {}".format(j, par, fh[j].name))
            dat = fh[j].data
            d[par] = np.array(dat, dtype=dat.dtype[0])
        self.param_values = d

        # main array of convolved model fluxes, do this last so the
        # dimensionality is checked by the setter
        self.fnujy_sr = fh[1].data
        header1 = fh[1].header
        if header1['BUNIT'] != (u.jansky/u.sr).to_string(format='fits'):
            fnuunit = header1['BUNIT']
            self.fnujy_sr = self.fnujy_sr * u.Unit(fnuunit).to('Jy / sr')

        # integer indices of filters/wavelengths
        if type_ == 'PhotModel':
            self.i = np.arange(len(self.filters))
        elif type_ == 'SpecModel':
            self.i = np.arange(len(self.wavelength))
        self.n_i = len(self.i)

        # create the hashed version
        self.fill_log_fnujy_sr_hashed()
        
        return self

    def write_file(self, file, overwrite=False):
        """Write model to a FITS file.
        
        Number of parameters and their names are stored in the primary
        HDU, flux density cube is in the first HDU, and arrays with the
        filter names and parameter ranges in subsequent HDUs.
        """

        # primary HDU (for metadata)
        hdu0 = fits.PrimaryHDU()
        hdu0.header['NAME'] = self.name
        hdu0.header['NPARAM'] = len(self.parameters)
        # keywords for parameters
        for i, par in enumerate(self.parameters):
            hdu0.header['PARAM'+str(i)] = par

        # fluxes
        hdu1 = fits.ImageHDU(self.fnujy_sr, name='MODELS')
        hdu1.header['BUNIT'] = (u.jansky/u.sr).to_string(format='fits')

        # see what type of model we're writing
        if isinstance(self, PhotModel):
        
            hdu0.header['SDFTYPE'] = 'PhotModel'
            # filter names
            t = Table()
            t.add_column(Table.Column(self.filters, name='FILTERS'))
            hdu2 = fits.BinTableHDU(np.array(t), name='FILTERS')
        
        elif isinstance(self, SpecModel):
        
            hdu0.header['SDFTYPE'] = 'SpecModel'
            # Wavelengths
            t = Table()
            t.add_column(Table.Column(self.wavelength, name='WAVLNTHS'))
            hdu2 = fits.BinTableHDU(np.array(t), name='WAVLNTHS')
        
        # parameter ranges
        hdup = []
        for par in self.parameters:
            t = Table()
            t.add_column(Table.Column(self.param_values[par],
                                      name=str.upper(par)))
            hdup.append(fits.BinTableHDU(np.array(t), name=str.upper(par)))

        hdus = [hdu0, hdu1, hdu2]
        hdus.extend(hdup)
        hdulist = fits.HDUList(hdus)
        hdulist.writeto(file, overwrite=overwrite)

    def fill_log_fnujy_sr_hashed(self):
        """Get a hashed copy of log10 fnujy_sr."""

        # log10 just the positive values, others have tiny value
        # using the out keyword avoids a bug in numpy 1.13.0
        pos = self.fnujy_sr > 0.0
        tmp = np.log10(self.fnujy_sr, where=pos,
                       out=np.zeros(self.fnujy_sr.shape)+np.log10(cfg.tiny))

        self.log_fnujy_sr_hashed = utils.hashable(tmp)

    def rginterpolator(self):
        """Return a regular grid interpolator.
            
        Memoizing doesn't appear to save any time.
        
        This was the chunk of code in fnujy below:
        
        # scipy.RegularGridInterpolator, save a bit of time since the
        # interpolator object is the same for each model. first we
        # create grid points we want, each row is just the different
        # filter numbers assigned above with the parameters appended
#        pargrid = np.tile(par, len(wave_arr)).reshape((len(wave_arr), len(par)))
#        pargrid = np.insert(pargrid, 0, wave_arr, axis=1)
#        f = self.rginterpolator()
#        fluxes = f(pargrid)

        """
    
        if isinstance(self, PhotModel):
            wave_arr = np.arange(len(self.filters))
        elif isinstance(self, SpecModel):
            wave_arr = np.arange(len(self.wavelength))
 
        points = (wave_arr, )
        for param in self.parameters:
            points += (self.param_values[param], )
        
        f = RegularGridInterpolator(points, self.fnujy_sr,
                                    bounds_error=False, fill_value=np.inf)
        return f

    def fnujy(self, param):
        """Return fluxes for a model with a specific solid angle

        Parameter in addition to those specificed for the model is the
        log10 of the area in steradian in units of Solar radii at 1pc,
        appended.
        
        This is spline interpolation. This doesn't matter too much since
        any high dynamic range parameters are already log spaced.

        .. todo:: this is the core of the sdf code in terms of execution
        time. Experiments so far find that map_coordinates is faster
        than RegularGridInterpolator, but is hindered somewhat by the
        need to do interpolation first to find the grid points needed by
        map_coordinates. np.interp seems to be faster than scipy
        UnivariateSpline or simple 1pt interpolation for this step.
            
        """

        # prepend this to the interpolation to return the results at all
        # filters/wavelengths
        wave_arr = self.i
        nwav = self.n_i

        # make sure par is a numpy array
        par_len = len(param)-1

        area_sr = cfg.ssr * 10**(param[-1])
        par = param[:par_len]

        # scipy.ndimage.map_coordinates, only real difference compared
        # to RegularGridInerpolator is that the coordinates are given
        # in pixels, so must be interpolated from the parameters first
        # using a homegrown 1pt interpolation linterp was no faster than
        # np.interp
        coords = []
        for i, p in enumerate(self.parameters):
            coords.append(np.interp(par[i], self.param_values[p],
                                    np.arange(len(self.param_values[p]))))

        pargrid = np.tile(np.array(coords), nwav).reshape((nwav, par_len))
        pargrid = np.insert(pargrid, 0, wave_arr, axis=1)

        # interpolation, sped up by doing spline_filter first and
        # memoizing the result, order must be the same in both calls
        ff = utils.spline_filter_mem(self.log_fnujy_sr_hashed, order=2)
        fluxes = map_coordinates(ff, pargrid.T, order=2, prefilter=False)
        # convert back to real fluxes
        fluxes = 10**fluxes

        # hack to avoid negative fluxes arising from ringing
        fluxes[fluxes < cfg.tiny] = cfg.tiny

        # per-filter normalisation for photometry (leave colours)
        if isinstance(self, PhotModel):
            filt = filter.iscolour(tuple(self.filters.tolist()))
            norm = np.zeros(len(self.filters)) + area_sr
            if np.any(filt):
                norm[filt] = 1.0
        else:
            norm = area_sr

        return norm * fluxes

    def copy(self):
        """Return a copy"""
        
        return copy.deepcopy(self)

    def param_shape(self):
        """Get the shape of the parameter values"""
        
        dim = ()
        for param in self.parameters:
            dim += (len(self.param_values[param]),)
        return dim

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = utils.validate_string(value)

    @property
    def parameters(self):
        return self._parameters

    @parameters.setter
    def parameters(self, value):
        self._parameters = utils.validate_1d(value, None, dtype=str)
    
    @property
    def param_values(self):
        return self._param_values

    @param_values.setter
    def param_values(self, value):
        self._param_values = utils.validate_dict(value)

    @property
    def filters(self):
        return self._filters

    @filters.setter
    def filters(self, value):
        self._filters = utils.validate_1d(value, None, dtype=str)

    @property
    def wavelength(self):
        return self._wavelength

    @wavelength.setter
    def wavelength(self, value):
        self._wavelength = utils.validate_1d(value, None)
    
    @property
    def fnujy_sr(self):
        return self._fnujy_sr

    @fnujy_sr.setter
    def fnujy_sr(self, value):
        if self.parameters is None:
            expected_dim = None
        else:
            expected_dim = len(self.parameters) + 1  # extra for area_sr
        value = utils.validate_nd(value, expected_dim)
        if self.param_values is not None:
            if isinstance(self, PhotModel):
                if value.shape[0] != len(self.filters):
                    raise utils.SdfError("expected {} elements in first dim, got {}".
                                         format(len(self.filters), value.shape[0]))
            elif isinstance(self, SpecModel):
                if value.shape[0] != len(self.wavelength):
                    raise utils.SdfError("expected {} elements in first dim, got {}".
                                         format(len(self.wavelength), value.shape[0]))
            for i, key in enumerate(self.parameters):
                if len(self.param_values[key]) != value.shape[i+1]:
                    raise utils.SdfError("expected dimension {} to have size {}\
                                    but got {}".format(i,
                                                       len(self.param_values[key]),
                                                       value.shape[i+1]))
        self._fnujy_sr = value


class PhotModel(Model):
    """Class to hold model grids with convolved fluxes
    
    The main array is a cube with n+1 dimensions, where n is the number
    of parameters for a given model, found in the ConvolvedModel class.
    The dimensions are [nf, p1, p2, ..., pn], where px is a parameter range
    and nf is the number of filters.
    
    Colour_bases provides info where a filter is a colour/index, in the
    form of an array of dicts, each of which holds an array for the
    *relative* indices (i.e. i_filter - i_colour) locating the base
    filters for that colour, and an array of their additive weights
    (when in mags). This info is not stored for now, but populated when
    the PhotModel is read in. """

    def __init__(self, name=None, parameters=None, param_values=None,
                 filters=None, colour_bases=None, fnujy_sr=None):
        self.name = name
        self.parameters = parameters
        self.param_values = param_values
        self.filters = filters
        self.colour_bases = colour_bases
        self.fnujy_sr = fnujy_sr

    def write_model(self, name, overwrite=False):
        """Write PhotModel as a FITS file.
            
        The location to write to is given by config. Directory is 
        created if it doesn't exist.
        
        Parameters
        ----------
        name : str
            The name of the model to write, this dictates location and
            name of the file as [name]_PhotModel.fits.
        overwrite : bool, optional
            Force overwrite of extant file.
        """
        
        dir_ = cfg.file['model_root']+'/'+name+'/'
        if not os.path.exists(dir_):
            os.mkdir(dir_)
        
        self.write_file(dir_+name+'_PhotModel.fits', overwrite=overwrite)
        
    def read_model(name):
        """Read a named model, location given in config"""
        
        self = Model.read_file(cfg.model_loc[name]+name+'_PhotModel.fits')
        self.fill_colour_bases()
        return self

    @classmethod
    def cmlist2model(cls, cm):
        """Turn a list of ConvolvedModel objects into a Model."""
        
        self = cls()
        
        # check for consistency as we go
        self.name = cm[0].name
        self.parameters = cm[0].parameters
        self.param_values = cm[0].param_values
        filters = np.array([])
        
        cubedim = np.append(len(cm), cm[0].param_shape())
        cube = np.ndarray(cubedim, dtype=float)
        for i in range(len(cm)):
            filters = np.append(filters, cm[i].filter)
            if not self.name == cm[i].name:
                raise utils.SdfError("name {} in not the same as {} ".
                                     format(cm[i].name, self.name))
            if not np.all(self.parameters == cm[i].parameters):
                raise utils.SdfError("parameters {} not the same as {}".
                                     format(cm[i].parameters, self.parameters))
            for param in cm[i].parameters:
                if not np.all(np.equal(self.param_values[param],
                                       cm[i].param_values[param])):
                    raise utils.SdfError("parameter {} values {} not the same as {}".
                                         format(param, cm[i].param_values[param],
                                                self.param_values[param]))
            cube[i] = cm[i].fnujy_sr
        
        self.filters = filters
        self.fnujy_sr = cube
        return self

    @classmethod
    def read_convolved_models(cls, name, filters='all'):
        """Load all model files for a given set of filters."""
        
        # full models, avoid reading if they exist
        models = glob.glob(cfg.model_loc[name]+name+'_*Model.fits')
        
        # array of models, one element for each filter
        files = glob.glob(cfg.model_loc[name]+'*fits')
        for model in models:
            if model in files:
                files.remove(model)
        cm = [convolve.ConvolvedModel() for i in range(len(files))]
        for i, f in enumerate(files):
            cm[i] = convolve.ConvolvedModel.read_file(f)
        
        self = PhotModel.cmlist2model(cm)
        if filters != 'all':
            self.keep_filters(filters)
        return self

    def fill_colour_bases(self):
        """Fill in colour_bases info.

        Fill an attribute called colour_bases which is a dict pointing
        to where base filters for colours/indices are. The indices are
        relative to the colour locations.
        
        See Also
        --------
        model.PhotModel.keep_filters
        """
    
        self.colour_bases = [[] for i in self.filters]
        for k, f in enumerate(self.filters):
            if filter.iscolour(f):
                col = filter.Colour.get(f)
                filteri = []
                for i, cf in enumerate(col.filters):
                    fi = np.where(cf == self.filters)[0][0]
                    filteri.append(fi - k)
                self.colour_bases[k] = {'filteri': filteri,
                                        'filterw': col.weights}

    def keep_filters(self, filternames, colour_bases=False):
        """Keep only desired filters from a PhotModel.
            
        Parameters
        ----------
        filternames : list
            A list of the filter names to keep from the model. Duplicate
            filters are kept, as is the order, as each slice in the
            model must line up with the corresponding one from the
            Photometry that was read in.
        
        colour_bases : bool, optional
            Keep the base filters that are used to compute
            colours/indices. These are added to the end of the model,
            beyond the filters that were asked for.
            
        See Also
        --------
        model.PhotModel.fill_colour_bases
        """
        
        keep = np.array([], dtype=bool)
        extras = np.array([])
        for f in filternames:
            
            # grab base filters for colours
            if filter.iscolour(f):
                col = filter.Colour.get(f)
                extras = np.append(extras, col.filters)
            
            if f in self.filters:
                fi = np.where(f == self.filters)[0]
                keep = np.append(keep, fi)
            else:
                raise utils.SdfError("filter {} not found in PhotModel. "
                                     "This probably means the PhotModel "
                                     "needs to be updated (using "
                                     "model_setup.setup_phot()).".format(f))

        # now add the base filters
        if len(extras) > 0 and colour_bases:
            for f in extras:
                if f in self.filters:
                    fi = np.where(f == self.filters)[0]
                    if fi not in keep:
                        keep = np.append(keep, fi)
        
        self.filters = self.filters[keep]
        self.fnujy_sr = self.fnujy_sr[keep]
        if colour_bases:
            self.fill_colour_bases()
        else:
            self.colour_bases = []

        self.i = np.arange(len(self.filters))
        self.n_i = len(self.i)

        # and update the hashed version
        self.fill_log_fnujy_sr_hashed()
            

class SpecModel(Model):
    """Class to hold grids of model spectra
        
    The main array is a cube with n+1 dimensions, where n is the number
    of parameters for a given model. The dimensions are
    [nw, p1, p2, ..., pn], where px is a parameter range and nw is the
    number of wavelengths.

    """
    
    def __init__(self, name=None, parameters=None, param_values=None,
                 wavelength=None, fnujy_sr=None):
        self.name = name
        self.parameters = parameters
        self.param_values = param_values
        self.wavelength = wavelength
        self.fnujy_sr = fnujy_sr

    def write_model(self, name, overwrite=False):
        """Write SpecModel as a FITS file.
            
        The location to write to is given by config. Directory is 
        created if it doesn't exist.
        
        Parameters
        ----------
        name : str
            The name of the model to write, this dictates location and
            name of the file as [name]_PhotModel.fits.
        overwrite : bool, optional
            Force overwrite of extant file.
        """
        
        dir_ = cfg.file['model_root']+'/'+name+'/'
        if not os.path.exists(dir_):
            os.mkdir(dir_)
        
        self.write_file(dir_ + name+'_SpecModel.fits', overwrite=overwrite)

    def read_model(name):
        """Read a named model, location given in config."""
        
        return Model.read_file(cfg.model_loc[name]+name+'_SpecModel.fits')
    
    @classmethod
    def read_kurucz(cls, file):
        """Read a Kurucz model grid file and return a SpecModel.
            
        The model grid will almost certainly not be filled out
        completely, so require some cropping before it can be used.
        
        """
        self = cls()
        
        # get the spectra, all have the same wavelength grid
        m, teff, logg, mh = spectrum.ModelSpectrum.read_kurucz(file)
        teffarr = np.unique(teff)
        loggarr = np.unique(logg)
        
        self.name = m[0].name
        self.wavelength = m[0].wavelength
        self.parameters = ['Teff', 'logg']
        self.param_values = {'Teff': teffarr,
                             'logg': loggarr}
        
        # put spectra in their place
        self.fnujy_sr = np.zeros((len(self.wavelength),
                                  len(teffarr),
                                  len(loggarr)), dtype=float)
        for i, mod in enumerate(m):
            j = np.where(teff[i] == teffarr)[0][0]
            k = np.where(logg[i] == loggarr)[0][0]
            self.fnujy_sr[:, j, k] = mod.fnujy_sr
        
        # see if the grid was filled (spectrum.read_kurucz sets any
        # zero values in the spectra to cfg.tiny)
        if np.min(self.fnujy_sr) < cfg.tiny:
            print("WARNING: model grid not filled, spectra with zeros exist")

        return self

    @classmethod
    def bb_disk_r(cls, name='bb_disk_r',
                  wavelengths=cfg.models['default_wave'],
                  temperatures=10**np.arange(0, 3.31, 0.1),
                  lam0=None, beta=None,
                  write=False, overwrite=False):
        """Generate a set of blackbody spectra."""
        
        # don't do the calculation if there will be a write error
        if write and overwrite is False:
            if os.path.exists(cfg.model_loc[name]+name+'.fits'):
                raise utils.SdfError("{} exists, will not overwrite".
                                     format(cfg.model_loc[name]+name+'.fits'))
    
        self = cls()

        m = [spectrum.ModelSpectrum.bnu_wave_micron(wavelengths, t,
                                                    lam0=lam0,
                                                    beta=beta)
             for t in temperatures]

        self.name = m[0].name
        self.wavelength = m[0].wavelength
        if 'star' in name:
            self.parameters = ['Teff']
            self.param_values = {'Teff': temperatures}
        else:
            self.parameters = ['log_Temp']
            self.param_values = {'log_Temp': np.log10(temperatures)}

        # put spectra in their place
        self.fnujy_sr = np.zeros((len(self.wavelength),
                                  len(temperatures)), dtype=float)
        for i, mod in enumerate(m):
            self.fnujy_sr[:, i] = mod.fnujy_sr

        if write:
            self.write_model(name, overwrite=overwrite)

        return self

    @classmethod
    def modbb_disk_r(cls, name='modbb_disk_r',
                     wavelengths=cfg.models['default_wave'],
                     temperatures=10**np.arange(0, 3.01, 0.1),
                     lam0=10**np.arange(1, 3.01, 0.1),
                     beta=np.arange(0, 3.01, 0.1),
                     write=False, overwrite=False):
        """Generate a set of modified blackbody spectra"""
        
        # don't do the calculation if there will be a write error
        if write and overwrite is False:
            if os.path.exists(cfg.model_loc[name]+name+'.fits'):
                raise utils.SdfError("{} exists, will not overwrite".
                                     format(cfg.model_loc[name]+name+'.fits'))
    
        self = cls()

        self.fnujy_sr = np.zeros((len(wavelengths),
                                  len(temperatures),
                                  len(lam0),
                                  len(beta)), dtype=float)
        for i, temp in enumerate(temperatures):
            for j, l0 in enumerate(lam0):
                for k, b in enumerate(beta):
                    m = spectrum.ModelSpectrum.bnu_wave_micron(wavelengths, temp,
                                                               lam0=l0, beta=b)
                    self.fnujy_sr[:, i, j, k] = m.fnujy_sr

        self.name = m.name
        self.wavelength = m.wavelength
        self.parameters = ['log_Temp', 'log_lam0', 'beta']
        self.param_values = {'log_Temp': np.log10(temperatures),
                             'log_lam0': np.log10(lam0),
                             'beta': beta}

        if write:
            self.write_model(name, overwrite=overwrite)

        return self

    @classmethod
    def modbb_disk_dr(cls, name='modbb_disk_dr',
                      wavelengths=cfg.models['default_wave'],
                      t_in_min=200.0, t_in_max=2000.0, n_t_in=20,
                      t_out_min=10.0, t_out_max=100.0, n_t_out=10,
                      alpha=np.linspace(-2, 2, 17),
                      beta=np.linspace(0, 2, 9),
                      write=False, overwrite=False):
        """Generate a set of wide-disk modified blackbody spectra.

        The disk area is normalised so that the full area is one. Thus
        very wide disks (those with cool outer edges) are much fainter
        in terms of model flux than smaller ones. 

        Parameters
        ----------
        name : str, optional
            Name of the model.
        wavelengths : array, optional
            Array of wavelengths for the model.
        t_in_min : float, optional
            Minimum inner edge temperature.
        t_in_max : float, optional
            Maximum inner edge temperature.
        n_t_in : int, optional
            Number of inner temperatures.
        t_out_min : float, optional
            Minimum outer edge temperature.
        t_out_max : float, optional
            Maximum inner edge temperature.
        n_t_out : int, optional
            Number of outer temperatures.
        alpha : array, optional
            Array of power law indices for optical depth.
        beta : array, optional
            Array of betas, lambda_0's fixed near blackbody peak.
        write : bool, optional
            Write the model to disk.
        overwrite : bool, optional
            Overwrite any existing model.
        """
            
        # don't do the calculation if there will be a write error
        if write and overwrite is False:
            if os.path.exists(cfg.model_loc[name]+name+'.fits'):
                raise utils.SdfError("{} exists, will not overwrite".
                                     format(cfg.model_loc[name]+name+'.fits'))
    
        self = cls()

        # set up temperature arrays
        log_t_in = np.linspace(np.log10(t_in_min), np.log10(t_in_max), n_t_in)
        t_in = 10**log_t_in
        log_t_out = np.linspace(np.log10(t_out_min), np.log10(t_out_max), n_t_out)
        t_out = 10**log_t_out

        self.fnujy_sr = np.zeros((len(wavelengths),
                                  n_t_in, n_t_out, len(alpha),
                                  len(beta)), dtype=float)

        # loop to fill model
        n_r = 100
        for i, t1 in enumerate(t_in):
            for j, t2 in enumerate(t_out):
                for k, a in enumerate(alpha):
                    for l, b in enumerate(beta):
                    
                        tau_tot = 0.0
                        # generate temps and a range of pseudo radii
                        t_edges = 10**np.linspace(np.log10(t1),
                                                  np.log10(t2),
                                                  n_r+1, dtype=float)
                        temps = (t_edges[1:]+t_edges[:-1])/2.
                        radius = lambda x_: (278.3/x_)**2
                        r = radius(temps)
                        r_edges = radius(t_edges)
                        dr = np.diff(r_edges)
                        tau = r**a
                        
                        for x in range(n_r):
                            m = spectrum.ModelSpectrum.bnu_wave_micron(
                                                      wavelengths,
                                                      temps[x],
                                                      lam0=3*2900.0/temps[x],
                                                      beta=b
                                                      )
                            annulus_area = tau[x] * 2.0*np.pi*r[x]*dr[x]
                            spec = m.fnujy_sr * annulus_area
                            self.fnujy_sr[:, i, j, k, l] += spec
                            tau_tot += annulus_area

                        # normalise area for all models
                        self.fnujy_sr[:, i, j, k, l] /= tau_tot

        self.name = m.name
        self.wavelength = m.wavelength
        self.parameters = ['log_T_in', 'log_T_out', 'alpha', 'beta']
        self.param_values = {'log_T_in': log_t_in, 'log_T_out': log_t_out,
                             'alpha': alpha, 'beta': beta}

        if write:
            self.write_model(name, overwrite=overwrite)

        return self

    @classmethod
    def sd_spectra(cls, name='sd_disk_r',
                   wavelengths=cfg.models['default_wave'],
                   temperatures=10**np.arange(0, 3.01, 0.1),
                   smin=10**np.arange(-1, 2.01, 0.1),
                   q=np.arange(1.67, 2.001, 0.02),
                   smax=100000, nsz=100,
                   write=False, overwrite=False):
        """Generate a set of size distribution models.

        Simple analytic grain model after Backman & Paresce, assumes
        that grains have Qabs that is 1 for s<pi.lambda, and Qabs
        decreasing as lambda^-n beyond. What n is depends on the dust
        properties, but it appears to be >1 because some disks have
        (sub)mm slopes steeper than Fnu oc nu^3. Models look like they
        have n~2, e.g. Draine astrosilicate.

        The sub-mm slopes are not as steep as can be obtained from the
        real grain models. This is something to do with the details of
        the absorption/emission efficiencies.

        This models assumes that the peak wavelength of the stellar
        spectrum peaks at lambda shorter than the grain size to estimate
        the temperatures, which is a bit suspect but necessary for an
        analytic solution. The results are a weak function [T^(n/(4+n)]
        of the stellar temperature anyway.
        """

        # constants
        xt = np.pi  # turnover in Qabs, pi is like "real" temperatures
        n = 2.0     # slope of Qabs beyond xt, 2 is like "real" dust
        cw = 5100.  # peak of blackbody emission in micron/K
        ts = 6000.  # assumed stellar temperature

        self = cls()

        self.fnujy_sr = np.zeros((len(wavelengths),
                                  len(temperatures),
                                  len(smin),
                                  len(q)), dtype=float)
        for i, tbb in enumerate(temperatures):
            for j, smini in enumerate(smin):
                for k, qi in enumerate(q):

                    # sizes
                    s = 10**np.linspace(np.log10(smini),
                                        np.log10(smax), nsz)
                    logs = np.log10(s)

                    # calculate grain temperatures
                    dbb = cw / tbb / xt
                    dsm = cw / ts / xt
                    sm = s < dsm
                    bb = s > dbb
                    s_temp = tbb * (dbb/s)**(n/(4+n))
                    s_temp[sm] = tbb**(4/(4+n)) * ts**(n/(4+n))
                    s_temp[bb] = tbb

                    # compute bnu for each size
                    bnu = np.ones((len(s), len(wavelengths)))
                    for l, st in enumerate(s_temp):
                        bnu[l, :] = utils.bnu_wav_micron(wavelengths, st)

                    # qabs
                    qabs = np.ones((len(s), len(wavelengths)))
                    for l, si in enumerate(s):
                        x = wavelengths / si
                        gtx = x > xt
                        qabs[l, gtx] = (xt/x[gtx])**n

                    # add up size distribution, this is taken straight
                    # from Wyatt's IDL sigmadbar
                    qfact = 5 - 3*qi
                    sigmadbar = qfact * np.log(10) * \
                        (10**(logs*qfact)) / \
                        (smax**qfact - smini**qfact)

                    for l in range(len(wavelengths)):
                        self.fnujy_sr[l, i, j, k] = \
                            utils.sdf_int(qabs[:, l]*bnu[:, l]*sigmadbar, logs)

#                    return wavelengths, s, s_temp, bnu, qabs, sigmadbar, self.fnujy_sr[:]

        self.name = 'sd'
        self.wavelength = wavelengths
        self.parameters = ['log_Temp', 'log_Dmin', 'q']
        self.param_values = {'log_Temp': np.log10(temperatures),
                             'log_Dmin': np.log10(smin), 'q': q}

        if write:
            self.write_model(name, overwrite=overwrite)

        return self

    def interp_to_wavelengths(self, wavelength, log=True):
        """Interpolate the model to the given wavelengths.

        .. todo:: this only needs to be run at the beginning of a fit but
        is very slow, especially when there are spectra, speed it up!

        .. todo:: this is straight linear/log interpolation, but could
        smooth the spectra first since the given wavelength grid will
        almost certainly be near the spectral resolution of whatever
        instrument it came from. Probably use resample.
        """

        # check we need to do something
        if np.all(self.wavelength == wavelength):
            return

        if log:
            neg = self.fnujy_sr <= 0.0
            self.fnujy_sr[neg] = cfg.tiny
            cube = np.log10(self.fnujy_sr)
            wave = np.log10(self.wavelength)
            wave_interp = np.log10(wavelength)
        else:
            cube = self.fnujy_sr
            wave = self.wavelength
            wave_interp = wavelength
        
        # get a function that will return interpolated values, ensure
        # error if extrapolation in wavelength requested (i.e. model
        # doesn't cover as wide as was requested)
        f = interp1d(wave, cube, axis=0, kind='linear',
                     bounds_error=True)
        cube_interp = f(wave_interp)

        self.wavelength = wavelength
        if log:
            self.fnujy_sr = np.power(10, cube_interp)
        else:
            self.fnujy_sr = cube_interp

        self.i = np.arange(len(self.wavelength))
        self.n_i = len(self.i)

        # and update the hashed version
        self.fill_log_fnujy_sr_hashed()


def model_fluxes(m, param, obs_nel, phot_only=False):
    """Get model fluxes and put in arrays.
    
    all_fnu is everything added up, with colours/indices added properly,
    comp_fnu[i] contains fluxes from the i-th model component and
    comp_fnu_col[i] contains these with colours computed.
    """
    
    comp_fnu = []
    all_fnu = []
    i0 = 0
    # loop over model components
    for comp in m:
        # loop over phot/spectra for this component if they exist
        if not isinstance(comp, tuple):
            comp = (comp, )
        flux = np.array([])
        # params same for all in each component
        nparam = len(comp[0].parameters)+1
        for mod in comp:
            if phot_only:
                if not isinstance(mod, PhotModel):
                    continue
            fnu = mod.fnujy(param[i0:i0+nparam])
            flux = np.append(flux, fnu)
        
        # since we don't know how long all_fnu will be, make sure
        # comp_fnu has first dimension equal number of components 
        if len(all_fnu) == 0:
            all_fnu = flux
            comp_fnu = np.array([flux])
        else:
            all_fnu += flux
            comp_fnu = np.vstack((comp_fnu, flux))
        i0 += nparam
    
    # fill colours, for total and components
    mod_fnu = fill_colours(m[0], all_fnu, obs_nel)
    comp_fnu_col = np.zeros((len(comp_fnu), len(mod_fnu)))
    for i, fnu in enumerate(comp_fnu):
        comp_fnu_col[i] = fill_colours(m[0], fnu, obs_nel)

    return mod_fnu, comp_fnu_col


def crop(m, param, range_):
    """Crop a model to ranges specified for a parameter."""

    out = m.copy()
    
    if param not in out.parameters:
        raise utils.SdfError("parameter {} not in model (has {})".format(param, out.parameters))

    # get axis to cut, and locations
    ax = np.where(out.parameters == param)[0][0] + 1
    locs = np.searchsorted(out.param_values[param], range_)

    print("cutting parameters along axis {} to indices {}".format(ax, locs))

    out.param_values[param] = out.param_values[param][locs[0]:locs[1]]
    arrin = out.fnujy_sr
    arr = np.rollaxis(out.fnujy_sr, ax)
    arr = arr[locs[0]:locs[1]]
    out.fnujy_sr = np.rollaxis(arr, 0, ax+1)
    print("cropped model from {} to {}".
          format(arrin.shape, out.fnujy_sr.shape))

    return out


def reduce_zerod(m, parameters):
    """Reduce a model to zero dimensions (i.e. a spectrum).
        
    Parameters
    ----------
    m : model object
        The model to reduce.
    parameters : list
        Parameter values, excluding the last (normalisation) parameter.
    """

    out = m.copy()

    # modify the attributes
    out.parameters = []
    out.param_values = {}

    # get the new spectrum and update the hashed version
    out.fnujy_sr = m.fnujy(np.append(parameters, -np.log10(cfg.ssr)))
    out.fill_log_fnujy_sr_hashed()
    
    return out


def reduce_squeeze(m):
    """Reduce model by removing length=one dimensions."""

    out = m.copy()
    print('input model has shape:{}'.format(m.param_shape()))

    keep = np.array(m.param_shape()) != 1

    for k, p in zip(keep, m.parameters):
        if not k:
            del out.param_values[p]
    out.parameters = m.parameters[keep]
    out.fnujy_sr = np.squeeze(m.fnujy_sr)

    print('output model has shape:{}'.format(out.param_shape()))
    return out


def append_parameter(m, name, value):
    """Append a single parameter to a model.
        
    Purpose is to prepare a model for addition models via concat.
    
    """
    out = m.copy()
    out.parameters = np.append(out.parameters, name)
    out.param_values[name] = np.array([value])
    out.fnujy_sr = np.reshape(out.fnujy_sr, out.fnujy_sr.shape+(1, ))
    return out


def concat(m0, m):
    """Add a model to the one we have.
        
    So far can only add models when both have the same sets of
    parameters, so for example joining two models with different
    metallicities.
        
    """

    out = m0.copy()

    # check types, parameters, wavelengths, filters are the same, allow
    # for small differences in wavelengths, which can apparently occur
    # when the arrays are calculated on different machines
    for i, par in enumerate(out.parameters):
        if par != m.parameters[i]:
            raise utils.SdfError("parameters {} and {} different".
                                 format(out.parameters, m.parameters))
    
    if type(out) != type(m):
        raise utils.SdfError("can't join models of type {} and {}".
                             format(type(out), type(m)))
    
    if isinstance(out, PhotModel):
        for i, filt in enumerate(out.filters):
            if filt != m.filters[i]:
                raise utils.SdfError("filters {} and {} different".
                                     format(out.filters, m.filters))
    
    if isinstance(out, SpecModel):
        if not np.allclose(out.wavelength, m.wavelength, rtol=1e-12, atol=1e-12):
            raise utils.SdfError("wavelengths {} and {} different".
                                 format(out.wavelength, m.wavelength))

    # parameters to add and their locations in out
    padd = []
    pax = []
    ploc = []
    for i, p in enumerate(out.parameters):
        if not np.all(np.equal(out.param_values[p], m.param_values[p])):
            pax.append(i+1)  # +1 since first dim is wav/filters
            padd.append(p)
            arrs = np.split(m.fnujy_sr, len(m.param_values[p]), axis=i)
            for val in m.param_values[p]:
                ploc.append(np.searchsorted(out.param_values[p], val))

    if len(padd) != 1:
        raise utils.SdfError("model parameters can't be joined (padd={})".
                             format(padd))
    else:
        padd = padd[0]

    print("Adding parameter {} at location(s) {} along axis {}".
          format(padd, ploc, pax))

    for i, loc in enumerate(ploc):
        
        if m.param_values[padd][i] in out.param_values[padd]:
            raise utils.SdfError("model already has {}={}".
                                 format(padd, m.param_values[padd][i]))

        print("  adding {}={} (dims {} to {}) at {}".
              format(padd, m.param_values[padd][i],
                     arrs[i].squeeze().shape, out.fnujy_sr.shape, loc))

        out.param_values[padd] = np.insert(out.param_values[padd],
                                           loc, m.param_values[padd][i])
        out.fnujy_sr = np.insert(out.fnujy_sr, loc,
                                 arrs[i].squeeze(), axis=pax[i])
            
    print("  new {} array is {}".format(padd, out.param_values[padd]))

    return out


def fill_colours(comp, mod_fnu, obs_nel):
    """Replace locations in mod_fnu with colours
        
    Colours cannot simply be added without knowing the absolute
    measurements. The model components have already been added so we
    only need to figure out where the filters associated with the
    colours are, and calculate the colours from these.
    
    By using obs_nel, the number of observed fluxes/colours per
    Phot/SpecModel component, the extra columns containing the base
    filters for colours are not included in the returned result.

    """
    final_fnu = np.array([])
    if not isinstance(comp, tuple):
        comp = (comp, )
    i0comp = 0  # zeroth index of the current Phot/SpecModel
    for k, mod in enumerate(comp):
        if isinstance(mod, PhotModel):
            # loop over each filter/index in the model, skip
            # if there are no colour_bases (i.e. a filter)
            for i, cb in enumerate(mod.colour_bases):
                if len(cb) > 0:
                    mags = 0.0
                    # loop over the filters in this colour/index
                    # and add the magnitude with the correct weight
                    # for this colour/index
                    for j, filteri in enumerate(cb['filteri']):
                        fname = mod.filters[i+filteri]
                        irel = i0comp+i+cb['filteri'][j]
                        filt = filter.Filter.get(fname)
                        mag = filt.flux2mag(mod_fnu[irel])
                        mags += mag * cb['filterw'][j]
                    mod_fnu[i] = mags
            final_fnu = np.append(final_fnu,
                                  mod_fnu[i0comp:i0comp+obs_nel[k]])
            i0comp += len(mod.filters)
        else:
            final_fnu = np.append(final_fnu,
                                  mod_fnu[i0comp:i0comp+obs_nel[k]])
            i0comp += len(mod.wavelength)

    return final_fnu


def get_models(obs, names):
    """Get tuples of models for given observations.
        
    Returns two tuples of models, the first is for fitting and contains
    extra filters beyond the set in the observations if there are
    colours in the photometry. The second exludes the extra filters, and
    has a max of one spectrum (for plotting purposes).

    """
    allmod = ()
    fullmod = ()
    for name in names:
       
        omod = ()
        fmod = ()
        
        # load models we will need
        ph = PhotModel.read_model(name)
        sp = SpecModel.read_model(name)

        for o in obs:
        
            if isinstance(o, photometry.Photometry):
                phmod = ph.copy()
                phmod.keep_filters(o.filters, colour_bases=True)
                omod = omod + (phmod, )
                phmod = ph.copy()
                phmod.keep_filters(o.filters, colour_bases=False)
                fmod = fmod + (phmod, )
            
            if isinstance(o, spectrum.ObsSpectrum):
                spmod = sp.copy()
                spmod.interp_to_wavelengths(o.wavelength)
                omod = omod + (spmod, )
                
        # always want a full spectrum in full model, but only one
        fmod = fmod + (sp, )

        allmod = allmod + (omod, )
        fullmod = fullmod + (fmod, )

    return allmod, fullmod


def models_info(m):
    """Return some info for a tuple of model(s).
        
    Also do some basic sanity checking.
    
    Here is where the ranges for model and spectra normalisation is set.
    
    """
    info = {'name': '',
            'ndim': 0,
            'ncomp': [],
            'nspec': [],
            'type': [],
            'p_rng': [],
            'parameters': [],
            'nmodels': len(m)}
    for comp in m:
        nspec = 0
        if not isinstance(comp, tuple):
            comp = (comp, )
        info['ncomp'].append(len(comp))
        if info['name'] == '':
            info['name'] = comp[0].name
        else:
            info['name'] += cfg.fitting['model_join']+comp[0].name
        info['ndim'] += len(comp[0].parameters)+1
        for par in comp[0].parameters:
            info['p_rng'].append((comp[0].param_values[par][0],
                                  comp[0].param_values[par][-1]))
            info['parameters'].append(par)
        # this is the range of allowed solid angles
        info['p_rng'].append(cfg.fitting['model_om_range'])
        info['parameters'].append('norm')
        for mod in comp:
            if isinstance(mod, SpecModel):
                nspec += 1
            info['type'].append(type(mod))
        info['nspec'].append(nspec)
    info['ndim'] += info['nspec'][0]

    # spectra normalisations last (nspec same for each comp)
    for i in range(nspec):
        info['p_rng'].append(cfg.fitting['spectra_norm_range'])
        info['parameters'].append('spec_norm')

    # check structure looks OK
    if len(np.unique(info['ncomp'])) > 1:
        raise utils.SdfError("model structure {} should have same number of\
                             subcomponents in each component, not {}".
                             format(m, info['ncomp']))
    if len(np.unique(info['nspec'])) > 1:
        raise utils.SdfError("model structure {} should have same number of spectra\
                             in each component, not {}".format(m, info['nspec']))

    return info
