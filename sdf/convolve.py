import numpy as np
import astropy.units as u
from . import filter
from . import spectrum
from . import utils
from . import config as cfg

c_micron = u.micron.to(u.Hz,equivalencies=u.spectral())

class ConvolvedModel(object):
    """Class holding fluxes in one filter for a set of models

    The models have n parameters, each of which spans some range
    given by arrays in the param_values dict. The fluxes are in an
    n-dimensional array. Filter is a string with the name.
    """

    def __init__(self,name=None,parameters=None,param_values=None,
                 filter=None,fnujy_sr=None):
        self.name = name
        self.parameters = parameters
        self.param_values = param_values
        self.filter = filter
        self.fnujy_sr = fnujy_sr


    @classmethod
    def read_file(cls,file):
        """Read in a specific model flux file"""
        self = cls()

        from astropy.io import fits
        from astropy.table import Table

        self = cls()

        # parameter names
        fh = fits.open(file)
        keywords = fh[0].header
        type = keywords['SDFTYPE']
        if type != 'convolvedmodel':
            raise utils.SdfError("file {} not a convolvedmodel, is a {}".format(file,type))

        self.name = keywords['NAME']
        self.filter = keywords['FILTER']
        nparam = keywords['NPARAM']
        self.parameters = np.array([keywords['PARAM'+str(i)] for i in range(nparam)])

        # parameter ranges, assume that all values have the same dtype
        # so it's OK if the resulting ndarray does
        d = {}
        for i,par in enumerate(self.parameters):
            j = i+2
            if str.upper(par) != fh[j].name:
                raise utils.SdfError("{}th parameter {} not equal to HDU with name {}".format(j,par,fh[j].name))
            dat = fh[j].data
            d[par] = np.array(dat,dtype=dat.dtype[0])
        self.param_values = d

        # main array of convolved model fluxes, do this last so the
        # dimensionality is checked by the setter
        self.fnujy_sr = fh[1].data
        header1 = fh[1].header
        if header1['BUNIT'] != (u.jansky/u.sr).to_string(format='fits'):
            fnuunit = header1['BUNIT']
            self.fnujy_sr = self.fnujy_sr * u.Unit(fnuunit).to('Jy / sr')

        return self
    

    def write_file(self,file,overwrite=False):
        """Write fluxes to a FITS file
        
        Number of parameters and their names are stored in the primary 
        HDU, flux density cube is in the first HDU, and arrays with the
        parameter ranges in subsequend HDUs.
        """

        from astropy.io import fits
        from astropy.table import Table

        # primary HDU (for metadata)
        hdu0 = fits.PrimaryHDU()
        hdu0.header['SDFTYPE'] = 'convolvedmodel'
        hdu0.header['NAME'] = self.name
        hdu0.header['FILTER'] = self.filter
        hdu0.header['NPARAM'] = len(self.parameters)
        # keywords for parameters
        for i,par in enumerate(self.parameters):
            hdu0.header['PARAM'+str(i)] = par

        # convolved fluxes
        hdu1 = fits.ImageHDU(self.fnujy_sr, name='CONVOLVED MODEL')
        hdu1.header['BUNIT'] = (u.jansky/u.sr).to_string(format='fits')

        # parameter ranges
        hdup = []
        for par in self.parameters:
            t = Table()
            t.add_column( Table.Column(self.param_values[par],name=str.upper(par)) )
            hdup.append( fits.BinTableHDU(np.array(t),name=str.upper(par)) )

        hdus = [hdu0,hdu1]
        hdus.extend(hdup)
        hdulist = fits.HDUList(hdus)
        hdulist.writeto(file, overwrite=overwrite)

        
    @classmethod
    def bb(cls,filtername,temperatures):
        """Return a ConvolvedModel object for black bodies
            
        Set temperature range on a log scale to make it easier for 
        fitting methods to find low values.
        """
        self = cls()
        self.name = 'bb'
        self.parameters = ['log_Temp']
        self.param_values = {'log_Temp':np.log10(temperatures)}
        
        self.filter = filtername
        if not filter.iscolour(filtername):
            f = filter.Filter.get(filtername)
            wav = c_micron / f.nu_hz
        else:
            col = filter.Colour.get(filtername)
            wav = np.array([])
            for f in col.filters:
                filt = filter.Filter.get(f)
                wav = np.append(wav,c_micron / filt.nu_hz)
            wav = np.unique(wav)

        fnujy_sr = np.array([],dtype=float)
        for temp in temperatures:
            spec = spectrum.ModelSpectrum.bnu_wave_micron(wav,temp)
            conv,cc = spec.synthphot(filtername)
            fnujy_sr = np.append(fnujy_sr,conv)

        self.fnujy_sr = fnujy_sr
        return self


    @classmethod
    def modbb(cls,filtername,temperatures,lam0,beta):
        """Return a ConvolvedModel object for modified black bodies
            
        Set temperature and lam0 range on a log scale to make it easier
        for fitting methods to find low values.
        """
        
        self = cls()
        self.name = 'modbb'
        self.parameters = ['log_Temp','log_lam0','beta']
        self.param_values = {'log_Temp':np.log10(temperatures)}
        self.param_values['log_lam0'] = np.log10(lam0)
        self.param_values['beta'] = beta

        self.filter = filtername
        if not filter.iscolour(filtername):
            f = filter.Filter.get(filtername)
            wav = c_micron / f.nu_hz
        else:
            col = filter.Colour.get(filtername)
            wav = np.array([])
            for f in col.filters:
                filt = filter.Filter.get(f)
                wav = np.append(wav,c_micron / filt.nu_hz)
            wav = np.unique(wav)
        
        fnujy_sr = np.zeros((len(temperatures),
                             len(lam0),
                             len(beta)),dtype=float)
        for i,temp in enumerate(temperatures):
            for j,l0 in enumerate(lam0):
                for k,b in enumerate(beta):
                    spec = spectrum.ModelSpectrum.bnu_wave_micron(wav,temp,
                                                                  lam0=l0,beta=b)
                    conv,cc = spec.synthphot(filtername)
                    fnujy_sr[i,j,k] = conv

        self.fnujy_sr = fnujy_sr
        return self


    def param_shape(self):
        """Get the shape of the parameter values
        """
        dim = ()
        for param in self.parameters:
            dim += ( len(self.param_values[param]), )
        return dim


    @property
    def name(self):
        return self._name
    @name.setter
    def name(self,value):
        self._name = utils.validate_string(value)

    @property
    def parameters(self):
        return self._parameters
    @parameters.setter
    def parameters(self,value):
        self._parameters = utils.validate_1d(value,None,dtype=str)

    @property
    def param_values(self):
        return self._param_values
    @param_values.setter
    def param_values(self,value):
        self._param_values = utils.validate_dict(value)

    @property
    def filter(self):
        return self._filter
    @filter.setter
    def filter(self,value):
        self._filter = utils.validate_string(value)

    @property
    def fnujy_sr(self):
        return self._fnujy_sr
    @fnujy_sr.setter
    def fnujy_sr(self,value):
        if self.parameters is None:
            expected_dim = None
        else:
            expected_dim = len(self.parameters)
        value = utils.validate_nd(value,expected_dim)
        if self.param_values is not None and value is not None:
            for i,key in enumerate(self.parameters):
                if len(self.param_values[key]) != value.shape[i]:
                    raise utils.SdfError("expected dimension {} to have size {} but got {}".format(i,len(self.param_values[key]),value.shape[i]))
        self._fnujy_sr = value
