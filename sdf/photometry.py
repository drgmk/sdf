import numpy as np
import astropy.units as u
from astropy.table import Table

from . import filter
from . import utils
from .utils import SdfError
from . import config as cfg

class Photometry(object):
    """Photometry class.

    This class is used for observations and synthetic photometry. It contains
    some meta info and lists of photometry.
    """
    def __init__(self,filters=None,measurement=None,e_measurement=None,unit=None,
                 s_measurement=None,bibcode=None,upperlim=None,ignore=None,
                 fnujy=None,e_fnujy=None):
        self.filters = filters
        self.measurement = measurement
        self.e_measurement = e_measurement
        self.s_measurement = s_measurement
        self.unit = unit
        self.bibcode = bibcode
        self.upperlim = upperlim
        self.ignore = ignore
        self.fnujy = fnujy
        self.e_fnujy = e_fnujy

    def addto(self,p):
        """Add Photometry object to another."""

        if self.nphot is None:
            self.filters = p.filters
            self.measurement = p.measurement
            self.e_measurement = p.e_measurement
            self.s_measurement = p.s_measurement
            self.unit = p.unit
            self.bibcode = p.bibcode
            self.upperlim =  p.upperlim
            self.ignore = p.ignore
            self.fnujy = p.fnujy
            self.e_fnujy = p.e_fnujy
        elif p.nphot is None:
            pass
        else:
            self.filters = np.append( self.filters, p.filters)
            self.measurement = np.append( self.measurement, p.measurement)
            self.e_measurement = np.append(self.e_measurement,p.e_measurement)
            self.s_measurement = np.append(self.s_measurement,p.s_measurement)
            self.unit = np.append( self.unit, p.unit)
            self.bibcode = np.append( self.bibcode, p.bibcode)
            self.upperlim = np.append( self.upperlim, p.upperlim)
            self.ignore = np.append( self.ignore, p.ignore)
            self.fnujy = np.append( self.fnujy, p.fnujy)
            self.e_fnujy = np.append( self.e_fnujy, p.e_fnujy)

                    
    @classmethod
    def read_sdb_file(cls,file,keep_filters=None):
        """ Load photometry from a file and return a Photometry object.

        The file format is set by sdb_getphot.py, and is ascii.ipac. To
        get any spectra use spectrum.ObsSpectrum.get_sdb_file, and for
        keywords use utils.get_sdb_keywords.

        """

        self = cls()

        # get the photometry, checking there is any
        phot = Table.read(file,format='ascii.ipac')
        if len(phot) == 0:
            return None

        # loop over the rows and fill the object if we're
        # keeping this filter
        for i in range(len(phot)):
            
            if keep_filters is not None:
                if phot[i]['Band'] not in keep_filters:
                    continue

            row = phot[i]

            # quick sanity checking
            for par in ['Phot','Err','Sys']:
                if not np.isfinite(row[par]):
                    SdfError("non-finite {} value {} in {}".format(par,row[par],file))
            
            p = Photometry()
            p.filters = np.array([row['Band']])
            p.measurement   = np.array([row['Phot']])
            p.e_measurement = np.abs( np.array([row['Err']]) )
            p.s_measurement = np.abs( np.array([row['Sys']]) )
            p.unit = np.array([u.Unit(row['Unit'])])
            p.bibcode = np.array([row['bibcode']])
            p.upperlim = np.array([ row['Lim'] == 1 ])
            p.ignore = np.array([ row['exclude'] == 1 ])
            
            if row['Band'] in cfg.fitting['exclude_filters']:
                p.ignore = np.array([True])

            self.addto(p)

        self.fill_fnujy()
        self.sort()
        return self
    

    @property
    def nphot(self):
        if self.filters is not None:
            return len(self.filters)
        else:
            return None


    @property
    def nused(self):
        return np.sum(self.ignore == False)


    def fill_fnujy(self):
        """Convert measurements to Jy and fill fnujy/e_fnujy arrays.
            
        Where no uncertainties are given 10% (or 0.1mag) is assumed.
        Colours/indices are left as is.
        
        """
        fnu = np.zeros(self.nphot)
        efnu = np.zeros(self.nphot)
        strom = [-1,-1,-1] # indices of b-y, m1, c1
        for i in range(self.nphot):

            # attempt to combine uncertainties
            if np.isfinite(self.e_measurement[i]) and np.isfinite(self.s_measurement[i]):
                etot = np.sqrt(self.e_measurement[i]**2 + self.s_measurement[i]**2)
            elif np.isfinite(self.e_measurement[i]) and not np.isfinite(self.s_measurement[i]):
                etot = self.e_measurement[i]
            elif not np.isfinite(self.e_measurement[i]) and np.isfinite(self.s_measurement[i]):
                etot = self.s_measurement[i]
            else:
                print("WARNING no uncertainties given, assuming nan")
                etot = np.nan

            # convert flux and uncertainty
            if self.unit[i] != 'mag':
                fnu[i] = (self.measurement[i]*self.unit[i]).to('Jy').value
                if np.isfinite(etot) and etot > 0.:
                    efnu[i] = (etot * self.unit[i]).to('Jy').value
                else:
                    efnu[i] = 0.1 * fnu[i] # assume 10%
            else:
                # use zero point to convert
                if not filter.iscolour(self.filters[i]):
                    filt = filter.Filter.get(self.filters[i])
                    fnu[i] = filt.mag2flux(self.measurement[i])
                    if np.isfinite(etot) and etot > 0.:
                        efnu[i] = fnu[i] * etot / 1.09 # small uncertainties trick
                    else:
                        efnu[i] = 0.1 * fnu[i] # assume 10%
                # leave colours/indices as is
                else:
                    fnu[i] = self.measurement[i]
                    if np.isfinite(etot) and etot > 0.:
                        efnu[i] = etot
                    else:
                        efnu[i] = 0.1 # assume 10%
                    
                    # note stromgren
                    if self.filters[i] == 'BS_YS' and strom[0] == -1:
                         strom[0] = i
                    if self.filters[i] == 'STROMM1' and strom[1] == -1:
                         strom[1] = i
                    if self.filters[i] == 'STROMC1' and strom[2] == -1:
                         strom[2] = i

        # convert uvby as a group (lowest possible sum of indices = 0+1+2)
        if np.sum(strom) > 2:
            fnu[strom[0]],fnu[strom[1]],fnu[strom[2]] = \
                    utils.uvby_convert(fnu[strom[0]],fnu[strom[1]],fnu[strom[2]])

        self.fnujy = fnu
        self.e_fnujy = efnu

        
    def mean_wavelength(self):
        """Return the mean wavelength of the filters."""
        
        mw = np.array([])
        for f in self.filters:
            if filter.iscolour(f):
                col = filter.Colour.get(f)
                mw = np.append( mw, col.mean_wavelength )
            else:
                filt = filter.Filter.get(f)
                mw = np.append( mw, filt.mean_wavelength )
        return mw

    
    def sort(self):
        """Sort arrays in increasing wavelength order."""
        
        _,srt = np.unique( self.mean_wavelength(), return_index=True )
        self.filters = self.filters[srt]
        self.measurement = self.measurement[srt]
        self.e_measurement = self.e_measurement[srt]
        self.s_measurement = self.s_measurement[srt]
        self.unit = self.unit[srt]
        self.bibcode = self.bibcode[srt]
        self.upperlim =  self.upperlim[srt]
        self.ignore = self.ignore[srt]
        self.fnujy = self.fnujy[srt]
        self.e_fnujy = self.e_fnujy[srt]

        
    @property
    def measurement(self):
        return self._measurement
    @measurement.setter
    def measurement(self, value):
        self._measurement = utils.validate_1d(value,self.nphot)

    @property
    def e_measurement(self):
        return self._e_measurement
    @e_measurement.setter
    def e_measurement(self, value):
        self._e_measurement = utils.validate_1d(value,self.nphot)

    @property
    def s_measurement(self):
        return self._s_measurement
    @s_measurement.setter
    def s_measurement(self, value):
        self._s_measurement = utils.validate_1d(value,self.nphot)

    @property
    def unit(self):
        return self._unit
    @unit.setter
    def unit(self, value):
        self._unit = utils.validate_1d(value,self.nphot,dtype=u.Unit)

    @property
    def bibcode(self):
        return self._bibcode
    @bibcode.setter
    def bibcode(self, value):
        self._bibcode = utils.validate_1d(value,self.nphot,dtype=str)


