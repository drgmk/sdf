"""Filter bandpassses, zero points, and offsets.

Overview
--------

Bandpasses
^^^^^^^^^^

Many bandpasses are from Mann & von Braun (2015, MvB,
2015PASP..127..102M), which sought a self-consistent set of bandpasses
and zero points, mainly for "heritage" photometric systems that include
observations of bright stars. These filters are all in units of energy,
so are integrated directly (also commonly called relative system
response, or RSRs).

Synthetic photometry
^^^^^^^^^^^^^^^^^^^^

The Spitzer colour correction tests pass in comparison with those on the
SSC pages, with precision of a few percent or better (worst for IRAC),
so the correction to the fluxes is different by a few 10^-4 and good
enough. The differences appear to be in the details of the interpolation
and integration, but exactly where isn't clear (e.g.  np.trapz and
scipy.integrate.simps give the same results). Anyway, the success of
these tests suggests that the synthetic photometry works.

Absolute calibration
^^^^^^^^^^^^^^^^^^^^

To ensure a consistent absolute calibration all zero points are re-
derived using the CALSPEC Vega spectrum, as this has been used by many
authors to derive zero point magnitudes (i.e. the magnitude of Vega in a
given photometric system, also known as the zero point offset). In
converting magnitudes to flux densities these offsets are subtracted, so
a positive offset yields a brighter flux for a given magnitude.

Optical
~~~~~~~

In the optical (shortward of 1 micron) the absolute calibration finds
Vega to be slightly fainter than zero magnitude. Systems still use Vega
for their zero point, but include an offset, which is the magnitude of
Vega in each filter, and is 0.027 in Johnson V.

The zero points are derived from the CALSPEC Vega spectrum, and offsets
have been derived by various authors e.g. Maiz-Appelaniz (2006), Bessell
& Murphy (2012), Mann & von Braun (2015).

In MvB there are no zero point offsets, so the zero points can be used 
to derive ZPOs relative to the CALSPEC Vega spectrum. Their zero points
tend to be larger numbers, meaning that the ZPO (i.e. Vega) is a more
positive magnitude.

In most cases the zero point offsets have been re-derived with a
minimisation process that attempts to reach a happy medium across all
filters.


Infrared
~~~~~~~~

In the IR zero magnitude is defined by Vega, but again small offsets may
be needed. The CALSPEC spectrum appears to have the same IR flux as the
one in Cohen et al 1992 (i.e. gives same fluxes for 2MASS in Cohen et al
2003 and IRAS in Cohen et al 1992). The Vega spectrum in the 1992 paper
underpins most IR photometric systems (apart from the one proposed by
Rieke et al 2008, which is ~2% brighter). Thus, zero point offsets in
the IR shold be minimal. So while the K magnitude for Vega was found to
-0.36 by Rieke et al 2008, the K magnitude of "Vega" used for
calibration purposes is zero by definition.

Filters
-------

This file contains all the information specific to all bandpasses in a
dictionary called filters. Each entry has the name used by that bandpass
(e.g. UJ for Johnson U), and is itself a dictionary containing the
necessary subset of infomration for filter.Filter.get to construct that
filter.

Most of the information comes from the Spanish VO filter service at
http://svo2.cab.inta-csic.es/theory/fps3/, but is supplemented by extra
information to avoid various errors (e.g. mainly if a given bandpass is
in photon or energy units, but sometimes zero points and always zero
point offsets). The specific treatment of a given filter depends on
provenance and use.

An additional list of "dummy" filters is given, containing i) generic
bandpasses centered on a specific wavelength, and ii) colours/indices
that are combinations of other filters. These contain little actual
information aside from the name, which uses the convention of an
underscore to form a colour (e.g. "UJ_BJ" is Johnson U-B). Practically
this is done using the filter.Colour object, which contains a list of
filters and weights for a given colour/index that says how to add the
magnitudes.

These are roughly in order of wavelength.

.. todo:: all to be updated when ZPOs are sorted

GALEX
  Filters are from SVO, and are effective area so can be considered
  photon counting. Catalogue supplies fluxes and magnitudes, so for now
  the former is used and no zero point needed. As-yet no evidence that
  this system is or isn't consistent with "Vega". Magnitude system is 
  AB. Generally don't include in fitting due to the possibility of a UV
  excess.

Johnson
  Filters are from MvB. Photometry in this system comes from Mermilliod
  (UBV means), Bessell (UBV, also Cousins RI), Koen (UBV).

Stromgren
  Filters are from Bessell 2011 (MvB exist), used because conversions
  from that paper are also used. Using these corrections didn't make any
  obvious differences.
  
  Vega has non-zero magnitudes in this system (e.g. Gray 1998), and 
  measurements in this system are only in colours/indices so the 
  absolute calibration is not important.
  
Hipparcos/Tycho
  Filters are from MvB. Photometry in this system comes from the 
  original Hipparcos catalogue, and the Tycho 2 catalogue. Tycho-2 is
  taken to have the "correct" absolute calibration in deriving others.

Gaia
  Filters are from SVO.

Cousins
  Filters are from MvB. Photometry in this system comes from Bessell.
  Haven't been able to find which system most Gliese (CNS3) RI photmetry
  is in, so not included.

Kepler
  Filter from SVO. No photometry from this system used, so filter exists
  to make predictions.
  
DDO
  Filter from SVO, KPNO/Mosaic.D51 assumed to be similar to that used
  for the survey of the Kepler field.
  
2MASS
  Filters from Cohen (via SVO). Photometry from the 2MASS PSC, separated
  into read1 (e.g. 2MR1H) and read2 (e.g. 2MR2H) after suggestion in 
  Rieke (2008) that there is a few percent offset for stars observed in 
  both modes. No evidence for this offset seen in SED fitting (yet).
  This system is used as one of those with the "correct" absolute
  calibration.

DENIS
  Filters from Fouque et al 2000. Calibration assumes Vega is zero mag.
  No zero point offsets.

Spitzer/IRAC
  Filters from SSC (via SVO), need to be converted to energy. Photometry
  is from SEIP. No zero point offsets used.

Spitzer/IRS PUI
  Filters from SSC (via SVO), need to be converted to energy. No
  photometry in this system currently used.

Spitzer/MIPS
  Filters from SSC (via SVO), in photon units. Zero points derived from
  "Vega", where "Vega" has been slightly rescaled (as returned by
  spectrum.ObsSpectrum.vega() to ensure 7.17Jy as 23.675um. Assume
  70 micron calibration used in published papers good enough (i.e.
  photometric uncertainties much larger than calibration unceratinty).
 
WISE
  Filters from SVO in energy units. Photometry from ALLWISE. Zero point
  offset needed for W3.

AKARI/IRC
  Filters from SVO in energy units. Photometry from IRC PSC. Catalogue
  has fluxes so zero points not needed. Calibration uses Cohen network
  so should have similar offsets to IRAC. Looks to be a systematic
  change in S18/S09 as a function of S18 (or S09) for stars in narrow
  B-V range, which could be the cause of apparent warm excesses arising
  from only S18 photometry.

MSX
  Filters from SVO as RSRs. No photometry in this system currently
  included.

LBTI/NOMIC
  Filter from Denis Defrere, assumed QE-based so needs to be converted
  to photon units. No photometry, used for predictions.

IRAS
  Filters from SVO in energy units. Photometry from PSC and FSC. At 12
  micron Rieke conclude that the fluxes need to be reduced by a factor
  of 0.992, and by 0.98 at 25 micron. Original calibration used.

Herschel/PACS
  Filters from SVO in energy units. Use fluxes as published or in GMK's
  personal catalogue, assume calibration uncertainty much smaller than
  photometric uncertainty.

Herschel/SPIRE
  Filters from SVO in energy units. Use fluxes as published or in GMK's
  personal catalogue, assume calibration uncertainty much smaller than
  photometric uncertainty.

"""

import os
import glob
import numpy as np
import astropy.units as u
from . import utils

c_micron = u.micron.to(u.Hz,equivalencies=u.spectral())


filters = {}

# dummy filters for colours/indices or artificial bandpasses
# that we want in the "all" list
extras = ['BS_YS','STROMM1','STROMC1',
          'UJ_BJ','BJ_VJ','VJ_IC','VJ_RC','RC_IC',
          'WAV350','WAV450','WAV610',
          'WAV800','WAV850','WAV860','WAV870','WAV880',
          'WAV1100','WAV1200','WAV1240','WAV1250',
          'WAV1300','WAV1330','WAV1350',
          'WAV2000','WAV2700',
          'WAV3000','WAV3190','WAV3300',
          'WAV6800',
          'WAV9000']

for f in extras:
    filters[f] = {
                  'magnitude_system':     None,
                  'zero_point':           None,
                  'zero_point_ref':       None,
                  'ref_wavelength':       None,
                  'ref_spectrum':         None,
                  'response_type':        None,
                  'response_ref':         None,
                  'wav_micron':           None,
                  'response':             None
                  }

# many filters can come from SVO, but photon/energy and probably zero
# points suspect, so do manually here
# http://svo2.cab.inta-csic.es/theory/fps3/

# GALEX responses are effective area, so equivalent
# to photon counting. Need to get these in AB system
filters['GALFUV'] = {'svo_name': 'PhotCalID=GALEX/GALEX.FUV/AB',
                     'response_type': 'photon'}
filters['GALNUV'] = {'svo_name': 'PhotCalID=GALEX/GALEX.NUV/AB',
                     'response_type': 'photon'}

# SDSS, zpos untested
filters['USDSS'] = {'svo_name': 'SLOAN/SDSS.u',
                    'magnitude_system': 'AB',
                    'zero_point_offset': 0.04,
                    'response_type': 'energy'}
filters['GSDSS'] = {'svo_name': 'SLOAN/SDSS.g',
                    'magnitude_system': 'AB',
                    'zero_point_offset': 0.01,
                    'response_type': 'energy'}
filters['RSDSS'] = {'svo_name': 'SLOAN/SDSS.r',
                    'magnitude_system': 'AB',
                    'zero_point_offset': 0.01,
                    'response_type': 'energy'}
filters['ISDSS'] = {'svo_name': 'SLOAN/SDSS.i',
                    'magnitude_system': 'AB',
                    'zero_point_offset': 0.01,
                    'response_type': 'energy'}
filters['ZSDSS'] = {'svo_name': 'SLOAN/SDSS.z',
                    'magnitude_system': 'AB',
                    'zero_point_offset': 0.0,
                    'response_type': 'energy'}

# APASS the same filters as SDSS, slight changes in zero point offsets
filters['GAPASS'] = {'svo_name': 'SLOAN/SDSS.g',
                     'magnitude_system': 'AB',
                     'zero_point_offset': 0.005,
                     'response_type': 'energy'}
filters['RAPASS'] = {'svo_name': 'SLOAN/SDSS.r',
                     'magnitude_system': 'AB',
                     'zero_point_offset': 0.005,
                     'response_type': 'energy'}
filters['IAPASS'] = {'svo_name': 'SLOAN/SDSS.i',
                     'magnitude_system': 'AB',
                     'zero_point_offset': -0.01,
                     'response_type': 'energy'}

# zero point offsets from Bessel & Murphy 2012 are 0.04, 0.022, 0.027.
# from MvB are 0.0188, 0.0185, 0.027 (latter fixed to 0.027)
# consider Johnson V magnitude of Vega to be immutable at 0.027
filters['UJ'] = {'svo_name': 'GCPD/Johnson.U',
                 'zero_point_offset': -0.04,
                 'response_type': 'energy'}
filters['BJ'] = {'svo_name': 'GCPD/Johnson.B',
                 'zero_point_offset': 0.058,
                 'response_type': 'energy'}
filters['VJ'] = {'svo_name': 'GCPD/Johnson.V',
                 'zero_point_offset': 0.027,
                 'response_type': 'energy'}

# APASS, uses Landolt standards
filters['BL'] = {'svo_name': 'GCPD/Johnson.B_Landolt',
                     'zero_point_offset': 0.008,
                     'response_type': 'energy'}
filters['VL'] = {'svo_name': 'GCPD/Johnson.V_Landolt',
                     'zero_point_offset': 0.017,
                     'response_type': 'energy'}
filters['BAPASS'] = filters['BL']
filters['VAPASS'] = filters['VL']

# zero point offsets from Bessel & Murphy 2012 are 0.027,0.028
# from MvB are 0.0212, 0.0091. Note that their R bandpasses look a bit
# different, so perhaps expect different ZPOs
filters['RC'] = {'svo_name': 'GCPD/Cousins.R',
                 'zero_point_offset': 0.047,
                 'response_type': 'energy'}
filters['IC'] = {'svo_name': 'GCPD/Cousins.I',
                 'zero_point_offset': 0.035,
                 'response_type': 'energy',
                 'zero_point': 2510.0}             # wrong on SVO

# Stromgren (uvby) from Maiz-Appelaniz 2006 is 1.435, 0.182, 0.021,
# 0.014. Comparing Mann & von Braun and CALSPEC zero points gives 1.401,
# 0.175, 0.0256, 0.031. Gray 1998 gives 1.445, 0.195, 0.034, 0.03 and
# GCPD is almost exactly the same. Bessel gives coefficients to convert
# between observed and synthetic photometry, having assumed the GCPD
# zero points, which are a function of b-y, so needs to be implemented
# elsewhere (i.e. file read time)
filters['US'] = {#'svo_name': 'GCPD/Stromgren.u',
                 'magnitude_system': 'Vega',
                 'zero_point_offset': 1.435, #1.257,
                 'response_type': 'photon',
                 # Bessell 2011 responses, photon counting
                 'wav_micron': [0.3150, 0.3175, 0.3200, 0.3225, 0.3250,
0.3275, 0.3300, 0.3325, 0.3350, 0.3375, 0.3400, 0.3425, 0.3450, 0.3475,
0.3500, 0.3525, 0.3550, 0.3575, 0.3600, 0.3625, 0.3650, 0.3675, 0.3700,
0.3725, 0.3750, 0.3775, 0.3800, 0.3825, 0.3850],
                 'response': [0.000, 0.004, 0.050, 0.122, 0.219, 0.341,
0.479, 0.604, 0.710, 0.809, 0.886, 0.939, 0.976, 1.000, 0.995, 0.981,
0.943, 0.880, 0.782, 0.659, 0.525, 0.370, 0.246, 0.151, 0.071, 0.030,
0.014, 0.000, 0.000]
                 }
filters['VS'] = {#'svo_name': 'GCPD/Stromgren.v',
                 'magnitude_system': 'Vega',
                 'zero_point_offset': 0.182, #0.272,
                 'response_type': 'photon',
                 'wav_micron': [0.3750, 0.3775, 0.3800, 0.3825,
0.3850, 0.3875, 0.3900, 0.3925, 0.3950, 0.3975, 0.4000, 0.4025, 0.4050,
0.4075, 0.4100, 0.4125, 0.4150, 0.4175, 0.4200, 0.4225, 0.4250, 0.4275,
0.4300, 0.4325, 0.4350, 0.4375, 0.4400, 0.4425, 0.4450],
                 'response': [0.000, 0.003, 0.006, 0.016, 0.029, 0.044,
0.060, 0.096, 0.157, 0.262, 0.404, 0.605, 0.810, 0.958, 1.000, 0.973,
0.882, 0.755, 0.571, 0.366, 0.224, 0.134, 0.079, 0.053, 0.039, 0.027,
0.014, 0.006, 0.000]
                 }
filters['BS'] = {#'svo_name': 'GCPD/Stromgren.b',
                 'magnitude_system': 'Vega',
                 'zero_point_offset': 0.021, #0.055,
                 'response_type': 'photon',
                 'wav_micron': [0.4350, 0.4375, 0.4400, 0.4425, 0.4450,
0.4475, 0.4500, 0.4525, 0.4550, 0.4575, 0.4600, 0.4625, 0.4650, 0.4675,
0.4700, 0.4725, 0.4750, 0.4775, 0.4800, 0.4825, 0.4850, 0.4875, 0.4900,
0.4925, 0.4950, 0.4975, 0.5000, 0.5025, 0.5050],
                 'response': [0.000, 0.010, 0.023, 0.039, 0.056, 0.086,
0.118, 0.188, 0.287, 0.457, 0.681,  0.896, 0.998, 1.000, 0.942, 0.783,
0.558, 0.342, 0.211, 0.130, 0.072, 0.045,  0.027, 0.021, 0.015, 0.011,
0.007, 0.003, 0.000]
                 }
filters['YS'] = {#'svo_name': 'GCPD/Stromgren.y',
                 'magnitude_system': 'Vega',
                 'zero_point_offset': 0.014, #0.03,
                 'response_type': 'photon',
                 'wav_micron': [0.5150, 0.5175, 0.5200, 0.5225, 0.5250,
0.5275, 0.5300, 0.5325, 0.5350, 0.5375, 0.5400, 0.5425, 0.5450, 0.5475,
0.5500, 0.5525, 0.5550, 0.5575, 0.5600, 0.5625, 0.5650, 0.5675, 0.5700,
0.5725, 0.5750, 0.5775, 0.5800, 0.5825, 0.5850],
                 'response': [0.000, 0.022, 0.053, 0.082, 0.116, 0.194,
0.274, 0.393, 0.579, 0.782, 0.928, 0.985, 0.999, 1.000, 0.997, 0.938,
0.789, 0.574, 0.388, 0.232, 0.143, 0.090, 0.054, 0.031, 0.016, 0.010,
0.009, 0.004, 0.000]
                 }

# Skymapper
filters['SkyMapper.u'] = {'svo_name': 'SkyMapper/SkyMapper.u',
                          'zero_point_offset': 0.0,
                          'magnitude_system': 'AB',
                          'response_type': 'energy'}
filters['SkyMapper.v'] = {'svo_name': 'SkyMapper/SkyMapper.v',
                          'zero_point_offset': 0.0,
                          'magnitude_system': 'AB',
                          'response_type': 'energy'}
filters['SkyMapper.g'] = {'svo_name': 'SkyMapper/SkyMapper.g',
                          'zero_point_offset': 0.0,
                          'magnitude_system': 'AB',
                          'response_type': 'energy'}
filters['SkyMapper.r'] = {'svo_name': 'SkyMapper/SkyMapper.r',
                          'zero_point_offset': 0.0,
                          'magnitude_system': 'AB',
                          'response_type': 'energy'}
filters['SkyMapper.i'] = {'svo_name': 'SkyMapper/SkyMapper.i',
                          'zero_point_offset': 0.0,
                          'magnitude_system': 'AB',
                          'response_type': 'energy'}
filters['SkyMapper.z'] = {'svo_name': 'SkyMapper/SkyMapper.z',
                          'zero_point_offset': 0.0,
                          'magnitude_system': 'AB',
                          'response_type': 'energy'}


# zero point offsets from Bessel & Murphy 2012 (0.03,0.023,0.038)
# from MvB 0.0232, 0.0118, 0.0196
filters['BT'] = {'svo_name': 'TYCHO/TYCHO.B_MvB',
                 'zero_point_offset': 0.03,
                 'response_type': 'energy'}
filters['VT'] = {'svo_name': 'TYCHO/TYCHO.V_MvB',
                 'zero_point_offset': 0.023,
                 'response_type': 'energy'}
filters['HP'] = {'svo_name': 'Hipparcos/Hipparcos.Hp_MvB',
                 'zero_point_offset': 0.032,
                 'response_type': 'energy'}

# Gaia. Bandpasses appear to be photon counting according to Evans+2012,
# but are multiplied by lambda in SVO to convert to energy counting.
# equation 2. DR2 assumes Vega has flux 3.66e-9 at 550nm, but CALSPEC
# spectrum has 3.54e-9, so DR2 magnitudes ~3% too bright
filters['GAIA.G'] = {'svo_name': 'GAIA/GAIA2r.G',
                     'zero_point_offset': 0.0,
                     'response_type': 'energy'}
filters['GAIA.BP'] = {'svo_name': 'GAIA/GAIA2r.Gbp',
                      'zero_point_offset': 0.0,
                      'response_type': 'energy'}
filters['GAIA.RP'] = {'svo_name': 'GAIA/GAIA2r.Grp',
                      'zero_point_offset': 0.0,
                      'response_type': 'energy'}

# Kepler, assume photon and same zero point offset as V (0.027)
filters['KP'] = {'svo_name': 'Kepler/Kepler.K',
                 'zero_point_offset': 0.027,
                 'response_type': 'photon'}
filters['D51'] = {'svo_name': 'KPNO/Mosaic.D51',
                  'zero_point_offset': 0.027,
                  'response_type': 'photon'}

# 2MASS, assume no difference between R1/R2 until evidence otherwise,
# zero point offsets from Cohen+2003 (who assume Vega=0)
filters['2MJ'] = {'svo_name': '2MASS/2MASS.J',
                 'zero_point_offset': -0.001 - 0.03,
                  'response_type': 'energy'}
filters['2MH'] = {'svo_name': '2MASS/2MASS.H',
                 'zero_point_offset': 0.019 + 0.005,
                  'response_type': 'energy'}
filters['2MKS'] = {'svo_name': '2MASS/2MASS.Ks',
                   'zero_point_offset': -0.017 + 0.01,
                   'response_type': 'energy'}
filters['2MR1J'] = filters['2MJ']
filters['2MR1H'] = filters['2MH']
filters['2MR1KS'] = filters['2MKS']
filters['2MR2J'] = filters['2MJ']
filters['2MR2H'] = filters['2MH']
filters['2MR2KS'] = filters['2MKS']

# DENIS, assume energy counting
filters['IDENIS'] = {'svo_name': 'DENIS/DENIS.I',
                     'zero_point_offset': 0.0,
                     'response_type': 'energy'}
filters['JDENIS'] = {'svo_name': 'DENIS/DENIS.J',
                     'zero_point_offset': -0.02,
                     'response_type': 'energy'}
filters['KSDENIS'] = {'svo_name': 'DENIS/DENIS.Ks',
                     'zero_point_offset': -0.01,
                     'response_type': 'energy'}


# WISE, RSRs already converted to energy, ref spectrum
# is F_nu oc 1/nu^2. residuals from seds suggest small changes
# in zero point offsets, perhaps because "Vega" is fainter than Vega
# Patel+2014 find Ks-W1=0.031, and give saturated calibrations for W1/2
# empirically find W3 needs to go up a bit
# 3.3% shift in W4 bandpass recommended by 2014PASA...31...49B
def w1_cal_func(x):
    """Calibration fit for W1 from Patel+2014."""
    if x < 8.0:
        return -0.1359+0.0396*x-0.0023*x**2
    else:
        return 0.0

def w2_cal_func(x):
    """Calibration fit for W2 from Patel+2014."""
    if x < 5.3:
        return 1.5777 - 0.3495 * x + 0.016 * x**2
    elif 5.3 <= x < 6.7:
        return -0.353 + 0.8826 * x - 0.238 * x**2 + 0.017 * x**3
    else:
        return 0.0

filters['WISE3P4'] = {'svo_name': 'WISE/WISE.W1',
                      'response_type': 'energy',
                      'zero_point_offset': -0.015,
                      'measurement_calibration': lambda x: x + w1_cal_func(x),
                      'ref_wavelength': 3.3526,
                      'ref_spectrum': lambda nu: 1.0/nu/nu}
filters['WISE4P6'] = {'svo_name': 'WISE/WISE.W2',
                      'response_type': 'energy',
                      'zero_point_offset': 0.01,
                      'measurement_calibration': lambda x: x + w2_cal_func(x),
                      'ref_wavelength': 4.6028,
                      'ref_spectrum': lambda nu: 1.0/nu/nu}
filters['WISE12'] = {'svo_name': 'WISE/WISE.W3',
                     'response_type': 'energy',
                     'zero_point_offset': 0.03,
                     'ref_wavelength': 11.5608,
                     'ref_spectrum': lambda nu: 1.0/nu/nu}
filters['WISE22'] = {'svo_name': 'WISE/WISE.W4',
                     'response_type': 'energy',
                     'zero_point_offset': 0.0,
                     'ref_wavelength': 22.0883,
                     'ref_spectrum': lambda nu: 1.0/nu/nu}
#                     'magnitude_system': 'Vega',
#                     'wav_micron':[19.71997,19.7303,19.74063,19.75096,19.76129,19.77162,19.78195,19.79228,19.80261,19.81294,19.96789,19.97822,19.98855,19.99888,20.00921,20.01954,20.02987,20.0402,20.05053,20.06086,20.07119,20.08152,20.09185,20.10218,20.11251,20.12284,20.13317,20.1435,20.15383,20.16416,20.17449,20.18482,20.19515,20.20548,20.21581,20.22614,20.23647,20.2468,20.25713,20.26746,20.27779,20.28812,20.29845,20.30878,20.31911,20.32944,20.33977,20.3501,20.36043,20.37076,20.38109,20.39142,20.40175,20.41208,20.42241,20.43274,20.44307,20.4534,20.46373,20.47406,20.48439,20.49472,20.50505,20.51538,20.52571,20.53604,20.54637,20.5567,20.56703,20.57736,20.58769,20.59802,20.60835,20.61868,20.62901,20.63934,20.64967,20.66,20.67033,20.68066,20.69099,20.70132,20.71165,20.72198,20.73231,20.74264,20.75297,20.7633,20.77363,20.78396,20.79429,20.80462,20.81495,20.82528,20.83561,20.84594,20.85627,20.8666,20.87693,20.88726,20.89759,20.90792,20.91825,20.92858,20.93891,20.94924,20.95957,20.9699,20.98023,20.99056,21.00089,21.01122,21.02155,21.03188,21.04221,21.05254,21.06287,21.0732,21.08353,21.09386,21.10419,21.11452,21.12485,21.13518,21.14551,21.15584,21.16617,21.1765,21.18683,21.19716,21.20749,21.21782,21.22815,21.23848,21.24881,21.25914,21.26947,21.2798,21.29013,21.30046,21.31079,21.32112,21.33145,21.34178,21.35211,21.36244,21.37277,21.3831,21.39343,21.40376,21.41409,21.42442,21.43475,21.44508,21.45541,21.46574,21.47607,21.4864,21.49673,21.50706,21.51739,21.52772,21.53805,21.54838,21.55871,21.56904,21.57937,21.5897,21.60003,21.61036,21.62069,21.63102,21.64135,21.65168,21.66201,21.67234,21.68267,21.693,21.70333,21.71366,21.72399,21.73432,21.74465,21.75498,21.76531,21.77564,21.78597,21.7963,21.80663,21.81696,21.82729,21.83762,21.84795,21.85828,21.86861,21.87894,21.88927,21.8996,21.90993,21.92026,21.93059,21.94092,21.95125,21.96158,21.97191,21.98224,21.99257,22.0029,22.01323,22.02356,22.03389,22.04422,22.05455,22.06488,22.07521,22.08554,22.09587,22.1062,22.11653,22.12686,22.13719,22.14752,22.15785,22.16818,22.17851,22.18884,22.19917,22.2095,22.21983,22.23016,22.24049,22.25082,22.26115,22.27148,22.28181,22.29214,22.30247,22.3128,22.32313,22.33346,22.34379,22.35412,22.36445,22.37478,22.38511,22.39544,22.40577,22.4161,22.42643,22.43676,22.44709,22.45742,22.46775,22.47808,22.48841,22.49874,22.50907,22.5194,22.52973,22.54006,22.55039,22.56072,22.57105,22.58138,22.59171,22.60204,22.61237,22.6227,22.63303,22.64336,22.65369,22.66402,22.67435,22.68468,22.69501,22.70534,22.71567,22.726,22.73633,22.74666,22.75699,22.76732,22.77765,22.78798,22.79831,22.80864,22.81897,22.8293,22.83963,22.84996,22.86029,22.87062,22.88095,22.89128,22.90161,22.91194,22.92227,22.9326,22.94293,22.95326,22.96359,22.97392,22.98425,22.99458,23.00491,23.01524,23.02557,23.0359,23.04623,23.05656,23.06689,23.07722,23.08755,23.09788,23.10821,23.11854,23.12887,23.1392,23.14953,23.15986,23.17019,23.18052,23.19085,23.20118,23.21151,23.22184,23.23217,23.2425,23.25283,23.26316,23.27349,23.28382,23.29415,23.30448,23.31481,23.32514,23.33547,23.3458,23.35613,23.36646,23.37679,23.38712,23.39745,23.40778,23.41811,23.42844,23.43877,23.4491,23.45943,23.46976,23.48009,23.49042,23.50075,23.51108,23.52141,23.53174,23.54207,23.5524,23.56273,23.57306,23.58339,23.59372,23.60405,23.61438,23.62471,23.63504,23.64537,23.6557,23.66603,23.67636,23.68669,23.69702,23.70735,23.71768,23.72801,23.73834,23.74867,23.759,23.76933,23.77966,23.78999,23.80032,23.81065,23.82098,23.83131,23.84164,23.85197,23.8623,23.87263,23.88296,23.89329,23.90362,23.91395,23.92428,23.93461,23.94494,23.95527,23.9656,23.97593,23.98626,23.99659,24.00692,24.01725,24.02758,24.03791,24.04824,24.05857,24.0689,24.07923,24.08956,24.09989,24.11022,24.12055,24.13088,24.14121,24.15154,24.16187,24.1722,24.18253,24.19286,24.20319,24.21352,24.22385,24.23418,24.24451,24.25484,24.26517,24.2755,24.28583,24.29616,24.30649,24.31682,24.32715,24.33748,24.34781,24.35814,24.36847,24.3788,24.38913,24.39946,24.40979,24.42012,24.43045,24.44078,24.45111,24.46144,24.47177,24.4821,24.49243,24.50276,24.51309,24.52342,24.53375,24.54408,24.55441,24.56474,24.57507,24.5854,24.59573,24.60606,24.61639,24.62672,24.63705,24.64738,24.65771,24.66804,24.67837,24.6887,24.69903,24.70936,24.71969,24.73002,24.74035,24.75068,24.76101,24.77134,24.78167,24.792,24.80233,24.81266,24.82299,24.83332,24.84365,24.85398,24.86431,24.87464,24.88497,24.8953,24.90563,24.91596,24.92629,24.93662,24.94695,24.95728,24.96761,24.97794,24.98827,24.9986,25.00893,25.01926,25.02959,25.03992,25.05025,25.06058,25.07091,25.08124,25.09157,25.1019,25.11223,25.12256,25.13289,25.14322,25.15355,25.16388,25.17421,25.18454,25.19487,25.2052,25.21553,25.22586,25.23619,25.24652,25.25685,25.26718,25.27751,25.28784,25.29817,25.3085,25.31883,25.32916,25.33949,25.34982,25.36015,25.37048,25.38081,25.39114,25.40147,25.4118,25.42213,25.43246,25.44279,25.45312,25.46345,25.47378,25.48411,25.49444,25.50477,25.5151,25.52543,25.53576,25.54609,25.55642,25.56675,25.57708,25.58741,25.59774,25.60807,25.6184,25.62873,25.63906,25.64939,25.65972,25.67005,25.68038,25.69071,25.70104,25.71137,25.7217,25.73203,25.74236,25.75269,25.76302,25.77335,25.78368,25.79401,25.80434,25.81467,25.825,25.83533,25.84566,25.85599,25.86632,25.87665,25.88698,25.89731,25.90764,25.91797,25.9283,25.93863,25.94896,25.95929,25.96962,25.97995,25.99028,26.00061,26.01094,26.02127,26.0316,26.04193,26.05226,26.06259,26.07292,26.08325,26.09358,26.10391,26.11424,26.12457,26.1349,26.14523,26.15556,26.16589,26.17622,26.18655,26.19688,26.20721,26.21754,26.22787,26.2382,26.24853,26.25886,26.26919,26.27952,26.28985,26.30018,26.31051,26.32084,26.33117,26.3415,26.35183,26.36216,26.37249,26.38282,26.39315,26.40348,26.41381,26.42414,26.43447,26.4448,26.45513,26.46546,26.47579,26.48612,26.49645,26.50678,26.51711,26.52744,26.53777,26.5481,26.55843,26.56876,26.57909,26.58942,26.59975,26.61008,26.62041,26.63074,26.64107,26.6514,26.66173,26.67206,26.68239,26.69272,26.70305,26.71338,26.72371,26.73404,26.74437,26.7547,26.76503,26.77536,26.78569,26.79602,26.80635,26.81668,26.82701,26.83734,26.84767,26.858,26.86833,26.87866,26.88899,26.89932,26.90965,26.91998,26.93031,26.94064,26.95097,26.9613,26.97163,26.98196,26.99229,27.00262,27.01295,27.02328,27.03361,27.04394,27.05427,27.0646,27.07493,27.08526,27.09559,27.10592,27.11625,27.12658,27.13691,27.14724,27.15757,27.1679,27.17823,27.18856,27.19889,27.20922,27.21955,27.22988,27.24021,27.25054,27.26087,27.2712,27.28153,27.29186,27.30219,27.31252,27.32285,27.33318,27.34351,27.35384,27.36417,27.3745,27.38483,27.39516,27.40549,27.41582,27.42615,27.43648,27.44681,27.45714,27.46747,27.4778,27.48813,27.49846,27.50879,27.51912,27.52945,27.53978,27.55011,27.56044,27.57077,27.5811,27.59143,27.60176,27.61209,27.62242,27.63275,27.64308,27.65341,27.66374,27.67407,27.6844,27.69473,27.70506,27.71539,27.72572,27.73605,27.74638,27.75671,27.76704,27.77737,27.7877,27.79803,27.80836,27.81869,27.82902,27.83935,27.84968,27.86001,27.87034,27.88067,27.891,27.90133,27.91166,27.92199,27.93232,27.94265,27.95298,27.96331,27.97364,27.98397,27.9943,28.00463,28.01496,28.02529,28.03562,28.04595,28.05628,28.06661,28.07694,28.08727,28.0976,28.10793,28.11826,28.12859,28.13892,28.14925,28.15958,28.16991,28.18024,28.19057,28.2009,28.21123,28.22156,28.23189,28.24222,28.25255,28.26288,28.67608,28.68641,28.69674,28.70707,28.7174,28.72773,28.73806,28.74839,28.75872,28.76905,28.77938,28.78971,28.80004,28.81037,28.8207,28.83103,28.84136,28.85169,28.86202,28.87235,28.88268,28.89301,28.90334],
#                         'response':[0.00167933,0.00231167,0.00280867,0.00295900,0.00301867,0.00289700,0.00271100,0.00239167,0.00196367,0.00144667,0.00102167,0.00116533,0.00135667,0.00155100,0.00196267,0.00243967,0.00327200,0.00426667,0.00528333,0.00621000,0.00698333,0.00754333,0.00805000,0.00841333,0.00877000,0.00904667,0.00932000,0.00972000,0.0102100,0.0109567,0.0119633,0.0132333,0.0149433,0.0168667,0.0195633,0.0224100,0.0262700,0.0303333,0.0354667,0.0407667,0.0467333,0.0533667,0.0610000,0.0693000,0.0782667,0.0885667,0.0995333,0.113133,0.127067,0.143933,0.161500,0.182000,0.202800,0.225933,0.250433,0.276567,0.305033,0.334333,0.371000,0.407333,0.447333,0.487000,0.527000,0.567000,0.600333,0.633667,0.664000,0.691000,0.714667,0.738333,0.758333,0.778667,0.795333,0.809000,0.815667,0.822667,0.822667,0.823000,0.819667,0.819667,0.820000,0.820000,0.823333,0.827000,0.830333,0.833667,0.837000,0.837000,0.837333,0.840667,0.847333,0.857333,0.867333,0.880667,0.890667,0.900333,0.907000,0.907000,0.910333,0.907000,0.907000,0.903667,0.900333,0.897333,0.897333,0.890667,0.884000,0.877333,0.870667,0.864000,0.860667,0.860667,0.860667,0.867000,0.873667,0.883667,0.890333,0.897000,0.900333,0.903667,0.907000,0.910333,0.917000,0.923667,0.933667,0.940333,0.950333,0.950333,0.953667,0.950333,0.947000,0.940333,0.937000,0.933667,0.930333,0.930333,0.930333,0.927000,0.923667,0.920333,0.913667,0.910333,0.903667,0.903667,0.907000,0.913667,0.923667,0.930333,0.937000,0.940333,0.940333,0.937000,0.930333,0.923667,0.917000,0.910333,0.910333,0.910333,0.913667,0.913667,0.917000,0.920333,0.920333,0.917000,0.917000,0.917000,0.917000,0.920333,0.920333,0.923667,0.923667,0.923667,0.917000,0.910333,0.903667,0.893667,0.883667,0.880333,0.873667,0.877000,0.877000,0.883667,0.890333,0.893667,0.900333,0.907000,0.910333,0.917000,0.923667,0.930333,0.940333,0.950333,0.960000,0.970000,0.976667,0.983333,0.983333,0.983333,0.983333,0.983333,0.986667,0.986667,0.990000,0.996667,1.00000,1.00000,1.00000,0.993333,0.986667,0.976667,0.963333,0.956667,0.946667,0.940000,0.936667,0.933333,0.930000,0.926667,0.923333,0.920000,0.913333,0.910000,0.906667,0.906667,0.913333,0.920000,0.930000,0.943333,0.956667,0.966667,0.976667,0.983333,0.986667,0.990000,0.990000,0.990000,0.990000,0.990000,0.990000,0.990000,0.983333,0.980000,0.970000,0.956667,0.943333,0.926667,0.910000,0.896667,0.886667,0.880000,0.873333,0.873333,0.873333,0.873333,0.876667,0.876667,0.876667,0.876667,0.880000,0.886667,0.893333,0.906667,0.916667,0.933000,0.949667,0.959667,0.969667,0.976333,0.973000,0.969667,0.963000,0.953000,0.946333,0.939667,0.933000,0.929667,0.926333,0.923000,0.923000,0.919667,0.916333,0.909667,0.906333,0.899667,0.893000,0.890000,0.886667,0.883333,0.883333,0.880000,0.880000,0.876667,0.870000,0.866667,0.860000,0.853333,0.850000,0.843333,0.840000,0.840000,0.840000,0.840000,0.840000,0.836667,0.833333,0.830000,0.823333,0.816667,0.810000,0.800000,0.796667,0.790000,0.786667,0.786667,0.783333,0.780000,0.776667,0.770000,0.763333,0.753333,0.740000,0.726667,0.713333,0.700000,0.690000,0.683333,0.676667,0.673333,0.673333,0.670000,0.666667,0.663333,0.660000,0.653333,0.650000,0.646667,0.643333,0.643333,0.643333,0.650000,0.653333,0.660000,0.663333,0.670000,0.670000,0.670000,0.670000,0.666667,0.663333,0.663333,0.663333,0.663333,0.666667,0.670000,0.673333,0.676667,0.676667,0.673333,0.673333,0.666667,0.660000,0.650000,0.643333,0.636667,0.633333,0.630000,0.630000,0.630000,0.633333,0.633333,0.636667,0.640000,0.639667,0.639667,0.639667,0.643000,0.643000,0.649667,0.653000,0.659667,0.666333,0.673000,0.676333,0.679667,0.679667,0.676333,0.676333,0.673000,0.666333,0.663000,0.659667,0.656667,0.653333,0.653333,0.650000,0.643333,0.640000,0.633333,0.623333,0.616667,0.603333,0.590000,0.580000,0.566667,0.556667,0.546667,0.540000,0.533333,0.530000,0.523333,0.516667,0.513333,0.506667,0.500000,0.493333,0.486667,0.480000,0.473333,0.470000,0.466667,0.463333,0.460000,0.460000,0.456667,0.453333,0.450000,0.443333,0.436667,0.430000,0.423333,0.413333,0.406667,0.403333,0.396667,0.393333,0.393333,0.390000,0.390000,0.393333,0.393333,0.393333,0.393333,0.393333,0.393333,0.393333,0.393333,0.393333,0.396667,0.396667,0.400000,0.403333,0.406667,0.406667,0.406667,0.406667,0.406667,0.403333,0.403333,0.400000,0.400000,0.396667,0.396667,0.400000,0.400000,0.403333,0.410000,0.413333,0.420000,0.423333,0.426667,0.430000,0.433333,0.436667,0.440000,0.440000,0.443333,0.446667,0.450000,0.453333,0.460000,0.463333,0.470000,0.476667,0.480000,0.483333,0.483333,0.486667,0.486667,0.486667,0.486667,0.486667,0.486667,0.486667,0.490000,0.493333,0.496667,0.500000,0.503333,0.503333,0.503333,0.503333,0.500000,0.500000,0.496667,0.493333,0.490000,0.490000,0.486667,0.486667,0.490000,0.493333,0.493333,0.496667,0.496667,0.500000,0.500000,0.500000,0.500000,0.496667,0.496667,0.493333,0.493333,0.493333,0.490000,0.490000,0.490000,0.486667,0.486667,0.483333,0.480000,0.473667,0.470333,0.463667,0.457000,0.450333,0.443667,0.437000,0.430333,0.427000,0.423667,0.420333,0.420333,0.420333,0.417000,0.417000,0.417000,0.413667,0.410333,0.410333,0.407000,0.400333,0.397000,0.393667,0.390333,0.390333,0.387000,0.387000,0.387000,0.387000,0.387000,0.383667,0.383667,0.380333,0.377000,0.373667,0.367000,0.363667,0.360333,0.357000,0.353667,0.350333,0.350333,0.347000,0.343667,0.343667,0.340333,0.337000,0.333333,0.329733,0.324400,0.318767,0.313100,0.307133,0.301167,0.295167,0.290200,0.285900,0.281600,0.278267,0.275600,0.272933,0.270567,0.268533,0.266500,0.264467,0.262800,0.260800,0.259133,0.257500,0.256200,0.254900,0.254267,0.253633,0.253333,0.253667,0.254000,0.254000,0.254000,0.254000,0.253667,0.252000,0.250633,0.248633,0.245967,0.243267,0.240600,0.237933,0.235233,0.232900,0.230900,0.228900,0.226900,0.224900,0.222267,0.219933,0.216300,0.212000,0.208000,0.203033,0.197400,0.192100,0.186800,0.182167,0.177533,0.173533,0.170900,0.168567,0.166267,0.164600,0.162933,0.161300,0.158633,0.155333,0.152367,0.148067,0.143133,0.137867,0.132567,0.126933,0.121333,0.116033,0.112033,0.107700,0.104033,0.102433,0.100800,0.0992000,0.0992667,0.0993667,0.0997667,0.100500,0.101533,0.102933,0.103967,0.105600,0.106933,0.108567,0.110167,0.111767,0.113400,0.114367,0.115667,0.116667,0.117000,0.117300,0.117333,0.117000,0.116033,0.114733,0.113467,0.111500,0.109567,0.107633,0.105400,0.103133,0.100900,0.0990000,0.0971000,0.0952333,0.0933667,0.0915333,0.0899667,0.0880667,0.0861667,0.0842333,0.0820333,0.0794333,0.0767667,0.0741667,0.0715000,0.0691000,0.0667000,0.0646000,0.0634667,0.0626000,0.0614667,0.0620000,0.0622667,0.0628333,0.0635000,0.0642000,0.0649000,0.0652667,0.0656000,0.0656333,0.0656667,0.0656667,0.0653667,0.0650333,0.0650667,0.0650667,0.0650667,0.0651000,0.0651667,0.0652333,0.0652667,0.0643667,0.0635000,0.0626000,0.0611333,0.0589667,0.0568000,0.0546667,0.0524667,0.0499667,0.0477333,0.0452000,0.0431667,0.0409000,0.0386667,0.0364000,0.0341667,0.0319000,0.0297400,0.0275833,0.0253667,0.0237300,0.0234867,0.0229767,0.0220867,0.0799000,0.0889333,0.0887333,0.0896667,0.0867000,0.0837333,0.0810667,0.0784000,0.0757333,0.0731333,0.0701667,0.0670667,0.0641000,0.0612667,0.0585333,0.0560000,0.0535333,0.0511667,0.0489667,0.0468333,0.0448000,0.0428667,0.0409333,0.0391000,0.0373667,0.0356667,0.0340000,0.0319900,0.0296733,0.0273733,0.0239900,0.0146633,0.0157133,0.0169200,0.0183433,0.0193867,0.0204467,0.0215200,0.0222533,0.0226367,0.0229467,0.0232767,0.0231633,0.0230633,0.0229467,0.0227167,0.0224200,0.0221600,0.0219167,0.0214500,0.0208700,0.0203533,0.0196467,0.0182533,0.0168700,0.0155700,0.0137000,0.0116500,0.00979333,0.00798000,0.00476000,0.00282900,0.00171433,0.00132567,0.00267667,0.00356333,0.00467000,0.00599667,0.00701667,0.00728333,0.00747667,0.00774333,0.00782333,0.00784000,0.00786333,0.00789667,0.00785000,0.00784000,0.00783667,0.00771333,0.00709667,0.00656667,0.00605333,0.00506000,0.00371000,0.00269700,0.00197233]}

# AKARI, already converted to energy. Based on flux ratio plot against
# flux (where factor 0.98 was already applied), add 18um calibration
# conversion
def s18_cal_func(x):
    """Calibration to non-constant part, by eye."""
    if x < 5.0:
        return (x/5)**0.02
    else:
        return 1.0


filters['AKARI9'] = {'svo_name': 'AKARI/IRC.S9W',
                     'response_type': 'energy',
                     'ref_wavelength': 9.0,
                     'ref_spectrum': lambda nu: 1.0/nu}
filters['AKARI18'] = {'svo_name': 'AKARI/IRC.L18W',
                      'response_type': 'energy',
                      'ref_wavelength': 18.0,
                      'measurement_calibration': lambda x: x*0.98*s18_cal_func(x),
                      'ref_spectrum': lambda nu: 1.0/nu}

# Spitzer, IRAC/IRS PUI is photon counting, MIPS energy
# IRAC/IRS ref spectrum F_nu oc 1/nu, MIPS 10k BB
filters['IRAC3P6'] = {'svo_name': 'Spitzer/IRAC.I1',
                      'response_type': 'photon',
                      'ref_wavelength': 3.550,
                      'ref_spectrum': lambda nu: 1.0/nu}
filters['IRAC4P5'] = {'svo_name': 'Spitzer/IRAC.I2',
                      'response_type': 'photon',
                      'ref_wavelength': 4.493,
                      'ref_spectrum': lambda nu: 1.0/nu}
filters['IRAC5P8'] = {'svo_name': 'Spitzer/IRAC.I3',
                      'response_type': 'photon',
                      'ref_wavelength': 5.731,
                      'ref_spectrum': lambda nu: 1.0/nu}
filters['IRAC8'] = {'svo_name': 'Spitzer/IRAC.I4',
                    'response_type': 'photon',
                    'ref_wavelength': 7.872,
                    'ref_spectrum': lambda nu: 1.0/nu}

filters['IRSPUB'] = {
    'magnitude_system':    'Vega',
    'zero_point':           0.00000,
    'zero_point_ref':      'bla',
    'ref_wavelength':       15.8,
    'ref_spectrum':         lambda nu: 1.0/nu,
    'response_type':       'photon',
    'response_ref':        'http://irsa.ipac.caltech.edu/data/SPITZER/docs/\
                            dataanalysistools/cookbook/14/',
    'wav_micron':[12.288, 12.317, 12.346, 12.376, 12.405, 12.435,
12.465, 12.495, 12.525, 12.555, 12.586, 12.616, 12.647, 12.678, 12.709,
12.741, 12.772, 12.803, 12.835, 12.867, 12.899, 12.931, 12.964, 12.996,
13.029, 13.061, 13.094, 13.128, 13.161, 13.194, 13.228, 13.262, 13.296,
13.330, 13.364, 13.399, 13.434, 13.469, 13.504, 13.539, 13.574, 13.610,
13.646, 13.682, 13.718, 13.754, 13.791, 13.828, 13.865, 13.902, 13.939,
13.977, 14.015, 14.052, 14.091, 14.129, 14.168, 14.207, 14.246, 14.285,
14.324, 14.364, 14.404, 14.444, 14.484, 14.525, 14.566, 14.607, 14.648,
14.689, 14.731, 14.773, 14.815, 14.858, 14.901, 14.943, 14.987, 15.030,
15.074, 15.118, 15.162, 15.206, 15.251, 15.296, 15.341, 15.387, 15.433,
15.479, 15.525, 15.572, 15.619, 15.666, 15.713, 15.761, 15.809, 15.857,
15.906, 15.955, 16.004, 16.054, 16.104, 16.154, 16.204, 16.255, 16.306,
16.358, 16.409, 16.462, 16.514, 16.567, 16.620, 16.673, 16.727, 16.781,
16.836, 16.890, 16.946, 17.001, 17.057, 17.114, 17.170, 17.227, 17.285,
17.342, 17.401, 17.459, 17.518, 17.578, 17.637, 17.698, 17.758, 17.819,
17.881, 17.942, 18.005, 18.067, 18.131, 18.194, 18.258, 18.323, 18.388,
18.453, 18.519, 18.586, 18.653, 18.720, 18.788, 18.856, 18.925, 18.994,
19.064, 19.134, 19.205, 19.277, 19.348, 19.421],
    'response':[0.001, 0.002, 0.003, 0.002, 0.002, 0.002, 0.002, 0.001,
0.002, 0.002, 0.002, 0.004, 0.006, 0.009, 0.011, 0.013, 0.018, 0.026,
0.038, 0.060, 0.077, 0.093, 0.107, 0.132, 0.171, 0.219, 0.290, 0.355,
0.432, 0.503, 0.578, 0.639, 0.676, 0.706, 0.723, 0.720, 0.709, 0.723,
0.736, 0.764, 0.791, 0.861, 0.909, 0.939, 0.946, 0.963, 0.936, 0.932,
0.936, 0.949, 0.953, 0.963, 0.963, 0.959, 0.953, 0.943, 0.929, 0.929,
0.929, 0.936, 0.936, 0.939, 0.949, 0.953, 0.939, 0.949, 0.959, 0.970,
0.976, 0.990, 1.000, 0.993, 0.993, 0.983, 0.966, 0.936, 0.916, 0.895,
0.882, 0.865, 0.845, 0.828, 0.807, 0.794, 0.784, 0.770, 0.764, 0.760,
0.767, 0.797, 0.828, 0.848, 0.872, 0.902, 0.902, 0.892, 0.882, 0.875,
0.858, 0.841, 0.831, 0.814, 0.797, 0.784, 0.774, 0.770, 0.777, 0.784,
0.787, 0.787, 0.791, 0.791, 0.791, 0.784, 0.774, 0.770, 0.767, 0.764,
0.770, 0.787, 0.794, 0.780, 0.794, 0.818, 0.821, 0.841, 0.845, 0.845,
0.841, 0.804, 0.770, 0.757, 0.750, 0.753, 0.747, 0.736, 0.740, 0.757,
0.747, 0.740, 0.760, 0.747, 0.703, 0.635, 0.530, 0.399, 0.271, 0.169,
0.095, 0.056, 0.033, 0.020, 0.013, 0.009, 0.005, 0.002]}

filters['IRSPUR'] = {
    'magnitude_system':    'Vega',
    'zero_point':           0.00000,
    'zero_point_ref':      'bla',
    'ref_wavelength':       22.3,
    'ref_spectrum':         lambda nu: 1.0/nu,
    'response_type':       'photon',
    'response_ref':        'http://irsa.ipac.caltech.edu/data/SPITZER/docs/\
                            dataanalysistools/cookbook/14/',
    'wav_micron':[17.170, 17.227, 17.285, 17.342, 17.401, 17.459,
17.518, 17.578, 17.637, 17.698, 17.758, 17.819, 17.881, 17.942, 18.005,
18.067, 18.131, 18.194, 18.258, 18.323, 18.388, 18.453, 18.519, 18.586,
18.653, 18.720, 18.788, 18.856, 18.925, 18.994, 19.064, 19.134, 19.205,
19.277, 19.348, 19.421, 19.494, 19.567, 19.642, 19.716, 19.792, 19.867,
19.944, 20.021, 20.098, 20.177, 20.255, 20.335, 20.415, 20.496, 20.577,
20.659, 20.742, 20.825, 20.909, 20.993, 21.079, 21.165, 21.252, 21.339,
21.427, 21.516, 21.606, 21.696, 21.787, 21.879, 21.972, 22.066, 22.160,
22.255, 22.351, 22.448, 22.545, 22.644, 22.743, 22.843, 22.944, 23.046,
23.149, 23.253, 23.358, 23.463, 23.570, 23.677, 23.786, 23.896, 24.006,
24.118, 24.231, 24.344, 24.459, 24.575, 24.692, 24.810, 24.930, 25.050,
25.172, 25.295, 25.419, 25.544, 25.670, 25.798, 25.927, 26.057, 26.189,
26.322, 26.456, 26.592, 26.729, 26.867, 27.007, 27.149, 27.292, 27.436,
27.582, 27.729, 27.878, 28.029, 28.181, 28.335, 28.491, 28.649, 28.808,
28.969, 29.131, 29.296, 29.462, 29.631, 29.801, 29.973, 30.147, 30.324,
30.502, 30.683, 30.865, 31.050, 31.237, 31.427, 31.618, 31.812, 32.008,
32.207, 32.409, 32.612, 32.819, 33.028, 33.240, 33.454, 33.671, 33.891,
34.114, 34.340, 34.569, 34.801, 35.036],
    'response':[ 0.001, 0.000, 0.002, 0.003, 0.004, 0.007, 0.003, 0.004,
0.006, 0.005, 0.001, 0.006, 0.003, 0.005, 0.008, 0.015, 0.025, 0.040,
0.059, 0.094, 0.149, 0.238, 0.351, 0.500, 0.651, 0.780, 0.862, 0.931,
0.954, 0.986, 1.000, 0.959, 0.913, 0.904, 0.940, 0.950, 0.945, 0.959,
1.000, 0.991, 0.995, 0.972, 0.986, 0.954, 0.945, 0.931, 0.922, 0.913,
0.908, 0.899, 0.885, 0.872, 0.858, 0.849, 0.849, 0.853, 0.862, 0.881,
0.894, 0.904, 0.913, 0.922, 0.922, 0.922, 0.917, 0.913, 0.904, 0.899,
0.890, 0.881, 0.862, 0.844, 0.830, 0.821, 0.812, 0.812, 0.817, 0.817,
0.803, 0.789, 0.771, 0.748, 0.729, 0.725, 0.739, 0.761, 0.784, 0.803,
0.817, 0.817, 0.798, 0.789, 0.766, 0.748, 0.725, 0.706, 0.679, 0.651,
0.619, 0.592, 0.555, 0.523, 0.491, 0.459, 0.430, 0.398, 0.378, 0.356,
0.336, 0.318, 0.304, 0.284, 0.266, 0.247, 0.226, 0.204, 0.185, 0.167,
0.154, 0.137, 0.124, 0.109, 0.098, 0.084, 0.077, 0.068, 0.061, 0.053,
0.046, 0.039, 0.032, 0.026, 0.020, 0.013, 0.008, 0.005, 0.004, 0.003,
0.002, 0.002, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000,-0.000,-0.000]}

filters['MIPS24'] = {'svo_name': 'Spitzer/MIPS.24mu',
                     'response_type': 'energy',
                     'ref_wavelength': 23.675,
                     'ref_spectrum': lambda nu: utils.bnu_nu_hz(nu,10000.0)}
filters['MIPS70'] = {'svo_name': 'Spitzer/MIPS.70mu',
                     'response_type': 'energy',
                     'ref_wavelength': 71.42,
                     'ref_spectrum': lambda nu: utils.bnu_nu_hz(nu,10000.0)}
filters['MIPS160'] = {'svo_name': 'Spitzer/MIPS.160mu',
                      'response_type': 'energy',
                      'ref_wavelength': 155.9,
                      'ref_spectrum': lambda nu: utils.bnu_nu_hz(nu,10000.0)}

# JWST NIRCAM, units of electrons/photon
nrc_filt_loc = os.path.dirname(os.path.abspath(__file__))+ \
                                    '/data/filters/nircam/'
for file in glob.glob(nrc_filt_loc+'*.txt'):
    filt_name = 'NIRCAM.'+os.path.basename(file).split('_')[0]
    filters[filt_name] = {
        'magnitude_system': 'Vega',
        'response_type':'photon',
        'ref_wavelength': None,
        'ref_spectrum': None,
        'response_ref': 'https://jwst-docs.stsci.edu/display/JTI/NIRCam+Filters',
        'wav_micron': np.loadtxt(file,skiprows=1,usecols=0),
        'response': np.loadtxt(file,skiprows=1,usecols=1)
                          }

# JWST MIRI, units of electrons/photon
miri_file = os.path.dirname(os.path.abspath(__file__))+ \
                        '/data/filters/ImPCE_TN-00072-ATC-Iss2.csv'
for i,filt in enumerate(['F560W','F770W','F1000W','F1280W','F1130W',
                         'F1500W','F1800W','F2100W','F2550W']):
    filt_name = 'MIRI.'+filt
    filters[filt_name] = {
        'magnitude_system': 'Vega',
        'response_type':'photon',
        'ref_wavelength': None,
        'ref_spectrum': None,
        'response_ref': 'https://jwst-docs.stsci.edu/display/JTI/MIRI+Filters+and+Dispersers',
        'wav_micron': np.loadtxt(miri_file,delimiter=',',skiprows=2,usecols=0),
        'response': np.loadtxt(miri_file,delimiter=',',skiprows=2,usecols=i+1)
                          }

# IRAS, RSRs, calibrations empirical, and see Rieke+2008
filters['IRAS12'] = {'svo_name': 'IRAS/IRAS.12mu',
                     'response_type': 'energy',
                     'ref_wavelength': 12.0,
                     'measurement_calibration': lambda x: 0.976*x,
                     'ref_spectrum': lambda nu: 1.0/nu}
filters['IRAS25'] = {'svo_name': 'IRAS/IRAS.25mu',
                     'response_type': 'energy',
                     'ref_wavelength': 25.0,
                     'measurement_calibration': lambda x: 0.94*x,
                     'ref_spectrum': lambda nu: 1.0/nu}
filters['IRAS60'] = {'svo_name': 'IRAS/IRAS.60mu',
                     'response_type': 'energy',
                     'ref_wavelength': 60.0,
                     'ref_spectrum': lambda nu: 1.0/nu}
filters['IRAS100'] = {'svo_name': 'IRAS/IRAS.100mu',
                      'response_type': 'energy',
                      'ref_wavelength': 100.0,
                      'ref_spectrum': lambda nu: 1.0/nu}

# MSX, RSRs, ref spectrum is F_nu oc 1/nu^2
filters['MSX8'] = {'svo_name': 'MSX/MSX.A',
                   'response_type': 'energy',
                   'ref_wavelength': 8.28,
                   'ref_spectrum': lambda nu: 1.0/nu}
filters['MSX12'] = {'svo_name': 'MSX/MSX.C',
                    'response_type': 'energy',
                    'ref_wavelength': 12.13,
                    'ref_spectrum': lambda nu: 1.0/nu}
filters['MSX15'] = {'svo_name': 'MSX/MSX.D',
                    'response_type': 'energy',
                    'ref_wavelength': 14.65,
                    'ref_spectrum': lambda nu: 1.0/nu}
filters['MSX21'] = {'svo_name': 'MSX/MSX.E',
                    'response_type': 'energy',
                    'ref_wavelength': 21.34,
                    'ref_spectrum': lambda nu: 1.0/nu}

# PACS/SPIRE, RSRs, ref spectrum F_nu oc 1/nu
filters['PACS70'] = {'svo_name': 'Herschel/Pacs.blue',
                     'response_type': 'energy',
                     'ref_wavelength': 70.0,
                     'ref_spectrum': lambda nu: 1.0/nu}
filters['PACS100'] = {'svo_name': 'Herschel/Pacs.green',
                     'response_type': 'energy',
                     'ref_wavelength': 100.0,
                     'ref_spectrum': lambda nu: 1.0/nu}
filters['PACS160'] = {'svo_name': 'Herschel/Pacs.red',
                     'response_type': 'energy',
                     'ref_wavelength': 160.0,
                     'ref_spectrum': lambda nu: 1.0/nu}
filters['SPIRE250'] = {'svo_name': 'Herschel/SPIRE.PSW',
                     'response_type': 'energy',
                     'ref_wavelength': 250.0,
                     'ref_spectrum': lambda nu: 1.0/nu}
filters['SPIRE350'] = {'svo_name': 'Herschel/SPIRE.PMW',
                     'response_type': 'energy',
                     'ref_wavelength': 350.0,
                     'ref_spectrum': lambda nu: 1.0/nu}
filters['SPIRE500'] = {'svo_name': 'Herschel/SPIRE.PLW',
                     'response_type': 'energy',
                     'ref_wavelength': 500.0,
                     'ref_spectrum': lambda nu: 1.0/nu}

# LBTI NOMIC N, assume QE based
filters['NOMICN'] = {
    'magnitude_system':    'Vega',
    'zero_point':           31.3207,
    'zero_point_ref':      'bla',
    'ref_wavelength':       None,
    'ref_spectrum':         None,
    'response_type':       'photon',
    'response_ref':        'bla',
    'wav_micron':[
9.000,  9.001,  9.002,  9.003,  9.004,  9.005,  9.006,  9.007,  9.008,
9.009,  9.010,  9.011,  9.012,  9.013,  9.014,  9.015,  9.016,  9.017,
9.018,  9.019,  9.020,  9.021,  9.022,  9.023,  9.024,  9.025,  9.026,
9.027,  9.028,  9.029,  9.030,  9.031,  9.032,  9.033,  9.034,  9.035,
9.036,  9.037,  9.038,  9.039,  9.040,  9.041,  9.042,  9.043,  9.044,
9.045,  9.046,  9.047,  9.048,  9.049,  9.050,  9.051,  9.052,  9.053,
9.054,  9.055,  9.056,  9.057,  9.058,  9.059,  9.060,  9.061,  9.062,
9.063,  9.064,  9.065,  9.066,  9.067,  9.068,  9.069,  9.070,  9.071,
9.072,  9.073,  9.074,  9.075,  9.076,  9.077,  9.078,  9.079,  9.080,
9.081,  9.082,  9.083,  9.084,  9.085,  9.086,  9.087,  9.088,  9.089,
9.090,  9.091,  9.092,  9.093,  9.094,  9.095,  9.096,  9.097,  9.098,
9.099,  9.100,  9.101,  9.102,  9.103,  9.104,  9.105,  9.106,  9.107,
9.108,  9.109,  9.110,  9.111,  9.112,  9.113,  9.114,  9.115,  9.116,
9.117,  9.118,  9.119,  9.120,  9.121,  9.122,  9.123,  9.124,  9.125,
9.126,  9.127,  9.128,  9.129,  9.130,  9.131,  9.132,  9.133,  9.134,
9.135,  9.136,  9.137,  9.138,  9.139,  9.140,  9.141,  9.142,  9.143,
9.144,  9.145,  9.146,  9.147,  9.148,  9.149,  9.150,  9.151,  9.152,
9.153,  9.154,  9.155,  9.156,  9.157,  9.158,  9.159,  9.160,  9.161,
9.162,  9.163,  9.164,  9.165,  9.166,  9.167,  9.168,  9.169,  9.170,
9.171,  9.172,  9.173,  9.174,  9.175,  9.176,  9.177,  9.178,  9.179,
9.180,  9.181,  9.182,  9.183,  9.184,  9.185,  9.186,  9.187,  9.188,
9.189,  9.190,  9.191,  9.192,  9.193,  9.194,  9.195,  9.196,  9.197,
9.198,  9.199,  9.200,  9.201,  9.202,  9.203,  9.204,  9.205,  9.206,
9.207,  9.208,  9.209,  9.210,  9.211,  9.212,  9.213,  9.214,  9.215,
9.216,  9.217,  9.218,  9.219,  9.220,  9.221,  9.222,  9.223,  9.224,
9.225,  9.226,  9.227,  9.228,  9.229,  9.230,  9.231,  9.232,  9.233,
9.234,  9.235,  9.236,  9.237,  9.238,  9.239,  9.240,  9.241,  9.242,
9.243,  9.244,  9.245,  9.246,  9.247,  9.248,  9.249,  9.250,  9.251,
9.252,  9.253,  9.254,  9.255,  9.256,  9.257,  9.258,  9.259,  9.260,
9.261,  9.262,  9.263,  9.264,  9.265,  9.266,  9.267,  9.268,  9.269,
9.270,  9.271,  9.272,  9.273,  9.274,  9.275,  9.276,  9.277,  9.278,
9.279,  9.280,  9.281,  9.282,  9.283,  9.284,  9.285,  9.286,  9.287,
9.288,  9.289,  9.290,  9.291,  9.292,  9.293,  9.294,  9.295,  9.296,
9.297,  9.298,  9.299,  9.300,  9.301,  9.302,  9.303,  9.304,  9.305,
9.306,  9.307,  9.308,  9.309,  9.310,  9.311,  9.312,  9.313,  9.314,
9.315,  9.316,  9.317,  9.318,  9.319,  9.320,  9.321,  9.322,  9.323,
9.324,  9.325,  9.326,  9.327,  9.328,  9.329,  9.330,  9.331,  9.332,
9.333,  9.334,  9.335,  9.336,  9.337,  9.338,  9.339,  9.340,  9.341,
9.342,  9.343,  9.344,  9.345,  9.346,  9.347,  9.348,  9.349,  9.350,
9.351,  9.352,  9.353,  9.354,  9.355,  9.356,  9.357,  9.358,  9.359,
9.360,  9.361,  9.362,  9.363,  9.364,  9.365,  9.366,  9.367,  9.368,
9.369,  9.370,  9.371,  9.372,  9.373,  9.374,  9.375,  9.376,  9.377,
9.378,  9.379,  9.380,  9.381,  9.382,  9.383,  9.384,  9.385,  9.386,
9.387,  9.388,  9.389,  9.390,  9.391,  9.392,  9.393,  9.394,  9.395,
9.396,  9.397,  9.398,  9.399,  9.400,  9.401,  9.402,  9.403,  9.404,
9.405,  9.406,  9.407,  9.408,  9.409,  9.410,  9.411,  9.412,  9.413,
9.414,  9.415,  9.416,  9.417,  9.418,  9.419,  9.420,  9.421,  9.422,
9.423,  9.424,  9.425,  9.426,  9.427,  9.428,  9.429,  9.430,  9.431,
9.432,  9.433,  9.434,  9.435,  9.436,  9.437,  9.438,  9.439,  9.440,
9.441,  9.442,  9.443,  9.444,  9.445,  9.446,  9.447,  9.448,  9.449,
9.450,  9.451,  9.452,  9.453,  9.454,  9.455,  9.456,  9.457,  9.458,
9.459,  9.460,  9.461,  9.462,  9.463,  9.464,  9.465,  9.466,  9.467,
9.468,  9.469,  9.470,  9.471,  9.472,  9.473,  9.474,  9.475,  9.476,
9.477,  9.478,  9.479,  9.480,  9.481,  9.482,  9.483,  9.484,  9.485,
9.486,  9.487,  9.488,  9.489,  9.490,  9.491,  9.492,  9.493,  9.494,
9.495,  9.496,  9.497,  9.498,  9.499,  9.500,  9.501,  9.502,  9.503,
9.504,  9.505,  9.506,  9.507,  9.508,  9.509,  9.510,  9.511,  9.512,
9.513,  9.514,  9.515,  9.516,  9.517,  9.518,  9.519,  9.520,  9.521,
9.522,  9.523,  9.524,  9.525,  9.526,  9.527,  9.528,  9.529,  9.530,
9.531,  9.532,  9.533,  9.534,  9.535,  9.536,  9.537,  9.538,  9.539,
9.540,  9.541,  9.542,  9.543,  9.544,  9.545,  9.546,  9.547,  9.548,
9.549,  9.550,  9.551,  9.552,  9.553,  9.554,  9.555,  9.556,  9.557,
9.558,  9.559,  9.560,  9.561,  9.562,  9.563,  9.564,  9.565,  9.566,
9.567,  9.568,  9.569,  9.570,  9.571,  9.572,  9.573,  9.574,  9.575,
9.576,  9.577,  9.578,  9.579,  9.580,  9.581,  9.582,  9.583,  9.584,
9.585,  9.586,  9.587,  9.588,  9.589,  9.590,  9.591,  9.592,  9.593,
9.594,  9.595,  9.596,  9.597,  9.598,  9.599,  9.600,  9.601,  9.602,
9.603,  9.604,  9.605,  9.606,  9.607,  9.608,  9.609,  9.610,  9.611,
9.612,  9.613,  9.614,  9.615,  9.616,  9.617,  9.618,  9.619,  9.620,
9.621,  9.622,  9.623,  9.624,  9.625,  9.626,  9.627,  9.628,  9.629,
9.630,  9.631,  9.632,  9.633,  9.634,  9.635,  9.636,  9.637,  9.638,
9.639,  9.640,  9.641,  9.642,  9.643,  9.644,  9.645,  9.646,  9.647,
9.648,  9.649,  9.650,  9.651,  9.652,  9.653,  9.654,  9.655,  9.656,
9.657,  9.658,  9.659,  9.660,  9.661,  9.662,  9.663,  9.664,  9.665,
9.666,  9.667,  9.668,  9.669,  9.670,  9.671,  9.672,  9.673,  9.674,
9.675,  9.676,  9.677,  9.678,  9.679,  9.680,  9.681,  9.682,  9.683,
9.684,  9.685,  9.686,  9.687,  9.688,  9.689,  9.690,  9.691,  9.692,
9.693,  9.694,  9.695,  9.696,  9.697,  9.698,  9.699,  9.700,  9.701,
9.702,  9.703,  9.704,  9.705,  9.706,  9.707,  9.708,  9.709,  9.710,
9.711,  9.712,  9.713,  9.714,  9.715,  9.716,  9.717,  9.718,  9.719,
9.720,  9.721,  9.722,  9.723,  9.724,  9.725,  9.726,  9.727,  9.728,
9.729,  9.730,  9.731,  9.732,  9.733,  9.734,  9.735,  9.736,  9.737,
9.738,  9.739,  9.740,  9.741,  9.742,  9.743,  9.744,  9.745,  9.746,
9.747,  9.748,  9.749,  9.750,  9.751,  9.752,  9.753,  9.754,  9.755,
9.756,  9.757,  9.758,  9.759,  9.760,  9.761,  9.762,  9.763,  9.764,
9.765,  9.766,  9.767,  9.768,  9.769,  9.770,  9.771,  9.772,  9.773,
9.774,  9.775,  9.776,  9.777,  9.778,  9.779,  9.780,  9.781,  9.782,
9.783,  9.784,  9.785,  9.786,  9.787,  9.788,  9.789,  9.790,  9.791,
9.792,  9.793,  9.794,  9.795,  9.796,  9.797,  9.798,  9.799,  9.800,
9.801,  9.802,  9.803,  9.804,  9.805,  9.806,  9.807,  9.808,  9.809,
9.810,  9.811,  9.812,  9.813,  9.814,  9.815,  9.816,  9.817,  9.818,
9.819,  9.820,  9.821,  9.822,  9.823,  9.824,  9.825,  9.826,  9.827,
9.828,  9.829,  9.830,  9.831,  9.832,  9.833,  9.834,  9.835,  9.836,
9.837,  9.838,  9.839,  9.840,  9.841,  9.842,  9.843,  9.844,  9.845,
9.846,  9.847,  9.848,  9.849,  9.850,  9.851,  9.852,  9.853,  9.854,
9.855,  9.856,  9.857,  9.858,  9.859,  9.860,  9.861,  9.862,  9.863,
9.864,  9.865,  9.866,  9.867,  9.868,  9.869,  9.870,  9.871,  9.872,
9.873,  9.874,  9.875,  9.876,  9.877,  9.878,  9.879,  9.880,  9.881,
9.882,  9.883,  9.884,  9.885,  9.886,  9.887,  9.888,  9.889,  9.890,
9.891,  9.892,  9.893,  9.894,  9.895,  9.896,  9.897,  9.898,  9.899,
9.900,  9.901,  9.902,  9.903,  9.904,  9.905,  9.906,  9.907,  9.908,
9.909,  9.910,  9.911,  9.912,  9.913,  9.914,  9.915,  9.916,  9.917,
9.918,  9.919,  9.920,  9.921,  9.922,  9.923,  9.924,  9.925,  9.926,
9.927,  9.928,  9.929,  9.930,  9.931,  9.932,  9.933,  9.934,  9.935,
9.936,  9.937,  9.938,  9.939,  9.940,  9.941,  9.942,  9.943,  9.944,
9.945,  9.946,  9.947,  9.948,  9.949,  9.950,  9.951,  9.952,  9.953,
9.954,  9.955,  9.956,  9.957,  9.958,  9.959,  9.960,  9.961,  9.962,
9.963,  9.964,  9.965,  9.966,  9.967,  9.968,  9.969,  9.970,  9.971,
9.972,  9.973,  9.974,  9.975,  9.976,  9.977,  9.978,  9.979,  9.980,
9.981,  9.982,  9.983,  9.984,  9.985,  9.986,  9.987,  9.988,  9.989,
9.990,  9.991,  9.992,  9.993,  9.994,  9.995,  9.996,  9.997,  9.998,
9.999, 10.000, 10.001, 10.002, 10.003, 10.004, 10.005, 10.006, 10.007,
10.008, 10.009, 10.010, 10.011, 10.012, 10.013, 10.014, 10.015, 10.016,
10.017, 10.018, 10.019, 10.020, 10.021, 10.022, 10.023, 10.024, 10.025,
10.026, 10.027, 10.028, 10.029, 10.030, 10.031, 10.032, 10.033, 10.034,
10.035, 10.036, 10.037, 10.038, 10.039, 10.040, 10.041, 10.042, 10.043,
10.044, 10.045, 10.046, 10.047, 10.048, 10.049, 10.050, 10.051, 10.052,
10.053, 10.054, 10.055, 10.056, 10.057, 10.058, 10.059, 10.060, 10.061,
10.062, 10.063, 10.064, 10.065, 10.066, 10.067, 10.068, 10.069, 10.070,
10.071, 10.072, 10.073, 10.074, 10.075, 10.076, 10.077, 10.078, 10.079,
10.080, 10.081, 10.082, 10.083, 10.084, 10.085, 10.086, 10.087, 10.088,
10.089, 10.090, 10.091, 10.092, 10.093, 10.094, 10.095, 10.096, 10.097,
10.098, 10.099, 10.100, 10.101, 10.102, 10.103, 10.104, 10.105, 10.106,
10.107, 10.108, 10.109, 10.110, 10.111, 10.112, 10.113, 10.114, 10.115,
10.116, 10.117, 10.118, 10.119, 10.120, 10.121, 10.122, 10.123, 10.124,
10.125, 10.126, 10.127, 10.128, 10.129, 10.130, 10.131, 10.132, 10.133,
10.134, 10.135, 10.136, 10.137, 10.138, 10.139, 10.140, 10.141, 10.142,
10.143, 10.144, 10.145, 10.146, 10.147, 10.148, 10.149, 10.150, 10.151,
10.152, 10.153, 10.154, 10.155, 10.156, 10.157, 10.158, 10.159, 10.160,
10.161, 10.162, 10.163, 10.164, 10.165, 10.166, 10.167, 10.168, 10.169,
10.170, 10.171, 10.172, 10.173, 10.174, 10.175, 10.176, 10.177, 10.178,
10.179, 10.180, 10.181, 10.182, 10.183, 10.184, 10.185, 10.186, 10.187,
10.188, 10.189, 10.190, 10.191, 10.192, 10.193, 10.194, 10.195, 10.196,
10.197, 10.198, 10.199, 10.200, 10.201, 10.202, 10.203, 10.204, 10.205,
10.206, 10.207, 10.208, 10.209, 10.210, 10.211, 10.212, 10.213, 10.214,
10.215, 10.216, 10.217, 10.218, 10.219, 10.220, 10.221, 10.222, 10.223,
10.224, 10.225, 10.226, 10.227, 10.228, 10.229, 10.230, 10.231, 10.232,
10.233, 10.234, 10.235, 10.236, 10.237, 10.238, 10.239, 10.240, 10.241,
10.242, 10.243, 10.244, 10.245, 10.246, 10.247, 10.248, 10.249, 10.250,
10.251, 10.252, 10.253, 10.254, 10.255, 10.256, 10.257, 10.258, 10.259,
10.260, 10.261, 10.262, 10.263, 10.264, 10.265, 10.266, 10.267, 10.268,
10.269, 10.270, 10.271, 10.272, 10.273, 10.274, 10.275, 10.276, 10.277,
10.278, 10.279, 10.280, 10.281, 10.282, 10.283, 10.284, 10.285, 10.286,
10.287, 10.288, 10.289, 10.290, 10.291, 10.292, 10.293, 10.294, 10.295,
10.296, 10.297, 10.298, 10.299, 10.300, 10.301, 10.302, 10.303, 10.304,
10.305, 10.306, 10.307, 10.308, 10.309, 10.310, 10.311, 10.312, 10.313,
10.314, 10.315, 10.316, 10.317, 10.318, 10.319, 10.320, 10.321, 10.322,
10.323, 10.324, 10.325, 10.326, 10.327, 10.328, 10.329, 10.330, 10.331,
10.332, 10.333, 10.334, 10.335, 10.336, 10.337, 10.338, 10.339, 10.340,
10.341, 10.342, 10.343, 10.344, 10.345, 10.346, 10.347, 10.348, 10.349,
10.350, 10.351, 10.352, 10.353, 10.354, 10.355, 10.356, 10.357, 10.358,
10.359, 10.360, 10.361, 10.362, 10.363, 10.364, 10.365, 10.366, 10.367,
10.368, 10.369, 10.370, 10.371, 10.372, 10.373, 10.374, 10.375, 10.376,
10.377, 10.378, 10.379, 10.380, 10.381, 10.382, 10.383, 10.384, 10.385,
10.386, 10.387, 10.388, 10.389, 10.390, 10.391, 10.392, 10.393, 10.394,
10.395, 10.396, 10.397, 10.398, 10.399, 10.400, 10.401, 10.402, 10.403,
10.404, 10.405, 10.406, 10.407, 10.408, 10.409, 10.410, 10.411, 10.412,
10.413, 10.414, 10.415, 10.416, 10.417, 10.418, 10.419, 10.420, 10.421,
10.422, 10.423, 10.424, 10.425, 10.426, 10.427, 10.428, 10.429, 10.430,
10.431, 10.432, 10.433, 10.434, 10.435, 10.436, 10.437, 10.438, 10.439,
10.440, 10.441, 10.442, 10.443, 10.444, 10.445, 10.446, 10.447, 10.448,
10.449, 10.450, 10.451, 10.452, 10.453, 10.454, 10.455, 10.456, 10.457,
10.458, 10.459, 10.460, 10.461, 10.462, 10.463, 10.464, 10.465, 10.466,
10.467, 10.468, 10.469, 10.470, 10.471, 10.472, 10.473, 10.474, 10.475,
10.476, 10.477, 10.478, 10.479, 10.480, 10.481, 10.482, 10.483, 10.484,
10.485, 10.486, 10.487, 10.488, 10.489, 10.490, 10.491, 10.492, 10.493,
10.494, 10.495, 10.496, 10.497, 10.498, 10.499, 10.500, 10.501, 10.502,
10.503, 10.504, 10.505, 10.506, 10.507, 10.508, 10.509, 10.510, 10.511,
10.512, 10.513, 10.514, 10.515, 10.516, 10.517, 10.518, 10.519, 10.520,
10.521, 10.522, 10.523, 10.524, 10.525, 10.526, 10.527, 10.528, 10.529,
10.530, 10.531, 10.532, 10.533, 10.534, 10.535, 10.536, 10.537, 10.538,
10.539, 10.540, 10.541, 10.542, 10.543, 10.544, 10.545, 10.546, 10.547,
10.548, 10.549, 10.550, 10.551, 10.552, 10.553, 10.554, 10.555, 10.556,
10.557, 10.558, 10.559, 10.560, 10.561, 10.562, 10.563, 10.564, 10.565,
10.566, 10.567, 10.568, 10.569, 10.570, 10.571, 10.572, 10.573, 10.574,
10.575, 10.576, 10.577, 10.578, 10.579, 10.580, 10.581, 10.582, 10.583,
10.584, 10.585, 10.586, 10.587, 10.588, 10.589, 10.590, 10.591, 10.592,
10.593, 10.594, 10.595, 10.596, 10.597, 10.598, 10.599, 10.600, 10.601,
10.602, 10.603, 10.604, 10.605, 10.606, 10.607, 10.608, 10.609, 10.610,
10.611, 10.612, 10.613, 10.614, 10.615, 10.616, 10.617, 10.618, 10.619,
10.620, 10.621, 10.622, 10.623, 10.624, 10.625, 10.626, 10.627, 10.628,
10.629, 10.630, 10.631, 10.632, 10.633, 10.634, 10.635, 10.636, 10.637,
10.638, 10.639, 10.640, 10.641, 10.642, 10.643, 10.644, 10.645, 10.646,
10.647, 10.648, 10.649, 10.650, 10.651, 10.652, 10.653, 10.654, 10.655,
10.656, 10.657, 10.658, 10.659, 10.660, 10.661, 10.662, 10.663, 10.664,
10.665, 10.666, 10.667, 10.668, 10.669, 10.670, 10.671, 10.672, 10.673,
10.674, 10.675, 10.676, 10.677, 10.678, 10.679, 10.680, 10.681, 10.682,
10.683, 10.684, 10.685, 10.686, 10.687, 10.688, 10.689, 10.690, 10.691,
10.692, 10.693, 10.694, 10.695, 10.696, 10.697, 10.698, 10.699, 10.700,
10.701, 10.702, 10.703, 10.704, 10.705, 10.706, 10.707, 10.708, 10.709,
10.710, 10.711, 10.712, 10.713, 10.714, 10.715, 10.716, 10.717, 10.718,
10.719, 10.720, 10.721, 10.722, 10.723, 10.724, 10.725, 10.726, 10.727,
10.728, 10.729, 10.730, 10.731, 10.732, 10.733, 10.734, 10.735, 10.736,
10.737, 10.738, 10.739, 10.740, 10.741, 10.742, 10.743, 10.744, 10.745,
10.746, 10.747, 10.748, 10.749, 10.750, 10.751, 10.752, 10.753, 10.754,
10.755, 10.756, 10.757, 10.758, 10.759, 10.760, 10.761, 10.762, 10.763,
10.764, 10.765, 10.766, 10.767, 10.768, 10.769, 10.770, 10.771, 10.772,
10.773, 10.774, 10.775, 10.776, 10.777, 10.778, 10.779, 10.780, 10.781,
10.782, 10.783, 10.784, 10.785, 10.786, 10.787, 10.788, 10.789, 10.790,
10.791, 10.792, 10.793, 10.794, 10.795, 10.796, 10.797, 10.798, 10.799,
10.800, 10.801, 10.802, 10.803, 10.804, 10.805, 10.806, 10.807, 10.808,
10.809, 10.810, 10.811, 10.812, 10.813, 10.814, 10.815, 10.816, 10.817,
10.818, 10.819, 10.820, 10.821, 10.822, 10.823, 10.824, 10.825, 10.826,
10.827, 10.828, 10.829, 10.830, 10.831, 10.832, 10.833, 10.834, 10.835,
10.836, 10.837, 10.838, 10.839, 10.840, 10.841, 10.842, 10.843, 10.844,
10.845, 10.846, 10.847, 10.848, 10.849, 10.850, 10.851, 10.852, 10.853,
10.854, 10.855, 10.856, 10.857, 10.858, 10.859, 10.860, 10.861, 10.862,
10.863, 10.864, 10.865, 10.866, 10.867, 10.868, 10.869, 10.870, 10.871,
10.872, 10.873, 10.874, 10.875, 10.876, 10.877, 10.878, 10.879, 10.880,
10.881, 10.882, 10.883, 10.884, 10.885, 10.886, 10.887, 10.888, 10.889,
10.890, 10.891, 10.892, 10.893, 10.894, 10.895, 10.896, 10.897, 10.898,
10.899, 10.900, 10.901, 10.902, 10.903, 10.904, 10.905, 10.906, 10.907,
10.908, 10.909, 10.910, 10.911, 10.912, 10.913, 10.914, 10.915, 10.916,
10.917, 10.918, 10.919, 10.920, 10.921, 10.922, 10.923, 10.924, 10.925,
10.926, 10.927, 10.928, 10.929, 10.930, 10.931, 10.932, 10.933, 10.934,
10.935, 10.936, 10.937, 10.938, 10.939, 10.940, 10.941, 10.942, 10.943,
10.944, 10.945, 10.946, 10.947, 10.948, 10.949, 10.950, 10.951, 10.952,
10.953, 10.954, 10.955, 10.956, 10.957, 10.958, 10.959, 10.960, 10.961,
10.962, 10.963, 10.964, 10.965, 10.966, 10.967, 10.968, 10.969, 10.970,
10.971, 10.972, 10.973, 10.974, 10.975, 10.976, 10.977, 10.978, 10.979,
10.980, 10.981, 10.982, 10.983, 10.984, 10.985, 10.986, 10.987, 10.988,
10.989, 10.990, 10.991, 10.992, 10.993, 10.994, 10.995, 10.996, 10.997,
10.998, 10.999, 11.000, 11.001, 11.002, 11.003, 11.004, 11.005, 11.006,
11.007, 11.008, 11.009, 11.010, 11.011, 11.012, 11.013, 11.014, 11.015,
11.016, 11.017, 11.018, 11.019, 11.020, 11.021, 11.022, 11.023, 11.024,
11.025, 11.026, 11.027, 11.028, 11.029, 11.030, 11.031, 11.032, 11.033,
11.034, 11.035, 11.036, 11.037, 11.038, 11.039, 11.040, 11.041, 11.042,
11.043, 11.044, 11.045, 11.046, 11.047, 11.048, 11.049, 11.050, 11.051,
11.052, 11.053, 11.054, 11.055, 11.056, 11.057, 11.058, 11.059, 11.060,
11.061, 11.062, 11.063, 11.064, 11.065, 11.066, 11.067, 11.068, 11.069,
11.070, 11.071, 11.072, 11.073, 11.074, 11.075, 11.076, 11.077, 11.078,
11.079, 11.080, 11.081, 11.082, 11.083, 11.084, 11.085, 11.086, 11.087,
11.088, 11.089, 11.090, 11.091, 11.092, 11.093, 11.094, 11.095, 11.096,
11.097, 11.098, 11.099, 11.100, 11.101, 11.102, 11.103, 11.104, 11.105,
11.106, 11.107, 11.108, 11.109, 11.110, 11.111, 11.112, 11.113, 11.114,
11.115, 11.116, 11.117, 11.118, 11.119, 11.120, 11.121, 11.122, 11.123,
11.124, 11.125, 11.126, 11.127, 11.128, 11.129, 11.130, 11.131, 11.132,
11.133, 11.134, 11.135, 11.136, 11.137, 11.138, 11.139, 11.140, 11.141,
11.142, 11.143, 11.144, 11.145, 11.146, 11.147, 11.148, 11.149, 11.150,
11.151, 11.152, 11.153, 11.154, 11.155, 11.156, 11.157, 11.158, 11.159,
11.160, 11.161, 11.162, 11.163, 11.164, 11.165, 11.166, 11.167, 11.168,
11.169, 11.170, 11.171, 11.172, 11.173, 11.174, 11.175, 11.176, 11.177,
11.178, 11.179, 11.180, 11.181, 11.182, 11.183, 11.184, 11.185, 11.186,
11.187, 11.188, 11.189, 11.190, 11.191, 11.192, 11.193, 11.194, 11.195,
11.196, 11.197, 11.198, 11.199, 11.200, 11.201, 11.202, 11.203, 11.204,
11.205, 11.206, 11.207, 11.208, 11.209, 11.210, 11.211, 11.212, 11.213,
11.214, 11.215, 11.216, 11.217, 11.218, 11.219, 11.220, 11.221, 11.222,
11.223, 11.224, 11.225, 11.226, 11.227, 11.228, 11.229, 11.230, 11.231,
11.232, 11.233, 11.234, 11.235, 11.236, 11.237, 11.238, 11.239, 11.240,
11.241, 11.242, 11.243, 11.244, 11.245, 11.246, 11.247, 11.248, 11.249,
11.250, 11.251, 11.252, 11.253, 11.254, 11.255, 11.256, 11.257, 11.258,
11.259, 11.260, 11.261, 11.262, 11.263, 11.264, 11.265, 11.266, 11.267,
11.268, 11.269, 11.270, 11.271, 11.272, 11.273, 11.274, 11.275, 11.276,
11.277, 11.278, 11.279, 11.280, 11.281, 11.282, 11.283, 11.284, 11.285,
11.286, 11.287, 11.288, 11.289, 11.290, 11.291, 11.292, 11.293, 11.294,
11.295, 11.296, 11.297, 11.298, 11.299, 11.300, 11.301, 11.302, 11.303,
11.304, 11.305, 11.306, 11.307, 11.308, 11.309, 11.310, 11.311, 11.312,
11.313, 11.314, 11.315, 11.316, 11.317, 11.318, 11.319, 11.320, 11.321,
11.322, 11.323, 11.324, 11.325, 11.326, 11.327, 11.328, 11.329, 11.330,
11.331, 11.332, 11.333, 11.334, 11.335, 11.336, 11.337, 11.338, 11.339,
11.340, 11.341, 11.342, 11.343, 11.344, 11.345, 11.346, 11.347, 11.348,
11.349, 11.350, 11.351, 11.352, 11.353, 11.354, 11.355, 11.356, 11.357,
11.358, 11.359, 11.360, 11.361, 11.362, 11.363, 11.364, 11.365, 11.366,
11.367, 11.368, 11.369, 11.370, 11.371, 11.372, 11.373, 11.374, 11.375,
11.376, 11.377, 11.378, 11.379, 11.380, 11.381, 11.382, 11.383, 11.384,
11.385, 11.386, 11.387, 11.388, 11.389, 11.390, 11.391, 11.392, 11.393,
11.394, 11.395, 11.396, 11.397, 11.398, 11.399, 11.400, 11.401, 11.402,
11.403, 11.404, 11.405, 11.406, 11.407, 11.408, 11.409, 11.410, 11.411,
11.412, 11.413, 11.414, 11.415, 11.416, 11.417, 11.418, 11.419, 11.420,
11.421, 11.422, 11.423, 11.424, 11.425, 11.426, 11.427, 11.428, 11.429,
11.430, 11.431, 11.432, 11.433, 11.434, 11.435, 11.436, 11.437, 11.438,
11.439, 11.440, 11.441, 11.442, 11.443, 11.444, 11.445, 11.446, 11.447,
11.448, 11.449, 11.450, 11.451, 11.452, 11.453, 11.454, 11.455, 11.456,
11.457, 11.458, 11.459, 11.460, 11.461, 11.462, 11.463, 11.464, 11.465,
11.466, 11.467, 11.468, 11.469, 11.470, 11.471, 11.472, 11.473, 11.474,
11.475, 11.476, 11.477, 11.478, 11.479, 11.480, 11.481, 11.482, 11.483,
11.484, 11.485, 11.486, 11.487, 11.488, 11.489, 11.490, 11.491, 11.492,
11.493, 11.494, 11.495, 11.496, 11.497, 11.498, 11.499, 11.500, 11.501,
11.502, 11.503, 11.504, 11.505, 11.506, 11.507, 11.508, 11.509, 11.510,
11.511, 11.512, 11.513, 11.514, 11.515, 11.516, 11.517, 11.518, 11.519,
11.520, 11.521, 11.522, 11.523, 11.524, 11.525, 11.526, 11.527, 11.528,
11.529, 11.530, 11.531, 11.532, 11.533, 11.534, 11.535, 11.536, 11.537,
11.538, 11.539, 11.540, 11.541, 11.542, 11.543, 11.544, 11.545, 11.546,
11.547, 11.548, 11.549, 11.550, 11.551, 11.552, 11.553, 11.554, 11.555,
11.556, 11.557, 11.558, 11.559, 11.560, 11.561, 11.562, 11.563, 11.564,
11.565, 11.566, 11.567, 11.568, 11.569, 11.570, 11.571, 11.572, 11.573,
11.574, 11.575, 11.576, 11.577, 11.578, 11.579, 11.580, 11.581, 11.582,
11.583, 11.584, 11.585, 11.586, 11.587, 11.588, 11.589, 11.590, 11.591,
11.592, 11.593, 11.594, 11.595, 11.596, 11.597, 11.598, 11.599, 11.600,
11.601, 11.602, 11.603, 11.604, 11.605, 11.606, 11.607, 11.608, 11.609,
11.610, 11.611, 11.612, 11.613, 11.614, 11.615, 11.616, 11.617, 11.618,
11.619, 11.620, 11.621, 11.622, 11.623, 11.624, 11.625, 11.626, 11.627,
11.628, 11.629, 11.630, 11.631, 11.632, 11.633, 11.634, 11.635, 11.636,
11.637, 11.638, 11.639, 11.640, 11.641, 11.642, 11.643, 11.644, 11.645,
11.646, 11.647, 11.648, 11.649, 11.650, 11.651, 11.652, 11.653, 11.654,
11.655, 11.656, 11.657, 11.658, 11.659, 11.660, 11.661, 11.662, 11.663,
11.664, 11.665, 11.666, 11.667, 11.668, 11.669, 11.670, 11.671, 11.672,
11.673, 11.674, 11.675, 11.676, 11.677, 11.678, 11.679, 11.680, 11.681,
11.682, 11.683, 11.684, 11.685, 11.686, 11.687, 11.688, 11.689, 11.690,
11.691, 11.692, 11.693, 11.694, 11.695, 11.696, 11.697, 11.698, 11.699,
11.700, 11.701, 11.702, 11.703, 11.704, 11.705, 11.706, 11.707, 11.708,
11.709, 11.710, 11.711, 11.712, 11.713, 11.714, 11.715, 11.716, 11.717,
11.718, 11.719, 11.720, 11.721, 11.722, 11.723, 11.724, 11.725, 11.726,
11.727, 11.728, 11.729, 11.730, 11.731, 11.732, 11.733, 11.734, 11.735,
11.736, 11.737, 11.738, 11.739, 11.740, 11.741, 11.742, 11.743, 11.744,
11.745, 11.746, 11.747, 11.748, 11.749, 11.750, 11.751, 11.752, 11.753,
11.754, 11.755, 11.756, 11.757, 11.758, 11.759, 11.760, 11.761, 11.762,
11.763, 11.764, 11.765, 11.766, 11.767, 11.768, 11.769, 11.770, 11.771,
11.772, 11.773, 11.774, 11.775, 11.776, 11.777, 11.778, 11.779, 11.780,
11.781, 11.782, 11.783, 11.784, 11.785, 11.786, 11.787, 11.788, 11.789,
11.790, 11.791, 11.792, 11.793, 11.794, 11.795, 11.796, 11.797, 11.798,
11.799, 11.800, 11.801, 11.802, 11.803, 11.804, 11.805, 11.806, 11.807,
11.808, 11.809, 11.810, 11.811, 11.812, 11.813, 11.814, 11.815, 11.816,
11.817, 11.818, 11.819, 11.820, 11.821, 11.822, 11.823, 11.824, 11.825,
11.826, 11.827, 11.828, 11.829, 11.830, 11.831, 11.832, 11.833, 11.834,
11.835, 11.836, 11.837, 11.838, 11.839, 11.840, 11.841, 11.842, 11.843,
11.844, 11.845, 11.846, 11.847, 11.848, 11.849, 11.850, 11.851, 11.852,
11.853, 11.854, 11.855, 11.856, 11.857, 11.858, 11.859, 11.860, 11.861,
11.862, 11.863, 11.864, 11.865, 11.866, 11.867, 11.868, 11.869, 11.870,
11.871, 11.872, 11.873, 11.874, 11.875, 11.876, 11.877, 11.878, 11.879,
11.880, 11.881, 11.882, 11.883, 11.884, 11.885, 11.886, 11.887, 11.888,
11.889, 11.890, 11.891, 11.892, 11.893, 11.894, 11.895, 11.896, 11.897,
11.898, 11.899, 11.900, 11.901, 11.902, 11.903, 11.904, 11.905, 11.906,
11.907, 11.908, 11.909, 11.910, 11.911, 11.912, 11.913, 11.914, 11.915,
11.916, 11.917, 11.918, 11.919, 11.920, 11.921, 11.922, 11.923, 11.924,
11.925, 11.926, 11.927, 11.928, 11.929, 11.930, 11.931, 11.932, 11.933,
11.934, 11.935, 11.936, 11.937, 11.938, 11.939, 11.940, 11.941, 11.942,
11.943, 11.944, 11.945, 11.946, 11.947, 11.948, 11.949, 11.950, 11.951,
11.952, 11.953, 11.954, 11.955, 11.956, 11.957, 11.958, 11.959, 11.960,
11.961, 11.962, 11.963, 11.964, 11.965, 11.966, 11.967, 11.968, 11.969,
11.970, 11.971, 11.972, 11.973, 11.974, 11.975, 11.976, 11.977, 11.978,
11.979, 11.980, 11.981, 11.982, 11.983, 11.984, 11.985, 11.986, 11.987,
11.988, 11.989, 11.990, 11.991, 11.992, 11.993, 11.994, 11.995, 11.996,
11.997, 11.998, 11.999, 12.000, 12.001, 12.002, 12.003, 12.004, 12.005,
12.006, 12.007, 12.008, 12.009, 12.010, 12.011, 12.012, 12.013, 12.014,
12.015, 12.016, 12.017, 12.018, 12.019, 12.020, 12.021, 12.022, 12.023,
12.024, 12.025, 12.026, 12.027, 12.028, 12.029, 12.030, 12.031, 12.032,
12.033, 12.034, 12.035, 12.036, 12.037, 12.038, 12.039, 12.040, 12.041,
12.042, 12.043, 12.044, 12.045, 12.046, 12.047, 12.048, 12.049, 12.050,
12.051, 12.052, 12.053, 12.054, 12.055, 12.056, 12.057, 12.058, 12.059,
12.060, 12.061, 12.062, 12.063, 12.064, 12.065, 12.066, 12.067, 12.068,
12.069, 12.070, 12.071, 12.072, 12.073, 12.074, 12.075, 12.076, 12.077,
12.078, 12.079, 12.080, 12.081, 12.082, 12.083, 12.084, 12.085, 12.086,
12.087, 12.088, 12.089, 12.090, 12.091, 12.092, 12.093, 12.094, 12.095,
12.096, 12.097, 12.098, 12.099, 12.100, 12.101, 12.102, 12.103, 12.104,
12.105, 12.106, 12.107, 12.108, 12.109, 12.110, 12.111, 12.112, 12.113,
12.114, 12.115, 12.116, 12.117, 12.118, 12.119, 12.120, 12.121, 12.122,
12.123, 12.124, 12.125, 12.126, 12.127, 12.128, 12.129, 12.130, 12.131,
12.132, 12.133, 12.134, 12.135, 12.136, 12.137, 12.138, 12.139, 12.140,
12.141, 12.142, 12.143, 12.144, 12.145, 12.146, 12.147, 12.148, 12.149,
12.150, 12.151, 12.152, 12.153, 12.154, 12.155, 12.156, 12.157, 12.158,
12.159, 12.160, 12.161, 12.162, 12.163, 12.164, 12.165, 12.166, 12.167,
12.168, 12.169, 12.170, 12.171, 12.172, 12.173, 12.174, 12.175, 12.176,
12.177, 12.178, 12.179, 12.180, 12.181, 12.182, 12.183, 12.184, 12.185,
12.186, 12.187, 12.188, 12.189, 12.190, 12.191, 12.192, 12.193, 12.194,
12.195, 12.196, 12.197, 12.198, 12.199, 12.200, 12.201, 12.202, 12.203,
12.204, 12.205, 12.206, 12.207, 12.208, 12.209, 12.210, 12.211, 12.212,
12.213, 12.214, 12.215, 12.216, 12.217, 12.218, 12.219, 12.220, 12.221,
12.222, 12.223, 12.224, 12.225, 12.226, 12.227, 12.228, 12.229, 12.230,
12.231, 12.232, 12.233, 12.234, 12.235, 12.236, 12.237, 12.238, 12.239,
12.240, 12.241, 12.242, 12.243, 12.244, 12.245, 12.246, 12.247, 12.248,
12.249, 12.250, 12.251, 12.252, 12.253, 12.254, 12.255, 12.256, 12.257,
12.258, 12.259, 12.260, 12.261, 12.262, 12.263, 12.264, 12.265, 12.266,
12.267, 12.268, 12.269, 12.270, 12.271, 12.272, 12.273, 12.274, 12.275,
12.276, 12.277, 12.278, 12.279, 12.280, 12.281, 12.282, 12.283, 12.284,
12.285, 12.286, 12.287, 12.288, 12.289, 12.290, 12.291, 12.292, 12.293,
12.294, 12.295, 12.296, 12.297, 12.298, 12.299, 12.300, 12.301, 12.302,
12.303, 12.304, 12.305, 12.306, 12.307, 12.308, 12.309, 12.310, 12.311,
12.312, 12.313, 12.314, 12.315, 12.316, 12.317, 12.318, 12.319, 12.320,
12.321, 12.322, 12.323, 12.324, 12.325, 12.326, 12.327, 12.328, 12.329,
12.330, 12.331, 12.332, 12.333, 12.334, 12.335, 12.336, 12.337, 12.338,
12.339, 12.340, 12.341, 12.342, 12.343, 12.344, 12.345, 12.346, 12.347,
12.348, 12.349, 12.350, 12.351, 12.352, 12.353, 12.354, 12.355, 12.356,
12.357, 12.358, 12.359, 12.360, 12.361, 12.362, 12.363, 12.364, 12.365,
12.366, 12.367, 12.368, 12.369, 12.370, 12.371, 12.372, 12.373, 12.374,
12.375, 12.376, 12.377, 12.378, 12.379, 12.380, 12.381, 12.382, 12.383,
12.384, 12.385, 12.386, 12.387, 12.388, 12.389, 12.390, 12.391, 12.392,
12.393, 12.394, 12.395, 12.396, 12.397, 12.398, 12.399, 12.400, 12.401,
12.402, 12.403, 12.404, 12.405, 12.406, 12.407, 12.408, 12.409, 12.410,
12.411, 12.412, 12.413, 12.414, 12.415, 12.416, 12.417, 12.418, 12.419,
12.420, 12.421, 12.422, 12.423, 12.424, 12.425, 12.426, 12.427, 12.428,
12.429, 12.430, 12.431, 12.432, 12.433, 12.434, 12.435, 12.436, 12.437,
12.438, 12.439, 12.440, 12.441, 12.442, 12.443, 12.444, 12.445, 12.446,
12.447, 12.448, 12.449, 12.450, 12.451, 12.452, 12.453, 12.454, 12.455,
12.456, 12.457, 12.458, 12.459, 12.460, 12.461, 12.462, 12.463, 12.464,
12.465, 12.466, 12.467, 12.468, 12.469, 12.470, 12.471, 12.472, 12.473,
12.474, 12.475, 12.476, 12.477, 12.478, 12.479, 12.480, 12.481, 12.482,
12.483, 12.484, 12.485, 12.486, 12.487, 12.488, 12.489, 12.490, 12.491,
12.492, 12.493, 12.494, 12.495, 12.496, 12.497, 12.498, 12.499, 12.500,
12.501, 12.502, 12.503, 12.504, 12.505, 12.506, 12.507, 12.508, 12.509,
12.510, 12.511, 12.512, 12.513, 12.514, 12.515, 12.516, 12.517, 12.518,
12.519, 12.520, 12.521, 12.522, 12.523, 12.524, 12.525, 12.526, 12.527,
12.528, 12.529, 12.530, 12.531, 12.532, 12.533, 12.534, 12.535, 12.536,
12.537, 12.538, 12.539, 12.540, 12.541, 12.542, 12.543, 12.544, 12.545,
12.546, 12.547, 12.548, 12.549, 12.550, 12.551, 12.552, 12.553, 12.554,
12.555, 12.556, 12.557, 12.558, 12.559, 12.560, 12.561, 12.562, 12.563,
12.564, 12.565, 12.566, 12.567, 12.568, 12.569, 12.570, 12.571, 12.572,
12.573, 12.574, 12.575, 12.576, 12.577, 12.578, 12.579, 12.580, 12.581,
12.582, 12.583, 12.584, 12.585, 12.586, 12.587, 12.588, 12.589, 12.590,
12.591, 12.592, 12.593, 12.594, 12.595, 12.596, 12.597, 12.598, 12.599,
12.600, 12.601, 12.602, 12.603, 12.604, 12.605, 12.606, 12.607, 12.608,
12.609, 12.610, 12.611, 12.612, 12.613, 12.614, 12.615, 12.616, 12.617,
12.618, 12.619, 12.620, 12.621, 12.622, 12.623, 12.624, 12.625, 12.626,
12.627, 12.628, 12.629, 12.630, 12.631, 12.632, 12.633, 12.634, 12.635,
12.636, 12.637, 12.638, 12.639, 12.640, 12.641, 12.642, 12.643, 12.644,
12.645, 12.646, 12.647, 12.648, 12.649, 12.650, 12.651, 12.652, 12.653,
12.654, 12.655, 12.656, 12.657, 12.658, 12.659, 12.660, 12.661, 12.662,
12.663, 12.664, 12.665, 12.666, 12.667, 12.668, 12.669, 12.670, 12.671,
12.672, 12.673, 12.674, 12.675, 12.676, 12.677, 12.678, 12.679, 12.680,
12.681, 12.682, 12.683, 12.684, 12.685, 12.686, 12.687, 12.688, 12.689,
12.690, 12.691, 12.692, 12.693, 12.694, 12.695, 12.696, 12.697, 12.698,
12.699, 12.700, 12.701, 12.702, 12.703, 12.704, 12.705, 12.706, 12.707,
12.708, 12.709, 12.710, 12.711, 12.712, 12.713, 12.714, 12.715, 12.716,
12.717, 12.718, 12.719, 12.720, 12.721, 12.722, 12.723, 12.724, 12.725,
12.726, 12.727, 12.728, 12.729, 12.730, 12.731, 12.732, 12.733, 12.734,
12.735, 12.736, 12.737, 12.738, 12.739, 12.740, 12.741, 12.742, 12.743,
12.744, 12.745, 12.746, 12.747, 12.748, 12.749, 12.750, 12.751, 12.752,
12.753, 12.754, 12.755, 12.756, 12.757, 12.758, 12.759, 12.760, 12.761,
12.762, 12.763, 12.764, 12.765, 12.766, 12.767, 12.768, 12.769, 12.770,
12.771, 12.772, 12.773, 12.774, 12.775, 12.776, 12.777, 12.778, 12.779,
12.780, 12.781, 12.782, 12.783, 12.784, 12.785, 12.786, 12.787, 12.788,
12.789, 12.790, 12.791, 12.792, 12.793, 12.794, 12.795, 12.796, 12.797,
12.798, 12.799, 12.800, 12.801, 12.802, 12.803, 12.804, 12.805, 12.806,
12.807, 12.808, 12.809, 12.810, 12.811, 12.812, 12.813, 12.814, 12.815,
12.816, 12.817, 12.818, 12.819, 12.820, 12.821, 12.822, 12.823, 12.824,
12.825, 12.826, 12.827, 12.828, 12.829, 12.830, 12.831, 12.832, 12.833,
12.834, 12.835, 12.836, 12.837, 12.838, 12.839, 12.840, 12.841, 12.842,
12.843, 12.844, 12.845, 12.846, 12.847, 12.848, 12.849, 12.850, 12.851,
12.852, 12.853, 12.854, 12.855, 12.856, 12.857, 12.858, 12.859, 12.860,
12.861, 12.862, 12.863, 12.864, 12.865, 12.866, 12.867, 12.868, 12.869,
12.870, 12.871, 12.872, 12.873, 12.874, 12.875, 12.876, 12.877, 12.878,
12.879, 12.880, 12.881, 12.882, 12.883, 12.884, 12.885, 12.886, 12.887,
12.888, 12.889, 12.890, 12.891, 12.892, 12.893, 12.894, 12.895, 12.896,
12.897, 12.898, 12.899, 12.900, 12.901, 12.902, 12.903, 12.904, 12.905,
12.906, 12.907, 12.908, 12.909, 12.910, 12.911, 12.912, 12.913, 12.914,
12.915, 12.916, 12.917, 12.918, 12.919, 12.920, 12.921, 12.922, 12.923,
12.924, 12.925, 12.926, 12.927, 12.928, 12.929, 12.930, 12.931, 12.932,
12.933, 12.934, 12.935, 12.936, 12.937, 12.938, 12.939, 12.940, 12.941,
12.942, 12.943, 12.944, 12.945, 12.946, 12.947, 12.948, 12.949, 12.950,
12.951, 12.952, 12.953, 12.954, 12.955, 12.956, 12.957, 12.958, 12.959,
12.960, 12.961, 12.962, 12.963, 12.964, 12.965, 12.966, 12.967, 12.968,
12.969, 12.970, 12.971, 12.972, 12.973, 12.974, 12.975, 12.976, 12.977,
12.978, 12.979, 12.980, 12.981, 12.982, 12.983, 12.984, 12.985, 12.986,
12.987, 12.988, 12.989, 12.990, 12.991, 12.992, 12.993, 12.994, 12.995,
12.996, 12.997, 12.998, 12.999, 13.000, 13.001, 13.002, 13.003, 13.004,
13.005, 13.006, 13.007, 13.008, 13.009, 13.010, 13.011, 13.012, 13.013,
13.014, 13.015, 13.016, 13.017, 13.018, 13.019, 13.020, 13.021, 13.022,
13.023, 13.024, 13.025, 13.026, 13.027, 13.028, 13.029, 13.030, 13.031,
13.032, 13.033, 13.034, 13.035, 13.036, 13.037, 13.038, 13.039, 13.040,
13.041, 13.042, 13.043, 13.044, 13.045, 13.046, 13.047, 13.048, 13.049,
13.050, 13.051, 13.052, 13.053, 13.054, 13.055, 13.056, 13.057, 13.058,
13.059, 13.060, 13.061, 13.062, 13.063, 13.064, 13.065, 13.066, 13.067,
13.068, 13.069, 13.070, 13.071, 13.072, 13.073, 13.074, 13.075, 13.076,
13.077, 13.078, 13.079, 13.080, 13.081, 13.082, 13.083, 13.084, 13.085,
13.086, 13.087, 13.088, 13.089, 13.090, 13.091, 13.092, 13.093, 13.094,
13.095, 13.096, 13.097, 13.098, 13.099, 13.100, 13.101, 13.102, 13.103,
13.104, 13.105, 13.106, 13.107, 13.108, 13.109, 13.110, 13.111, 13.112,
13.113, 13.114, 13.115, 13.116, 13.117, 13.118, 13.119, 13.120, 13.121,
13.122, 13.123, 13.124, 13.125, 13.126, 13.127, 13.128, 13.129, 13.130,
13.131, 13.132, 13.133, 13.134, 13.135, 13.136, 13.137, 13.138, 13.139,
13.140, 13.141, 13.142, 13.143, 13.144, 13.145, 13.146, 13.147, 13.148,
13.149, 13.150, 13.151, 13.152, 13.153, 13.154, 13.155, 13.156, 13.157,
13.158, 13.159, 13.160, 13.161, 13.162, 13.163, 13.164, 13.165, 13.166,
13.167, 13.168, 13.169, 13.170, 13.171, 13.172, 13.173, 13.174, 13.175,
13.176, 13.177, 13.178, 13.179, 13.180, 13.181, 13.182, 13.183, 13.184,
13.185, 13.186, 13.187, 13.188, 13.189, 13.190, 13.191, 13.192, 13.193,
13.194, 13.195, 13.196, 13.197, 13.198, 13.199, 13.200, 13.201, 13.202,
13.203, 13.204, 13.205, 13.206, 13.207, 13.208, 13.209, 13.210, 13.211,
13.212, 13.213, 13.214, 13.215, 13.216, 13.217, 13.218, 13.219, 13.220,
13.221, 13.222, 13.223, 13.224, 13.225, 13.226, 13.227, 13.228, 13.229,
13.230, 13.231, 13.232, 13.233, 13.234, 13.235, 13.236, 13.237, 13.238,
13.239, 13.240, 13.241, 13.242, 13.243, 13.244, 13.245, 13.246, 13.247,
13.248, 13.249, 13.250, 13.251, 13.252, 13.253, 13.254, 13.255, 13.256,
13.257, 13.258, 13.259, 13.260, 13.261, 13.262, 13.263, 13.264, 13.265,
13.266, 13.267, 13.268, 13.269, 13.270, 13.271, 13.272, 13.273, 13.274,
13.275, 13.276, 13.277, 13.278, 13.279, 13.280, 13.281, 13.282, 13.283,
13.284, 13.285, 13.286, 13.287, 13.288, 13.289, 13.290, 13.291, 13.292,
13.293, 13.294, 13.295, 13.296, 13.297, 13.298, 13.299, 13.300, 13.301,
13.302, 13.303, 13.304, 13.305, 13.306, 13.307, 13.308, 13.309, 13.310,
13.311, 13.312, 13.313, 13.314, 13.315, 13.316, 13.317, 13.318, 13.319,
13.320, 13.321, 13.322, 13.323, 13.324, 13.325, 13.326, 13.327, 13.328,
13.329, 13.330, 13.331, 13.332, 13.333, 13.334, 13.335, 13.336, 13.337,
13.338, 13.339, 13.340],
    'response':[
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002,
0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002,
0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002,
0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002,
0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002,
0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002,
0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002,
0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.003, 0.003,
0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003,
0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003,
0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003,
0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003,
0.003, 0.003, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004,
0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004,
0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004,
0.004, 0.004, 0.004, 0.004, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005,
0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005,
0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.006, 0.006, 0.006, 0.006,
0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006,
0.006, 0.006, 0.006, 0.006, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007,
0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007,
0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008,
0.008, 0.008, 0.008, 0.008, 0.008, 0.009, 0.009, 0.009, 0.009, 0.009,
0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.010,
0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010,
0.010, 0.010, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011,
0.011, 0.011, 0.012, 0.012, 0.012, 0.012, 0.012, 0.012, 0.012, 0.012,
0.012, 0.012, 0.013, 0.013, 0.013, 0.013, 0.013, 0.013, 0.013, 0.013,
0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.015, 0.015,
0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.016, 0.016, 0.016, 0.016,
0.016, 0.016, 0.016, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.018,
0.018, 0.018, 0.018, 0.018, 0.018, 0.019, 0.019, 0.019, 0.019, 0.019,
0.019, 0.020, 0.020, 0.020, 0.020, 0.020, 0.021, 0.021, 0.021, 0.021,
0.021, 0.022, 0.022, 0.022, 0.022, 0.022, 0.023, 0.023, 0.023, 0.023,
0.024, 0.024, 0.024, 0.024, 0.025, 0.025, 0.025, 0.025, 0.025, 0.026,
0.026, 0.026, 0.027, 0.027, 0.027, 0.027, 0.028, 0.028, 0.028, 0.028,
0.029, 0.029, 0.029, 0.030, 0.030, 0.030, 0.030, 0.031, 0.031, 0.031,
0.032, 0.032, 0.032, 0.033, 0.033, 0.033, 0.034, 0.034, 0.034, 0.035,
0.035, 0.035, 0.036, 0.036, 0.036, 0.037, 0.037, 0.037, 0.038, 0.038,
0.038, 0.039, 0.039, 0.040, 0.040, 0.040, 0.041, 0.041, 0.042, 0.042,
0.042, 0.043, 0.043, 0.044, 0.044, 0.045, 0.045, 0.046, 0.046, 0.046,
0.047, 0.047, 0.048, 0.048, 0.049, 0.049, 0.050, 0.050, 0.051, 0.051,
0.052, 0.052, 0.053, 0.053, 0.054, 0.055, 0.055, 0.056, 0.056, 0.057,
0.057, 0.058, 0.059, 0.059, 0.060, 0.060, 0.061, 0.062, 0.062, 0.063,
0.064, 0.064, 0.065, 0.066, 0.066, 0.067, 0.068, 0.068, 0.069, 0.070,
0.070, 0.071, 0.072, 0.073, 0.073, 0.074, 0.075, 0.076, 0.076, 0.077,
0.078, 0.079, 0.079, 0.080, 0.081, 0.082, 0.083, 0.084, 0.084, 0.085,
0.086, 0.087, 0.088, 0.089, 0.090, 0.091, 0.092, 0.093, 0.094, 0.095,
0.096, 0.097, 0.098, 0.099, 0.100, 0.101, 0.102, 0.103, 0.104, 0.105,
0.106, 0.107, 0.108, 0.109, 0.111, 0.112, 0.113, 0.114, 0.115, 0.117,
0.118, 0.119, 0.120, 0.121, 0.123, 0.124, 0.125, 0.127, 0.128, 0.129,
0.131, 0.132, 0.133, 0.135, 0.136, 0.138, 0.139, 0.140, 0.142, 0.143,
0.145, 0.146, 0.148, 0.149, 0.151, 0.153, 0.154, 0.156, 0.157, 0.159,
0.161, 0.162, 0.164, 0.166, 0.167, 0.169, 0.171, 0.173, 0.174, 0.176,
0.178, 0.180, 0.182, 0.183, 0.185, 0.187, 0.189, 0.191, 0.193, 0.195,
0.197, 0.199, 0.201, 0.203, 0.205, 0.207, 0.209, 0.211, 0.214, 0.216,
0.218, 0.220, 0.222, 0.224, 0.227, 0.229, 0.231, 0.234, 0.236, 0.238,
0.241, 0.243, 0.245, 0.248, 0.250, 0.253, 0.255, 0.258, 0.260, 0.263,
0.265, 0.268, 0.270, 0.273, 0.276, 0.278, 0.281, 0.284, 0.286, 0.289,
0.292, 0.295, 0.297, 0.300, 0.303, 0.306, 0.309, 0.312, 0.315, 0.318,
0.320, 0.323, 0.326, 0.329, 0.332, 0.335, 0.339, 0.342, 0.345, 0.348,
0.351, 0.354, 0.357, 0.360, 0.364, 0.367, 0.370, 0.373, 0.376, 0.380,
0.383, 0.386, 0.390, 0.393, 0.396, 0.400, 0.403, 0.407, 0.410, 0.413,
0.417, 0.420, 0.424, 0.427, 0.431, 0.434, 0.438, 0.441, 0.445, 0.448,
0.452, 0.456, 0.459, 0.463, 0.466, 0.470, 0.474, 0.477, 0.481, 0.485,
0.488, 0.492, 0.496, 0.499, 0.503, 0.507, 0.511, 0.514, 0.518, 0.522,
0.525, 0.529, 0.533, 0.537, 0.540, 0.544, 0.548, 0.552, 0.555, 0.559,
0.563, 0.567, 0.570, 0.574, 0.578, 0.582, 0.585, 0.589, 0.593, 0.597,
0.600, 0.604, 0.608, 0.612, 0.615, 0.619, 0.623, 0.626, 0.630, 0.634,
0.637, 0.641, 0.645, 0.648, 0.652, 0.656, 0.659, 0.663, 0.666, 0.670,
0.674, 0.677, 0.681, 0.684, 0.688, 0.691, 0.695, 0.698, 0.702, 0.705,
0.709, 0.712, 0.715, 0.719, 0.722, 0.726, 0.729, 0.732, 0.736, 0.739,
0.742, 0.745, 0.749, 0.752, 0.755, 0.758, 0.761, 0.765, 0.768, 0.771,
0.774, 0.777, 0.780, 0.783, 0.786, 0.789, 0.792, 0.795, 0.798, 0.801,
0.804, 0.807, 0.810, 0.812, 0.815, 0.818, 0.821, 0.823, 0.826, 0.829,
0.832, 0.834, 0.837, 0.839, 0.842, 0.845, 0.847, 0.850, 0.852, 0.854,
0.857, 0.859, 0.862, 0.864, 0.866, 0.869, 0.871, 0.873, 0.875, 0.878,
0.880, 0.882, 0.884, 0.886, 0.888, 0.890, 0.892, 0.894, 0.896, 0.898,
0.900, 0.902, 0.904, 0.906, 0.908, 0.910, 0.912, 0.913, 0.915, 0.917,
0.919, 0.920, 0.922, 0.924, 0.925, 0.927, 0.928, 0.930, 0.931, 0.933,
0.934, 0.936, 0.937, 0.939, 0.940, 0.941, 0.943, 0.944, 0.945, 0.947,
0.948, 0.949, 0.951, 0.952, 0.953, 0.954, 0.955, 0.956, 0.957, 0.959,
0.960, 0.961, 0.962, 0.963, 0.964, 0.965, 0.966, 0.967, 0.968, 0.968,
0.969, 0.970, 0.971, 0.972, 0.973, 0.974, 0.974, 0.975, 0.976, 0.977,
0.977, 0.978, 0.979, 0.980, 0.980, 0.981, 0.982, 0.982, 0.983, 0.983,
0.984, 0.985, 0.985, 0.986, 0.986, 0.987, 0.987, 0.988, 0.988, 0.989,
0.989, 0.990, 0.990, 0.991, 0.991, 0.992, 0.992, 0.992, 0.993, 0.993,
0.993, 0.994, 0.994, 0.994, 0.995, 0.995, 0.995, 0.996, 0.996, 0.996,
0.996, 0.997, 0.997, 0.997, 0.997, 0.998, 0.998, 0.998, 0.998, 0.998,
0.998, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999,
1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
1.000, 1.000, 1.000, 1.000, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999,
0.999, 0.999, 0.999, 0.999, 0.999, 0.998, 0.998, 0.998, 0.998, 0.998,
0.998, 0.998, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.996, 0.996,
0.996, 0.996, 0.996, 0.996, 0.995, 0.995, 0.995, 0.995, 0.995, 0.994,
0.994, 0.994, 0.994, 0.994, 0.993, 0.993, 0.993, 0.993, 0.993, 0.992,
0.992, 0.992, 0.992, 0.991, 0.991, 0.991, 0.991, 0.990, 0.990, 0.990,
0.990, 0.990, 0.989, 0.989, 0.989, 0.989, 0.988, 0.988, 0.988, 0.988,
0.987, 0.987, 0.987, 0.987, 0.986, 0.986, 0.986, 0.986, 0.985, 0.985,
0.985, 0.985, 0.984, 0.984, 0.984, 0.983, 0.983, 0.983, 0.983, 0.982,
0.982, 0.982, 0.982, 0.981, 0.981, 0.981, 0.981, 0.980, 0.980, 0.980,
0.980, 0.979, 0.979, 0.979, 0.979, 0.978, 0.978, 0.978, 0.978, 0.977,
0.977, 0.977, 0.977, 0.976, 0.976, 0.976, 0.976, 0.975, 0.975, 0.975,
0.975, 0.974, 0.974, 0.974, 0.974, 0.973, 0.973, 0.973, 0.973, 0.972,
0.972, 0.972, 0.972, 0.971, 0.971, 0.971, 0.971, 0.970, 0.970, 0.970,
0.970, 0.969, 0.969, 0.969, 0.969, 0.968, 0.968, 0.968, 0.967, 0.967,
0.967, 0.967, 0.966, 0.966, 0.966, 0.966, 0.965, 0.965, 0.965, 0.965,
0.964, 0.964, 0.964, 0.963, 0.963, 0.963, 0.963, 0.962, 0.962, 0.962,
0.962, 0.961, 0.961, 0.961, 0.960, 0.960, 0.960, 0.960, 0.959, 0.959,
0.959, 0.959, 0.958, 0.958, 0.958, 0.958, 0.957, 0.957, 0.957, 0.956,
0.956, 0.956, 0.956, 0.956, 0.955, 0.955, 0.955, 0.955, 0.954, 0.954,
0.954, 0.954, 0.954, 0.953, 0.953, 0.953, 0.953, 0.953, 0.953, 0.953,
0.952, 0.952, 0.952, 0.952, 0.952, 0.952, 0.952, 0.952, 0.952, 0.951,
0.951, 0.951, 0.951, 0.951, 0.951, 0.951, 0.951, 0.951, 0.951, 0.951,
0.951, 0.951, 0.951, 0.951, 0.951, 0.951, 0.951, 0.951, 0.951, 0.951,
0.951, 0.952, 0.952, 0.952, 0.952, 0.952, 0.952, 0.952, 0.952, 0.952,
0.952, 0.952, 0.952, 0.952, 0.953, 0.953, 0.953, 0.953, 0.953, 0.953,
0.953, 0.953, 0.953, 0.953, 0.954, 0.954, 0.954, 0.954, 0.954, 0.954,
0.954, 0.954, 0.954, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.955,
0.955, 0.955, 0.956, 0.956, 0.956, 0.956, 0.956, 0.956, 0.956, 0.956,
0.956, 0.956, 0.957, 0.957, 0.957, 0.957, 0.957, 0.957, 0.957, 0.957,
0.957, 0.957, 0.958, 0.958, 0.958, 0.958, 0.958, 0.958, 0.958, 0.958,
0.958, 0.958, 0.958, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959,
0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.960, 0.960, 0.960, 0.960,
0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.960,
0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.961, 0.961,
0.961, 0.961, 0.961, 0.961, 0.961, 0.961, 0.961, 0.961, 0.961, 0.961,
0.961, 0.961, 0.961, 0.961, 0.961, 0.961, 0.961, 0.961, 0.960, 0.960,
0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.960,
0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.959,
0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959,
0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959,
0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959,
0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959,
0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959,
0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959,
0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.960,
0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.959,
0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959,
0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.959, 0.958, 0.958, 0.958,
0.958, 0.958, 0.958, 0.958, 0.958, 0.958, 0.958, 0.958, 0.957, 0.957,
0.957, 0.957, 0.957, 0.957, 0.957, 0.957, 0.957, 0.956, 0.956, 0.956,
0.956, 0.956, 0.956, 0.956, 0.956, 0.956, 0.956, 0.955, 0.955, 0.955,
0.955, 0.955, 0.955, 0.955, 0.955, 0.955, 0.954, 0.954, 0.954, 0.954,
0.954, 0.954, 0.954, 0.954, 0.954, 0.953, 0.953, 0.953, 0.953, 0.953,
0.953, 0.953, 0.953, 0.952, 0.952, 0.952, 0.952, 0.952, 0.952, 0.952,
0.952, 0.952, 0.951, 0.951, 0.951, 0.951, 0.951, 0.951, 0.951, 0.951,
0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.949, 0.949,
0.949, 0.949, 0.949, 0.949, 0.949, 0.949, 0.948, 0.948, 0.948, 0.948,
0.948, 0.948, 0.948, 0.948, 0.947, 0.947, 0.947, 0.947, 0.947, 0.947,
0.947, 0.947, 0.947, 0.946, 0.946, 0.946, 0.946, 0.946, 0.946, 0.946,
0.946, 0.946, 0.946, 0.945, 0.945, 0.945, 0.945, 0.945, 0.945, 0.945,
0.945, 0.945, 0.945, 0.945, 0.945, 0.945, 0.945, 0.944, 0.944, 0.944,
0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944,
0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944,
0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944,
0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944,
0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944,
0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944,
0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944,
0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944,
0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944, 0.944,
0.944, 0.944, 0.944, 0.944, 0.945, 0.945, 0.945, 0.945, 0.945, 0.945,
0.945, 0.945, 0.945, 0.945, 0.945, 0.945, 0.945, 0.945, 0.945, 0.945,
0.945, 0.945, 0.945, 0.945, 0.945, 0.945, 0.945, 0.945, 0.945, 0.946,
0.946, 0.946, 0.946, 0.946, 0.946, 0.946, 0.946, 0.946, 0.946, 0.946,
0.946, 0.946, 0.946, 0.946, 0.946, 0.946, 0.946, 0.946, 0.946, 0.946,
0.946, 0.946, 0.947, 0.947, 0.947, 0.947, 0.947, 0.947, 0.947, 0.947,
0.947, 0.947, 0.947, 0.947, 0.947, 0.947, 0.947, 0.947, 0.947, 0.947,
0.947, 0.947, 0.947, 0.947, 0.947, 0.947, 0.948, 0.948, 0.948, 0.948,
0.948, 0.948, 0.948, 0.948, 0.948, 0.948, 0.948, 0.948, 0.948, 0.948,
0.948, 0.948, 0.948, 0.948, 0.948, 0.948, 0.948, 0.948, 0.948, 0.948,
0.948, 0.948, 0.948, 0.948, 0.949, 0.949, 0.949, 0.949, 0.949, 0.949,
0.949, 0.949, 0.949, 0.949, 0.949, 0.949, 0.949, 0.949, 0.949, 0.949,
0.949, 0.949, 0.949, 0.949, 0.949, 0.949, 0.949, 0.949, 0.949, 0.949,
0.949, 0.949, 0.949, 0.949, 0.949, 0.949, 0.949, 0.949, 0.949, 0.949,
0.949, 0.949, 0.949, 0.949, 0.949, 0.949, 0.950, 0.950, 0.950, 0.950,
0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950,
0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950,
0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950,
0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950, 0.950,
0.950, 0.950, 0.950, 0.950, 0.949, 0.949, 0.949, 0.949, 0.949, 0.949,
0.949, 0.949, 0.949, 0.949, 0.948, 0.948, 0.948, 0.948, 0.948, 0.948,
0.947, 0.947, 0.947, 0.947, 0.947, 0.946, 0.946, 0.946, 0.946, 0.945,
0.945, 0.945, 0.945, 0.944, 0.944, 0.944, 0.943, 0.943, 0.943, 0.942,
0.942, 0.942, 0.941, 0.941, 0.940, 0.940, 0.939, 0.939, 0.939, 0.938,
0.938, 0.937, 0.937, 0.936, 0.936, 0.935, 0.935, 0.934, 0.934, 0.933,
0.932, 0.932, 0.931, 0.931, 0.930, 0.930, 0.929, 0.928, 0.928, 0.927,
0.926, 0.926, 0.925, 0.925, 0.924, 0.923, 0.923, 0.922, 0.921, 0.921,
0.920, 0.919, 0.919, 0.918, 0.917, 0.917, 0.916, 0.915, 0.915, 0.914,
0.914, 0.913, 0.912, 0.912, 0.911, 0.911, 0.910, 0.909, 0.909, 0.908,
0.908, 0.907, 0.907, 0.906, 0.906, 0.905, 0.905, 0.904, 0.904, 0.904,
0.903, 0.903, 0.902, 0.902, 0.902, 0.901, 0.901, 0.901, 0.901, 0.900,
0.900, 0.900, 0.900, 0.899, 0.899, 0.899, 0.899, 0.899, 0.899, 0.899,
0.899, 0.899, 0.899, 0.899, 0.899, 0.899, 0.899, 0.899, 0.899, 0.899,
0.899, 0.899, 0.899, 0.899, 0.900, 0.900, 0.900, 0.900, 0.900, 0.901,
0.901, 0.901, 0.901, 0.902, 0.902, 0.902, 0.902, 0.903, 0.903, 0.903,
0.904, 0.904, 0.905, 0.905, 0.905, 0.906, 0.906, 0.907, 0.907, 0.907,
0.908, 0.908, 0.909, 0.909, 0.910, 0.910, 0.911, 0.911, 0.912, 0.912,
0.913, 0.913, 0.914, 0.914, 0.915, 0.915, 0.916, 0.916, 0.917, 0.917,
0.918, 0.919, 0.919, 0.920, 0.920, 0.921, 0.921, 0.922, 0.923, 0.923,
0.924, 0.924, 0.925, 0.926, 0.926, 0.927, 0.928, 0.928, 0.929, 0.929,
0.930, 0.931, 0.931, 0.932, 0.933, 0.933, 0.934, 0.934, 0.935, 0.936,
0.936, 0.937, 0.938, 0.938, 0.939, 0.939, 0.940, 0.941, 0.941, 0.942,
0.942, 0.943, 0.944, 0.944, 0.945, 0.945, 0.946, 0.946, 0.947, 0.947,
0.948, 0.948, 0.949, 0.949, 0.950, 0.950, 0.951, 0.951, 0.952, 0.952,
0.953, 0.953, 0.954, 0.954, 0.955, 0.955, 0.955, 0.956, 0.956, 0.957,
0.957, 0.957, 0.958, 0.958, 0.958, 0.959, 0.959, 0.959, 0.960, 0.960,
0.960, 0.960, 0.961, 0.961, 0.961, 0.961, 0.962, 0.962, 0.962, 0.962,
0.962, 0.962, 0.963, 0.963, 0.963, 0.963, 0.963, 0.963, 0.963, 0.964,
0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964,
0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964,
0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964,
0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964,
0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964,
0.964, 0.964, 0.964, 0.964, 0.963, 0.963, 0.963, 0.963, 0.963, 0.963,
0.963, 0.963, 0.963, 0.963, 0.963, 0.963, 0.963, 0.963, 0.963, 0.963,
0.963, 0.963, 0.963, 0.963, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964,
0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964,
0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964,
0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964,
0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964,
0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964,
0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.964, 0.963, 0.963,
0.963, 0.963, 0.963, 0.963, 0.963, 0.963, 0.963, 0.963, 0.962, 0.962,
0.962, 0.962, 0.962, 0.962, 0.962, 0.962, 0.961, 0.961, 0.961, 0.961,
0.961, 0.961, 0.960, 0.960, 0.960, 0.960, 0.960, 0.959, 0.959, 0.959,
0.959, 0.959, 0.958, 0.958, 0.958, 0.958, 0.957, 0.957, 0.957, 0.957,
0.956, 0.956, 0.956, 0.955, 0.955, 0.955, 0.955, 0.954, 0.954, 0.954,
0.953, 0.953, 0.953, 0.952, 0.952, 0.952, 0.951, 0.951, 0.951, 0.950,
0.950, 0.950, 0.949, 0.949, 0.949, 0.948, 0.948, 0.948, 0.947, 0.947,
0.947, 0.947, 0.946, 0.946, 0.946, 0.945, 0.945, 0.945, 0.944, 0.944,
0.944, 0.943, 0.943, 0.943, 0.943, 0.942, 0.942, 0.942, 0.941, 0.941,
0.941, 0.941, 0.940, 0.940, 0.940, 0.940, 0.939, 0.939, 0.939, 0.939,
0.939, 0.938, 0.938, 0.938, 0.938, 0.938, 0.937, 0.937, 0.937, 0.937,
0.937, 0.937, 0.936, 0.936, 0.936, 0.936, 0.936, 0.936, 0.936, 0.935,
0.935, 0.935, 0.935, 0.935, 0.935, 0.935, 0.935, 0.934, 0.934, 0.934,
0.934, 0.934, 0.934, 0.934, 0.934, 0.933, 0.933, 0.933, 0.933, 0.933,
0.933, 0.933, 0.933, 0.932, 0.932, 0.932, 0.932, 0.932, 0.932, 0.932,
0.932, 0.931, 0.931, 0.931, 0.931, 0.931, 0.931, 0.931, 0.930, 0.930,
0.930, 0.930, 0.930, 0.930, 0.929, 0.929, 0.929, 0.929, 0.929, 0.928,
0.928, 0.928, 0.928, 0.928, 0.928, 0.927, 0.927, 0.927, 0.927, 0.927,
0.926, 0.926, 0.926, 0.926, 0.925, 0.925, 0.925, 0.925, 0.925, 0.924,
0.924, 0.924, 0.924, 0.923, 0.923, 0.923, 0.923, 0.923, 0.922, 0.922,
0.922, 0.922, 0.921, 0.921, 0.921, 0.921, 0.921, 0.920, 0.920, 0.920,
0.920, 0.919, 0.919, 0.919, 0.919, 0.918, 0.918, 0.918, 0.918, 0.917,
0.917, 0.917, 0.917, 0.916, 0.916, 0.916, 0.916, 0.916, 0.915, 0.915,
0.915, 0.915, 0.914, 0.914, 0.914, 0.914, 0.913, 0.913, 0.913, 0.913,
0.912, 0.912, 0.912, 0.912, 0.912, 0.911, 0.911, 0.911, 0.911, 0.910,
0.910, 0.910, 0.910, 0.910, 0.909, 0.909, 0.909, 0.909, 0.909, 0.908,
0.908, 0.908, 0.908, 0.907, 0.907, 0.907, 0.907, 0.907, 0.906, 0.906,
0.906, 0.906, 0.906, 0.906, 0.905, 0.905, 0.905, 0.905, 0.905, 0.904,
0.904, 0.904, 0.904, 0.904, 0.904, 0.903, 0.903, 0.903, 0.903, 0.903,
0.902, 0.902, 0.902, 0.902, 0.902, 0.902, 0.902, 0.901, 0.901, 0.901,
0.901, 0.901, 0.901, 0.901, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900,
0.900, 0.899, 0.899, 0.899, 0.899, 0.899, 0.899, 0.899, 0.899, 0.898,
0.898, 0.898, 0.898, 0.898, 0.898, 0.898, 0.898, 0.898, 0.897, 0.897,
0.897, 0.897, 0.897, 0.897, 0.897, 0.897, 0.897, 0.897, 0.897, 0.896,
0.896, 0.896, 0.896, 0.896, 0.896, 0.896, 0.896, 0.896, 0.896, 0.896,
0.896, 0.896, 0.896, 0.896, 0.895, 0.895, 0.895, 0.895, 0.895, 0.895,
0.895, 0.895, 0.895, 0.895, 0.895, 0.895, 0.895, 0.895, 0.895, 0.895,
0.895, 0.895, 0.895, 0.895, 0.895, 0.895, 0.895, 0.895, 0.895, 0.895,
0.895, 0.895, 0.895, 0.895, 0.895, 0.895, 0.895, 0.895, 0.895, 0.895,
0.895, 0.895, 0.896, 0.896, 0.896, 0.896, 0.896, 0.896, 0.896, 0.896,
0.896, 0.896, 0.896, 0.896, 0.896, 0.896, 0.896, 0.896, 0.897, 0.897,
0.897, 0.897, 0.897, 0.897, 0.897, 0.897, 0.897, 0.897, 0.898, 0.898,
0.898, 0.898, 0.898, 0.898, 0.898, 0.898, 0.898, 0.898, 0.899, 0.899,
0.899, 0.899, 0.899, 0.899, 0.899, 0.899, 0.899, 0.899, 0.899, 0.900,
0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900,
0.900, 0.900, 0.900, 0.901, 0.901, 0.901, 0.901, 0.901, 0.901, 0.901,
0.901, 0.901, 0.901, 0.901, 0.901, 0.901, 0.901, 0.901, 0.901, 0.901,
0.901, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900, 0.900,
0.900, 0.900, 0.899, 0.899, 0.899, 0.899, 0.899, 0.899, 0.898, 0.898,
0.898, 0.898, 0.897, 0.897, 0.897, 0.897, 0.896, 0.896, 0.896, 0.895,
0.895, 0.895, 0.894, 0.894, 0.894, 0.893, 0.893, 0.892, 0.892, 0.891,
0.891, 0.890, 0.890, 0.889, 0.889, 0.888, 0.888, 0.887, 0.886, 0.886,
0.885, 0.884, 0.884, 0.883, 0.882, 0.881, 0.880, 0.880, 0.879, 0.878,
0.877, 0.876, 0.875, 0.874, 0.873, 0.872, 0.871, 0.870, 0.869, 0.868,
0.867, 0.866, 0.864, 0.863, 0.862, 0.861, 0.860, 0.858, 0.857, 0.856,
0.854, 0.853, 0.851, 0.850, 0.848, 0.847, 0.845, 0.844, 0.842, 0.841,
0.839, 0.837, 0.835, 0.834, 0.832, 0.830, 0.828, 0.827, 0.825, 0.823,
0.821, 0.819, 0.817, 0.815, 0.813, 0.811, 0.809, 0.806, 0.804, 0.802,
0.800, 0.798, 0.795, 0.793, 0.791, 0.788, 0.786, 0.784, 0.781, 0.779,
0.776, 0.774, 0.771, 0.768, 0.766, 0.763, 0.761, 0.758, 0.755, 0.752,
0.750, 0.747, 0.744, 0.741, 0.738, 0.735, 0.732, 0.729, 0.726, 0.723,
0.720, 0.717, 0.714, 0.711, 0.708, 0.705, 0.702, 0.698, 0.695, 0.692,
0.689, 0.685, 0.682, 0.679, 0.675, 0.672, 0.669, 0.665, 0.662, 0.658,
0.655, 0.651, 0.648, 0.644, 0.641, 0.637, 0.634, 0.630, 0.626, 0.623,
0.619, 0.616, 0.612, 0.608, 0.605, 0.601, 0.597, 0.593, 0.590, 0.586,
0.582, 0.579, 0.575, 0.571, 0.567, 0.564, 0.560, 0.556, 0.552, 0.548,
0.545, 0.541, 0.537, 0.533, 0.530, 0.526, 0.522, 0.518, 0.514, 0.511,
0.507, 0.503, 0.499, 0.496, 0.492, 0.488, 0.484, 0.480, 0.477, 0.473,
0.469, 0.466, 0.462, 0.458, 0.454, 0.451, 0.447, 0.443, 0.440, 0.436,
0.432, 0.429, 0.425, 0.422, 0.418, 0.414, 0.411, 0.407, 0.404, 0.400,
0.397, 0.393, 0.390, 0.386, 0.383, 0.379, 0.376, 0.373, 0.369, 0.366,
0.363, 0.359, 0.356, 0.353, 0.349, 0.346, 0.343, 0.340, 0.337, 0.333,
0.330, 0.327, 0.324, 0.321, 0.318, 0.315, 0.312, 0.309, 0.306, 0.303,
0.300, 0.297, 0.294, 0.291, 0.288, 0.286, 0.283, 0.280, 0.277, 0.274,
0.272, 0.269, 0.266, 0.264, 0.261, 0.258, 0.256, 0.253, 0.251, 0.248,
0.246, 0.243, 0.241, 0.238, 0.236, 0.233, 0.231, 0.229, 0.226, 0.224,
0.222, 0.219, 0.217, 0.215, 0.213, 0.211, 0.208, 0.206, 0.204, 0.202,
0.200, 0.198, 0.196, 0.194, 0.192, 0.190, 0.188, 0.186, 0.184, 0.182,
0.180, 0.178, 0.176, 0.174, 0.173, 0.171, 0.169, 0.167, 0.166, 0.164,
0.162, 0.160, 0.159, 0.157, 0.155, 0.154, 0.152, 0.151, 0.149, 0.147,
0.146, 0.144, 0.143, 0.141, 0.140, 0.138, 0.137, 0.135, 0.134, 0.133,
0.131, 0.130, 0.128, 0.127, 0.126, 0.124, 0.123, 0.122, 0.121, 0.119,
0.118, 0.117, 0.116, 0.114, 0.113, 0.112, 0.111, 0.110, 0.109, 0.107,
0.106, 0.105, 0.104, 0.103, 0.102, 0.101, 0.100, 0.099, 0.098, 0.097,
0.096, 0.095, 0.094, 0.093, 0.092, 0.091, 0.090, 0.089, 0.088, 0.087,
0.086, 0.085, 0.084, 0.084, 0.083, 0.082, 0.081, 0.080, 0.079, 0.078,
0.078, 0.077, 0.076, 0.075, 0.075, 0.074, 0.073, 0.072, 0.072, 0.071,
0.070, 0.069, 0.069, 0.068, 0.067, 0.067, 0.066, 0.065, 0.065, 0.064,
0.063, 0.063, 0.062, 0.061, 0.061, 0.060, 0.060, 0.059, 0.058, 0.058,
0.057, 0.057, 0.056, 0.055, 0.055, 0.054, 0.054, 0.053, 0.053, 0.052,
0.052, 0.051, 0.051, 0.050, 0.050, 0.049, 0.049, 0.048, 0.048, 0.047,
0.047, 0.046, 0.046, 0.045, 0.045, 0.045, 0.044, 0.044, 0.043, 0.043,
0.042, 0.042, 0.042, 0.041, 0.041, 0.040, 0.040, 0.040, 0.039, 0.039,
0.038, 0.038, 0.038, 0.037, 0.037, 0.037, 0.036, 0.036, 0.036, 0.035,
0.035, 0.035, 0.034, 0.034, 0.034, 0.033, 0.033, 0.033, 0.032, 0.032,
0.032, 0.031, 0.031, 0.031, 0.031, 0.030, 0.030, 0.030, 0.029, 0.029,
0.029, 0.029, 0.028, 0.028, 0.028, 0.027, 0.027, 0.027, 0.027, 0.026,
0.026, 0.026, 0.026, 0.026, 0.025, 0.025, 0.025, 0.025, 0.024, 0.024,
0.024, 0.024, 0.024, 0.023, 0.023, 0.023, 0.023, 0.022, 0.022, 0.022,
0.022, 0.022, 0.022, 0.021, 0.021, 0.021, 0.021, 0.021, 0.020, 0.020,
0.020, 0.020, 0.020, 0.019, 0.019, 0.019, 0.019, 0.019, 0.019, 0.019,
0.018, 0.018, 0.018, 0.018, 0.018, 0.018, 0.017, 0.017, 0.017, 0.017,
0.017, 0.017, 0.017, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016,
0.016, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.014,
0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.013, 0.013,
0.013, 0.013, 0.013, 0.013, 0.013, 0.013, 0.013, 0.013, 0.012, 0.012,
0.012, 0.012, 0.012, 0.012, 0.012, 0.012, 0.012, 0.012, 0.011, 0.011,
0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011,
0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010,
0.010, 0.010, 0.010, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009,
0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.008, 0.008,
0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008,
0.008, 0.008, 0.008, 0.008, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007,
0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007,
0.007, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006,
0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006,
0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005,
0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005,
0.005, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004,
0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004,
0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.003, 0.003, 0.003,
0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003,
0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003,
0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003,
0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003,
0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.002, 0.002, 0.002,
0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002,
0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002,
0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002,
0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002,
0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002,
0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002,
0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002,
0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002,
0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002,
0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002,
0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002,
0.002, 0.002, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
0.000]}

