import pytest
import numpy as np

from .context import sdf

def test_read_phoenix():
    s = sdf.spectrum.ModelSpectrum.read_phoenix(
                    'models/lte009.5-3.5-0.0.BT-Settl.7.bz2'
                                                )
    s = sdf.spectrum.ModelSpectrum.read_phoenix(
                    'models/lte058-4.5-0.0a+0.0.BT-Settl.7.bz2'
                                                )
