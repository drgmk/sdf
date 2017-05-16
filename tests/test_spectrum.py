import os
import pytest
import numpy as np

from .context import sdf

def test_read_phoenix():
    dir = sdf.config.file['model_root']
    s = sdf.spectrum.ModelSpectrum.read_phoenix(
                    dir+'lte009.5-3.5-0.0.BT-Settl.7.bz2'
                                                )
    os.unlink(dir+'lte009.5-3.5-0.0.BT-Settl.7.bz2.npy')
    
    s = sdf.spectrum.ModelSpectrum.read_phoenix(
                    dir+'lte058-4.5-0.0a+0.0.BT-Settl.7.bz2'
                                                )
    os.unlink(dir+'lte058-4.5-0.0a+0.0.BT-Settl.7.bz2.npy')
