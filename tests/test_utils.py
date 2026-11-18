import numpy as np
import pytest

from .context import sdf

def test_utils_plot_join_line():
    t = {'id':np.array(['a','b','c','b']),
         'x' :np.array([1  ,2  ,3  ,4  ]),
         'y' :np.array([5  ,6  ,7  ,8  ])}
    x,y = sdf.utils.plot_join_line(t,'id','x','y')
    assert(np.all(np.equal(x,[2,4])))
    assert(np.all(np.equal(y,[6,8])))


def test_spectral_resample_matrix_preserves_constant_flux():
    wave_in = np.geomspace(5.0, 15.0, 2000)
    wave_out = np.geomspace(6.0, 14.0, 100)

    operator = sdf.utils.spectral_resample_matrix(
        wave_in, wave_out, resolving_power=100.0,
        native_resolving_power=1000.0,
    )

    assert operator.shape == (len(wave_out), len(wave_in))
    assert np.allclose(operator @ np.ones(len(wave_in)), 1.0)


def test_spectral_resample_matrix_refuses_to_sharpen_model():
    with pytest.raises(sdf.utils.SdfError, match='must exceed'):
        sdf.utils.spectral_resample_matrix(
            np.geomspace(5.0, 15.0, 100),
            np.geomspace(6.0, 14.0, 50),
            resolving_power=600.0,
            native_resolving_power=100.0,
        )
