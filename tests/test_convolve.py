from .context import sdf

def test_convolve_model_wave_range_large_enuf():
    fs = sdf.filter.Filter.all
    temps = [1,1e4]
    for f in fs:
        cm = sdf.convolve.ConvolvedModel.bb(f,temps)
