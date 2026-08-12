import shutil

from .context import sdf


def test_read_phoenix(tmp_path):
    model_root = sdf.config.file['model_root']
    for name in ('lte009.5-3.5-0.0.BT-Settl.7.bz2',
                 'lte058-4.5-0.0a+0.0.BT-Settl.7.bz2'):
        model_file = tmp_path / name
        shutil.copyfile(model_root + name, model_file)

        spectrum = sdf.spectrum.ModelSpectrum.read_phoenix(str(model_file))

        assert spectrum.name == 'phoenix'
        assert (tmp_path / (name + '.npy')).is_file()
