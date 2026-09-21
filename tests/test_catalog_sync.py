from .test_catalog_server import _catalog
from .test_dataset import _write_rawphot
import sdf.catalog_sync


def test_sync_indexes_and_queues_new_exported_package(tmp_path):
    service = _catalog(tmp_path)
    synchronizer = sdf.catalog_sync.CatalogSynchronizer(service)
    package = tmp_path / 'new-package'
    package.mkdir()
    _write_rawphot(
        package / 'new-rawphot.txt', sdbid='sdb-new', main_id='HD 999'
    )

    result = synchronizer.sync()

    assert result['added'] == ('new-package',)
    assert len(result['queued']) == 1
    assert service.package('new-package')['status'] == 'unfit'
    assert service.state.jobs()[0].actor == 'sdf watcher'
    assert service.state.last_actor() is None


def test_sync_refreshes_changed_input_without_duplicate_active_job(tmp_path):
    service = _catalog(tmp_path)
    synchronizer = sdf.catalog_sync.CatalogSynchronizer(service)
    service.queue_package('package', actor='Grant')
    rawphot = tmp_path / 'package' / 'target-rawphot.txt'
    rawphot.write_text(rawphot.read_text() + '\n', encoding='utf-8')

    result = synchronizer.sync()

    assert result['changed'] == ('package',)
    assert result['queued'] == ()
    assert len(service.state.jobs()) == 1


def test_input_signature_tracks_spectral_manifest_and_asset(tmp_path):
    package = tmp_path / 'package'
    spectra = package / 'spectra'
    spectra.mkdir(parents=True)
    _write_rawphot(
        package / 'target-rawphot.txt', sdbid='sdb-test', main_id='HD 123',
    )
    spectrum = spectra / 'irs.csv'
    spectrum.write_text('wavelength,flux,error\n8,1,0.1\n', encoding='utf-8')
    manifest = package / 'spectra.yml'
    manifest.write_text('version: 1\nspectra: []\n', encoding='utf-8')

    first = sdf.catalog_sync._input_signature(package)
    spectrum.write_text(
        spectrum.read_text(encoding='utf-8') + '9,1.1,0.1\n',
        encoding='utf-8',
    )
    second = sdf.catalog_sync._input_signature(package)
    manifest.write_text(
        manifest.read_text(encoding='utf-8') + '# changed\n',
        encoding='utf-8',
    )
    third = sdf.catalog_sync._input_signature(package)

    assert second != first
    assert third != second
