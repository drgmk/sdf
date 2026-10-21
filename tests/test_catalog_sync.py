import fcntl
from threading import Event, Lock, Thread

from .test_catalog_server import _catalog
from .test_dataset import _write_rawphot
import sdf.catalog_sync


def test_sync_takes_process_lock_before_export_lock(tmp_path):
    service = _catalog(tmp_path)
    synchronizer = sdf.catalog_sync.CatalogSynchronizer(service)
    gate = Lock()
    attempted = Event()
    finished = Event()
    gate.acquire()

    class ObservedLock:
        def __enter__(self):
            attempted.set()
            gate.acquire()

        def __exit__(self, exc_type, exc_value, traceback):
            gate.release()

    service.export_sync_lock = ObservedLock()

    def sync():
        synchronizer.sync()
        finished.set()

    thread = Thread(target=sync, daemon=True)
    thread.start()
    assert attempted.wait(1)
    try:
        with (tmp_path / '.sdb-export.lock').open('a+') as stream:
            fcntl.flock(
                stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB
            )
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
    finally:
        gate.release()
    thread.join(1)
    assert finished.is_set()


def test_sync_waits_while_sdb_export_lock_is_held(tmp_path):
    service = _catalog(tmp_path)
    synchronizer = sdf.catalog_sync.CatalogSynchronizer(service)
    rawphot = tmp_path / 'package' / 'target-rawphot.txt'
    rawphot.write_text(rawphot.read_text() + '\n', encoding='utf-8')

    with (tmp_path / '.sdb-export.lock').open('a+') as stream:
        fcntl.flock(stream.fileno(), fcntl.LOCK_EX)
        result = synchronizer.sync()

    assert result['export_busy'] is True
    assert result['changed'] == ()
    assert service.state.jobs() == ()


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
