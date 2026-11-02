from sdf.external_services import ExternalServices, simbad_identifier_route


def test_external_services_read_sdb_mirrors(tmp_path):
    config = tmp_path / 'sdb.toml'
    config.write_text(
        '[mirrors]\n'
        'simbad = "https://simbad.idia.ac.za/simbad"\n'
        'vizier = "https://vizier.idia.ac.za/viz-bin"\n',
        encoding='utf-8',
    )

    services = ExternalServices.from_sdb_config(config)

    assert services.simbad_identifier_url('HD 1326') == (
        'https://simbad.idia.ac.za/simbad/sim-id?'
        'submit=submit+id&Ident=HD+1326'
    )
    assert services.vizier_catalog_url('II/328/allwise') == (
        'https://vizier.idia.ac.za/viz-bin/VizieR?'
        '-source=II%2F328%2Fallwise'
    )
    assert services.vizier_row_url(
        'II/328/allwise', 'AllWISE', 'J025959.07+000020.3',
    ) == (
        'https://vizier.idia.ac.za/viz-bin/VizieR-5?'
        '-out.add=.&-source=II%2F328%2Fallwise'
        '&AllWISE===J025959.07%2B000020.3'
    )


def test_simbad_route_is_stable_and_relative():
    assert simbad_identifier_route('HD 1326') == (
        '/external/simbad/id?name=HD+1326'
    )


def test_external_services_use_first_ordered_vizier_mirror(tmp_path):
    config = tmp_path / 'sdb.toml'
    config.write_text(
        '[mirrors]\n'
        'vizier = ["vizier.idia.ac.za", "vizier.cds.unistra.fr"]\n',
        encoding='utf-8',
    )

    services = ExternalServices.from_sdb_config(config)

    assert services.vizier_server == 'vizier.idia.ac.za'


def test_external_services_honour_ordered_vizier_environment(monkeypatch):
    monkeypatch.setenv(
        'SDB_VIZIER_SERVERS',
        'vizier.idia.ac.za,vizier.cds.unistra.fr',
    )

    services = ExternalServices.from_sdb_config()

    assert services.vizier_server == 'vizier.idia.ac.za'
