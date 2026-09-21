from .context import sdf
from .test_coordinator import FakeSampler, _fit_dataset
from .test_dataset import _write_rawphot
from urllib.parse import parse_qs, urlparse
import sdf.coordinator
import sdf.joint_plotting
import sdf.model_selection
import sdf.sed_page


def test_results_page_loads_ordinary_single_file_package(tmp_path):
    package = tmp_path / 'single'
    package.mkdir()
    _write_rawphot(
        package / 'target-rawphot.txt',
        sdbid='sdb-single', main_id='HD 123',
    )
    dataset = sdf.dataset.FitDataset.from_path(package)
    candidate = sdf.model_selection.ModelCandidate.for_single_target(
        dataset.targets[0].key, ('phoenix_sol',)
    )
    sdf.coordinator.FitCoordinator(
        package, sampler=FakeSampler()
    ).fit_candidate(dataset, candidate)

    output = sdf.joint_plotting.write_joint_results_page(package)

    assert output == package / 'index.html'
    assert 'HD 123' in output.read_text(encoding='utf-8')


def test_joint_results_page_contains_candidate_and_scoped_plot(tmp_path):
    dataset = _fit_dataset(tmp_path)
    candidate = sdf.model_selection.ModelCandidate.for_single_target(
        dataset.targets[0].key, ('kurucz-0.0',)
    )
    coordinator = sdf.coordinator.FitCoordinator(
        tmp_path, sampler=FakeSampler()
    )
    coordinator.fit_candidate(dataset, candidate)

    output = sdf.joint_plotting.write_joint_results_page(
        tmp_path, dataset=dataset
    )
    html = output.read_text(encoding='utf-8')

    assert output == tmp_path / 'index.html'
    assert 'Bokeh' in html
    assert 'aria-label="Default fit"' in html
    assert '>★</span>' in html
    assert 'kurucz-0.0' in html
    assert 'ln Z 0.00' in html
    assert 'log(g)' in html
    assert 'Residual' in html
    assert 'current input' in html
    assert '<h2>System</h2>' in html
    assert '>Results</a>' not in html
    assert '>Data</a>' not in html
    assert 'result.json' in html
    assert 'posterior-' in html
    assert 'fit-controls' not in html
    assert 'aria-label="Search SDF catalog"' in html
    assert 'aria-label="Footer navigation"' in html
    assert '>sdf home</a>' in html
    assert 'toolbar_location":"below"' in html


def test_operator_results_page_extends_public_page_with_fit_drawer(tmp_path):
    dataset = _fit_dataset(tmp_path)
    candidate = sdf.model_selection.ModelCandidate.for_single_target(
        dataset.targets[0].key, ('kurucz-0.0',)
    )
    sdf.coordinator.FitCoordinator(
        tmp_path, sampler=FakeSampler()
    ).fit_candidate(dataset, candidate)

    output = sdf.joint_plotting.write_joint_results_page(
        tmp_path, dataset=dataset, operator=True,
    )
    html = output.read_text(encoding='utf-8')

    assert output == tmp_path / 'operator.html'
    assert 'class="sed-page operator-mode"' in html
    assert 'id="fit-controls"' in html
    assert 'id="open-fit-controls"' in html
    assert 'data-target-key="simbad:Test"' in html
    assert 'data-fit-reference=' in html


def test_scope_colours_reuse_owner_colour_only_for_resolved_observations():
    scopes = (
        {'observation_id': 'A', 'contributors': ('a',)},
        {'observation_id': 'AB', 'contributors': ('a', 'b')},
        {'observation_id': 'B', 'contributors': ('b',)},
    )

    colours = sdf.joint_plotting._scope_colours(
        scopes, {'a': 'blue', 'b': 'orange'},
    )

    assert colours['A'] == 'blue'
    assert colours['B'] == 'orange'
    assert colours['AB'] not in {'blue', 'orange'}


def test_system_projection_builds_compact_component_hierarchy():
    targets = [
        {
            'key': 'a', 'label': 'HD 123A', 'colour': 'blue',
            'spectral_type': 'F5V', 'distance_pc': 10.0,
            'input_file': 'a.txt', 'external_links': [],
        },
        {
            'key': 'b', 'label': 'HD 123B', 'colour': 'orange',
            'spectral_type': 'K0V', 'distance_pc': 10.0,
            'input_file': 'b.txt', 'external_links': [],
        },
    ]
    observations = [{
        'label': 'HD 123', 'contributor_keys': ['a', 'b'],
        'spectral_type': None, 'distance_pc': 10.0,
        'input_file': 'ab.txt', 'external_links': [],
    }]

    systems, labels = sdf.sed_page.system_page_rows(
        targets, observations
    )

    assert len(systems) == 1
    assert systems[0]['component'] == 'AB'
    assert [child['component'] for child in systems[0]['children']] == [
        'A', 'B'
    ]
    assert labels == {'a': 'A', 'b': 'B'}


def test_external_search_radii_expand_for_wide_systems():
    radius = sdf.sed_page._search_radius_arcsec(
        10.0, -20.0, ((10.0, -20.0), (10.1, -20.0))
    )
    links = {
        link['label']: link['url']
        for link in sdf.sed_page._external_links(
            'HD 123', 10.0, -20.0, radius_arcsec=radius
        )
    }

    cassis = parse_qs(urlparse(links['CASSIS']).query)
    spitzer = parse_qs(urlparse(links['Spitzer']).query)
    esasky = parse_qs(urlparse(links['ESASky']).query)
    assert float(cassis['radius'][0]) > 340.0
    assert float(spitzer['sr'][0][:-1]) > 340.0
    assert float(esasky['fov'][0]) > radius * 2 / 3600.0
