from .context import sdf
from .test_coordinator import FakeSampler, _fit_dataset
from .test_dataset import _write_rawphot
from urllib.parse import parse_qs, urlparse
from types import SimpleNamespace
import numpy as np
from bokeh.plotting import figure
import sdf.analytics
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


def test_results_page_links_exported_observation_products(tmp_path):
    package = tmp_path / 'single'
    package.mkdir()
    _write_rawphot(
        package / 'target-rawphot.txt',
        sdbid='sdb-single', main_id='HD 123',
    )
    (package / 'resources.yml').write_text(
        'version: 1\n'
        'resources:\n'
        '  - provider: tldr_pacs\n'
        '    label: Herschel/PACS images\n'
        '    url: https://example.test/herschel/sources/source.html\n'
        '    observation_ids: [1342223848]\n'
        '    observation_sdbids: [sdb-single]\n'
        '    target_sdbids: [sdb-single]\n',
        encoding='utf-8',
    )
    dataset = sdf.dataset.FitDataset.from_path(package)
    candidate = sdf.model_selection.ModelCandidate.for_single_target(
        dataset.targets[0].key, ('phoenix_sol',)
    )
    sdf.coordinator.FitCoordinator(
        package, sampler=FakeSampler()
    ).fit_candidate(dataset, candidate)

    output = sdf.joint_plotting.write_joint_results_page(package)
    html = output.read_text(encoding='utf-8')

    assert 'href="https://example.test/herschel/sources/source.html"' in html
    assert '>Herschel/PACS images</a>' in html


def test_joint_results_page_contains_candidate_and_scoped_plot(tmp_path):
    dataset = _fit_dataset(tmp_path)
    candidate = sdf.model_selection.ModelCandidate.for_single_target(
        dataset.targets[0].key, ('kurucz-0.0',)
    )
    coordinator = sdf.coordinator.FitCoordinator(
        tmp_path, sampler=FakeSampler()
    )
    coordinator.fit_candidate(dataset, candidate)
    stored = coordinator.library.find(candidate.canonical_id)
    (stored.path / 'corner.png').write_bytes(b'plot')
    (stored.path / 'distributions.png').write_bytes(b'plot')

    output = sdf.joint_plotting.write_joint_results_page(
        tmp_path, dataset=dataset,
        navigation={
            'admin_url': '/admin/system/test/',
            'review_url': '/review/target/test',
            'admin_access': True,
        },
    )
    html = output.read_text(encoding='utf-8')

    assert output == tmp_path / 'index.html'
    assert 'Bokeh' in html
    assert 'aria-label="Default fit"' in html
    assert '>★</span>' in html
    assert 'kurucz-0.0' in html
    assert 'ln Z 0.00' in html
    assert 'log(g)' not in html
    assert 'Residual' in html
    assert 'current input' not in html
    assert 'fit input has changed' not in html
    assert html.index('Fit results') < html.index('metadata (JSON)')
    assert html.index('metadata (JSON)') < html.index('ln Z 0.00')
    assert '<h2>System</h2>' in html
    assert '>Results</a>' not in html
    assert '>Data</a>' not in html
    assert 'result.json' in html
    assert 'posterior-' in html
    assert 'parameters (PNG)' in html
    assert 'derived (PNG)' in html
    assert r'metadata (JSON)&lt;/a&gt; \\u00b7 &lt;a' in html
    assert 'data-distribution-plot' in html
    assert 'id="distribution-dialog"' in html
    assert '>Download PNG</a>' in html
    assert 'fit-controls' not in html
    assert 'href="/admin/system/test/">Admin SED</a>' in html
    assert 'href="/review/target/test">Review data</a>' in html
    assert 'aria-label="Search SEDs"' in html
    assert 'aria-label="Site navigation"' in html
    assert '>Search</a>' not in html
    assert 'aria-label="Footer navigation"' in html
    assert '>home</a>' in html
    assert 'href="/docs">docs</a>' in html
    assert '<title>Test</title>' in html
    assert 'toolbar_location":"below"' in html
    assert '"title":"Limits"' in html
    assert 'excess range' not in html
    assert '.bk-header .bk-tab:last-child { margin-left: auto; }' in html

    with open(dataset.observations[0].source, 'a', encoding='utf-8') as handle:
        handle.write('\n\\ input changed after fitting\n')
    stale_html = sdf.joint_plotting.write_joint_results_page(
        tmp_path, dataset=dataset,
    ).read_text(encoding='utf-8')
    assert 'fit input has changed' in stale_html


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
    assert 'action="/admin/search"' in html
    assert 'href="/admin/jobs"' in html


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


def test_display_residuals_include_ignored_data_and_limits():
    arrays = {
        'residual': np.array([1.0, 0.0, -5.0, 0.0]),
        'observed_flux': np.array([12.0, 12.0, 12.0, 12.0]),
        'observed_uncertainty': np.array([2.0, 2.0, 0.0, 0.0]),
        'prediction_median': np.array([10.0, 6.0, 20.0, 6.0]),
        'observed_upper_limit': np.array([False, False, True, True]),
        'observed_ignore': np.array([False, True, False, True]),
    }

    residuals = sdf.joint_plotting._display_residual_values(arrays)

    assert np.allclose(residuals, [1.0, 3.0, -5.0, -1.5])


def test_residual_bounds_keep_minimum_and_do_not_clip_data():
    assert sdf.joint_plotting._residual_bounds([np.nan, 1.0]) == (-4.0, 4.0)

    low, high = sdf.joint_plotting._residual_bounds([-12.0, 7.0])

    assert low < -12.0
    assert high > 7.0


def test_detection_limits_match_noise_and_catalogue_upper_limits():
    arrays = {
        'wavelength_micron': np.array([12.0, 60.0, 100.0]),
        'filter': np.array(['WISE12', 'IRAS60', 'IRAS100']),
        'observed_flux': np.array([1.0, 0.5, 0.2]),
        'observed_uncertainty': np.array([0.01, 0.05, 0.02]),
        'observed_upper_limit': np.array([False, True, False]),
        'observed_ignore': np.array([False, False, True]),
        'component_prediction_median': np.array([[0.5, 0.1, 0.1]]),
        'component_prediction_p16': np.array([[0.49, 0.09, 0.09]]),
        'component_prediction_p84': np.array([[0.51, 0.11, 0.11]]),
    }
    target = {
        'parallax_arcsec': {'median': 0.1},
        'stellar_components': [{
            'luminosity_lsun': {'median': 1.0},
            'luminosity_at_1pc_lsun': {'median': 0.01},
        }],
    }
    candidate = SimpleNamespace(
        analysis_arrays=arrays,
        metadata={
            'components': [{
                'owner_key': 'simbad:A', 'role': 'stellar',
            }],
            'analysis': {
                'components': [{
                    'owner_key': 'simbad:A',
                    'parameters': {'Teff': {'median': 5800.0}},
                }],
                'targets': {'simbad:A': target},
                'data_scopes': [{
                    'kind': 'Photometry',
                    'contributors': ['simbad:A'],
                    'element_start': 0,
                    'element_stop': 3,
                }],
                'derived_parameters': [],
            },
        },
    )

    limits = sdf.joint_plotting._detection_limits(candidate, 'simbad:A')

    assert limits['filters'].tolist() == ['WISE12', 'IRAS60']
    assert limits['basis'].tolist() == [
        '3σ uncertainty', 'upper limit',
    ]
    assert limits['curves'].shape == (1000, 2)
    assert np.all(limits['envelope'] == np.min(limits['curves'], axis=1))
    expected_flux_limits = np.array([
        3.0 * np.sqrt(0.01**2 + 0.01**2),
        0.5,
    ])
    expected = sdf.analytics.BB_Disk(
        lstar=1.0, distance=10.0,
    ).f_limits(
        np.array([12.0, 60.0]), flux_limits=expected_flux_limits,
    )
    assert np.allclose(limits['curves'], expected)
    trial_fractional_luminosity = limits['envelope'] * 1.01
    assert np.all(np.any(
        trial_fractional_luminosity[:, None] >= limits['curves'], axis=1,
    ))
    assert limits['radius_au'] is not None


def test_fit_summaries_include_scientific_luminosities_and_modbb_shape():
    stellar = sdf.joint_plotting._stellar_summary(
        {'parameters': {
            'Teff': {'median': 6000.0, 'minus': 50.0, 'plus': 50.0},
            'logg': {'median': 4.5, 'minus': 0.1, 'plus': 0.1},
        }},
        {
            'radius_rsun': {'median': 1.1, 'minus': 0.1, 'plus': 0.1},
            'luminosity_lsun': {
                'median': 2.5, 'minus': 0.2, 'plus': 0.4,
            },
        },
    )
    disk = sdf.joint_plotting._disk_summary(
        {'parameters': {
            'beta': {'median': 1.2, 'minus': 0.1, 'plus': 0.2},
            'log_lam0': {'median': 2.0, 'minus': 0.1, 'plus': 0.1},
        }},
        {
            'temperature_k': {
                'median': 80.0, 'minus': 2.0, 'plus': 4.0,
            },
            'fractional_luminosity': {
                'median': 1e-3, 'minus': 1e-4, 'plus': 3e-4,
            },
        },
    )

    assert all('log(g)' not in value for value in stellar)
    assert 'L<sub>★</sub> 2.50e+00 ± 3.00e-01 L<sub>☉</sub>' in stellar
    assert 'β 1.20 ± 0.15' in disk
    assert any(value.startswith('λ<sub>0</sub> 100.00 ± ') for value in disk)
    assert 'L<sub>d</sub>/L<sub>★</sub> 1.00e-03 ± 2.00e-04' in disk


def test_photometry_symbols_cover_state_and_match_in_both_panels():
    arrays = {
        'wavelength_micron': np.array([1.0, 2.0, 3.0, 4.0]),
        'observed_flux': np.array([10.0, 9.0, 8.0, 7.0]),
        'observed_uncertainty': np.ones(4),
        'prediction_median': np.array([10.0, 8.0, 6.0, 5.0]),
        'residual': np.array([0.0, 0.0, -2.25, 0.0]),
        'observed_upper_limit': np.array([False, False, True, True]),
        'observed_ignore': np.array([False, True, False, True]),
        'filter': np.array(['A', 'B', 'C', 'D']),
    }
    candidate = type('Candidate', (), {'analysis_arrays': arrays})()
    scope = {'element_start': 0, 'element_stop': 4, 'kind': 'Photometry'}
    sed = figure()
    residual = figure()

    sdf.joint_plotting._add_scope_data(
        sed, residual, candidate, scope, 'AB', '#123456'
    )
    sdf.joint_plotting._add_symbol_legend(sed)

    legend_items = {
        item.label.value: item.renderers[0]
        for item in sed.legend[0].items
    }
    assert set(legend_items) == {
        'data', 'upper limit', 'excluded', 'model',
    }
    data = legend_items['data'].glyph
    model = legend_items['model'].glyph
    assert model.size / 2 - model.line_width / 2 > data.size / 2
    assert model.size == 14
    model_source = legend_items['model'].data_source
    assert len(model_source.data['wavelength']) == 1
    assert np.isnan(model_source.data['wavelength'][0])

    residual_markers = [
        renderer.glyph.marker for renderer in residual.renderers
    ]
    assert residual_markers == [
        'circle', 'circle', 'inverted_triangle', 'inverted_triangle'
    ]
    assert residual.renderers[1].glyph.fill_color == \
        sdf.joint_plotting.EXCLUDED_COLOUR
    assert residual.renderers[3].glyph.fill_color == \
        sdf.joint_plotting.EXCLUDED_COLOUR


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
