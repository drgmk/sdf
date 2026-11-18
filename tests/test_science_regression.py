import importlib.util
from pathlib import Path

import pytest


PATH = (
    Path(__file__).resolve().parents[1] /
    'examples' / 'fitting' / 'science_regression.py'
)
SPEC = importlib.util.spec_from_file_location('science_regression', PATH)
science_regression = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(science_regression)


def _summary(offset=0.0):
    measurement = {
        'median': 6000.0 + offset,
        'minus': 10.0,
        'plus': 12.0,
    }
    return {
        'schema': science_regression.SCHEMA,
        'kind': 'test',
        'model': 'phoenix_sol',
        'parameters': {'Teff': measurement},
        'stellar': {'Teff': measurement},
        'evidence': -10.0 + offset,
        'evidence_error': 0.2,
        'evidence_sampling_error': 2.0,
        'chi_squared': 4.0 + offset,
    }


def test_compare_summaries_normalizes_parameter_differences():
    comparison = science_regression.compare_summaries(
        _summary(), _summary(offset=11.0)
    )

    assert comparison['quantities'][0]['quantity'] == 'parameters.Teff'
    assert comparison['quantities'][0]['difference'] == 11.0
    assert comparison['quantities'][0]['z'] == pytest.approx(1 / 2**0.5)
    assert comparison['evidence']['difference'] == 11.0
    assert comparison['evidence']['z'] == pytest.approx(
        11.0 / (2.0 * 2**0.5)
    )
    assert comparison['chi_squared']['difference'] == 11.0


def test_compare_summaries_rejects_different_models():
    second = _summary()
    second['model'] = 'koester_wd'

    with pytest.raises(ValueError, match='cannot compare'):
        science_regression.compare_summaries(_summary(), second)


def test_parameter_only_comparison_omits_incomparable_fit_statistics():
    comparison = science_regression.compare_summaries(
        _summary(), _summary(), include_fit_statistics=False
    )

    assert 'evidence' not in comparison
    assert 'chi_squared' not in comparison
