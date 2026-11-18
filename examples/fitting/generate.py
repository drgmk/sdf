"""Regenerate deterministic rawphot inputs for the fitting examples."""

from __future__ import annotations

import json
from hashlib import sha256
from pathlib import Path

from astropy.table import Table
import numpy as np
import yaml

from sdf.model import PhotModel, SpecModel


ROOT = Path(__file__).resolve().parent
MODEL_NAME = 'phoenix_sol'
FILTERS = (
    'UJ', 'BJ', 'VJ', 'RC', 'IC', 'GAIA.G', '2MJ', '2MH', '2MKS',
    'WISE12', 'WISE22', 'PACS70',
)
TEFF = 6000.0
LOGG = 4.5
BRIGHT_NORM = -2.6
FAINT_NORM = -5.6


def model_flux(parameters):
    model = PhotModel.read_model(MODEL_NAME)
    model.keep_filters(FILTERS)
    return np.asarray(model.fnujy(parameters), dtype=float)


def write_rawphot(path, *, identifier, main_id, flux):
    flux = np.asarray(flux, dtype=float)
    table = Table(
        rows=[
            (
                band, value, 0.02 * value, 0.01 * value, 0, 'Jy',
                'synthetic', 'exact model photometry', 0,
            )
            for band, value in zip(FILTERS, flux)
        ],
        names=(
            'Band', 'Phot', 'Err', 'Sys', 'Lim', 'Unit', 'bibcode',
            'Note1', 'exclude',
        ),
    )
    table.meta['keywords'] = {
        'id': {'value': identifier},
        'main_id': {'value': main_id},
        'sp_type': {'value': 'G2V'},
        'plx_value': {'value': 50.0},
        'plx_err': {'value': 0.1},
    }
    table.write(path, format='ascii.ipac', overwrite=True)
    # Astropy pads IPAC rows to the declared column widths. Keep the fixture
    # representation stable without committing line-ending whitespace.
    path.write_text(
        ''.join(line.rstrip() + '\n' for line in path.read_text().splitlines()),
        encoding='utf-8',
    )


def write_policy(path, names):
    payload = {
        'targets': {
            name: {
                'stellar_models': [MODEL_NAME],
                'disk_models': ['modbb_disk_r'],
                'max_disks': 0,
                'search_mode': 'evidence-guided',
                'reddening': False,
            }
            for name in names
        }
    }
    path.write_text(
        yaml.safe_dump(payload, sort_keys=False), encoding='utf-8'
    )


def write_truth(path, components):
    path.write_text(
        json.dumps({
            'model': MODEL_NAME,
            'filters': FILTERS,
            'components': components,
        }, indent=2) + '\n',
        encoding='utf-8',
    )


def write_spectrum(path, parameters):
    wavelength = np.geomspace(0.8, 2.5, 60)
    model = SpecModel.read_model(MODEL_NAME).copy()
    model.prepare_for_spectrum(wavelength)
    flux = np.asarray(model.fnujy(parameters), dtype=float)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        ''.join(
            '{:.17g},{:.17g},{:.17g}\n'.format(
                wave, value, 0.02 * value
            )
            for wave, value in zip(wavelength, flux)
        ),
        encoding='utf-8',
    )


def generate_single():
    directory = ROOT / 'single-star'
    directory.mkdir(parents=True, exist_ok=True)
    parameters = [TEFF, LOGG, BRIGHT_NORM]
    write_rawphot(
        directory / 'single-star-rawphot.txt',
        identifier='synthetic-single', main_id='Synthetic Single',
        flux=model_flux(parameters),
    )
    write_policy(directory / 'sdf-fit.yml', ('Synthetic Single',))
    write_truth(directory / 'truth.json', {
        'Synthetic Single': parameters,
    })


def generate_binary(name, prefix, norm_b):
    directory = ROOT / name
    directory.mkdir(parents=True, exist_ok=True)
    parameters_a = [TEFF, LOGG, BRIGHT_NORM]
    parameters_b = [TEFF, LOGG, norm_b]
    flux_a = model_flux(parameters_a)
    flux_b = model_flux(parameters_b)
    names = (
        'Synthetic {} A'.format(prefix),
        'Synthetic {} B'.format(prefix),
    )
    write_rawphot(
        directory / 'A-rawphot.txt', identifier='A', main_id=names[0],
        flux=flux_a,
    )
    write_rawphot(
        directory / 'B-rawphot.txt', identifier='B', main_id=names[1],
        flux=flux_b,
    )
    write_rawphot(
        directory / 'AB-rawphot.txt', identifier='AB',
        main_id='Synthetic {} AB'.format(prefix), flux=flux_a + flux_b,
    )
    (directory / 'joint-fit.yml').write_text(
        'version: 1\nobservations:\n  AB:\n    - A\n    - B\n',
        encoding='utf-8',
    )
    write_policy(directory / 'sdf-fit.yml', names)
    write_truth(directory / 'truth.json', {
        names[0]: parameters_a,
        names[1]: parameters_b,
    })


def generate_spectrum():
    directory = ROOT / 'spectrum-star'
    directory.mkdir(parents=True, exist_ok=True)
    parameters = [TEFF, LOGG, BRIGHT_NORM]
    write_rawphot(
        directory / 'spectrum-star-rawphot.txt',
        identifier='S', main_id='Synthetic Spectrum',
        flux=model_flux(parameters),
    )
    spectrum_path = directory / 'spectra' / 'spectrum.csv'
    write_spectrum(spectrum_path, parameters)
    digest = sha256(spectrum_path.read_bytes()).hexdigest()
    (directory / 'spectra.yml').write_text(
        'version: 1\n'
        'spectra:\n'
        '  - id: synthetic:spectrum\n'
        '    file: spectra/spectrum.csv\n'
        '    format: csv\n'
        '    contributors: [S]\n'
        '    sha256: {}\n'.format(digest),
        encoding='utf-8',
    )
    write_policy(directory / 'sdf-fit.yml', ('Synthetic Spectrum',))
    truth = {
        'model': MODEL_NAME,
        'filters': FILTERS,
        'components': {'Synthetic Spectrum': parameters},
        'spectrum_norms': [1.0],
    }
    (directory / 'truth.json').write_text(
        json.dumps(truth, indent=2) + '\n', encoding='utf-8'
    )


def main():
    generate_single()
    generate_binary('faint-binary', 'Faint', FAINT_NORM)
    generate_binary('equal-binary', 'Equal', BRIGHT_NORM)
    generate_spectrum()


if __name__ == '__main__':
    main()
