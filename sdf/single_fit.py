"""Standalone package preparation for one-file fits."""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
import re
import shutil
from uuid import uuid4

from . import utils


@dataclass(frozen=True)
class SingleFitPackage:
    """The input and staged paths for a standalone fit."""

    input_path: Path
    directory: Path
    rawphot_path: Path


def prepare_single_fit_package(input_path, output=None):
    """Copy one IPAC table into a safe, self-contained fitting package.

    Existing results in the output package are retained. Re-running the
    command refreshes only the package's one rawphot input, which makes the
    retained candidates stale through their normal input fingerprint.
    """

    input_path = Path(input_path).expanduser().resolve()
    if not input_path.is_file():
        raise utils.SdfError(
            'single-fit input is not a file: {}'.format(input_path)
        )
    if _is_package_member(input_path):
        raise utils.SdfError(
            '{} belongs to an SDF package; fit that package with sdf-fits '
            'instead of detaching one member'.format(input_path)
        )

    stem = _input_stem(input_path.name)
    if output is None:
        directory = input_path.parent / '{}-sdf'.format(stem)
    else:
        directory = Path(output).expanduser().resolve()
    if directory.exists() and not directory.is_dir():
        raise utils.SdfError(
            'single-fit output is not a directory: {}'.format(directory)
        )
    directory.mkdir(parents=True, exist_ok=True)

    rawphot_path = directory / '{}-rawphot.txt'.format(stem)
    other_inputs = tuple(
        path for path in directory.glob('*-rawphot.txt')
        if path.resolve() != rawphot_path.resolve()
    )
    if other_inputs:
        raise utils.SdfError(
            'single-fit output already contains another rawphot input: {}'.
            format(', '.join(path.name for path in other_inputs))
        )

    if input_path != rawphot_path.resolve():
        temporary = directory / '.{}.{}.tmp'.format(
            rawphot_path.name, uuid4().hex
        )
        try:
            shutil.copy2(input_path, temporary)
            os.replace(temporary, rawphot_path)
        finally:
            if temporary.exists():
                temporary.unlink()

    return SingleFitPackage(
        input_path=input_path,
        directory=directory.resolve(),
        rawphot_path=rawphot_path.resolve(),
    )


def _input_stem(filename):
    suffix = '-rawphot.txt'
    value = filename[:-len(suffix)] if filename.endswith(suffix) \
        else Path(filename).stem
    value = re.sub(r'[^A-Za-z0-9_.-]+', '-', value).strip('.-')
    if not value:
        raise utils.SdfError(
            'could not derive an output name from {}'.format(filename)
        )
    return value


def _is_package_member(path):
    """Return whether a rawphot header explicitly belongs to a package."""

    with Path(path).open('r', encoding='utf-8') as handle:
        for line in handle:
            if not line.startswith('\\'):
                return False
            name, separator, _ = line[1:].partition('=')
            if separator and name.strip() in {
                'sdf_package', 'spectra_manifest'
            }:
                return True
    return False
