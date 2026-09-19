"""Portable, hand-authorable spectral inputs for an SDF fit package."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path

import yaml

from . import utils


SPECTRA_MANIFEST_VERSION = 1


@dataclass(frozen=True)
class SegmentOverride:
    """Non-default contributor or eligibility state for one segment."""

    contributors: tuple[str, ...] | None = None
    include: bool | None = None


@dataclass(frozen=True)
class SpectrumDefinition:
    """One provider-native spectrum referenced by a fit package."""

    spectrum_id: str
    path: Path
    format: str
    contributors: tuple[str, ...]
    include: bool = True
    sha256: str | None = None
    segment_overrides: dict[str, SegmentOverride] | None = None


@dataclass(frozen=True)
class SpectraManifest:
    """Validated contents of one ``spectra.yml`` file."""

    path: Path
    spectra: tuple[SpectrumDefinition, ...]
    version: int = SPECTRA_MANIFEST_VERSION


def read_spectra_manifest(path: str | Path) -> SpectraManifest:
    """Read and validate a portable spectral-input manifest."""

    path = Path(path).expanduser().resolve()
    try:
        payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as error:
        raise utils.SdfError(
            "could not read {}: {}".format(path, error)
        ) from error
    if not isinstance(payload, dict):
        raise utils.SdfError("spectra.yml must contain a mapping")
    unknown = set(payload) - {"version", "spectra"}
    if unknown:
        raise utils.SdfError(
            "spectra.yml contains unknown keys: {}".format(
                ", ".join(sorted(unknown))
            )
        )
    if payload.get("version") != SPECTRA_MANIFEST_VERSION:
        raise utils.SdfError(
            "unsupported spectra.yml version {}".format(
                payload.get("version")
            )
        )
    rows = payload.get("spectra")
    if not isinstance(rows, list):
        raise utils.SdfError("spectra.yml spectra must be a list")

    definitions = []
    seen_ids = set()
    for index, row in enumerate(rows, start=1):
        if not isinstance(row, dict):
            raise utils.SdfError(
                "spectra.yml spectrum {} must be a mapping".format(index)
            )
        unknown = set(row) - {
            "id", "file", "format", "contributors", "include", "sha256",
            "segment_overrides",
        }
        if unknown:
            raise utils.SdfError(
                "spectrum {} contains unknown keys: {}".format(
                    index, ", ".join(sorted(unknown))
                )
            )
        spectrum_id = _nonempty_text(row.get("id"), "spectrum id")
        if spectrum_id in seen_ids:
            raise utils.SdfError(
                "duplicate spectrum id {}".format(spectrum_id)
            )
        seen_ids.add(spectrum_id)
        relative = Path(_nonempty_text(row.get("file"), "spectrum file"))
        if relative.is_absolute() or ".." in relative.parts:
            raise utils.SdfError(
                "spectrum {} file must be a package-relative path".format(
                    spectrum_id
                )
            )
        asset = (path.parent / relative).resolve()
        try:
            asset.relative_to(path.parent)
        except ValueError as error:
            raise utils.SdfError(
                "spectrum {} file escapes the package".format(spectrum_id)
            ) from error
        if not asset.is_file():
            raise utils.SdfError(
                "spectrum {} file not found: {}".format(
                    spectrum_id, relative
                )
            )
        format_ = _nonempty_text(row.get("format"), "spectrum format")
        contributors = _contributors(
            row.get("contributors"),
            "contributors for {}".format(spectrum_id),
        )
        include = row.get("include", True)
        if not isinstance(include, bool):
            raise utils.SdfError(
                "include for {} must be true or false".format(spectrum_id)
            )
        digest = row.get("sha256")
        if digest is not None:
            digest = _nonempty_text(digest, "sha256")
            if len(digest) != 64 or any(
                    value not in "0123456789abcdefABCDEF" for value in digest):
                raise utils.SdfError(
                    "sha256 for {} must contain 64 hexadecimal characters".
                    format(spectrum_id)
                )
            actual = _sha256(asset)
            if actual.lower() != digest.lower():
                raise utils.SdfError(
                    "spectrum {} checksum mismatch".format(spectrum_id)
                )
            digest = digest.lower()
        overrides = _segment_overrides(
            row.get("segment_overrides", {}), spectrum_id
        )
        definitions.append(SpectrumDefinition(
            spectrum_id=spectrum_id,
            path=asset,
            format=format_,
            contributors=contributors,
            include=include,
            sha256=digest,
            segment_overrides=overrides,
        ))
    return SpectraManifest(path=path, spectra=tuple(definitions))


def _segment_overrides(value, spectrum_id):
    if not isinstance(value, dict):
        raise utils.SdfError(
            "segment_overrides for {} must be a mapping".format(spectrum_id)
        )
    overrides = {}
    for segment, row in value.items():
        segment = _nonempty_text(segment, "segment name")
        if not isinstance(row, dict):
            raise utils.SdfError(
                "segment override {}:{} must be a mapping".format(
                    spectrum_id, segment
                )
            )
        unknown = set(row) - {"contributors", "include"}
        if unknown:
            raise utils.SdfError(
                "segment override {}:{} contains unknown keys: {}".format(
                    spectrum_id, segment, ", ".join(sorted(unknown))
                )
            )
        contributors = (
            None if "contributors" not in row else
            _contributors(
                row["contributors"],
                "contributors for {}:{}".format(spectrum_id, segment),
            )
        )
        include = row.get("include")
        if include is not None and not isinstance(include, bool):
            raise utils.SdfError(
                "include for {}:{} must be true or false".format(
                    spectrum_id, segment
                )
            )
        overrides[segment] = SegmentOverride(
            contributors=contributors, include=include,
        )
    return overrides


def _contributors(value, label):
    if (
        not isinstance(value, list) or not value or
        any(not isinstance(item, str) or not item.strip() for item in value)
    ):
        raise utils.SdfError(
            "{} must be a non-empty list of strings".format(label)
        )
    normalized = tuple(item.strip() for item in value)
    if len(set(normalized)) != len(normalized):
        raise utils.SdfError("{} repeats a contributor".format(label))
    return normalized


def _nonempty_text(value, label):
    if not isinstance(value, str) or not value.strip():
        raise utils.SdfError("{} must be a non-empty string".format(label))
    return value.strip()


def _sha256(path):
    digest = sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
