"""Worker boundary for slow operator-requested fits."""

from __future__ import annotations

import time

from .coordinator import fit_package
from .fit_management import FitManager


def run_next_job(state, *, sampler=None, raise_errors=False):
    """Run one queued job, returning it or ``None`` when the queue is empty."""

    job = state.claim_next()
    if job is None:
        return None
    try:
        fit_package(
            job.package_path, sampler=sampler, policy_store=state,
        )
        FitManager(job.package_path).render()
    except Exception as error:
        state.finish(
            job.id,
            error='{}: {}'.format(type(error).__name__, error),
        )
        if raise_errors:
            raise
    else:
        state.finish(job.id)
    return job


def run_worker(state, *, poll_seconds=5.0):
    """Continuously drain the queue without coupling work to HTTP requests."""

    while True:
        if run_next_job(state) is None:
            time.sleep(poll_seconds)
