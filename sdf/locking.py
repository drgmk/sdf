"""Small advisory locks for work shared by local processes."""

from __future__ import annotations

import fcntl
from pathlib import Path


class InterProcessLock:
    """Hold a non-blocking advisory lock until explicitly released."""

    def __init__(self, path, *, description='work'):
        self.path = Path(path)
        self.description = description
        self._handle = None

    def acquire(self):
        if self._handle is not None:
            return self
        self.path.parent.mkdir(parents=True, exist_ok=True)
        handle = self.path.open('a+', encoding='utf-8')
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            handle.close()
            raise RuntimeError(
                '{} is already active (lock: {})'.format(
                    self.description, self.path
                )
            ) from error
        self._handle = handle
        return self

    def release(self):
        if self._handle is None:
            return
        try:
            fcntl.flock(self._handle.fileno(), fcntl.LOCK_UN)
        finally:
            self._handle.close()
            self._handle = None

    def __enter__(self):
        return self.acquire()

    def __exit__(self, exc_type, exc_value, traceback):
        self.release()

