"""Simple context-manager timer."""

import time


class Timer:
    def __init__(self, verbose=False):
        self.verbose = verbose
        self.start = 0.0
        self.end = 0.0
        self.secs = 0.0
        self.msecs = 0.0

    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args):
        self.end = time.time()
        self.secs = self.end - self.start
        self.msecs = self.secs * 1000
        if self.verbose:
            print(f"elapsed time: {self.msecs} ms")
