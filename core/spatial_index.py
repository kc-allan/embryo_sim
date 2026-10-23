"""Deterministic uniform spatial hash for local cell queries."""

from collections import defaultdict
from math import floor


class SpatialIndex:
    """Index a fixed cell arrangement; rebuild after any position/list change.

    Buckets are square with the requested positive ``bucket_size``. Queries
    inspect every bucket touched by the query radius, then apply exact distance
    filtering in the caller. Results are stable indices into ``cells``.
    """

    def __init__(self, cells, bucket_size: float):
        if bucket_size <= 0:
            raise ValueError("bucket_size must be positive")
        self.cells = list(cells)
        self.bucket_size = float(bucket_size)
        self.buckets = defaultdict(list)
        for index, cell in enumerate(self.cells):
            self.buckets[self._key(cell.position)].append(index)

    def _key(self, position):
        return (floor(float(position[0]) / self.bucket_size),
                floor(float(position[1]) / self.bucket_size))

    def query(self, position, radius: float) -> list[int]:
        if radius < 0:
            raise ValueError("radius must be nonnegative")
        cx, cy = self._key(position)
        reach = int(radius // self.bucket_size) + 1
        result = []
        for bx in range(cx - reach, cx + reach + 1):
            for by in range(cy - reach, cy + reach + 1):
                result.extend(self.buckets.get((bx, by), ()))
        return sorted(result)

    def candidate_pairs(self, cutoff: float):
        """Yield each unordered candidate pair once in deterministic order."""
        pairs = set()
        for i, cell in enumerate(self.cells):
            for j in self.query(cell.position, cutoff):
                if j > i:
                    pairs.add((i, j))
        yield from sorted(pairs)
