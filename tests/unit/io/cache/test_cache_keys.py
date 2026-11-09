# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from pathlib import Path

import pytest

from pysurveyor.io.cache.cache_keys import (
    FileCacheKey,
    _make_hashable,
    create_file_cache_key,
)


# ======================================================================================
# Test for hashing the cache key values
# ======================================================================================
class TestMakeHashable:
    """Test conversion of cache-key values into hashable representations."""

    # ==================================================================================
    @pytest.mark.parametrize(
        ("value", "expected"),
        [
            pytest.param(1, 1, id="scalar"),
            pytest.param(
                {"b": [2, 3], "a": {4, 5}},
                (
                    "mapping",
                    (
                        ("a", ("set", frozenset({4, 5}))),
                        ("b", ("sequence", (2, 3))),
                    ),
                ),
                id="nested-mapping",
            ),
            pytest.param(
                [1, (2, 3)],
                ("sequence", (1, ("sequence", (2, 3)))),
                id="sequence",
            ),
            pytest.param(
                {1, 2, 3},
                ("set", frozenset({1, 2, 3})),
                id="set",
            ),
        ],
    )
    def test_convert_values(self, value, expected):
        """Convert supported values recursively and preserve their categories."""
        result = _make_hashable(value)

        # ------------------------------------------------------------------------------

        assert result == expected
        assert hash(result) is not None

    # ==================================================================================
    def test_mapping_order(self):
        """Normalize mappings independently of insertion order."""
        first = _make_hashable({"first": 1, "second": 2})
        second = _make_hashable({"second": 2, "first": 1})

        # ------------------------------------------------------------------------------

        assert first == second

    # ==================================================================================
    def test_convert_callable(self):
        """Represent callables by their importable identity."""
        result = _make_hashable(Path.exists)

        # ------------------------------------------------------------------------------

        assert result == (
            "callable",
            Path.exists.__module__,
            Path.exists.__qualname__,
        )

    # ==================================================================================
    def test_unhashable_value(self):
        """Raise an error for unsupported unhashable values."""
        unhashable = object.__new__(type("Unhashable", (), {"__hash__": None}))

        with pytest.raises(TypeError, match=f"{unhashable!r}"):
            _make_hashable(unhashable)


# ======================================================================================
# Test for creating file cache keys
# ======================================================================================
class TestFileCacheKey:
    """Test construction of file cache keys."""

    # ==================================================================================
    def test_file_metadata_and_options(self, tmp_text_file):
        """Include the resolved path, normalized options, and file metadata."""
        filepath = tmp_text_file("data.csv", "value\n1\n")

        key = create_file_cache_key(
            filepath,
            {"nrows": 1, "skiprows": [0, 2]},
        )

        # ------------------------------------------------------------------------------

        assert isinstance(key, FileCacheKey)
        assert key.resolved_path == str(filepath.resolve())
        assert key.reader_options == (
            ("nrows", 1),
            ("skiprows", ("sequence", (0, 2))),
        )
        assert key.modified_time_ns == filepath.stat().st_mtime_ns
        assert key.file_size == filepath.stat().st_size

    # ==================================================================================
    def test_equivalent_paths_and_options(self, tmp_text_file):
        """Create equal keys for equivalent paths and option mappings."""
        filepath = tmp_text_file("data.csv", "value\n1\n")
        relative_path = filepath.parent / "." / filepath.name

        first = create_file_cache_key(filepath, {"delimiter": ",", "nrows": 1})
        second = create_file_cache_key(relative_path, {"nrows": 1, "delimiter": ","})

        # ------------------------------------------------------------------------------

        assert first == second
        assert hash(first) == hash(second)

    # ==================================================================================
    def test_key_changes(self, tmp_text_file):
        """Invalidate a key when the file size changes."""
        filepath = tmp_text_file("data.csv", "value\n1\n")
        first = create_file_cache_key(filepath, {})

        filepath.write_text("value\n100\n")
        second = create_file_cache_key(filepath, {})

        # ------------------------------------------------------------------------------

        assert first != second
        assert first.file_size != second.file_size

    # ==================================================================================
    def test_missing_file_raises(self, tmp_path):
        """Propagate the missing-file error from metadata lookup."""
        with pytest.raises(FileNotFoundError):
            create_file_cache_key(tmp_path / "missing.csv", {})
