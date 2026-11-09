# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

import pandas as pd
import pytest

from pysurveyor.io.cache.cache_keys import FileCacheKey
from pysurveyor.io.cache.loader_cache import LoaderCache

"""Test loading, invalidation, and scoping of cached data."""


# ======================================================================================
# Automaically reset the shared LoaderCache state before each test
# ======================================================================================
@pytest.fixture(autouse=True)
def reset_loader_cache():
    """Reset LoaderCache state before and after each test."""
    original_enabled = LoaderCache.enabled

    LoaderCache.clear()
    LoaderCache.enabled = True

    yield

    LoaderCache.clear()
    LoaderCache.enabled = original_enabled


# ======================================================================================
# Tests for the LoaderCache class
# ======================================================================================
def test_shared_cache_instance():
    """Expose the shared cache instance via class-level access."""

    instance = LoaderCache.shared()

    # ----------------------------------------------------------------------------------

    assert instance is LoaderCache.shared()


# ======================================================================================
def test_explicit_instance_is_separate_from_shared_cache():
    """Keep explicitly created caches separate from class-level operations."""

    separate_cache = LoaderCache(enabled=True)
    calls = 0

    # ----------------------------------------------------------------------------------
    # Define a loader function that increments a counter on each call
    # ----------------------------------------------------------------------------------
    def loader():
        nonlocal calls
        calls += 1
        return pd.DataFrame({"value": [calls]})

    # ----------------------------------------------------------------------------------
    # Load data using the shared cache and the separate cache
    # ----------------------------------------------------------------------------------
    shared_data = LoaderCache.get_data("key", loader)
    separate_data = separate_cache.get_data("key", loader)

    # ----------------------------------------------------------------------------------

    assert calls == 2
    assert shared_data.iloc[0, 0] == 1
    assert separate_data.iloc[0, 0] == 2


# ======================================================================================
def test_disabled_cache():
    """Do not retain data when caching is disabled."""

    LoaderCache.enabled = False

    calls = 0

    # ----------------------------------------------------------------------------------
    # Define a loader function that increments a counter on each call
    # ----------------------------------------------------------------------------------
    def loader():
        nonlocal calls
        calls += 1
        return pd.DataFrame({"value": [calls]})

    # ----------------------------------------------------------------------------------

    first = LoaderCache.get_data("key", loader)
    second = LoaderCache.get_data("key", loader)

    # ----------------------------------------------------------------------------------

    assert calls == 2
    assert first.iloc[0, 0] == 1
    assert second.iloc[0, 0] == 2


# ======================================================================================
def test_enabled_cache_and_copy():
    """Cache one result and protect it from mutations to returned copies."""

    calls = 0

    # ----------------------------------------------------------------------------------
    # Define a loader function that increments a counter on each call
    # ----------------------------------------------------------------------------------
    def loader():
        nonlocal calls
        calls += 1
        return pd.DataFrame({"value": [1]})

    first = LoaderCache.get_data("key", loader)
    second = LoaderCache.get_data("key", loader)

    # ----------------------------------------------------------------------------------
    # Modify the first retained data
    # ----------------------------------------------------------------------------------
    first.loc[0, "value"] = 2

    # ----------------------------------------------------------------------------------

    assert calls == 1
    assert second.loc[0, "value"] == 1
    assert second is not first


# ======================================================================================
def test_enabled_cache_and_disabled_copy():
    """Return the cached object itself when copy is disabled."""

    data = pd.DataFrame({"value": [1]})

    result = LoaderCache.get_data("key", lambda: data, copy=False)

    # ----------------------------------------------------------------------------------

    assert result is data


# ======================================================================================
def test_clear_and_remove():
    """Clear all entries or remove one selected entry."""

    calls = {"first": 0, "second": 0}

    # ----------------------------------------------------------------------------------
    # Define loader functions that increment a counter on each call
    # ----------------------------------------------------------------------------------
    def load_first():
        calls["first"] += 1
        return pd.DataFrame({"value": [calls["first"]]})

    def load_second():
        calls["second"] += 1
        return pd.DataFrame({"value": [calls["second"]]})

    # ----------------------------------------------------------------------------------

    LoaderCache.get_data("first", load_first)
    LoaderCache.get_data("second", load_second)

    # Both loaders should have been called once
    assert calls == {"first": 1, "second": 1}

    # ----------------------------------------------------------------------------------

    LoaderCache.remove("first")
    LoaderCache.remove("missing")
    LoaderCache.get_data("first", load_first)
    LoaderCache.get_data("second", load_second)

    # ----------------------------------------------------------------------------------
    # Only the first loader should have been called again,
    # the second should not as it was still cached
    # ----------------------------------------------------------------------------------
    assert calls == {"first": 2, "second": 1}

    # ----------------------------------------------------------------------------------

    LoaderCache.clear()
    LoaderCache.get_data("first", load_first)
    LoaderCache.get_data("second", load_second)

    # ----------------------------------------------------------------------------------
    # Both loaders should have been called now again, as the cache was cleared
    # ----------------------------------------------------------------------------------
    assert calls == {"first": 3, "second": 2}


# ======================================================================================
def test_remove_file(tmp_text_file):
    """Remove all entries for one file while retaining entries for other files."""
    """Remove all entries for one file while retaining entries for other files."""

    first_path = tmp_text_file("first.csv", "first").resolve()
    second_path = tmp_text_file("second.csv", "second").resolve()

    first_key = FileCacheKey(str(first_path), (), 1, 1)
    second_key = FileCacheKey(str(second_path), (), 1, 1)
    third_key = FileCacheKey(str(first_path), (("delimiter", ";"),), 1, 1)

    calls = {"first": 0, "second": 0, "third": 0}

    # ----------------------------------------------------------------------------------
    # Define loader functions that increment a counter on each call
    # ----------------------------------------------------------------------------------
    def load_first():
        calls["first"] += 1
        return pd.DataFrame({"value": [calls["first"]]})

    def load_second():
        calls["second"] += 1
        return pd.DataFrame({"value": [calls["second"]]})

    def load_third():
        calls["third"] += 1
        return pd.DataFrame({"value": [calls["third"]]})

    LoaderCache.get_data(first_key, load_first)
    LoaderCache.get_data(second_key, load_second)
    LoaderCache.get_data(third_key, load_third)

    assert calls == {"first": 1, "second": 1, "third": 1}

    # ----------------------------------------------------------------------------------
    # Remove all entries for the first file
    # ----------------------------------------------------------------------------------
    LoaderCache.remove_file(first_path)

    LoaderCache.get_data(first_key, load_first)
    LoaderCache.get_data(second_key, load_second)
    LoaderCache.get_data(third_key, load_third)

    # ----------------------------------------------------------------------------------
    # Both entries for the first file should have been removed, while the second file
    # should remain cached so the second loader is not called
    # ----------------------------------------------------------------------------------

    assert calls == {"first": 2, "second": 1, "third": 2}


# ======================================================================================
def test_enabled_property_uses_instance_override(monkeypatch):
    """Prefer an instance setting over the class-wide switch."""

    monkeypatch.setattr(LoaderCache, "enabled", True)

    # ----------------------------------------------------------------------------------

    assert LoaderCache().is_enabled
    assert not LoaderCache(enabled=False).is_enabled


# ======================================================================================
def test_shared_and_current_cache():
    """Use the shared cache unless a scoped cache is active."""

    LoaderCache.enabled = False

    first = LoaderCache.shared()
    second = LoaderCache.current()

    # ----------------------------------------------------------------------------------

    assert first is second
    assert not LoaderCache.enabled
    assert not first.is_enabled
    assert not second.is_enabled

    # ----------------------------------------------------------------------------------

    with LoaderCache.session() as session:
        assert LoaderCache.current() is session
        assert session.is_enabled

    # ----------------------------------------------------------------------------------

    assert LoaderCache.current() is first
    assert not first.is_enabled


# ======================================================================================
def test_session_clears_and_restores_cache():
    """Clear session data and restore the previous active cache on exit."""

    # ----------------------------------------------------------------------------------
    # Cache one entry in the shared cache to verify that it is retained
    # ----------------------------------------------------------------------------------

    data = LoaderCache.get_data("key", lambda: pd.DataFrame({"value": [1]}))

    assert len(LoaderCache.current()._cached_data) == 1
    assert data.iloc[0, 0] == 1

    # ----------------------------------------------------------------------------------
    # Create a new session and cache one entry in it
    # ----------------------------------------------------------------------------------
    with LoaderCache.session() as session:
        data = session.get_data("key", lambda: pd.DataFrame({"value": [2]}))

        assert len(session._cached_data) == 1
        assert data.iloc[0, 0] == 2

    # ----------------------------------------------------------------------------------
    # Create another session and cache one entry in it
    # ----------------------------------------------------------------------------------
    with LoaderCache.session() as another_session:
        data = another_session.get_data("key", lambda: pd.DataFrame({"value": [3]}))

        assert len(another_session._cached_data) == 1
        assert data.iloc[0, 0] == 3

    # ----------------------------------------------------------------------------------
    # Ensure that the previous cache is restored and the loader is not called again
    # ----------------------------------------------------------------------------------
    data = LoaderCache.get_data("key", lambda: pd.DataFrame({"value": [4]}))

    assert len(LoaderCache.current()._cached_data) == 1
    assert data.iloc[0, 0] == 1


# ======================================================================================
def test_session_restores_context():
    """Restore the previous cache even when a session body raises."""

    previous = LoaderCache.current()

    # ----------------------------------------------------------------------------------

    with pytest.raises(ValueError, match="failure"):
        # ------------------------------------------------------------------------------
        # Enter a new session
        # ------------------------------------------------------------------------------
        with LoaderCache.session():
            assert LoaderCache.current() is not previous

            # --------------------------------------------------------------------------
            # Simulate a failure in the session body
            # --------------------------------------------------------------------------
            raise ValueError("failure")

    # ----------------------------------------------------------------------------------
    # Ensure that the previous cache is restored on exit, even after an exception
    # ----------------------------------------------------------------------------------
    assert LoaderCache.current() is previous
