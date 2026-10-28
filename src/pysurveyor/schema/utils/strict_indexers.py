# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from typing import Any

from pysurveyor.schema.utils.model_validation import assert_valid_model


# ======================================================================================
# Base class for strict model indexers
# ======================================================================================
class _StrictIndexer:
    """Wrap a DataFrame indexer with assignment validation.

    Parameters
    ----------
    schema : StrictSchema
        The schema whose model is accessed.
    indexer_name : str
        The underlying DataFrame indexer name.
    """

    __slots__ = ("_schema", "_indexer_name")

    # ==================================================================================
    # Dunder methods
    # ==================================================================================
    def __init__(self, schema: Any, indexer_name: str) -> None:
        """Initialize the indexer wrapper."""
        self._schema = schema
        self._indexer_name = indexer_name

    # ==================================================================================
    def __getitem__(self, key: Any) -> Any:
        """Return data selected by the underlying DataFrame indexer."""
        indexer = getattr(self._schema._model, self._indexer_name)
        return indexer[key]

    # ==================================================================================
    def __setitem__(self, key: Any, value: Any) -> None:
        """Apply assignment, then validate and commit."""
        schema = self._schema
        schema_name = schema.__class__.__name__

        # Copy the model to avoid modifying the original until validation passes -------
        candidate = schema._model.copy(deep=True)
        indexer = getattr(candidate, self._indexer_name)
        indexer[key] = value

        # ------------------------------------------------------------------------------
        assert_valid_model(
            candidate,
            schema._expected_columns,
            schema._expected_index,
            schema_name,
        )

        # Commit only after every validation succeeds ----------------------------------
        object.__setattr__(schema, "_model", candidate)


# ======================================================================================
# Label-based indexer
# ======================================================================================
class StrictLocIndexer(_StrictIndexer):
    """Provide label-based model access with validated assignments."""

    __slots__ = ()

    def __init__(self, schema: Any) -> None:
        """Initialize the label-based indexer."""
        super().__init__(schema, "loc")


# ======================================================================================
# Position-based indexer
# ======================================================================================
class StrictIlocIndexer(_StrictIndexer):
    """Provide position-based model access with validated assignments."""

    __slots__ = ()

    def __init__(self, schema: Any) -> None:
        """Initialize the position-based indexer."""
        super().__init__(schema, "iloc")


# ======================================================================================
# Scalar label-based indexer
# ======================================================================================
class StrictAtIndexer(_StrictIndexer):
    """Provide scalar label-based access with validated assignments."""

    __slots__ = ()

    def __init__(self, schema: Any) -> None:
        """Initialize the scalar label-based indexer."""
        super().__init__(schema, "at")


# ======================================================================================
# Scalar position-based indexer
# ======================================================================================
class StrictIatIndexer(_StrictIndexer):
    """Provide scalar position-based access with validated assignments."""

    __slots__ = ()

    def __init__(self, schema: Any) -> None:
        """Initialize the scalar position-based indexer."""
        super().__init__(schema, "iat")
