# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from collections.abc import Set

from pysurveyor.backend.backend import DataFrame, Index, RangeIndex, isna


# ======================================================================================
# Model validation for PySurveyor schemas
# ======================================================================================
def validate_model(
    model: DataFrame,
    expected_columns: dict[str, type],
    expected_index: Index,
    schema_name: str,
) -> DataFrame:
    """Standardize and validate a model's index, columns, and non-null column values."""

    expected_index_name = expected_index.name
    # ----------------------------------------------------------------------------------
    if expected_index_name in model.columns:
        if not isinstance(model.index, RangeIndex):
            raise ValueError(
                f"Cannot set index for {schema_name}: column {expected_index_name!r} "
                f"exists, but the model already has an index {model.index.name!r}."
            )
        model = model.set_index(expected_index_name)

    # ----------------------------------------------------------------------------------
    assert_valid_model(model, expected_columns, expected_index, schema_name)
    # ----------------------------------------------------------------------------------
    return model


# ======================================================================================
def assert_valid_model(
    model: DataFrame,
    expected_columns: dict[str, type],
    expected_index: Index,
    schema_name: str,
) -> None:
    """Validate a model's index, columns, and non-null column values."""
    assert_valid_index(model.index, expected_index, schema_name)
    assert_valid_columns(model, expected_columns, schema_name)


# ======================================================================================
def assert_valid_index(index: Index, expected_index: Index, schema_name: str) -> None:
    """Validate the index name and labels, ignoring label order."""
    if (
        index.name == expected_index.name
        and index.is_unique
        and expected_index.is_unique
        and len(index) == len(expected_index)
        and index.isin(expected_index).all()
    ):
        return

    # ----------------------------------------------------------------------------------
    raise ValueError(
        f"Invalid index for {schema_name}. "
        f"Expected {expected_index.to_list()} with name "
        f"{expected_index.name!r}, got {index.to_list()} "
        f"with name {index.name!r}."
    )


# ======================================================================================
def assert_valid_columns(
    model: DataFrame, expected_columns: dict[str, type], schema_name: str
) -> None:
    """Validate column names and non-null values against expected types."""

    # ----------------------------------------------------------------------------------
    assert_valid_columns_existence(model.columns, expected_columns.keys(), schema_name)
    # ----------------------------------------------------------------------------------
    assert_valid_columns_values(model, expected_columns, schema_name)


# ======================================================================================
def assert_valid_columns_existence(
    columns: Index, expected_columns: Set[str], schema_name: str
) -> None:
    """Validate column names against expected types."""
    if (
        columns.is_unique
        and len(columns) == len(expected_columns)
        and columns.isin(expected_columns).all()
    ):
        return

    # ----------------------------------------------------------------------------------
    missing = [column for column in expected_columns if column not in columns]
    unexpected = [column for column in columns if column not in expected_columns]

    # ----------------------------------------------------------------------------------
    raise ValueError(
        f"Invalid columns for {schema_name}. "
        f"Missing: {missing!r}; unexpected: {unexpected!r}; "
        f"got: {columns!r}; expected: {expected_columns!r}."
    )


# ======================================================================================
def assert_valid_columns_values(
    model: DataFrame, expected_columns: dict[str, type], schema_name: str
) -> None:
    """Validate column names and non-null values against expected types."""

    # Validate column data type --------------------------------------------------------
    for col, dtype in expected_columns.items():
        # ------------------------------------------------------------------------------
        invalid_values = [
            value
            for value in model[col]
            if not isna(value) and not isinstance(value, dtype)
        ]
        # ------------------------------------------------------------------------------
        if invalid_values:
            raise ValueError(
                f"Invalid dtype in column {col!r} for {schema_name}. "
                f"Expected {dtype.__name__}; "
                f"invalid values: {invalid_values!r}."
            )
