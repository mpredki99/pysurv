# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from os import PathLike
from pathlib import Path

from pysurveyor.backend.backend import NA, Series, isna, read_csv
from pysurveyor.basic.units import parse_unit
from pysurveyor.io.cache.cache_keys import create_file_cache_key
from pysurveyor.io.cache.loader_cache import LoaderCache
from pysurveyor.io.config.file_configs import CSVConfig
from pysurveyor.io.loader.abc.table_loader import TableLoader
from pysurveyor.schema.utils.validation_mode import ValidationMode
from pysurveyor.typing.vars import SchemaType
from pysurveyor.validators.utils.parser import parse_validator


# ======================================================================================
# Generic loader for schema models
# ======================================================================================
class SchemaLoader(TableLoader[SchemaType]):
    """Load schema objects from a model sources.

    Parameters
    ----------
    constructor : type[SchemaType]
        Schema class used to construct the loaded schema object.
    """

    # ==================================================================================
    def from_csv(
        self,
        path: PathLike[str] | str,
        config: CSVConfig | None = None,
        *,
        validation_mode: ValidationMode | str | None = None,
    ) -> SchemaType:
        """Load a schema object from a schema-model CSV file.

        Parameters
        ----------
        path : PathLike[str] or str
            Path to the schema-model CSV file.
        config : CSVConfig, optional
            CSV reader configuration. Defaults to CSVConfig().
        validation_mode : ValidationMode or str, optional
            Validation behavior forwarded to the schema constructor.
            If None, the schema class default is used.

        Returns
        -------
        SchemaType
            Instance created by the supplied schema class.

        Raises
        ------
        FileNotFoundError
            If the path does not point to an existing file.
        """
        path = Path(path).expanduser().resolve()

        # ------------------------------------------------------------------------------
        if not path.is_file():
            raise FileNotFoundError(f"CSV file not found: {path}")

        # ------------------------------------------------------------------------------
        config = config or CSVConfig()

        key = create_file_cache_key(path, config.kwargs)
        model = LoaderCache.current().get_data(
            key, lambda: read_csv(path, **config.kwargs)
        )

        # ------------------------------------------------------------------------------
        # Parse units and validators in the schema model
        # ------------------------------------------------------------------------------
        if "unit" in model.columns:
            model["unit"] = model["unit"].map(
                lambda u: parse_unit(u.strip()) if not isna(u) else NA
            )

        if "validator" in model.columns:
            model["validator"] = model["validator"].map(
                lambda v: parse_validator(v.strip()) if not isna(v) else NA
            )

        # ------------------------------------------------------------------------------
        if validation_mode is None:
            return self._constructor(model)

        return self._constructor(model, validation_mode=validation_mode)
