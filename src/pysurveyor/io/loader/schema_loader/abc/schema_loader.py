# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

from pysurveyor.io.loader.abc.table_loader import TableLoader
from pysurveyor.typing.vars import SchemaType


# ======================================================================================
# Abstract base class for schema loaders
# ======================================================================================
class SchemaLoader(TableLoader[SchemaType]):
    """Base loader for schema objects.

    Type Parameters
    ---------------
    SchemaType
        Concrete schema type returned by the loader.
    """

    __pysurveyor_abstract__ = True
