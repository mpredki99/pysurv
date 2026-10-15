# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from pysurveyor.typing.vars import DomainClassType

from .pysurveyor_loader import PySurveyorLoader


# ======================================================================================
# Base class for loaders bound to a domain class instance
# ======================================================================================
class ObjectBoundLoader(PySurveyorLoader[DomainClassType]):
    """Base class for loaders bound to a domain class instance.

    ``ObjectBoundLoader`` stores the domain class instance used by concrete loader
    implementations to create domain objects.

    Attributes
    ----------
    _object : DomainClassType
        Domain class instance used to create domain objects.
    """

    __pysurveyor_abstract__ = True

    # ==================================================================================
    def __init__(self, domain_class_object: DomainClassType) -> None:
        super().__init__()
        self._object = domain_class_object
