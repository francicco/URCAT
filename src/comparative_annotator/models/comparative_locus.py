from __future__ import annotations

"""
Compatibility shim for ComparativeLocus.

The canonical implementation lives in:
    comparative_annotator.models.comparative

Do not define a second ComparativeLocus class here, otherwise different parts
of the pipeline may import incompatible objects with the same name.
"""

from comparative_annotator.models.comparative import ComparativeLocus

__all__ = ["ComparativeLocus"]
