"""
Modifications made to :mod:`datreant` classes on import of module.

"""

from .treants import Treant
from .collections import Bundle
from .limbs import Tags, Categories
from .agglimbs import AggTags, AggCategories

# Treants get tags and categories
for limb in (Tags, Categories):
    Treant._attach_limb_class(limb)

# Bundles get tags and categories aggregation limbs
for agglimb in (AggTags, AggCategories):
    Bundle._attach_agglimb_class(agglimb)
