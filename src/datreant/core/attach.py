"""
Modifications made to :mod:`datreant` classes on import of module.

"""

from .treants import Treant, Group
from .collections import Bundle
from .limbs import Tags, Categories, MemberBundle
from .agglimbs import AggTags, AggCategories

# Treants get tags and categories
for limb in (Tags, Categories):
    Treant._attach_limb_class(limb)

# Groups get members
Group._attach_limb_class(MemberBundle)

# Bundles get tags and categories aggregation limbs
for agglimb in (AggTags, AggCategories):
    Bundle._attach_agglimb_class(agglimb)
