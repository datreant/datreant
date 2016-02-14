"""
Modifications made to :mod:`datreant` classes on import of module.

"""

from .treants import Treant, Group
from .limbs import Tags, Categories, MemberBundle

# Treants get tags and categories
for limb in (Tags, Categories):
    Treant._attach_limb_class(limb)

# Groups get members
Group._attach_limb_class(MemberBundle)
