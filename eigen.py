from core import boundary
from core import eigen

dom=boundary.Domain([0,20])
eigen.helmholtz(100,dom)
