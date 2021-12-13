import pathlib

_ROOT = pathlib.Path(__file__).parent.absolute()

# create the ipopt-watertap solver and register
# it as the default IDAES solver
import watertap.core.plugins
