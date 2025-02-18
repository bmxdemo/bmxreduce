from astropy.utils import iers
iers.conf.iers_auto_url="file:///gpfs02/astro/workarea/bmxdata/almanac/finals2000A.all"
iers.conf.iers_data = "/gpfs02/astro/workarea/bmxdata/almanac/finals2000A.all"
from .datamanager import datamanager
#from .farmit import farmit
from .reduce import reduce
from .satellites import Satellites

#from .reduce_init import reduce
#from .mapmanager import mapmanager
#from .reduce_coadd import coaddbygroup



