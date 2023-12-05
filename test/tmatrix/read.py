import pyarrow.feather as pf
import pandas as pd

d = pf.read_feather('long.arrow')

d.head()
