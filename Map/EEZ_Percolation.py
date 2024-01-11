#%%
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from shapely.ops import unary_union
import distinctipy as dpy
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.colors as colors
from pyproj import CRS

#%%

# Load data 

trayectories_df = pd.read_csv("../data/full_data_inds.csv.gz", sep=";", compression="gzip")
