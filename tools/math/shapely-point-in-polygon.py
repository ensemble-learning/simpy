# point in polygon
import numpy as np
from math import sqrt

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

point = Point(13.0, 0.0)
pt = np.array([[-1, 0], [-0.24*sqrt(3), 0.9], [0.24*sqrt(3), 0.9],
               [1, 0], [0.24*sqrt(3), -0.9], [-0.24*sqrt(3), -0.9],])
ns = 0
for ns in range(1,20,2):
    polygon = Polygon(ns*pt)
    if polygon.contains(point) == True:
        break
print(ns)
