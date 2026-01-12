import math
import random
from typing import List , Tuple , Optional , Union
from functools import cmp_to_key
import matplotlib as plt
import numpy as np
import pandas as pd
from bitalg.visualizer.main import Visualizer
from time import time
import os
import matplotlib.pyplot as plt
import kod_loboda_werminski as kod

#%matplotlib tk

POINT = (20, 20)

vis = Visualizer()
Pol = kod.draw_planar_graph()
V, E = kod.convert_regions_to_graph(Pol)
step1 = vis.add_point(V, color="blue")
step2 = vis.add_line_segment([[V[i], V[j]] for i, j in E], color="blue")
vis.save("step1")
V_, E_ = kod.process_polygons_to_mesh(Pol)
vis.remove_figure(step2)
step3 = vis.add_line_segment([[V[i], V[j]] for i, j in E_], color="red")
step4 = vis.add_line_segment([[V[i], V[j]] for i, j in E], color="blue")
vis.save("step2")
vis = kod.animate_point_location(V_, E_, POINT)
vis.save("step3")
vis.save_gif("vis")



