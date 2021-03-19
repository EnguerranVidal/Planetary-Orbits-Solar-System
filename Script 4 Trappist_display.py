# DATA MANIPULATION AND PLOTTING MODULES
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# FILE HANDLING MODULES
import os
import sys
import time

from main import*


################################################################################
#--------- SCRIPT  ---------#

X=Planetary_System('TRAPPIST -1 _System.txt')

X.RUN(0.001,10**2,0.1,method='Euler_semi_implicit')
X.display_3D ()
X.energy_conservation ()
