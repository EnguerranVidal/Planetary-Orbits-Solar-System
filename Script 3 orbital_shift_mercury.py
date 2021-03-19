# DATA MANIPULATION AND PLOTTING MODULES
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# FILE HANDLING MODULES
import os
import sys
import time

from man import*


################################################################################
#--------- SCRIPT ---------#

X=Planetary_System('Solar_System.txt')

X.load_session("semi_implicit.txt")
X.apsidal_precession(displayed=['Mercury'])
