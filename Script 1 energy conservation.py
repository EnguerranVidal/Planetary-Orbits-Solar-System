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
#--------- SCRIPT ---------#

X=Planetary_System('Solar_System.txt')

X.RUN(0.1,10**4,5,method='Euler_explicit')
X.energy_conservation()

X.new_session('Solar_System.txt')
X.RUN(0.1,10**5,5,method='Euler_semi_implicit')
X.energy_conservation()

X.new_session('Solar_System.txt')
X.RUN(0.1,3*10**4,5,method='Euler_symplectic')
X.energy_conservation()

X.new_session('Solar_System.txt')
X.RUN(0.1,10**4,5,method='Heun')
X.energy_conservation()

X.new_session('Solar_System.txt')
X.RUN(0.1,3*10**4,5,method='Runge_Kutta')
X.energy_conservation()
