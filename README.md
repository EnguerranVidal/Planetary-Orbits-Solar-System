# Presentation

This code was realised by Jonathan Oers and Enguerran Vidal as a numerical physics project as students undergoing studies in our first year SUTS Masters Degree at the University of Paul Sabatier in Toulouse, France. The goal of this project was to model the planetary orbits of the Solar System as well as their perturbations from their interactions. We had to compare different integration schemes to solve the set of differential equations and their influence on the "wobblyness" of our results. We were then instructed to track the perihelion shifts (apsidal precessions) as well as see thechanges in the trajectories of our massive bodies. This was done using Numpy and Matplolib update animation system. For further explanation please check the [Project Report](https://github.com/EnguerranVidal/Planetary-Orbits-Solar-System/blob/main/Numerical_Physics_Project.pdf ).

# Project Content :

In this repository youwill be able to find : 

- **\ephem** folder : folder containing initial positions of the orbiting planets and central body in a .txt format. **Kepler-79_System.txt** : initial positions for the Kepler-79 exoplanetary system, **Solar_System.txt** for the Solar System's initial conditions and **TRAPPIST-1_System.txt** for the TRAPPIST-1 system.
- **\logs** folder : folder containing the save logs of past calculations, currently empty. The **main.py** Python file will automatically create it if not present in its home directory.
- **Numerical_Physics_Report.pdf** : pdf version of our project final report.
- **main.py** : Python file containing all classes and functions of our project. 

# PreRequesites :

In order to work, this code need to be runned in an environment containing :
- Numpy
- Matplotlib

If you encounter any issues with this code, please write us an issue report in order to fix it as quickly as possible, we are not planning any updates at the moment, seen as this project was mandatory for an Uniervsity module.
