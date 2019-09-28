from MMTK import *
from MMTK.Environment import PathIntegrals
from MMTK.Trajectory import Trajectory
from sys import argv
from Scientific import N
from Scientific.Statistics import mean, standardDeviation
from math import sqrt

#This is the conversion factor to Units of K
Kper1overcm=11604.505/8065.54445
conv=Kper1overcm/1.196e-2

traj = argv[1]
trajectory = Trajectory(None, traj)
npoints = len(trajectory)
universe = trajectory.universe
natoms = universe.numberOfAtoms()
time=trajectory.time
np = universe.numberOfPoints()
P = np/natoms

# Print averages of the quantu energy estimator
print "Number of Atoms:", natoms
print "Number of Beads:", P
print "Number of Steps:", npoints

Temp = trajectory.temperature
Eprim = trajectory.quantum_energy_primitive
Ecenv = trajectory.quantum_energy_centroid_virial
Epot = trajectory.potential_energy
Erot = trajectory.quantum_energy_rotation

print "Temperature:", mean(Temp), "+/-",  standardDeviation(Temp)/sqrt(npoints), " K"
print "Primitive estimator:", mean(Eprim/Units.k_B), "+/-",  standardDeviation(Eprim/Units.k_B)/sqrt(npoints), " K"
print "Centroid Virial estimator:", mean(Ecenv/Units.k_B), "+/-", standardDeviation(Ecenv/Units.k_B)/sqrt(npoints), " K"
print "Potential estimator:", mean(Epot/Units.k_B), "+/-",  standardDeviation(Epot/Units.k_B)/sqrt(npoints), " K"
print "Kinetic estimator:",mean(Eprim-Epot)/Units.k_B, "+/-", standardDeviation(Eprim/Units.k_B)/sqrt(npoints)-standardDeviation(Epot/Units.k_B)/sqrt(npoints), " K"
print "Rotational Energy estimator:", mean(Erot/Units.k_B), "+/-",  standardDeviation(Erot/Units.k_B)/sqrt(npoints), " K"
print "Total Energy:", mean(Eprim+Erot)/Units.k_B, "+/-", sqrt((standardDeviation(Eprim)/Units.k_B/sqrt(npoints))**2+(standardDeviation(Erot)/Units.k_B/sqrt(npoints))**2), " K"

trajectory.close()
