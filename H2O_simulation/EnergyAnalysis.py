from MMTK import *
from MMTK.Environment import PathIntegrals
from MMTK.Trajectory import Trajectory
from sys import argv
from Scientific import N
from Scientific.Statistics import mean, standardDeviation
from math import sqrt
from numpy import zeros, pi, cos, sin, sqrt
from mbpol import mbpolForceField
from MMTK.ForceFields.ForceField import CompoundForceField
from qTIP4pFF import HarmonicAngleForceField, HarmonicBondForceField, QuarticBondForceField, ElectrostaticForceField, LennardJonesForceField
#This is the conversion factor to Units of K
Kper1overcm=11604.505/8065.54445
conv=Kper1overcm/1.196e-2
traj = argv[1]
fftype = "r-tip4p"
##################################################
######Construct a pseudo universe for mbpol#######
##################################################
universe= InfiniteUniverse()

#set parameters
temperature = float(traj[traj.find("T")+1:traj.find("-R")])
nH2O = int(traj[0:traj.find("H2O")])
P = int(traj[traj.find("P")+1:traj.find("-T")])
lattice_spacing = float(traj[traj.find("R")+1:traj.find("-local")])

universe.addObject(PathIntegrals(temperature))

for i in range(nH2O):
    universe.addObject(Molecule('water', position = Vector(0., 0., i*lattice_spacing)))

for atom in universe.atomList():
    atom.setNumberOfBeads(P)

natoms = len(universe.atomList())

###############  WATER MOLECULE PARAMETER LIST  ################

if fftype == "q-tip4p" or "r-tip4p":
    #q-tip4p quartic bond
    D_r = 116.09 * 4.184            #kJ/mol
    a_r = 2.287 * 10.            #nm^-1
    c1 = a_r ** 2
    c2 = a_r ** 3
    c3 = (a_r ** 4) * 7. / 12.
    r_eq = 0.09419                          #nm

    #q-tip4p harmonic angle
    theta_rad = 107.4 * pi / 180.                   #radians
    amber_k_theta = 87.85                           #kcal/mol/(rad**2)
    k_theta = 4.184 * amber_k_theta                 #kj/mol/(rad**2)

    #q-tip4p electrostatics
    fraction = 0.73612
    o_charge = -1.1128

    #q-tip4p lennard jones
    epsilon = 0.1852 * 4.184        #kj/mol
    sigma   = 0.31589               #nm

elif fftype == "q-spc-fw":
    #spcfw-q harmonic bond
    r_eq = 0.1
    amber_k_r = 529.581       #kcal/mol(Ang**2)
    k_r = 4.184*amber_k_r*0.01       #kj/mol(nm**2)

    #spcfw-q harmonic angle
    theta_rad = 112.04 * pi / 180.                 #radians
    amber_k_theta = 37.95                          #kcal/mol/(rad**2)
    k_theta = 2. * 4.184 * amber_k_theta           #kj/mol/(rad**2)

    #spcfw-q electrostatics
    fraction = 1.0
    o_charge = -0.84

    #spcfw-q Lennard-Jones
    epsilon = 0.650194             #parameters from nonbonded.c
    sigma = 0.316536               #parameters from nonbonded.c
################# END WATER PARAMETERS ##############################
ff = []
################# CONSTRUCT FORCEFIELD ###############################
if fftype != "mbpol":
    for o in range(nH2O):
        a_list=universe.objectList()[o].atomList()
        if fftype == "q-tip4p":
            ff.append(QuarticBondForceField([a_list[0],a_list[2]],D_r,c1,c2,c3,r_eq))
            ff.append(QuarticBondForceField([a_list[1],a_list[2]],D_r,c1,c2,c3,r_eq))
            ff.append(HarmonicAngleForceField([a_list[0],a_list[2],a_list[1]],theta_rad,k_theta))
        elif fftype == "q-spc-fw":
            ff.append(HarmonicBondForceField([a_list[0],a_list[2]],r_eq,k_r))
            ff.append(HarmonicBondForceField([a_list[1],a_list[2]],r_eq,k_r))
            ff.append(HarmonicAngleForceField([a_list[0],a_list[2],a_list[1]],theta_rad,k_theta))

    es = ElectrostaticForceField(universe,fraction,o_charge)
    ff.append(es)
    lj = LennardJonesForceField(universe,epsilon,sigma)
    ff.append(lj)
else:
    ff.append(mbpolForceField(universe))


universe.setForceField(CompoundForceField(*ff))
trajectory = Trajectory(universe, traj)
###################################################

npoints = len(trajectory)
universe = trajectory.universe
natoms = universe.numberOfAtoms()
time=trajectory.time
np = universe.numberOfPoints()
P = np/natoms

kjtokcal = 0.239001
print Units.k_B
# Print averages of the quantu energy estimator
print "Number of Atoms:", natoms
print "Number of Beads:", P
print "Number of Steps:", npoints

Temp = trajectory.temperature
Eprim = trajectory.quantum_energy_primitive
Ecenv = trajectory.quantum_energy_centroid_virial
Epot = trajectory.potential_energy
Erot = trajectory.quantum_energy_rotation

#print "Temperature:", mean(Temp), "+/-",  standardDeviation(Temp)/sqrt(npoints), " K"
#print "Primitive estimator:", mean(Eprim/Units.k_B), "+/-",  standardDeviation(Eprim/Units.k_B)/sqrt(npoints), " K"
#print "Centroid Virial estimator:", mean(Ecenv/Units.k_B), "+/-", standardDeviation(Ecenv/Units.k_B)/sqrt(npoints), " K"

#print "Potential estimator:", mean(Eprim)/Units.k_B, "+/-",  standardDeviation(Eprim)/Units.k_B, " K"
##print "Potential estimator:", mean(Epot/Units.k_B), "+/-",  standardDeviation(Epot/Units.k_B)/sqrt(npoints), " K"
#print "Kinetic estimator:",mean(Ecenv)/Units.k_B, "+/-", standardDeviation(Ecenv)/Units.k_B, " K"
##print "Kinetic estimator:",mean(Eprim-Epot)/Units.k_B, "+/-", standardDeviation(Eprim/Units.k_B)/sqrt(npoints)-standardDeviation(Epot/Units.k_B)/sqrt(npoints), " K"
#print "Rotational Energy estimator:", mean(Erot)/Units.k_B, "+/-",  standardDeviation(Erot)/Units.k_B, " K"
#
#print "Total Energy:", mean(Eprim)/Units.k_B+mean(Erot)/Units.k_B, "+/-", sqrt((standardDeviation(Eprim)/Units.k_B)**2+(standardDeviation(Erot)/Units.k_B)**2), " K"

print "Potential estimator:", mean(Eprim), "+/-",  standardDeviation(Eprim), " kj/mol"
#print "Potential estimator:", mean(Epot/Units.k_B), "+/-",  standardDeviation(Epot/Units.k_B)/sqrt(npoints), " K"
print "Kinetic estimator:",mean(Ecenv), "+/-", standardDeviation(Ecenv), " kj/mol"
#print "Kinetic estimator:",mean(Eprim-Epot)/Units.k_B, "+/-", standardDeviation(Eprim/Units.k_B)/sqrt(npoints)-standardDeviation(Epot/Units.k_B)/sqrt(npoints), " K"
print "Rotational Energy estimator:", mean(Erot), "+/-",  standardDeviation(Erot), " kj/mol"

print "Total Energy:", mean(Eprim)+mean(Erot), "+/-", sqrt((standardDeviation(Eprim))**2+(standardDeviation(Erot))**2), " kj/mol"


trajectory.close()
