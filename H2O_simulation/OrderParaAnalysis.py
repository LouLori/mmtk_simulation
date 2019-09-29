from MMTK import *
from mbpol import mbpolForceField
from MMTK.ForceFields.ForceFieldTest import gradientTest, forceConstantTest
from MMTK.Dynamics import VelocityVerletIntegrator, Heater, \
                          TranslationRemover, RotationRemover
from MMTK import Features
#from MMTK_PINormalModeIntegrator import PINormalModeIntegrator, PILangevinNormalModeIntegrator
from MMTK.Environment import PathIntegrals
from MMTK.NormalModes import VibrationalModes
from MMTK.Trajectory import Trajectory, TrajectoryOutput, \
                            RestartTrajectoryOutput, StandardLogOutput, \
                            trajectoryInfo
from sys import argv,exit
from Scientific.Statistics import mean, standardDeviation
#from nMOLDYN.Mathematics.Analysis import correlation
from Scientific import N
from Scientific.Geometry import Vector
from numpy import *
from numpy.linalg import *


##################################################
######Construct a pseudo universe for mbpol#######
##################################################
universe= InfiniteUniverse()
traj = argv[1]

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

universe.setForceField(mbpolForceField(universe))
##################################################

trajectory = Trajectory(universe, traj)

filename = traj[0:traj.find(".nc")]

universe = trajectory.universe
natoms = universe.numberOfAtoms()
np = universe.numberOfPoints()
P = np/natoms

#oufile = open("oufile-"+filename,"w")
ozfile = open("ozfile-"+filename+".dat","w")

ozA = 0.0
stp = -1

for step in trajectory:
    stp+=1
    universe.setConfiguration(step['configuration'])
    ozN = 0.0
    for i in range(nH2O):
        ozP = 0.0
        for j in range(P):
            xH1 = universe.objectList()[i].atomList()[0].beadPositions()[j]
            xH2 = universe.objectList()[i].atomList()[1].beadPositions()[j]
            xO  = universe.objectList()[i].atomList()[2].beadPositions()[j]
            mH1 = universe.objectList()[i].atomList()[0].mass()
            mH2 = universe.objectList()[i].atomList()[1].mass()
            mO = universe.objectList()[i].atomList()[2].mass()
            xcom = (xH1*mH1+xH2*mH2+xO*mO)/(mH1+mH2+mO)
            rx = [1.,0.,0.]
            rz = [0.,0.,1.]
            dv = xcom-xO
            uv = dv/dv.length()
            uz = uv[2]
            ozP += uz/float(P)
        ozN += ozP/float(nH2O)
    ozfile.write(str(ozN)+"\n")
    ozA += abs(ozN)

#for step in trajectory:
#    stp+=1
#    universe.setConfiguration(step['configuration'])
#    ozN = 0.0
#    for i in range(nH2O):
#        ozP = 0.0
#        xH1 = universe.objectList()[i].atomList()[0].position()
#        xH2 = universe.objectList()[i].atomList()[1].position()
#        xO  = universe.objectList()[i].atomList()[2].position()
#        mH1 = universe.objectList()[i].atomList()[0].mass()
#        mH2 = universe.objectList()[i].atomList()[1].mass()
#        mO = universe.objectList()[i].atomList()[2].mass()
#        xcom = (xH1*mH1+xH2*mH2+xO*mO)/(mH1+mH2+mO)
#        rx = [1.,0.,0.]
#        rz = [0.,0.,1.]
#        dv = xcom-xO
#        uv = dv/dv.length()
#        uz = uv[2]
#        ozN += uz/float(nH2O)
#    ozfile.write(str(ozN)+"\n")
#    ozA += abs(ozN)

ozA = ozA/stp
print "Average Absolute Order Parameters:", ozA

trajectory.close()
#oufile.close()
ozfile.close()
