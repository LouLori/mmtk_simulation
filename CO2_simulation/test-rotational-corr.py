from MMTK import *
from HeCO2TransFF import HeCO2TransForceField
from HeHeFF import HeHeForceField
from NoPotFF import NoPotForceField
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
from numpy import zeros, correlate, asarray, real, linalg, dot
from numpy.linalg import norm

traj = argv[1]
temperature=argv[2]
nCO2=1
label="co2test-1dz"
trajectory = Trajectory(None, traj)
universe = trajectory.universe
natoms = universe.numberOfAtoms()

np = universe.numberOfPoints()
P = np/natoms

tau=1.0/(float(P)*float(temperature))

rfile=open("final-rot-"+str(P)+"-"+label,"w")

c=zeros(P+1,float)

r=zeros((2,P,3),float)
counter=0
for step in trajectory:
    counter+=1
    universe.setConfiguration(step['configuration'])

    r[0]=asarray(universe.atomList()[1].beadPositions())
    r[1]=asarray(universe.atomList()[2].beadPositions())

    #bond=zeros((P+1,3),float)
    bond=zeros(P+1,float)
    for i in range(P):
        bond[i]=(r[0][i]-r[1][i])[2]/linalg.norm(r[0][i]-r[1][i])
        #for j in range(3):
        #    bond[i][j]=(r[0][i]-r[1][i])[j]/linalg.norm(r[0][i]-r[1][i])
            

    bond[P]=bond[0]
    for deltat in range(P+1):
        for corrtime in range(P+1-deltat):
            #c[deltat]+=dot(bond[corrtime],bond[corrtime+deltat])/(P+1-deltat)
            c[deltat]+=(bond[corrtime]*bond[corrtime+deltat])/(P+1-deltat)

for i in range(P+1):
    c[i]=c[i]/counter

for i in range(P+1):
    rfile.write(str(i*tau)+" "+str(c[i])+"\n")

           
rfile.close()
#vfile.close()
trajectory.close()
