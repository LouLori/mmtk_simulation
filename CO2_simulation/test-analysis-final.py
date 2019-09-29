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
from numpy import *
from numpy.linalg import *

traj = argv[1]
nCO2=1
label="norotskip-"+argv[2]
trajectory = Trajectory(None, traj)
universe = trajectory.universe
natoms = universe.numberOfAtoms()

np = universe.numberOfPoints()
P = np/natoms

stepcount = 0

rval=zeros(len(trajectory)*P,float)
cval=zeros(len(trajectory)*P,float)
#vfile=open("final-pot-"+str(P)+"-"+label,"w")
#rfile=open("hist-r-"+str(P)+"-"+label,"w")
#ctfile=open("hist-cost-"+str(P)+"-"+label,"w")

rfile=open("data-r-"+str(P)+"-"+label,"w")
ctfile=open("data-cost-"+str(P)+"-"+label,"w")

print universe.atomList()[0]
print universe.atomList()[1]
print universe.atomList()[2]

stepno=-1
for step in trajectory:
    stepno+=1
    universe.setConfiguration(step['configuration'])

    #datarot = trajectory.quantum_energy_rotation
    #datapot = trajectory.quantum_energy_primitive
    
    #rfile.write(str(datarot)+"\n")
    #vfile.write(str(datapot)+"\n")

    for i in range(3*nCO2,natoms):
        for j in range(P):
            r=universe.atomList()[0].beadPositions()[j]-universe.atomList()[i].beadPositions()[j]
            bond=universe.atomList()[1].beadPositions()[j]-universe.atomList()[2].beadPositions()[j]
            rval[stepno*P+j]=r.length()
            cval[stepno*P+j]=dot(bond/bond.length(),r/r.length())
            
            rfile.write(str(r.length())+"\n")
            ctfile.write(str(dot(bond/bond.length(),r/r.length()))+"\n")
   #     bond1=(universe.atomList()[0].beadPositions()[j]-universe.atomList()[1].beadPositions()[j]).length()
   #     bond2=(universe.atomList()[0].beadPositions()[j]-universe.atomList()[2].beadPositions()[j]).length()
   #     rfile.write(str(bond1)+"\n")
   #     rfile.write(str(bond2)+"\n")
    #        efile.write(str(dot(bond1,bond2)/(bond1.length()*bond2.length()))+"\n")
    #        #efile.write(str(bond1[2])+"\n")
    #        #xfile.write(str(abs(bond1[0])/bond1.length())+"\n")


rfile.close()
ctfile.close()
trajectory.close()
raise()
print "Now constructing histogram"
rhist=histogram(rval,bins=50, density=True)
chist=histogram(cval,bins=50, density=True)

dr=(rhist[1][1]-rhist[1][0])
rmin=rhist[1][0]+dr/2.0

dc=(chist[1][1]-chist[1][0])
cmin=chist[1][0]+dc/2.0


for i in range(50):
    rfile.write(str(rmin+i*dr)+" "+str(rhist[0][i])+"\n")
    ctfile.write(str(cmin+i*dc)+" "+str(chist[0][i])+"\n")

rfile.close()
ctfile.close()
trajectory.close()
