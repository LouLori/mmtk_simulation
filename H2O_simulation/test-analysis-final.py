from MMTK import *
from CoulombFF import CoulombForceField
from LJFF import LJForceField
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
filename = traj[:traj.find(".nc")]
nH2O = int(traj[:1])

trajectory = Trajectory(None, traj)
universe = trajectory.universe
natoms = universe.numberOfAtoms()

np = universe.numberOfPoints()
P = np/natoms

afile=open("afile-"+filename,"w")
rfile=open("rfile-"+filename,"w")
cfile=open("cfile-"+filename,"w")
bfile=open("bfile-"+filename,"w")

print universe.atomList()[0]
print universe.atomList()[1]
print universe.atomList()[2]

for step in trajectory:
    universe.setConfiguration(step['configuration'])
    for i in range(nH2O):
        for j in range(P):
            xm1H1=universe.objectList()[0].atomList()[0].beadPositions()[i]
            xm1H2=universe.objectList()[0].atomList()[1].beadPositions()[i]
            xm1O=universe.objectList()[0].atomList()[2].beadPositions()[i]
            xm2H1=universe.objectList()[1].atomList()[0].beadPositions()[i]
            xm2H2=universe.objectList()[1].atomList()[1].beadPositions()[i]
            xm2O=universe.objectList()[1].atomList()[2].beadPositions()[i]
            mH1=universe.objectList()[0].atomList()[0].mass()
            mH2=universe.objectList()[0].atomList()[1].mass()
            mO=universe.objectList()[0].atomList()[2].mass()

            xcom1=(xm1H1*mH1+xm1H2*mH2+xm1O*mO)/(mH1+mH2+mO)
            xcom2=(xm2H1*mH1+xm2H2*mH2+xm2O*mO)/(mH1+mH2+mO)

            rcom=xcom1-xcom2
            rOO=xm1O-xm2O
            rOH1=xm2O-xm2H1
            rOH2=xm2O-xm2H2
            aHOH=arccos(dot(rOH1/rOH1.length(),rOH2/rOH2.length()))
            rv1=xm1O-xcom1
            rv2=xm2O-xcom2
            rz = [0.,0.,1.]
            cost=rv2[2]/rv2.length()
            #afile.write(str(aHOH)+"\n")
            rfile.write(str(rcom.length())+"\n")
            cfile.write(str(cost)+"\n")
            #bfile.write(str(rOH2.length())+"\n")

rfile.close()
cfile.close()
afile.close()
bfile.close()
trajectory.close()

#print "Now constructing histogram"
#rhist=histogram(rval,bins=50, density=True)
#chist=histogram(cval,bins=50, density=True)
#
#dr=(rhist[1][1]-rhist[1][0])
#rmin=rhist[1][0]+dr/2.0
#
#dc=(chist[1][1]-chist[1][0])
#cmin=chist[1][0]+dc/2.0
#
#
#for i in range(50):
#    rfile.write(str(rmin+i*dr)+" "+str(rhist[0][i])+"\n")
#    ctfile.write(str(cmin+i*dc)+" "+str(chist[0][i])+"\n")
#
#rfile.close()
#ctfile.close()
#trajectory.close()
