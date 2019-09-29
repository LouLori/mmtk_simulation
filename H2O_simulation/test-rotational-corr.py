from MMTK import *
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
filename = traj[:traj.find(".nc")]
nH2O = int(traj[:1])
trajectory = Trajectory(None, traj)
universe = trajectory.universe
natoms = universe.numberOfAtoms()
stp = len(trajectory)

np = universe.numberOfPoints()
P = np/natoms


#crfile=open("rotcorrcent"+"-"+label+"-"+str(P),"w")
#cpfile=open("potcorrcent"+"-"+label+"-"+str(P),"w")
#crefile=open("rotcorreng"+"-"+label+"-"+str(P),"w")
cafile=open("rotcorravg"+"-"+filename,"w")

#datap = trajectory.quantum_energy_primitive
#datar = trajectory.quantum_energy_rotation
#cim=zeros(P+1,float)
cre1=zeros(stp+1,float)
#pre = zeros(stp+1,float)
#rre = zeros(stp+1,float)
#cre2=zeros(stp+1,float)

r=zeros((stp,P,3),float)
#rcm=zeros((stp,P,3),float)
#rcm_cent=zeros((stp,3),float)
#r3_cent=zeros((stp,3),float)

counter=0
#print universe.atomList()[0]
#print universe.atomList()[1]
#print universe.atomList()[2]
#print universe.atomList()[3]
m0=universe.atomList()[0].mass()
m1=universe.atomList()[1].mass()
m2=universe.atomList()[2].mass()
#m3=universe.atomList()[3].mass()

#dist1=zeros((stp+1),float)
#dist2=zeros((stp+1,P),float)
#uv=zeros((stp+1,P,3),float)
zv=zeros((3), float)
coscent=zeros((stp+1),float)
#datapcent=zeros((stp+1),float)
#datarcent=zeros((stp+1),float)

zv[2]=1.

for step in trajectory:
    universe.setConfiguration(step['configuration'])
    #r[0][counter]=asarray(universe.atomList()[0].beadPositions()) #H1 position
    #r[1][counter]=asarray(universe.atomList()[1].beadPositions()) #H2 position
    r[counter]=asarray(universe.objectList()[0].atomList()[2].beadPositions()) #O position
    #r[3][counter]=asarray(universe.atomList()[3].beadPositions()) #He position
    #rcm[counter]=(r[0][counter]*m0+r[1][counter]*m1+r[2][counter]*m2)/(m0+m1+m2)
    #uv[counter]=(r[2][counter]-r[1][counter])/(linalg.norm(r[2][counter]-r[1][counter]))
    for j in range(P):
        coscent[counter] += 1./P*(r[counter][j][2])/(linalg.norm(r[counter][j]))
#rcm_cent[counter] += 1./P*(rcm[counter][j])
#r3_cent[counter] += 1./P*(r[3][counter][j])
#print N.dot((r[2][counter][j]-r[1][counter][j]),zv)/(linalg.norm(r[2][counter][j]-r[1][counter][j]))
#        coscent[counter] += 1./P*(N.dot((r[2][counter][j]-r[1][counter][j]),zv))/(linalg.norm(r[2][counter][j]-r[1][counter][j]))
#    datapcent[counter] = datap[counter]
    #datarcent[counter] = datar[counter]
#print coscent[counter]
#dist1[counter]=linalg.norm(r3_cent[counter]-rcm_cent[counter])
    counter+=1
#raise()
#dist1[stp]=dist1[0]
#dist2[stp]=dist2[0]
#uv[stp]=uv[0]
coscent[stp]=coscent[0]
#datapcent[stp]=datap[0]
#datarcent[stp]=datar[0]

print "Cos theta filling finished"

for deltat in range(stp+1):
    for corrtime in range(stp+1-deltat):
#cre1[deltat]+=(dist1[corrtime]*dist1[corrtime+deltat])/(stp+1-deltat)
        cre1[deltat]+=(coscent[corrtime]*coscent[corrtime+deltat])/(stp+1-deltat)
#        pre[deltat]+=(datapcent[corrtime]*datapcent[corrtime+deltat])/(stp+1-deltat)
#        rre[deltat]+=(datarcent[corrtime]*datarcent[corrtime+deltat])/(stp+1-deltat)
#        print cre1[deltat]
#    crfile.write(str(deltat)+" "+str(cre1[deltat])+"\n")
#    cpfile.write(str(deltat)+" "+str(pre[deltat])+"\n")
#    crefile.write(str(deltat)+" "+str(rre[deltat])+"\n")
    cafile.write(str(deltat)+" "+str(cre1[deltat])+"\n")
#        for i in range(P):
#cre2[deltat]+=(dist2[corrtime][i]*dist2[corrtime+deltat][i])/(stp+1-deltat)
#            cre2[deltat]+=N.dot(uv[corrtime][i],uv[corrtime+deltat][i])/(P*(stp+1-deltat))
#            r2file.write(str(deltat)+" "+str(cre2[deltat])+"\n")
print "centroid correlation function finished"
#print "avg correlation function finished"

#crfile.close()
#cpfile.close()
#crefile.close()
cafile.close()
#r2file.close()
#vfile.close()
trajectory.close()
