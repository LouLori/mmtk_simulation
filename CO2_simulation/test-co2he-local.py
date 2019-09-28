from MMTK import *
import sys
from MMTK.ForceFields.ForceField import CompoundForceField
from MMTK.ForceFields.ForceFieldTest import gradientTest
from MMTK import Features
from MMTK.Rigid2DRotor_PINormalModeIntegrator import Rigid2DRotor_PINormalModeIntegrator, Rigid2DRotor_PILangevinNormalModeIntegrator
#from MMTK.Rot2DOnly_PINormalModeIntegrator import Rot2DOnly_PINormalModeIntegrator, Rot2DOnly_PILangevinNormalModeIntegrator
#from MMTK.ForceFields.Restraints import HarmonicDistanceRestraint
#from CoulombFF import CoulombForceField
#from dipoleFF import dipoleForceField
#from NoPotFF import NoPotForceField
from HeCO2TransFF import HeCO2TransForceField
from ExactHeCO2TransFF import ExactHeCO2TransForceField
from HeHeFF import HeHeForceField
from MMTK.Environment import PathIntegrals
from MMTK.NormalModes import VibrationalModes
from MMTK.Trajectory import Trajectory, TrajectoryOutput, \
                            RestartTrajectoryOutput, StandardLogOutput, \
                            trajectoryInfo
from sys import argv,exit
from Scientific.Statistics import mean, standardDeviation
#from nMOLDYN.Mathematics.Analysis import correlation
from MMTK.Minimization import ConjugateGradientMinimizer
from Scientific import N
from numpy import zeros,cos,sin,sqrt,pi, dot, asarray, sign, arctan, random
import subprocess
import re

#################################
# set the number of quantum beads
#################################

outdir="/work/x3lou/MMTK/CO2Simulations/LocalUpdates/Trans_Rot/"
#outdir="/work/x3lou/MMTK/CO2Simulations/LocalUpdates/Trans/MAY282019/"
#outdir="/work/x3lou/MMTK/CO2Simulations/LocalUpdates/Rot/MAY232019/"
updatetype = "local"
testtype = "trans-rot"
#testtype = "rot"
#testtype = "trans"
label= updatetype + "-" + testtype
testnum = "test"

######### Parameters ##########
densname = argv[1]
temperature = float(densname[densname.find("T")+1:densname.find("P")])*Units.K  # temperature
P = int(densname[densname.find("P")+1:densname.find(".dat")])                    # number of beads

lattice_spacing=10.05*Units.Ang
ndensn = int(densname[densname.find("n")+1:densname.find("e")])
ntens = int(densname[densname.find("e")+1:densname.find("T")])
ndens=ndensn*(10**ntens)
print ndens

rho=zeros(ndens,float)
roteng=zeros(ndens,float)

universe = InfiniteUniverse()
# nanometers

universe.addObject(PathIntegrals(temperature))

nCO2=1
nhelium=int(argv[2])
nobjects=nCO2+nhelium

rotstepval=float(argv[3])
print "ROT_STEP:", rotstepval

## CO2 will initially be at origin point
universe.addObject(Molecule('co2v', position = Vector(0., 0., 0.)))

for i in range(nhelium):
    if i < 5:
        apos=Vector([random.uniform(0,1)/100.,.306*cos(i*2.*pi/5.),.306*sin(i*2.*pi/5.)])
    else:
        if i < 10:
            apos=Vector([0.35+random.uniform(0,1)/100.,.306*cos(i*2.*pi/5.),.306*sin(i*2.*pi/5.)])
        else:
            if i <15:
                apos=Vector([-0.35+random.uniform(0,1)/100.,.306*cos(i*2.*pi/5.),.306*sin(i*2.*pi/5.)])
            else:
                apos=Vector([random.uniform(0,1)/100.,2.*.306*cos(i*2.*pi/10.),2.*.306*sin(i*2.*pi/10.)])

    universe.addObject(Atom('he', position=apos))
    #universe.addObject(Atom('he', position=Vector([4.0,4.0,4.0])*Units.Ang))


for atom in universe.atomList():
	atom.setNumberOfBeads(P)

print "ATOMS"
natoms = len(universe.atomList())

for i in range(natoms):
    print universe.atomList()[i], universe.atomList()[i].mass()

ff=[]

################################i##########################
############## He-CO2 and He-He POTENTIAL #################
###########################################################
for i in range(3*nCO2,nhelium+3*nCO2):
        ## Construct He-CO2 force field on CO2 with every He
        ff.append(ExactHeCO2TransForceField(universe.atomList()[i],universe.atomList()[1],universe.atomList()[0]))
        print "He index:", i, ", CO2 index:",0
        ## Construct He-He force field for each He
        for j in range(i+1,nhelium+3*nCO2):
                ff.append(HeHeForceField(universe.atomList()[i],universe.atomList()[j]))
                print "Pair He index:", i,j

universe.setForceField(CompoundForceField(*ff))

universe.writeToFile("u.pdb")
universe.initializeVelocitiesToTemperature(temperature)

#raise()
##################################################
############ READ ROT DENSITY FILE ###############
##################################################

densfile=open(densname,"r")
for i in range(ndens):
        dummy=densfile.readline().split()
        rho[i]=float(dummy[1])
        roteng[i]=float(dummy[2])*0.0083144621 #K to kJ/mol [MMTK Units of Energy]
densfile.close()

## SET TIME STEP
dt = 5.0*Units.fs
print "Using a time step of ", dt*1000., " fs"

# Initialize velocities
universe.initializeVelocitiesToTemperature(temperature)

# USE THE FRICTION PARAMETER FROM BEFORE
friction = 0.0
integrator = Rigid2DRotor_PILangevinNormalModeIntegrator(universe, delta_t=dt, centroid_friction = friction, densmat=rho, rotengmat=roteng, rotstep=float(rotstepval))
#integrator = Rot2DOnly_PILangevinNormalModeIntegrator(universe, delta_t=dt, centroid_friction = friction, densmat=rho,rotengmat=roteng, rotstep=float(rotstepval))

#integrator(steps=1000, actions = [ TrajectoryOutput(None,('configuration','time'), 0, None, 100)])
#raise()

RunSteps = 500.0*Units.ps/dt
print "RunSteps:", RunSteps
SkipSteps = 50.0*Units.fs/dt
print "SkipSteps:", SkipSteps

#trajectory = Trajectory(universe, outdir+str(nCO2)+"CO2He-P"+str(P)+"-"+"T"+str(temperature)+"-"+label+"-"+testnum+".nc", "w", "A simple test case")
trajectory = Trajectory(universe, outdir+str(nCO2)+"CO2"+"-"+str(nhelium)+"He-P"+str(P)+"-"+"T"+str(temperature)+"-"+label+"-"+testnum+".nc", "w", "A simple test case")
Nblocks=1

################################################################################################
########################### BEGIN TRANSLATION/ROTATION SIMULATION ##############################
################################################################################################

# RUN PIMD WITH PIMC ROTATION INCLUDED
print "We're going to run the Langevin integrator for ", RunSteps/SkipSteps, "independent steps of PIMD"
integrator(steps=RunSteps,
           # Remove global translation every 50 steps.
	   actions = [
		   TrajectoryOutput(trajectory, ("time", "thermodynamic", "configuration","energy", "auxiliary"),
                                    0, None, SkipSteps)])

gradientTest(universe)
trajectory.close()
