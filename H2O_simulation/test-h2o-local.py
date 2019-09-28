from MMTK import *
import sys
from MMTK.ForceFields.Amber import *
from MMTK.ForceFields.ForceField import CompoundForceField
from MMTK.ForceFields.ForceFieldTest import gradientTest
from MMTK import Features
#from RigidRotor_PINormalModeIntegrator import RigidRotor_PINormalModeIntegrator, RigidRotor_PILangevinNormalModeIntegrator
#from MMTK.Rigid3DRotor_PINormalModeIntegrator import Rigid3DRotor_PINormalModeIntegrator, Rigid3DRotor_PILangevinNormalModeIntegrator
#from MMTK.Rot3DOnly_PINormalModeIntegrator_test import Rot3DOnly_PINormalModeIntegrator, Rot3DOnly_PILangevinNormalModeIntegrator
from MMTK.Rot3DOnly_PINormalModeIntegrator import Rot3DOnly_PINormalModeIntegrator, Rot3DOnly_PILangevinNormalModeIntegrator
#from RotOnly_PINormalModeIntegrator import RotOnly_PINormalModeIntegrator, RotOnly_PILangevinNormalModeIntegrator
#from MMTK.ForceFields.Restraints import HarmonicDistanceRestraint
from CoulombFF import CoulombForceField
from LJFF import LJForceField
from qTIP4pFF import HarmonicAngleForceField, HarmonicBondForceField, QuarticBondForceField, ElectrostaticForceField, LennardJonesForceField
#from dipoleFF import dipoleForceField
#from HarmonicWellFF import HarmonicWellForceField
from mbpol import mbpolForceField
from HarmonicOscillatorFF import HarmonicOscillatorForceField
from HOFF import HOForceField
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

############################
# set the number of quantum beads
############################
# directory and filename
outdir="/work/x3lou/MMTK/H2OSimulations/LocalUpdates/q-TIP4P/N2-R0.3/"
updatetype = "local"
doftype = "rot"
#doftype = "trans"
#fftype = "q-spc-fw"
fftype = "r-tip4p"
#fftype = "mbpol"
label= updatetype + "-" + doftype + "-" + fftype
testnum = "test"
#testnum = "fun"

# Parameters
rhoname = argv[1]
erotname = argv[2]
esqname = argv[3]
rotstepval=argv[4]
print "ROT STEP:", rotstepval

temperature = float(rhoname[rhoname.find("T")+1:rhoname.find("P")])*Units.K  # temperature
P = int(rhoname[rhoname.find("P")+1:rhoname.find(".rho")])                    # number of beads
lattice_spacing=3.0*Units.Ang
print "Lattice spacing:", lattice_spacing, " nm"

ndens=23588101
rho=zeros(ndens,float)
erot=zeros(ndens,float)
esq=zeros(ndens,float)
universe = InfiniteUniverse()
# nanometers

universe.addObject(PathIntegrals(temperature))

nH2O = int(argv[5])
print "Number of Water:", nH2O
nobjects = nH2O

#center = zeros( (nobjects,3) , float)

#for i in range(nobjects):
#    universe.addObject(Molecule('water', position = Vector(0.,0.,i*lattice_spacing)))

## They will initially be aligned along the z-axis
for i in range(nobjects):
    universe.addObject(Molecule('water', position = Vector(0.,0.,i*lattice_spacing)))
#print universe.objectList()[1].atomList()[2].beadPositions()

for atom in universe.atomList():
	atom.setNumberOfBeads(P)

#print "ATOMS"
print  universe.atomList()
natoms = len(universe.atomList())
print universe.objectList()[0].mass()
print universe.objectList()[1].mass()
print universe.atomList()[0].mass()
print universe.atomList()[1].mass()
print universe.atomList()[2].mass()
#print universe.atomList()[3].mass()


######################  FORCEFIELD SETUP  ######################
ff=[]

###############  WATER MOLECULE PARAMETER LIST  ################

if fftype == "q-tip4p" or fftype == "r-tip4p":
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

################# CONSTRUCT FORCEFIELD ###############################
if fftype != "mbpol":
    for o in range(nobjects):
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

##################################################
############## COULOMB POTENTIAL #################
##################################################
#ff.append(NoPotForceField(universe.atomList()[0],universe.atomList()[1]))

#for i in range(nobjects):
#    for ia in range(len(universe.objectList()[i].atomList())):
#        for j in range(i+1,nobjects):
#            for ja in range(len(universe.objectList()[j].atomList())):
#                ff.append(CoulombForceField(universe.objectList()[i].atomList()[ia],universe.objectList()[j].atomList()[ja],ia,ja))
#
#for i in range(nobjects):
#    for j in range(i+1,nobjects):
#        ff.append(LJForceField(universe.objectList()[i].atomList()[2],universe.objectList()[j].atomList()[2],2,2))


##################################################
############## DIPOLE POTENTIAL #################
##################################################
#for i in range(nobjects):
#        for j in range(i+1,nobjects):
#                ff.append(dipoleForceField(universe.objectList()[i].atomList()[0],universe.objectList()[i].atomList()[1], universe.objectList()[i].atomList()[2],
#                                           universe.objectList()[j].atomList()[0],universe.objectList()[j].atomList()[1],universe.objectList()[j].atomList()[2]))
#

##################################################
############## HARMONIC WELL #####################
##################################################
#for i in range(nobjects):
#    for j in range(len(universe.objectList()[i].atomList())):
#        ff.append(HarmonicWellForceField(universe.objectList()[i].atomList()[j], Vector(0.,0.,0.), 100.))
#for i in range(nobjects):
#    for j in range(len(universe.objectList()[i].atomList())):
#        ff.append(HarmonicOscillatorForceField(universe.objectList()[i].atomList()[j], Vector(i*lattice_spacing,0.,0.), 100.))

########################################################
##############     COM RIGID BODY      #################
############## HARMONIC TRAP POTENTIAL #################
########################################################
#forceconst = 100.0
#print "Force Constant:", forceconst
#
#for i in range(nobjects):
#        ff.append(HOForceField(universe.objectList()[i].atomList()[0], universe.objectList()[i].atomList()[1], universe.objectList()[i].atomList()[2], forceconst))

###################################################
################### MB-pol Model ##################
###################################################
#ff.append(Amber99ForceField(mod_file='glycam_06g_spcfw_q.dat'))


universe.setForceField(CompoundForceField(*ff))
###################################################


#print 'ATOM LIST'
#for i in range(natoms):
#    print universe.atomList()[i].beadPositions()
#print 'END ATOM LIST \n'

#gradientTest(universe)
print "Initial Universe Energy"
print universe.energy()/Units.k_B, "K"
#print "last energy values"
#print universe.energyEvaluator().CEvaluator().last_energy_values
#
#raise()

universe.writeToFile("u.pdb")
#universe.environmentObjectList(Environment.PathIntegrals)[0].include_spring_terms = False
#universe._changed(True)

rhofile=open(rhoname,"r")
erotfile=open(erotname, "r")
esqfile=open(esqname, "r")

for i in range(ndens):
        rho[i]=float(rhofile.readline())
        erot[i]=float(erotfile.readline())*0.0083144621 #K to kJ/mol [MMTK Units of Energy]
        esq[i]=float(esqfile.readline())*0.0083144621*0.0083144621
rhofile.close()
erotfile.close()
esqfile.close()

dt = 1.0*Units.fs
print "Using a timestep of ", dt, " fs"

# Initialize velocities
universe.environmentObjectList(Environment.PathIntegrals)[0].include_spring_terms = False
universe._changed(True)
universe.initializeVelocitiesToTemperature(temperature)


# USE THE FRICTION PARAMETER FROM BEFORE
friction = 0.0
#integrator = Rigid3DRotor_PILangevinNormalModeIntegrator(universe, delta_t=dt, centroid_friction = friction, denrho = rho, denerot = erot, denesq = esq, rotstep=float(rotstepval))
integrator = Rot3DOnly_PILangevinNormalModeIntegrator(universe, delta_t=dt, centroid_friction = friction, denrho = rho, denerot = erot, denesq = esq, rotstep=float(rotstepval))
#integrator = Rot3DOnly_PILangevinNormalModeIntegrator(universe, delta_t=dt, centroid_friction = friction, denrho = rho, denerot = erot, denesq = esq, rotstep=float(rotstepval))

integrator(steps=1000, actions = [ TrajectoryOutput(None,('configuration','time'), 0, None, 100)] )
#integrator(steps=100, actions = [ TrajectoryOutput(None,('configuration','time'), 0, None, 10)] )
print "Preskip Universe Energy"
print universe.energy()/Units.k_B, "K"

#raise()

#RunSteps = 1000.0
#SkipSteps = 100.0

RunSteps = 1000.0*Units.ps/dt
SkipSteps = 100.0*Units.fs/dt

#RunSteps = 1000.0*Units.ps/dt
#SkipSteps = 100.0*Units.fs/dt

trajectory = Trajectory(universe, outdir+str(nH2O)+"H2O-P"+str(P)+"-"+"T"+str(temperature)+"-"+"R"+str(lattice_spacing)+"-"+label+"-"+str(rotstepval)+"-"+testnum+".nc", "w", "A simple test case")
Nblocks=1


############################## BEGIN TRANSLATION/ROTATION SIMULATION ##############################
# RUN PIMD WITH PIMC ROTATION INCLUDED
print "We're going to run the Langevin integrator for ", RunSteps/SkipSteps, "independent steps of PIMD"
integrator(steps=RunSteps, actions = [TrajectoryOutput(trajectory, ("time", "thermodynamic","configuration","energy", "auxiliary"), 0, None, SkipSteps)])
####################################################################################################

#gradientTest(universe)
trajectory.close()

