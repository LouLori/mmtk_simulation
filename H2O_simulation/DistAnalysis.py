from MMTK import *
from mbpol import mbpolForceField
from qTIP4pFF import HarmonicAngleForceField, HarmonicBondForceField, QuarticBondForceField, ElectrostaticForceField, LennardJonesForceField
from MMTK.ForceFields.ForceField import CompoundForceField
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
filename = traj[0:-3]
nthmol = int(argv[2])
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


afile = open("afile-"+str(nthmol+1)+"thmol"+filename+".dat", "w")

for step in trajectory:
    universe.setConfiguration(step['configuration'])
    xH1 = universe.objectList()[nthmol].atomList()[0].position()
    xH2 = universe.objectList()[nthmol].atomList()[1].position()
    xO  = universe.objectList()[nthmol].atomList()[2].position()
    mH1 = universe.objectList()[nthmol].atomList()[0].mass()
    mH2 = universe.objectList()[nthmol].atomList()[1].mass()
    mO = universe.objectList()[nthmol].atomList()[2].mass()
    xcom = (xH1*mH1+xH2*mH2+xO*mO)/(mH1+mH2+mO)
    rx = [1.,0.,0.]
    rz = [0.,0.,1.]
    dv = xcom-xO
    uv = dv/dv.length()
    uz = uv[2]
    afile.write(str(uz)+"\n")

trajectory.close()
afile.close()




