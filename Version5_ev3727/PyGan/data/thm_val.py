## PyGan script to generate thm procedure
# Date : 12/05/2025
# Author : Clément HUET, Raphaël GUASCH
# Purpose : test and validate  the thermalhydraulics on a single BWR pincell
from iapws import IAPWS97
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simps
import os, shutil
import lifo
import lcm
import cle2000
from assertS import *

import time

######## End helper functions declaration ##########

##
########## User input ##########
########## Algorithm parameters ###########
start_time = time.time()
solveConduction = True

########## Mesh parameters ###########
zPlotting = [] #If empty, no plotting of the axial distribution of the fields, otherwise, list of the axial positions where the fields are plotted
## Meshing parameters:
If = 8
I1 = 3

# Sensitivity to the meshing parameters
Iz1 = 10 # number of control volumes in the axial direction, added 70 for comparison with GeN-Foam
# Iz1 = 10, 20, 40, 50, 70, 80 and 160 are supported for the DONJON solution

PFiss = 10e3 # W, total fission power in the fuel rod, 40kW, 35kW, 30kW, 25kW, 20kW
# Test powers 40kW, 35kW, 30kW, 25kW, 20kW.
# Suspected limitation of pressure drop correlation at high powers in GeN-Foam comparison (21/12/2024)

########## Choice of Thermalhydraulics correlation ##########
voidFractionCorrel = 'EPRIvoidModel' # 'modBestion', 'HEM1', 'GEramp', 'EPRIvoidModel'
frfaccorel = "blasius" # 'base', 'blasius', 'Churchill', 
P2Pcorel = "lockhartMartinelli" # 'base', 'HEM1', 'HEM2', 'MNmodel', 'lockhartMartinelli'
numericalMethod = "BiCG" # "FVM": Solves the system using matrix inversion with preconditioning.
                         # "GaussSiedel" : Applies the Gauss-Seidel iterative solver.
                         # "BiCG" : Uses the BiConjugate Gradient method for solving non-symmetric or indefinite matrices.
                         # "BiCGStab" : Applies the BiCGStab (BiConjugate Gradient Stabilized) method to ensure faster and more stable convergence.

########## Thermal hydraulics parameters ##########
## Geometric parameters
canalType = "square" # "square", "cylindrical"
pitch = 1.295e-2 # m : ATRIUM10 pincell pitch
fuelRadius = 0.4435e-2 # m : fuel rod radius
gapRadius = 0.4520e-2 # m : expansion gap radius : "void" between fuel and clad - equivalent to inner clad radius
cladRadius = 0.5140e-2 # m : clad external radius
height = 1.555 # m : height : 3.8 m : active core height in BWRX-300 SMR, 1.555 m : for GeNFoam comparison.

area = pitch**2 - (np.pi * cladRadius**2) # m2 : cross section area of the fuel assembly

## Fluid parameters

# T_inlet, T_outlet = 270, 287 Celcius
#tInlet = 270 + 273.15 # K, for BWRX-300 SMR core, try lowering the inlet temperature to set boiling point back and reduce the void fraction increase in the first few cm
tInlet = 270 + 273.15 # K, for BWRX-300 SMR core
#Nominal operating pressure = 7.2 MPa (abs)
pOutlet =  7.2e6 # Pa 
# Nominal coolant flow rate = 1530 kg/s
massFlow = 1530  / (200*91)  # kg/s

## Material parameters
kFuel = 4.18 # W/m.K, TECHNICAL REPORTS SERIES No. 59 : Thermal Conductivity of Uranium Dioxide, IAEA, VIENNA, 1966
Hgap = 10000 

#Hgap = 9000
kClad = 21.5 #21.5 # W/m.K, Thermal Conductivity of Zircaloy-2 (as used in BWRX-300) according to https://www.matweb.com/search/datasheet.aspx?MatGUID=eb1dad5ce1ad4a1f9e92f86d5b44740d
# k_Zircaloy-4 = 21.6 W/m.K too so check for ATRIUM-10 clad material but should have the same thermal conductivity

############ Nuclear Parameters ###########
# Number of fuel rods and assemblies for a small modular Boiling Water Reactor core

n_rods = 91 # This value is for the ATRIUM-10 fuel assembly (10*10-9), in a GNF2 fuel assembly (BWRX-300) there are 92 fuel rods per assembly
n_assmblies = 240 # This value is for the BWRX-300 SMR core, ref : "Status Report – BWRX-300 (GE Hitachi and Hitachi GE Nuclear Energy)"
full_core_power = 870e6 # W, full core thermal power of the BWRX-300 SMR core, ref : "Status Report – BWRX-300 (GE Hitachi and Hitachi GE Nuclear Energy)"

# Total power of the fuel rod
power = 0.038 #MW, 0.14 MW for a single fuel rod

print(f"$$ - BEGIN Thermohydraulic case : Iz1 = {Iz1}, total power = {power} W")

TeffIni = 900 # K, initial effective temperature
TwaterIni = 270 + 273.15 # K, initial coolant temperature
rhoIni = 1000 # kg/m3, initial coolant density

def run_ThmDonjon(pitch, fuelRadius, gapRadius, cladRadius, height, power, tInlet, pOutlet, massFlow, area, TeffIni, TwaterIni, rhoIni):

  # 1) Create LCM object to store the TH data to be used in the neutronics solution
  print("Creating initial THData object")
  THData = lcm.new('LCM','THData')
  THData['TFuelList']    = np.array(TeffIni, dtype='f')
  THData['TCoolList'] = np.array(TwaterIni, dtype='f')
  THData['DCoolList'] = np.array(rhoIni/1000, dtype='f')
  THData.close() # close without erasing

  # 2) construct the Lifo stack for thm_val
  ipLifo1=lifo.new()
  ipLifo1.pushEmpty("Fmap", "LCM") # Fuel Map
  ipLifo1.pushEmpty("Matex", "LCM") # Material Indexation
  ipLifo1.pushEmpty("Cpo", "LCM") # Compo
  ipLifo1.pushEmpty("Track", "LCM") # Tracking data for FEM
  ipLifo1.pushEmpty("Thm", "LCM") # THM data empty
  ipLifo1.push(THData) # THData
  #ipLifo1.push(compo_name) # Compo name
  #ipLifo1.push(Fuel_power) # Mass of the fuel
  ipLifo1.push(power) # Total fission power in the fuel rod
  ipLifo1.push(height) # Height of the fuel rod
  ipLifo1.push(pitch) # pitch of the fuel assembly/cell
  ipLifo1.push(fuelRadius) # Fuel rod radius
  ipLifo1.push(gapRadius) # Gap radius
  ipLifo1.push(cladRadius) # Clad radius
  ipLifo1.push(tInlet) # Inlet temperature
  ipLifo1.push(pOutlet) # Outlet pressure
  ipLifo1.push(massFlow) # Mass flow rate
  ipLifo1.push(area) # Cross section area

  # 3) call thm_val Cle-2000 procedure
  thm_val = cle2000.new('thm_val',ipLifo1,1)
  thm_val.exec()
  print("thm_validation execution completed")

  # 4) recover the output LCM objects
  Fmap = ipLifo1.node("Fmap")
  Matex = ipLifo1.node("Matex")
  Cpo = ipLifo1.node("Cpo")
  Track = ipLifo1.node("Track")
  stateVector = Fmap["STATE-VECTOR"]
  mylength = stateVector[0]*stateVector[1]
  npar = stateVector[7]

  # 5) recover thermo-hydraulics information
  # 5.1) recover the THM data
  Thm = ipLifo1.node("Thm")
  THData = ipLifo1.node("THData")

  stationary_info = Thm["HISTORY-DATA"]["TIMESTEP0000"]["CHANNEL"][0]
  print("stationary_info : ", stationary_info)
  stationary_info.lib();
  pinlet=stationary_info["PINLET"][0]
  print("stationary inlet pressure=", pinlet, "Pa")

  # 5.2) recover the usefull TH data
  PCOOL = stationary_info["PRESSURE"]
  VCOOL = stationary_info["VELOCITIES"]
  TCOOL = stationary_info["COOLANT-TEMP"]
  TFUEL = stationary_info["CENTER-TEMPS"]
  EPS = stationary_info["EPSILON"]
  XFL = stationary_info["XFL"]

  # empty the ipLifo1 Lifo stack
  while ipLifo1.getMax() > 0:
    ipLifo1.pop()

  return PCOOL, VCOOL, TCOOL, TFUEL, EPS, XFL, Fmap, Matex, Cpo, Track, Thm

PCOOL, VCOOL, TCOOL, TFUEL, EPS, XFL, Fmap, Matex, Cpo, Track, Thm = run_ThmDonjon(pitch, fuelRadius, gapRadius, cladRadius, height, power, tInlet, pOutlet, massFlow, area, TeffIni, TwaterIni, rhoIni)

print("PCOOL : ", PCOOL)
print("VCOOL : ", VCOOL)
print("TCOOL : ", TCOOL)
print("TFUEL : ", TFUEL)
print("EPS : ", EPS)
print("XFL : ", XFL)