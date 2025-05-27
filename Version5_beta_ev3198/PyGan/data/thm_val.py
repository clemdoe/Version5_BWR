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
import pandas as pd

import time

######## End helper functions declaration ##########

########## User input ##########

########## Thermal hydraulics parameters ##########
## Geometric parameters
canalType = "square"
pitch = 0.0126 #m0.0133409 #
fuelRadius = 0.0027115493728018247#0.00467325 #
gapRadius = fuelRadius + 0.00000001  # m : expansion gap radius : "void" between fuel and clad - equivalent to inner clad radius
cladRadius =  0.0094996/2 #0.00474998 # m : clad external radius
#cladRadius = 0.0133221/2 # m : clad external radius
height = 1.555 # m :  : active core height in BWRX-300 SMR
area = pitch**2 - np.pi * cladRadius**2


## Material parameters
kFuel = 4.18 # W/m.K, TECHNICAL REPORTS SERIES No. 59 : Thermal Conductivity of Uranium Dioxide, IAEA, VIENNA, 1966
Hgap = 10000
kClad = 21.5 # W/m.K, Thermal Conductivity of Zircaloy-2 (as used in BWRX-300) according to https://www.matweb.com/search/datasheet.aspx?MatGUID=eb1dad5ce1ad4a1f9e92f86d5b44740d
# k_Zircaloy-4 = 21.6 W/m.K too so check for ATRIUM-10 clad material but should have the same thermal conductivity
########## Algorithm parameters ###########
nIter = 1000
tol = 1e-4

# Charger le fichier Excel (remplace 'chemin_du_fichier.xlsx' par le chemin réel)
chemin_fichier = '/home/clhue/Version5_BWR/Version5_ev3727/PyGan/data/BFBT.xlsx'
def chargement_excel(chemin_fichier):
  df = pd.read_excel(chemin_fichier)

  # Extraire les 7 premières colonnes
  colonnes = df.columns[:7]  # Sélectionne les noms des 7 premières colonnes

  # Mettre chaque colonne dans une liste
  listes = [df[colonne].tolist() for colonne in colonnes]

  powerList = listes[3]
  pressureList = listes[1]
  temperatureList = listes[4]
  flowList = listes[2]
  testNumber = listes[0]
  densityResult = listes[5]
  voidResult = listes[6]

  return powerList, pressureList, temperatureList, flowList, testNumber, densityResult, voidResult


powerList, pressureList, temperatureList, flowList, testNumber, densityResult, voidResult = chargement_excel(chemin_fichier)
print(f'powerList: {powerList} kW')
print(f'pressureList: {pressureList} bar')
print(f'temperatureList: {temperatureList} °C')
print(f'flowList: {flowList} m3/h')
print(f'testNumber: {testNumber}')


""" global ipLifo1
ipLifo1=lifo.new()
# 0) Create the Lifo stack for the THM validation
global Fmap
global Matex
global Cpo
global Track
global Thm """

ipLifo1=lifo.new()

rhoOutNum = []
epsOutNum = []

for i in range(len(powerList)):
  testBFBT = testNumber[i]
  testN = i+1
  pOutlet = pressureList[i] * 98066.5 # Pa
  tInlet = temperatureList[i] + 273.15 # K
  massFlow = flowList[i] *0.000107098*1000000/3600 # kg/s
  powi = powerList[i] / 1000  #/ (1.555*np.pi*(0.00271154937280182)**2) #MW

  # Print the parameters for the current test
  print(f"Numero de test: {testN}, numero de test BFBT: {testBFBT}")
  print(f"Power: {powi} MW, Power: {powi*1000} kW, Pressure: {pOutlet} Pa, Inlet Temperature: {tInlet} K, Mass Flow: {massFlow} kg/s")

  print("powi = ", powi)
  print("height = ", height)
  print("pitch = ", pitch)
  print("fuelRadius = ", fuelRadius)
  print("gapRadius = ", gapRadius)
  print("cladRadius = ", cladRadius)
  print("tInlet = ", tInlet)
  print("pOutlet = ", pOutlet)
  print("massFlow = ", massFlow)
  print("area = ", area)


  if testN == 1:
    print("testN = 1")
    print("Creating Lifo stack for thm_val")
    #ipLifo1.pushEmpty("Fmap", "LCM") # Fuel Map
    #ipLifo1.pushEmpty("Matex", "LCM") # Material Indexation
    ipLifo1.pushEmpty("Cpo", "LCM") # Compo
    #ipLifo1.pushEmpty("Track", "LCM") # Tracking data for FEM
    ipLifo1.pushEmpty("Thm", "LCM") # THM data empty
    #ipLifo1.push(THData) # THData

  else:
    print("testN = ", testN)
    print("updating Lifo stack for thm_val")
    #ipLifo1.push(Fmap) # Fuel Map
    #ipLifo1.push(Matex) # Material Indexation
    ipLifo1.push(Cpo) # Compo
    #ipLifo1.push(Track) # Tracking data for FEM
    ipLifo1.push(Thm) # THM data empty
    #ipLifo1.push(THData) # THData

  #ipLifo1.push(compo_name) # Compo name
  #ipLifo1.push(Fuel_power) # Mass of the fuel

  ipLifo1.push(int(testN))
  ipLifo1.push(powi*0.6) # Total fission power in the fuel rod
  print("powi*0.8 = ", powi*0.5)
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
  #Fmap = ipLifo1.node("Fmap")
  #Matex = ipLifo1.node("Matex")
  Cpo = ipLifo1.node("Cpo")
  #Track = ipLifo1.node("Track")

  # 5) recover thermo-hydraulics information
  # 5.1) recover the THM data
  Thm = ipLifo1.node("Thm")
  #THData = ipLifo1.node("THData")

  stationary_info = Thm["HISTORY-DATA"]["TIMESTEP0000"]["CHANNEL"][0]
  print("stationary_info : ", stationary_info)
  print("stationary_info keys : ", stationary_info.keys())
  print("stationary_info lib : ", stationary_info.lib())
  stationary_info.lib();

  # 5.2) recover the usefull TH data
  PCOOL = stationary_info["PRESSURE"]
  VCOOL = stationary_info["VELOCITIES"]
  TCOOL = stationary_info["COOLANT-TEMP"]
  TFUEL = stationary_info["CENTER-TEMPS"]
  RHO = stationary_info["DENSITY"]
  EPS = stationary_info["EPSILON"]
  XFL = stationary_info["XFL"]

  rhoOutNum.append(RHO[-1])
  epsOutNum.append(EPS[-1])

  # empty the ipLifo1 Lifo stack
  while ipLifo1.getMax() > 0:
    ipLifo1.pop()

  # save the results in a dictionnary
  print("PCOOL : ", PCOOL)
  print("VCOOL : ", VCOOL)
  print("TCOOL : ", TCOOL)
  print("TFUEL : ", TFUEL)
  print("EPS : ", EPS)
  print("XFL : ", XFL)

print("rhoOutNum : ", rhoOutNum)
print("epsOutNum : ", epsOutNum)