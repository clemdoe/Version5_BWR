import numpy as np 
import matplotlib.pyplot as plt 


# Données brutes sous forme de chaîne de caractères
Chexal = """PCOOL:   11183776.0       11159890.0       11134740.0       11109219.0       11083284.0       11056887.0       11029974.0       11002793.0       10943532.0       10800000.0    
 VCOOL:   9.24485111       9.24525547       9.40664387       9.58006287       9.76711273       9.96990490       10.1909695       10.3948908       14.3582106       17.5827770    
 DCOOL:   757.193542       757.161194       744.170898       730.701599       716.706970       702.124023       686.896057       673.466370       487.471008       398.117676    
 TCOOL:   553.518188       553.515747       560.382385       567.064636       573.542419       579.795105       585.791565       590.651489       590.845032       589.856201    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000      0.303399742      0.451105356    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       4.63160388E-02   8.65294188E-02
 TSAT:   592.477844       592.316772       592.146851       591.974060       591.798218       591.618896       591.435730       591.250366       590.845032       589.856201    
 MUT:   9.49877649E-05   9.49810346E-05   9.21364044E-05   8.93829929E-05   8.67021954E-05   8.40749562E-05   8.14856103E-05   7.93053841E-05   7.92316278E-05   7.96065724E-05
"""

GEramp = """PCOOL:   11186030.0       11162144.0       11136994.0       11111473.0       11085539.0       11059143.0       11032230.0       11005046.0       10946355.0       10800000.0    
 VCOOL:   9.24481964       9.24522591       9.40661812       9.58003426       9.76709080       9.96985626       10.1908588       10.3950062       14.2867613       17.7505283    
 DCOOL:   757.196289       757.164246       744.173706       730.704651       716.709839       702.126648       686.898010       673.449219       489.982910       394.563080    
 TCOOL:   553.518616       553.515991       560.382874       567.065247       573.543335       579.796265       585.793274       590.659912       590.864380       589.856201    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000      0.299242347      0.456894219    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       4.62049134E-02   8.82447287E-02
 TSAT:   592.493042       592.331970       592.162048       591.989380       591.813538       591.634277       591.451111       591.265747       590.864380       589.856201    
 MUT:   9.49883179E-05   9.49816676E-05   9.21369719E-05   8.93835677E-05   8.67027047E-05   8.40754365E-05   8.14859377E-05   7.93027357E-05   7.92242790E-05   7.96065724E-05
"""

HEM1 = """PCOOL:   11196972.0       11173086.0       11147936.0       11122415.0       11096480.0       11070085.0       11043175.0       11015988.0       10950382.0       10800000.0    
 VCOOL:   9.24464512       9.24503899       9.40647697       9.57993984       9.76704407       9.96988201       10.1909018       10.3963280       15.1034250       18.7330685    
 DCOOL:   757.212341       757.178711       744.188232       730.721924       716.725952       702.136475       686.898804       673.354004       462.879333       378.703979    
 TCOOL:   553.519043       553.517273       560.385071       567.067017       573.546570       579.803284       585.804626       590.704834       590.891968       589.856201    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000      0.343578786      0.482721567    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       4.61407453E-02   7.78478980E-02
 TSAT:   592.566772       592.405823       592.236023       592.063477       591.887756       591.708618       591.525635       591.340393       590.891968       589.856201    
 MUT:   9.49916866E-05   9.49846799E-05   9.21398314E-05   8.93868491E-05   8.67056951E-05   8.40771900E-05   8.14861851E-05   7.92880164E-05   7.92138089E-05   7.96065724E-05
"""

EPRI = """PCOOL:   11177021.0       11153134.0       11127983.0       11102461.0       11076525.0       11050124.0       11023207.0       10996041.0       10940693.0       10800000.0    
 VCOOL:   9.24496555       9.24537086       9.40686035       9.58036327       9.76749992       9.97038746       10.1916122       10.3945007       13.8703117       16.8903637    
 DCOOL:   757.184448       757.152100       744.159851       730.691406       716.690857       702.098267       686.857422       673.498596       504.618378       416.974915    
 TCOOL:   553.517456       553.515015       560.382080       567.063416       573.543152       579.798950       585.799011       590.633057       590.825562       589.856201    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000      0.275364488      0.420395344    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       4.64751050E-02   9.08466727E-02
 TSAT:   592.432312       592.271118       592.101135       591.928284       591.752319       591.572937       591.389648       591.204285       590.825562       589.856201    
 MUT:   9.49858804E-05   9.49791502E-05   9.21342071E-05   8.93810138E-05   8.66992195E-05   8.40704088E-05   8.14791274E-05   7.93102372E-05   7.92390128E-05   7.96065724E-05
"""

ModBestion = """PCOOL:   11183207.0       11159320.0       11134171.0       11108650.0       11082715.0       11056318.0       11029405.0       11002225.0       10943162.0       10800000.0    
 VCOOL:   9.24486256       9.24526596       9.40668869       9.58010578       9.76714325       9.96992970       10.1909838       10.3947210       14.3316851       17.5287781    
 DCOOL:   757.192993       757.159790       744.169983       730.701050       716.704651       702.119385       686.891174       673.469849       488.430573       399.512909    
 TCOOL:   553.518005       553.516052       560.382324       567.064392       573.542847       579.796448       585.792847       590.649719       590.842529       589.856201    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000      0.301833928      0.448833168    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       4.63093556E-02   8.68076161E-02
 TSAT:   592.474060       592.312927       592.143005       591.970215       591.794373       591.615051       591.431885       591.246521       590.842529       589.856201    
 MUT:   9.49876558E-05   9.49807218E-05   9.21362298E-05   8.93828837E-05   8.67017952E-05   8.40741413E-05   8.14847808E-05   7.93058935E-05   7.92325736E-05   7.96065724E-05
"""



# Traitement de chaque ligne pour extraire les valeurs numériques
def extraire(var_name, line_nb, data_name):
    data_lines = data_name.strip().split("\n")
    values = data_lines[line_nb].split(":")[1]
    # Séparation des valeurs en une liste de flottants
    for value in values.split():
        var_name.append(float(value))
    return

# Affichage du tableau des valeurs de pression triées
PHEM1=[]
PMB=[]
PEPRI=[]
PGER=[]
PCHEXAL=[]
extraire(PHEM1, 0, HEM1)
extraire(PMB, 0, ModBestion)
extraire(PEPRI, 0, EPRI)
extraire(PGER, 0, GEramp)
extraire(PCHEXAL, 0, Chexal)

"""extraire(VCOOL, 1)
extraire(DCOOL, 2)
extraire(TCOOL, 3)
extraire(EPS, 4)
extraire(XFL, 5)
extraire(TSAT, 6)
extraire(MUT, 7)"""
Zaxis=np.linspace(0,2,len(PHEM1))

figp, axp = plt.subplots()
axp.plot(Zaxis,PCHEXAL, label='Chexal')
axp.plot(Zaxis,PGER, label='GEramp')
axp.plot(Zaxis,PEPRI, label='EPRI')
axp.plot(Zaxis,PMB, label='ModBestion')
axp.plot(Zaxis,PHEM1, label='HEM1')
axp.set_xlabel("Axial position in m")
axp.set_ylabel("Pressure in Pa")
axp.set_title("Comparaison of different correlations")
axp.grid()
axp.legend(loc="best")
plt.show()
figp.savefig(rf'/home/p122173/Version5_BWR/Version5_ev3727/Donjon/data/results_marie/PCorrel.png', bbox_inches='tight')

# Affichage du tableau des valeurs de titre triées
XHEM1=[]
XMB=[]
XEPRI=[]
XGER=[]
XCHEXAL=[]
extraire(XHEM1, 5, HEM1)
extraire(XMB, 5, ModBestion)
extraire(XEPRI, 5, EPRI)
extraire(XGER, 5, GEramp)
extraire(XCHEXAL, 5, Chexal)

Zaxis=np.linspace(0,2,len(XHEM1))

figx, axx = plt.subplots()
axx.plot(Zaxis,XCHEXAL, label='Chexal')  
axx.plot(Zaxis,XGER, label='GEramp')
axx.plot(Zaxis,XEPRI, label='EPRI')
axx.plot(Zaxis,XMB, label='ModBestion')
axx.plot(Zaxis,XHEM1, label='HEM1')
axx.set_xlabel("Axial position in m")
axx.set_ylabel("Quality")
axx.set_title("Comparaison of different correlations")
axx.grid()
axx.legend(loc="best")
plt.show()
figx.savefig(rf'/home/p122173/Version5_BWR/Version5_ev3727/Donjon/data/results_marie/XCorrel.png', bbox_inches='tight')

# Affichage du tableau des valeurs de taux de vide triées
EPSHEM1=[]
EPSMB=[]   
EPSPEPRI=[]
EPSGER=[]
EPSCHEXAL=[]
extraire(EPSHEM1, 4, HEM1)
extraire(EPSMB, 4, ModBestion)
extraire(EPSPEPRI, 4, EPRI) 
extraire(EPSGER, 4, GEramp)
extraire(EPSCHEXAL, 4, Chexal)
Zaxis=np.linspace(0,2,len(EPSHEM1))
figeps, axeps = plt.subplots()
axeps.plot(Zaxis,EPSCHEXAL, label='Chexal')
axeps.plot(Zaxis,EPSGER, label='GEramp')
axeps.plot(Zaxis,EPSPEPRI, label='EPRI')
axeps.plot(Zaxis,EPSMB, label='ModBestion')
axeps.plot(Zaxis,EPSHEM1, label='HEM1')
axeps.set_xlabel("Axial position in m")
axeps.set_ylabel("Void fraction")
axeps.set_title("Comparaison of different correlations")
axeps.grid()
axeps.legend(loc="best")
plt.show()
figeps.savefig(rf'/home/p122173/Version5_BWR/Version5_ev3727/Donjon/data/results_marie/EPSCorrel.png', bbox_inches='tight')


