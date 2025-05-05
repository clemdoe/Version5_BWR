import numpy as np 
import matplotlib.pyplot as plt 


# Données brutes sous forme de chaîne de caractères

Chexal = """PCOOL:   11029811.0       11005952.0       10981161.0       10956115.0       10930795.0       10905178.0       10879240.0       10852953.0       10826412.0       10800000.0    
 VCOOL:   9.22811985       9.22851086       9.34786415       9.47372913       9.60673141       9.74761868       9.89730167       10.0567102       10.2270584       10.3620052    
 DCOOL:   758.577026       758.544006       748.859314       738.913330       728.684814       718.151611       707.286499       696.070007       684.468445       675.559021    
 TCOOL:   552.637878       552.636047       557.819214       562.900940       567.874451       572.731506       577.463074       582.054321       586.494385       589.699463    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000    
 TSAT:   591.434631       591.271973       591.102600       590.931274       590.757690       590.581787       590.403381       590.222168       590.038879       589.856201    
 MUT:   9.53041599E-05   9.52972259E-05   9.31402537E-05   9.10387462E-05   8.89841249E-05   8.69682481E-05   8.49830540E-05   8.30225690E-05   8.10788624E-05   7.96323147E-05
"""

GEramp = """PCOOL:   11029811.0       11005952.0       10981161.0       10956115.0       10930795.0       10905178.0       10879240.0       10852953.0       10826412.0       10800000.0    
 VCOOL:   9.22811985       9.22851086       9.34786415       9.47372913       9.60673141       9.74761868       9.89730167       10.0567102       10.2270584       10.3620052    
 DCOOL:   758.577026       758.544006       748.859314       738.913330       728.684814       718.151611       707.286499       696.070007       684.468445       675.559021    
 TCOOL:   552.637878       552.636047       557.819214       562.900940       567.874451       572.731506       577.463074       582.054321       586.494385       589.699463    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000    
 TSAT:   591.434631       591.271973       591.102600       590.931274       590.757690       590.581787       590.403381       590.222168       590.038879       589.856201    
 MUT:   9.53041599E-05   9.52972259E-05   9.31402537E-05   9.10387462E-05   8.89841249E-05   8.69682481E-05   8.49830540E-05   8.30225690E-05   8.10788624E-05   7.96323147E-05
"""

HEM1 = """PCOOL:   11029811.0       11005952.0       10981161.0       10956115.0       10930795.0       10905178.0       10879240.0       10852953.0       10826412.0       10800000.0    
 VCOOL:   9.22811985       9.22851086       9.34786415       9.47372913       9.60673141       9.74761868       9.89730167       10.0567102       10.2270584       10.3620052    
 DCOOL:   758.577026       758.544006       748.859314       738.913330       728.684814       718.151611       707.286499       696.070007       684.468445       675.559021    
 TCOOL:   552.637878       552.636047       557.819214       562.900940       567.874451       572.731506       577.463074       582.054321       586.494385       589.699463    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000    
 TSAT:   591.434631       591.271973       591.102600       590.931274       590.757690       590.581787       590.403381       590.222168       590.038879       589.856201    
 MUT:   9.53041599E-05   9.52972259E-05   9.31402537E-05   9.10387462E-05   8.89841249E-05   8.69682481E-05   8.49830540E-05   8.30225690E-05   8.10788624E-05   7.96323147E-05
"""

EPRI = """PCOOL:   11029811.0       11005952.0       10981161.0       10956115.0       10930795.0       10905178.0       10879240.0       10852953.0       10826412.0       10800000.0    
 VCOOL:   9.22811985       9.22851086       9.34786415       9.47372913       9.60673141       9.74761868       9.89730167       10.0567102       10.2270584       10.3620052    
 DCOOL:   758.577026       758.544006       748.859314       738.913330       728.684814       718.151611       707.286499       696.070007       684.468445       675.559021    
 TCOOL:   552.637878       552.636047       557.819214       562.900940       567.874451       572.731506       577.463074       582.054321       586.494385       589.699463    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000    
 TSAT:   591.434631       591.271973       591.102600       590.931274       590.757690       590.581787       590.403381       590.222168       590.038879       589.856201    
 MUT:   9.53041599E-05   9.52972259E-05   9.31402537E-05   9.10387462E-05   8.89841249E-05   8.69682481E-05   8.49830540E-05   8.30225690E-05   8.10788624E-05   7.96323147E-05
"""

ModBestion = """PCOOL:   11029811.0       11005952.0       10981161.0       10956115.0       10930795.0       10905178.0       10879240.0       10852953.0       10826412.0       10800000.0    
 VCOOL:   9.22811985       9.22851086       9.34786415       9.47372913       9.60673141       9.74761868       9.89730167       10.0567102       10.2270584       10.3620052    
 DCOOL:   758.577026       758.544006       748.859314       738.913330       728.684814       718.151611       707.286499       696.070007       684.468445       675.559021    
 TCOOL:   552.637878       552.636047       557.819214       562.900940       567.874451       572.731506       577.463074       582.054321       586.494385       589.699463    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000    
 TSAT:   591.434631       591.271973       591.102600       590.931274       590.757690       590.581787       590.403381       590.222168       590.038879       589.856201    
 MUT:   9.53041599E-05   9.52972259E-05   9.31402537E-05   9.10387462E-05   8.89841249E-05   8.69682481E-05   8.49830540E-05   8.30225690E-05   8.10788624E-05   7.96323147E-05
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
figp.savefig(rf'/home/p122173/Version5_BWR/Version5_ev3727/Donjon/data/results_marie/PCorrel12.png', bbox_inches='tight')

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
figx.savefig(rf'/home/p122173/Version5_BWR/Version5_ev3727/Donjon/data/results_marie/XCorrel12.png', bbox_inches='tight')

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
figeps.savefig(rf'/home/p122173/Version5_BWR/Version5_ev3727/Donjon/data/results_marie/EPSCorrel12.png', bbox_inches='tight')


