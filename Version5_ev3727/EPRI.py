import numpy as np 
import matplotlib.pyplot as plt 


# Données brutes sous forme de chaîne de caractères
P12 = """PCOOL:   11029811.0       11005952.0       10981161.0       10956115.0       10930795.0       10905178.0       10879240.0       10852953.0       10826412.0       10800000.0    
 VCOOL:   9.22811985       9.22851086       9.34786415       9.47372913       9.60673141       9.74761868       9.89730167       10.0567102       10.2270584       10.3620052    
 DCOOL:   758.577026       758.544006       748.859314       738.913330       728.684814       718.151611       707.286499       696.070007       684.468445       675.559021    
 TCOOL:   552.637878       552.636047       557.819214       562.900940       567.874451       572.731506       577.463074       582.054321       586.494385       589.699463    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000    
 TSAT:   591.434631       591.271973       591.102600       590.931274       590.757690       590.581787       590.403381       590.222168       590.038879       589.856201    
 MUT:   9.53041599E-05   9.52972259E-05   9.31402537E-05   9.10387462E-05   8.89841249E-05   8.69682481E-05   8.49830540E-05   8.30225690E-05   8.10788624E-05   7.96323147E-05
"""

P13 = """PCOOL:   11051477.0       11027610.0       11002728.0       10977564.0       10952094.0       10926290.0       10900120.0       10873551.0       10846598.0       10800000.0    
 VCOOL:   9.23254490       9.23295784       9.36281776       9.50038719       9.64650440       9.80219746       9.96858120       10.1470337       10.3314219       12.8197107    
 DCOOL:   758.208679       758.175659       747.669556       736.853210       725.695740       714.167542       702.238831       689.879456       677.557007       546.339905    
 TCOOL:   552.856689       552.854797       558.459595       563.945190       569.304749       574.526917       579.596802       584.496765       589.034119       589.856201    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000      0.209717751    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       3.13557647E-02
 TSAT:   591.582092       591.419617       591.249939       591.078064       590.903748       590.726807       590.547058       590.364197       590.178345       589.856201    
 MUT:   9.52200353E-05   9.52131086E-05   9.28832669E-05   9.06170608E-05   8.84028341E-05   8.62305824E-05   8.40912980E-05   8.19767156E-05   7.99568661E-05   7.96065724E-05
"""

P14 = """PCOOL:   11112291.0       11088418.0       11063449.0       11038170.0       11012550.0       10986554.0       10960144.0       10933279.0       10893418.0       10800000.0    
 VCOOL:   9.23631096       9.23671818       9.37691021       9.52609730       9.68535233       9.85599804       10.0396614       10.2384386       11.9344473       14.2900839    
 DCOOL:   757.894165       757.861023       746.532410       734.843811       722.761719       710.247681       697.260376       683.758667       587.212036       487.920380    
 TCOOL:   553.078491       553.076660       559.103149       564.989746       570.725159       576.295593       581.683777       586.867798       590.500916       589.856201    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000      0.141141951      0.304856956    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       1.93976518E-02   5.30423708E-02
 TSAT:   591.994873       591.833069       591.663513       591.491577       591.316956       591.139465       590.958862       590.774719       590.500916       589.856201    
 MUT:   9.51473921E-05   9.51404363E-05   9.26385765E-05   9.02095198E-05   8.78399078E-05   8.55170729E-05   8.32292426E-05   8.09659396E-05   7.93621875E-05   7.96065724E-05
"""

P16 = """PCOOL:   11177021.0       11153134.0       11127983.0       11102461.0       11076525.0       11050124.0       11023207.0       10996041.0       10940693.0       10800000.0    
 VCOOL:   9.24496555       9.24537086       9.40686035       9.58036327       9.76749992       9.97038746       10.1916122       10.3945007       13.8703117       16.8903637    
 DCOOL:   757.184448       757.152100       744.159851       730.691406       716.690857       702.098267       686.857422       673.498596       504.618378       416.974915    
 TCOOL:   553.517456       553.515015       560.382080       567.063416       573.543152       579.798950       585.799011       590.633057       590.825562       589.856201    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000      0.275364488      0.420395344    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       4.64751050E-02   9.08466727E-02
 TSAT:   592.432312       592.271118       592.101135       591.928284       591.752319       591.572937       591.389648       591.204285       590.825562       589.856201    
 MUT:   9.49858804E-05   9.49791502E-05   9.21342071E-05   8.93810138E-05   8.66992195E-05   8.40704088E-05   8.14791274E-05   7.93102372E-05   7.92390128E-05   7.96065724E-05
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
PP12=[]
PP13=[]
PP14=[]
PP16=[]
extraire(PP12, 0, P12)      
extraire(PP13, 0, P13)
extraire(PP14, 0, P14)
extraire(PP16, 0, P16)
Zaxis=np.linspace(0,2,len(PP12))
figp, axp = plt.subplots()  
axp.plot(Zaxis,PP12, label='P12')
axp.plot(Zaxis,PP13, label='P13')
axp.plot(Zaxis,PP14, label='P14')
axp.plot(Zaxis,PP16, label='P16')
axp.set_xlabel("Axial position in m")   
axp.set_ylabel("Pressure in Pa")
axp.set_title("Comparaison of different powers for EPRI correlation")
axp.grid()
axp.legend(loc="best")
plt.show()
figp.savefig(rf'/home/p122173/Version5_BWR/Version5_ev3727/Donjon/data/results_marie/PEPRI.png', bbox_inches='tight')

# Affichage du tableau des valeurs de titre triées
XP12=[]
XP13=[]
XP14=[]
XP16=[]         
extraire(XP12, 1, P12)
extraire(XP13, 1, P13)
extraire(XP14, 1, P14)
extraire(XP16, 1, P16)
Zaxis=np.linspace(0,2,len(XP12))
figx, axx = plt.subplots()
axx.plot(Zaxis,XP12, label='P12')
axx.plot(Zaxis,XP13, label='P13')
axx.plot(Zaxis,XP14, label='P14')
axx.plot(Zaxis,XP16, label='P16')
axx.set_xlabel("Axial position in m")
axx.set_ylabel("Quality")
axx.set_title("Comparaison of different qualities for EPRI correlation")
axx.grid()
axx.legend(loc="best")
plt.show()
figx.savefig(rf'/home/p122173/Version5_BWR/Version5_ev3727/Donjon/data/results_marie/XEPRI.png', bbox_inches='tight')

# Affichage du tableau des valeurs de vide triées
EPSP12=[]
EPSP13=[]
EPSP14=[]
EPSP16=[]
extraire(EPSP12, 4, P12)
extraire(EPSP13, 4, P13)
extraire(EPSP14, 4, P14)
extraire(EPSP16, 4, P16)
Zaxis=np.linspace(0,2,len(EPSP12))
figeps, axeps = plt.subplots()
axeps.plot(Zaxis,EPSP12, label='P12')
axeps.plot(Zaxis,EPSP13, label='P13')
axeps.plot(Zaxis,EPSP14, label='P14')
axeps.plot(Zaxis,EPSP16, label='P16')
axeps.set_xlabel("Axial position in m")
axeps.set_ylabel("Void fraction")
axeps.set_title("Comparaison of different void fractions for EPRI correlation")
axeps.grid()
axeps.legend(loc="best")
plt.show()
figeps.savefig(rf'/home/p122173/Version5_BWR/Version5_ev3727/Donjon/data/results_marie/EPSPERI.png', bbox_inches='tight')
