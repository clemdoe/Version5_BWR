import numpy as np 
import matplotlib.pyplot as plt 


# Données brutes sous forme de chaîne de caractères
Chexal = """PCOOL:   11054273.0       11030407.0       11005525.0       10980361.0       10954891.0       10929086.0       10902917.0       10876347.0       10849405.0       10800000.0    
 VCOOL:   9.23251629       9.23290825       9.36278725       9.50038719       9.64653969       9.80220890       9.96856403       10.1469727       10.3314667       13.1335011    
 DCOOL:   758.212280       758.179504       747.672241       736.856689       725.699585       714.169312       702.235718       689.874023       677.547363       532.877197    
 TCOOL:   552.857117       552.855042       558.460571       563.945984       569.305542       574.528809       579.600952       584.501770       589.040527       589.856201    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000      0.231642544    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       3.13469954E-02
 TSAT:   591.601135       591.438660       591.269043       591.097168       590.922852       590.745972       590.566284       590.383423       590.197693       589.856201    
 MUT:   9.52207629E-05   9.52139162E-05   9.28838126E-05   9.06177302E-05   8.84035617E-05   8.62309171E-05   8.40907451E-05   8.19758425E-05   7.99554109E-05   7.96065724E-05
"""

GEramp = """PCOOL:   11053564.0       11029698.0       11004816.0       10979652.0       10954182.0       10928377.0       10902208.0       10875638.0       10848685.0       10800000.0    
 VCOOL:   9.23252201       9.23292542       9.36279202       9.50038719       9.64653206       9.80221748       9.96857738       10.1469860       10.3314371       13.0752125    
 DCOOL:   758.211365       758.178223       747.671936       736.855835       725.698181       714.168518       702.237122       689.875061       677.548584       535.760193    
 TCOOL:   552.857056       552.855164       558.460144       563.945740       569.305542       574.528442       579.599670       584.500610       589.039368       589.856201    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000      0.226947382    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       3.13471928E-02
 TSAT:   591.596313       591.433838       591.264221       591.092285       590.918030       590.741089       590.561401       590.378540       590.192749       589.856201    
 MUT:   9.52205592E-05   9.52136324E-05   9.28837471E-05   9.06175774E-05   8.84032997E-05   8.62307861E-05   8.40909852E-05   8.19760171E-05   7.99555637E-05   7.96065724E-05
"""

HEM1 = """PCOOL:   11058459.0       11034593.0       11009710.0       10984546.0       10959076.0       10933272.0       10907102.0       10880533.0       10853591.0       10800000.0    
 VCOOL:   9.23244858       9.23284721       9.36275291       9.50035572       9.64651966       9.80220890       9.96860123       10.1470222       10.3316631       13.6328506    
 DCOOL:   758.218018       758.185242       747.678223       736.862000       725.703613       714.172607       702.236877       689.872620       677.531799       513.213257    
 TCOOL:   552.857422       552.855408       558.461060       563.947083       569.307556       574.531494       579.604675       584.506714       589.050659       589.856201    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000      0.263666213    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       3.13766561E-02
 TSAT:   591.629578       591.467224       591.297607       591.125793       590.951538       590.774719       590.595032       590.412292       590.226562       589.856201    
 MUT:   9.52219852E-05   9.52151167E-05   9.28850277E-05   9.06187634E-05   8.84043111E-05   8.62314773E-05   8.40909779E-05   8.19756533E-05   7.99529807E-05   7.96065724E-05
"""

EPRI = """PCOOL:   11051477.0       11027610.0       11002728.0       10977564.0       10952094.0       10926290.0       10900120.0       10873551.0       10846598.0       10800000.0    
 VCOOL:   9.23254490       9.23295784       9.36281776       9.50038719       9.64650440       9.80219746       9.96858120       10.1470337       10.3314219       12.8197107    
 DCOOL:   758.208679       758.175659       747.669556       736.853210       725.695740       714.167542       702.238831       689.879456       677.557007       546.339905    
 TCOOL:   552.856689       552.854797       558.459595       563.945190       569.304749       574.526917       579.596802       584.496765       589.034119       589.856201    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000      0.209717751    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       3.13557647E-02
 TSAT:   591.582092       591.419617       591.249939       591.078064       590.903748       590.726807       590.547058       590.364197       590.178345       589.856201    
 MUT:   9.52200353E-05   9.52131086E-05   9.28832669E-05   9.06170608E-05   8.84028341E-05   8.62305824E-05   8.40912980E-05   8.19767156E-05   7.99568661E-05   7.96065724E-05
"""

ModBestion = """PCOOL:   11054125.0       11030258.0       11005376.0       10980212.0       10954742.0       10928938.0       10902768.0       10876199.0       10849257.0       10800000.0    
 VCOOL:   9.23251247       9.23291874       9.36280537       9.50039101       9.64651680       9.80217934       9.96855450       10.1469955       10.3314962       13.1159163    
 DCOOL:   758.211487       758.178711       747.672729       736.856384       725.698730       714.168518       702.237244       689.877380       677.550171       533.556885    
 TCOOL:   552.857422       552.855347       558.460144       563.946045       569.305786       574.528992       579.600159       584.500366       589.039429       589.856201    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000      0.230535567    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       3.13575380E-02
 TSAT:   591.600098       591.437683       591.268005       591.096130       590.921875       590.744995       590.565247       590.382446       590.196655       589.856201    
 MUT:   9.52205883E-05   9.52137343E-05   9.28839363E-05   9.06176501E-05   8.84034089E-05   8.62307788E-05   8.40910216E-05   8.19763882E-05   7.99558329E-05   7.96065724E-05
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
figp.savefig(rf'/home/p122173/Version5_BWR/Version5_ev3727/Donjon/data/results_marie/PCorrel13.png', bbox_inches='tight')

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
figx.savefig(rf'/home/p122173/Version5_BWR/Version5_ev3727/Donjon/data/results_marie/XCorrel13.png', bbox_inches='tight')

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
figeps.savefig(rf'/home/p122173/Version5_BWR/Version5_ev3727/Donjon/data/results_marie/EPSCorrel13.png', bbox_inches='tight')


