import numpy as np 
import matplotlib.pyplot as plt 


# Données brutes sous forme de chaîne de caractères
Chexal = """PCOOL:   11117765.0       11093893.0       11068924.0       11043646.0       11018025.0       10992029.0       10965620.0       10938756.0       10897429.0       10800000.0    
 VCOOL:   9.23620605       9.23661423       9.37680340       9.52600956       9.68525982       9.85587692       10.0395527       10.2383661       12.1082964       14.7495108    
 DCOOL:   757.901184       757.867981       746.541809       734.854309       722.772278       710.256042       697.266785       683.766052       578.694397       473.150879    
 TCOOL:   553.079285       553.077393       559.103027       564.989563       570.725464       576.297363       581.686829       586.870972       590.528503       589.856201    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000      0.154979259      0.328909844    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       1.92806460E-02   5.11834733E-02
 TSAT:   592.031982       591.870178       591.700745       591.528809       591.354309       591.176880       590.996338       590.812317       590.528503       589.856201    
 MUT:   9.51488328E-05   9.51418988E-05   9.26404900E-05   9.02115644E-05   8.78418796E-05   8.55185863E-05   8.32303776E-05   8.09672056E-05   7.93517174E-05   7.96065724E-05
"""

GEramp = """PCOOL:   11117832.0       11093960.0       11068991.0       11043712.0       11018092.0       10992095.0       10965685.0       10938820.0       10897975.0       10800000.0    
 VCOOL:   9.23620892       9.23660946       9.37678814       9.52597713       9.68524456       9.85589790       10.0396318       10.2384615       12.0521021       14.7431879    
 DCOOL:   757.901855       757.868652       746.541626       734.852966       722.772583       710.260254       697.271606       683.765442       581.361023       473.274689    
 TCOOL:   553.078979       553.077087       559.103210       564.990234       570.725342       576.295593       581.684937       586.871216       590.532288       589.856201    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000      0.150608450      0.328708231    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       1.92892831E-02   5.18886782E-02
 TSAT:   592.032410       591.870667       591.701172       591.529297       591.354736       591.177368       590.996765       590.812744       590.532288       589.856201    
 MUT:   9.51489783E-05   9.51420516E-05   9.26404318E-05   9.02113097E-05   8.78419523E-05   8.55193503E-05   8.32312144E-05   8.09671328E-05   7.93502841E-05   7.96065724E-05
"""

HEM1 = """PCOOL:   11123414.0       11099542.0       11074574.0       11049295.0       11023676.0       10997685.0       10971281.0       10944419.0       10900905.0       10800000.0    
 VCOOL:   9.23615170       9.23655605       9.37665844       9.52576447       9.68498135       9.85566807       10.0393181       10.2375498       12.3622313       15.4077978    
 DCOOL:   757.908691       757.875854       746.548706       734.859131       722.780640       710.275574       697.297058       683.794189       566.285034       457.291016    
 TCOOL:   553.079956       553.077881       559.104431       564.992310       570.726868       576.294495       581.680481       586.866638       590.552429       589.856201    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000      0.175194591      0.354738474    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       1.91202182E-02   4.73767892E-02
 TSAT:   592.070190       591.908508       591.739075       591.567261       591.392822       591.215515       591.035034       590.851135       590.552429       589.856201    
 MUT:   9.51503825E-05   9.51435213E-05   9.26418288E-05   9.02124666E-05   8.78434585E-05   8.55221078E-05   8.32356091E-05   8.09718840E-05   7.93426443E-05   7.96065724E-05
"""

EPRI = """PCOOL:   11112291.0       11088418.0       11063449.0       11038170.0       11012550.0       10986554.0       10960144.0       10933279.0       10893418.0       10800000.0    
 VCOOL:   9.23631096       9.23671818       9.37691021       9.52609730       9.68535233       9.85599804       10.0396614       10.2384386       11.9344473       14.2900839    
 DCOOL:   757.894165       757.861023       746.532410       734.843811       722.761719       710.247681       697.260376       683.758667       587.212036       487.920380    
 TCOOL:   553.078491       553.076660       559.103149       564.989746       570.725159       576.295593       581.683777       586.867798       590.500916       589.856201    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000      0.141141951      0.304856956    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       1.93976518E-02   5.30423708E-02
 TSAT:   591.994873       591.833069       591.663513       591.491577       591.316956       591.139465       590.958862       590.774719       590.500916       589.856201    
 MUT:   9.51473921E-05   9.51404363E-05   9.26385765E-05   9.02095198E-05   8.78399078E-05   8.55170729E-05   8.32292426E-05   8.09659396E-05   7.93621875E-05   7.96065724E-05
"""

ModBestion = """PCOOL:   11117343.0       11093471.0       11068502.0       11043223.0       11017603.0       10991607.0       10965197.0       10938333.0       10897066.0       10800000.0    
 VCOOL:   9.23622227       9.23662186       9.37681389       9.52601814       9.68525791       9.85585880       10.0395412       10.2384081       12.1015825       14.7144508    
 DCOOL:   757.900879       757.868286       746.540833       734.852112       722.770508       710.255859       697.268311       683.768311       579.005188       474.233826    
 TCOOL:   553.079102       553.076904       559.103210       564.990234       570.725830       576.296997       581.685730       586.869629       590.526001       589.856201    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000      0.154478893      0.327146232    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       1.92958489E-02   5.12216128E-02
 TSAT:   592.029114       591.867371       591.697876       591.525940       591.351440       591.174011       590.993408       590.809387       590.526001       589.856201    
 MUT:   9.51487746E-05   9.51419643E-05   9.26402645E-05   9.02111351E-05   8.78415667E-05   8.55185644E-05   8.32306541E-05   8.09675985E-05   7.93526706E-05   7.96065724E-05
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
figp.savefig(rf'/home/p122173/Version5_BWR/Version5_ev3727/Donjon/data/results_marie/PCorrel14.png', bbox_inches='tight')

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
figx.savefig(rf'/home/p122173/Version5_BWR/Version5_ev3727/Donjon/data/results_marie/XCorrel14.png', bbox_inches='tight')

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
figeps.savefig(rf'/home/p122173/Version5_BWR/Version5_ev3727/Donjon/data/results_marie/EPSCorrel14.png', bbox_inches='tight')


