import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd 

#Marie Bellier, 2025
#extract TPF data from csv files and plot 
#requires updating in file paths, and to copy/paste THM data from THM results
#physical test geometry is described in the internship report


def readTPF(file_path):
    # Charger les données
    df = pd.read_csv(file_path, delimiter=',', skipinitialspace=True)
    z = df['zPosition']
    # Convertir les colonnes de température en float
    TempVap = df['TempVap'].astype(float)
    TempLiq = df['TempLiq'].astype(float)
    TempFuelC = df['TempFuelC'].astype(float)
    TempStru = df['TempStru'].astype(float)
    Void = df['Void'].astype(float)
    Pressure = df['Pressure'].astype(float)
    VelzVap = df['VelzVap'].astype(float)
    VelzLiq = df['VelzLiq'].astype(float)

    for i in range(len(z)):
        TempVap[i] = TempVap[i] + 273.15
        TempLiq[i] = TempLiq[i] + 273.15
        TempFuelC[i] = TempFuelC[i] + 273.15
        TempStru[i] = TempStru[i] + 273.15
        
    return z, TempVap, TempLiq, TempFuelC, TempStru, Void, Pressure, VelzVap, VelzLiq


# Données brutes sous forme de chaîne de caractères
raw_data_0="""PCOOL:   7214342.50       7213741.00       7213140.00       7212539.50       7211939.00       7211339.00       7210739.50       7210140.00       7209541.00       7208942.00       7208343.50       7207745.50       7207148.00       7206550.50       7205953.00       7205356.50       7204760.00       7204164.00       7203568.00       7202972.50       7202377.50       7201782.50       7201188.00       7200594.00       7200000.00    
 VCOOL:   1.28856552       1.28962088       1.29067159       1.29173410       1.29279339       1.29386175       1.29492640       1.29600394       1.29707742       1.29816079       1.29924500       1.30033231       1.30142510       1.30251884       1.30361819       1.30472124       1.30583024       1.30693865       1.30805612       1.30917299       1.31029510       1.31142163       1.31255317       1.31368697       1.31482887    
 DCOOL:   769.617004       768.987183       768.361206       767.729187       767.100098       766.466675       765.836548       765.199768       764.566528       763.928406       763.290955       762.652588       762.012146       761.372192       760.729431       760.086853       759.442017       758.797424       758.148865       757.502136       756.853394       756.202759       755.551208       754.898987       754.245117    
 TCOOL:   543.507751       543.864380       544.217712       544.570618       544.926147       545.281250       545.632874       545.987244       546.338013       546.691589       547.041565       547.394287       547.743713       548.095398       548.444031       548.794800       549.142029       549.492432       549.839417       550.188416       550.534485       550.881409       551.227783       551.573792       551.919312    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000    
 TSAT:   561.028198       561.022522       561.016846       561.011169       561.005493       560.999817       560.994202       560.988525       560.982849       560.977173       560.971558       560.965881       560.960205       560.954590       560.948914       560.943298       560.937622       560.932007       560.926392       560.920715       560.915100       560.909485       560.903870       560.898254       560.892639    
 MUT:   9.79735851E-05   9.78181415E-05   9.76642987E-05   9.75107978E-05   9.73563074E-05   9.72021517E-05   9.70496549E-05   9.68961031E-05   9.67442538E-05   9.65913277E-05   9.64400824E-05   9.62877893E-05   9.61370315E-05   9.59854297E-05   9.58352539E-05   9.56842778E-05   9.55349242E-05   9.53843191E-05   9.52352857E-05   9.50854883E-05   9.49370369E-05   9.47883163E-05   9.46399159E-05   9.44917556E-05   9.43438936E-05
 HCOOL:   1186342.25       1188151.75       1189961.12       1191770.62       1193580.12       1195389.75       1197199.12       1199008.75       1200818.12       1202627.75       1204437.12       1206246.75       1208056.25       1209865.88       1211675.38       1213485.00       1215294.50       1217104.12       1218913.62       1220723.12       1222532.62       1224342.12       1226151.62       1227961.12       1229770.62    
"""

raw_data_1="""PCOOL:   7214656.50       7214054.00       7213452.50       7212851.50       7212251.50       7211652.50       7211054.00       7210456.50       7209859.50       7209263.00       7208668.00       7208073.50       7207479.50       7206886.50       7206294.50       7205703.00       7205112.50       7204522.50       7203905.00       7203282.50       7202655.50       7202019.00       7201378.00       7200697.00       7200000.00    
 VCOOL:   1.28931248       1.29186964       1.29444551       1.29704142       1.29966629       1.30230761       1.30497110       1.30765975       1.31036675       1.31310368       1.31585729       1.31864083       1.32144344       1.32427454       1.32712877       1.33000934       1.33291602       1.33584905       1.35374999       1.37184203       1.39012909       1.41511953       1.43932295       1.49148846       1.54262388    
 DCOOL:   769.171143       767.648621       766.120972       764.585449       763.043091       761.495911       759.941528       758.378967       756.812073       755.234985       753.653687       752.063599       750.468384       748.864014       747.253479       745.634888       744.008972       742.375549       732.559021       722.897949       713.388367       700.789246       689.005981       664.906982       642.868652    
 TCOOL:   544.011658       544.869934       545.726501       546.580383       547.431885       548.281006       549.124878       549.968384       550.809265       551.647583       552.483215       553.316223       554.146484       554.974060       555.798828       556.620789       557.439880       558.256104       559.069397       559.879761       560.687073       560.469421       560.905640       560.899231       560.892639    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       1.16290273E-02   2.30797548E-02   3.43572497E-02   5.06668091E-02   6.74902350E-02  0.102009259      0.133576751    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       7.42421311E-04   1.49498601E-03   2.25794106E-03   3.03153973E-03   3.69640626E-03   5.07469056E-03   8.03337246E-03   1.09923799E-02
 TSAT:   561.031189       561.025513       561.019836       561.014099       561.008484       561.002808       560.997131       560.991516       560.985840       560.980225       560.974609       560.968994       560.963379       560.957764       560.952148       560.946594       560.940979       560.935425       560.929565       560.923706       560.917725       560.911743       560.905640       560.899231       560.892639    
 MUT:   9.77543677E-05   9.73813658E-05   9.70099936E-05   9.66406151E-05   9.62730483E-05   9.59072349E-05   9.55443684E-05   9.51822803E-05   9.48218949E-05   9.44631320E-05   9.41059989E-05   9.37504083E-05   9.33963602E-05   9.30437818E-05   9.26926659E-05   9.23429616E-05   9.19946469E-05   9.16476565E-05   9.13019685E-05   9.09575319E-05   9.06143396E-05   9.07065769E-05   9.05634806E-05   9.05659545E-05   9.05685010E-05
 HCOOL:   1188913.75       1193294.75       1197675.88       1202057.00       1206438.00       1210818.88       1215199.88       1219581.00       1223962.00       1228343.00       1232724.00       1237105.00       1241486.00       1245867.12       1250248.12       1254629.25       1259010.12       1263391.12       1267772.00       1272153.00       1276534.00       1280915.00       1285296.00       1289677.00       1294058.00    
"""

raw_data_2 = """PCOOL:   7220422.50       7219818.50       7219216.00       7218614.50       7218015.50       7217417.50       7216821.00       7216172.50       7215515.50       7214850.00       7214167.00       7213454.50       7212676.00       7211864.00       7211015.00       7210126.50       7209196.50       7208222.00       7207201.50       7206132.00       7205013.00       7203841.50       7202617.00       7201336.50       7200000.00    
 VCOOL:   1.29071617       1.29615140       1.30168760       1.30731738       1.31305265       1.31889021       1.32483757       1.36099207       1.39794183       1.43572259       1.48594081       1.55033004       1.65418100       1.75417852       1.85081434       1.94433033       2.03507400       2.12332368       2.20931339       2.29324031       2.37530255       2.45566463       2.53448462       2.61190367       2.68805242    
 DCOOL:   768.336609       765.113037       761.859192       758.578125       755.265564       751.921509       748.546387       728.661255       709.401794       690.734192       667.392395       639.723572       599.558411       565.377625       535.854919       510.078033       487.330505       467.071655       448.888794       432.455719       417.512604       403.845337       391.282349       379.680908       368.921661    
 TCOOL:   544.963501       546.766846       548.555664       550.336914       552.102844       553.858887       555.602844       557.334412       559.053162       560.758789       560.353149       561.019836       561.012451       561.004822       560.996704       560.988403       560.979614       560.970398       560.960754       560.950623       560.940063       560.928955       560.917358       560.905273       560.892639    
 EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       2.32942086E-02   4.58799861E-02   6.77969009E-02   9.81650129E-02  0.137801126      0.195346206      0.244317517      0.286615461      0.323546678      0.356137782      0.385163277      0.411214530      0.434758842      0.456168443      0.475750118      0.493749440      0.510370970      0.525785804    
 XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       1.50241051E-03   3.04742251E-03   4.63719666E-03   6.27405848E-03   7.61104189E-03   1.14353159E-02   1.76554266E-02   2.38761026E-02   3.00976280E-02   3.63199823E-02   4.25431281E-02   4.87673581E-02   5.49924634E-02   6.12187982E-02   6.74456432E-02   7.36737847E-02   7.99027905E-02   8.61327201E-02   9.23637301E-02
 TSAT:   561.085632       561.079956       561.074219       561.068542       561.062927       561.057251       561.051636       561.045471       561.039307       561.033020       561.026550       561.019836       561.012451       561.004822       560.996765       560.988403       560.979614       560.970398       560.960754       560.950623       560.940063       560.928955       560.917358       560.905273       560.892639    
 MUT:   9.73428178E-05   9.65623549E-05   9.57916054E-05   9.50270696E-05   9.42715778E-05   9.35222997E-05   9.27797009E-05   9.20433813E-05   9.13130571E-05   9.05883426E-05   9.07604117E-05   9.05194029E-05   9.05222551E-05   9.05251945E-05   9.05283086E-05   9.05315319E-05   9.05349225E-05   9.05384804E-05   9.05422057E-05   9.05461202E-05   9.05501947E-05   9.05544803E-05   9.05589550E-05   9.05636261E-05   9.05685010E-05
 HCOOL:   1193770.88       1203009.00       1212247.38       1221485.38       1230723.50       1239961.50       1249199.62       1258438.00       1267676.12       1276914.12       1286152.25       1295390.50       1304628.50       1313866.75       1323104.75       1332343.00       1341581.00       1350819.25       1360057.38       1369295.50       1378533.75       1387771.88       1397009.88       1406247.88       1415485.88    
"""

raw_data_4="""PCOOL:   7238661.00       7238053.00       7237356.00       7236644.50       7235917.00       7235158.00       7234343.00       7233396.50       7232375.00       7231268.50       7230069.00       7228769.00       7227362.00       7225842.00       7224203.00       7222439.00       7220545.00       7218515.00       7216342.50       7214023.00       7211552.00       7208920.00       7206122.00       7203151.00       7200000.00    
 VCOOL:   1.29330111       1.30418754       1.36979759       1.43813944       1.50945008       1.60119259       1.72050786       1.90908229       2.08647251       2.25471878       2.41538167       2.56975627       2.71890211       2.86330247       3.00411487       3.14177394       3.27679372       3.40962815       3.54068565       3.67033195       3.79675221       3.92363334       4.04978704       4.17543030       4.30081129    
 DCOOL:   766.799194       760.398193       723.977051       689.571838       656.994751       619.349365       576.411316       519.473633       475.309418       439.841370       410.583862       385.919189       364.749451       346.354401       330.119354       315.654633       302.647888       290.856812       280.090637       270.196716       261.199097       252.752258       244.878265       237.508652       230.584106    
 TCOOL:   546.712402       550.228271       553.700745       557.123352       560.498779       559.929749       561.216980       561.208008       561.198425       561.187988       561.176636       561.164368       561.151123       561.136780       561.121338       561.104675       561.086792       561.067627       561.047119       561.025208       561.001831       560.976990       560.950500       560.922424       560.892639    
 EPS:   0.00000000       0.00000000       4.17665653E-02   8.13012794E-02  0.118821479      0.166577101      0.228117988      0.309744358      0.373058379      0.423905671      0.465849370      0.501208663      0.531557143      0.557927489      0.581200838      0.601935685      0.620579839      0.637480497      0.652911186      0.667090595      0.679984450      0.692087650      0.703368664      0.713925540      0.723843277    
 XFL:   0.00000000       2.73411046E-03   5.61538246E-03   8.65802635E-03   1.18786581E-02   1.41239613E-02   2.17932370E-02   3.40319388E-02   4.62716930E-02   5.85130788E-02   7.07565323E-02   8.30014795E-02   9.52483267E-02  0.107497014      0.119747445      0.131999865      0.144253984      0.156509817      0.168767318      0.181026444      0.193287015      0.205548912      0.217812330      0.230076954      0.242342427    
 TSAT:   561.257629       561.251892       561.245361       561.238647       561.231812       561.224609       561.216980       561.208008       561.198425       561.187988       561.176636       561.164368       561.151123       561.136780       561.121338       561.104675       561.086792       561.067627       561.047119       561.025208       561.001831       560.976990       560.950500       560.922424       560.892639    
 MUT:   9.65919608E-05   9.50800750E-05   9.35965218E-05   9.21403043E-05   9.07065114E-05   9.09478986E-05   9.04433400E-05   9.04467961E-05   9.04504923E-05   9.04545232E-05   9.04589033E-05   9.04636327E-05   9.04687404E-05   9.04742774E-05   9.04802364E-05   9.04866611E-05   9.04935659E-05   9.05009583E-05   9.05088746E-05   9.05173292E-05   9.05263514E-05   9.05359411E-05   9.05461638E-05   9.05570050E-05   9.05685010E-05
 HCOOL:   1202723.25       1220913.88       1239104.38       1257294.88       1275485.25       1293675.88       1311866.38       1330056.88       1348247.38       1366437.62       1384628.25       1402818.62       1421009.12       1439199.62       1457390.12       1475580.62       1493771.12       1511961.62       1530152.12       1548342.62       1566533.12       1584723.62       1602914.12       1621104.75       1639295.00    
"""

z, TempVap, TempLiq, TempFuelC, TempStru, Void, Pressure, VelzVap, VelzLiq, VelzMoy = [], [], [], [], [], [], [], [], [], []

# Traitement de chaque ligne pour extraire les valeurs numériques
def extraire(var_name, line_nb, data_name):
    data_lines = data_name.strip().split("\n")
    values = data_lines[line_nb].split(":")[1]
    # Séparation des valeurs en une liste de flottants
    for value in values.split():
        var_name.append(float(value))
    return

# Valeurs THM sans Chruchill
PCOOL0=[]
PCOOL1=[]
PCOOL2=[]
PCOOL4=[]
extraire(PCOOL0, 0, raw_data_0)
extraire(PCOOL1, 0, raw_data_1)
extraire(PCOOL2, 0, raw_data_2)
extraire(PCOOL4, 0, raw_data_4)

VCOOL0=[]
VCOOL1=[]
VCOOL2=[]
VCOOL4=[]
extraire(VCOOL0, 1, raw_data_0)
extraire(VCOOL1, 1, raw_data_1)
extraire(VCOOL2, 1, raw_data_2)
extraire(VCOOL4, 1, raw_data_4)

EPS0=[]
EPS1=[]
EPS2=[]
EPS4=[]
extraire(EPS0, 4, raw_data_0)
extraire(EPS1, 4, raw_data_1)
extraire(EPS2, 4, raw_data_2)
extraire(EPS4, 4, raw_data_4)

TCOOL0=[]
TCOOL1=[]
TCOOL2=[]
TCOOL4=[]
extraire(TCOOL0, 3, raw_data_0)
extraire(TCOOL1, 3, raw_data_1)
extraire(TCOOL2, 3, raw_data_2)
extraire(TCOOL4, 3, raw_data_4)

VCOOL=[VCOOL4, VCOOL2, VCOOL1, VCOOL0]
PCOOL=[PCOOL4, PCOOL2, PCOOL1, PCOOL0]
EPS=[EPS4, EPS2, EPS1, EPS0]
TCOOL=[TCOOL4, TCOOL2, TCOOL1, TCOOL0]


powers = [4,2,1,0]
powerskW = [38.4, 19.2, 9.6, 3.8]

for i in powers:
    X = readTPF(rf"THM_prototype/TPF/graph_axial_{i}e9W.txt")
    z.append(X[0])
    TempVap.append(X[1])
    TempLiq.append(X[2])
    TempFuelC.append(X[3])
    TempStru.append(X[4])
    Void.append(X[5])
    Pressure.append(X[6])
    VelzVap.append(X[7])
    VelzLiq.append(X[8])
    VelzMoy.append(0)

colors = ['blue', 'red', 'green', 'orange']

fig, ax = plt.subplots()
for i in range(len(z)):
    ax.step(z[i], TempLiq[i], label=f"TPF {powerskW[i]} kW", linestyle='-', color = colors[i])
    ax.step(z[i], TCOOL[i], label=f"THM {powerskW[i]} kW", linestyle=':', color = colors[i])
ax.set_xlabel("Axial position in m")
ax.set_ylabel("Temperature in K")
ax.grid()
plt.savefig(rf'/home/p122173/Version5_BWR/Version5_ev3727/Donjon/data/results_marie/TPF_THM_temp.png', bbox_inches='tight') #modify all access paths as needed !

fig1, ax1 = plt.subplots()
for i in range(len(z)):
    ax1.step(z[i], Void[i], label=f"TPF {powerskW[i]} kW", linestyle='-', color = colors[i])
    ax1.step(z[i], EPS[i], label=f"THM {powerskW[i]} kW", linestyle=':', color = colors[i])
ax1.set_xlabel("Axial position in m")
ax1.set_ylabel("Void fraction")
ax1.grid()
fig1.savefig(rf'/home/p122173/Version5_BWR/Version5_ev3727/Donjon/data/results_marie/TPF_THM_void.png', bbox_inches='tight')

fig2, ax2 = plt.subplots()
for i in range(len(z)):
    ax2.step(z[i], Pressure[i], label=f"TPF {powerskW[i]} kW", linestyle='-', color = colors[i])
    ax2.step(z[i], PCOOL[i], label=f"THM {powerskW[i]} kW", linestyle=':', color = colors[i])
ax2.set_xlabel("Axial position in m")
ax2.set_ylabel("Pressure in Pa")
ax2.grid()
fig2.savefig(rf'/home/p122173/Version5_BWR/Version5_ev3727/Donjon/data/results_marie/TPF_THM_pressure.png', bbox_inches='tight')

fig21, ax21 = plt.subplots()
for i in (0,1):
    ax21.step(z[i], Pressure[i], label=f"TPF {powerskW[i]} kW", linestyle='-', color = colors[i])
    ax21.step(z[i], PCOOLc[i], label=f"THM {powerskW[i]} kW", linestyle=':', color = colors[i])
ax21.set_xlabel("Axial position in m")
ax21.set_ylabel("Pressure in Pa")   
ax21.grid()
fig21.savefig(rf'/home/p122173/Version5_BWR/Version5_ev3727/Donjon/data/results_marie/TPF_THM_pressure_1.png', bbox_inches='tight')

fig22, ax22 = plt.subplots()
for i in (2,3):
    ax22.step(z[i], Pressure[i], label=f"TPF {powerskW[i]} kW", linestyle='-', color = colors[i])
    ax22.step(z[i], PCOOLc[i], label=f"THM {powerskW[i]} kW", linestyle=':', color = colors[i])
ax22.set_xlabel("Axial position in m")
ax22.set_ylabel("Pressure in Pa")
ax22.grid()
fig22.savefig(rf'/home/p122173/Version5_BWR/Version5_ev3727/Donjon/data/results_marie/TPF_THM_pressure_2.png', bbox_inches='tight')


fig3, ax3 = plt.subplots()
for i in range(len(z)):
    VelzMoy[i] = VelzVap[i]*(1-Void[i]) + VelzLiq[i]*Void[i]
    ax3.step(z[i], VelzMoy[i], label=f"TPF {powerskW[i]} kW", linestyle='-', color = colors[i])
    ax3.step(z[i], VCOOLc[i], label=f"THM {powerskW[i]} kW", linestyle=':', color = colors[i])
ax3.set_xlabel("Axial position in m")
ax3.set_ylabel("Coolant velocity in m/s")
ax3.grid()
fig3.savefig(rf'/home/p122173/Version5_BWR/Version5_ev3727/Donjon/data/results_marie/TPF_THM_velocity.png', bbox_inches='tight')



# Calcul des différences point par point
for i in range(len(PCOOL)):
    diffsP = [a - b for a, b in zip(PCOOL[i], Pressure[i])]
    diffsPc =[a - b for a, b in zip(PCOOLc[i], Pressure[i])]
    diffsT = [a - b for a, b in zip(TCOOL[i], TempLiq[i])]
    diffsTc =[a - b for a, b in zip(TCOOLc[i], TempLiq[i])]
    diffsV = [a - b for a, b in zip(VCOOL[i], VelzMoy[i])]
    diffsVc =[a - b for a, b in zip(VCOOLc[i], VelzMoy[i])]
    diffsE = [a - b for a, b in zip(EPS[i], Void[i])]  
    diffsEc =[a - b for a, b in zip(EPSc[i], Void[i])]

    # RMS deviation in %
    rms_deviation_P = np.sqrt(sum(d**2 for d in diffsP) / len(diffsP))/np.mean(Pressure[i])*100
    rms_deviation_Pc = np.sqrt(sum(d**2 for d in diffsPc) / len(diffsPc))/np.mean(Pressure[i])*100
    rms_deviation_T = np.sqrt(sum(d**2 for d in diffsT) / len(diffsT))/np.mean(TempLiq[i])*100
    rms_deviation_Tc = np.sqrt(sum(d**2 for d in diffsTc) / len(diffsTc))/np.mean(TempLiq[i])*100
    rms_deviation_V = np.sqrt(sum(d**2 for d in diffsV) / len(diffsV))/np.mean(VelzMoy[i])*100
    rms_deviation_Vc = np.sqrt(sum(d**2 for d in diffsVc) / len(diffsVc))/np.mean(VelzMoy[i])*100
    rms_deviation_E = np.sqrt(sum(d**2 for d in diffsE) / len(diffsE))*100
    rms_deviation_Ec = np.sqrt(sum(d**2 for d in diffsEc) / len(diffsEc))*100
    # Max deviation in %
    max_deviation_P = max(abs(d) for d in diffsP)/np.mean(Pressure[i])*100
    max_deviation_Pc = max(abs(d) for d in diffsPc)/np.mean(Pressure[i])*100
    max_deviation_T = max(abs(d) for d in diffsT)/np.mean(TempLiq[i])*100
    max_deviation_Tc = max(abs(d) for d in diffsTc)/np.mean(TempLiq[i])*100
    max_deviation_V = max(abs(d) for d in diffsV)/np.mean(VelzMoy[i])*100
    max_deviation_Vc = max(abs(d) for d in diffsVc)/np.mean(VelzMoy[i])*100
    max_deviation_E = max(abs(d) for d in diffsE)*100
    max_deviation_Ec = max(abs(d) for d in diffsEc)*100

    # Affichage des résultats
    print(powerskW[i], "kW") 
    print(f"RMS Deviation Pressure wo C : {rms_deviation_P:.4f}")
    print(f"Max Deviation Pressure wo C : {max_deviation_P:.4f}")
    print(f"RMS Deviation Pressure w C : {rms_deviation_Pc:.4f}")
    print(f"Max Deviation Pressure w C : {max_deviation_Pc:.4f}")
    print(f"RMS Deviation Temperature wo C : {rms_deviation_T:.4f}")
    print(f"Max Deviation Temperature wo C : {max_deviation_T:.4f}")
    print(f"RMS Deviation Temperature w C : {rms_deviation_Tc:.4f}")
    print(f"Max Deviation Temperature w C : {max_deviation_Tc:.4f}")
    print(f"RMS Deviation Velocity wo C : {rms_deviation_V:.4f}")
    print(f"Max Deviation Velocity wo C : {max_deviation_V:.4f}")
    print(f"RMS Deviation Velocity w C : {rms_deviation_Vc:.4f}")
    print(f"Max Deviation Velocity w C : {max_deviation_Vc:.4f}")
    print(f"RMS Deviation Void wo C : {rms_deviation_E:.4f}")
    print(f"Max Deviation Void wo C : {max_deviation_E:.4f}")
    print(f"RMS Deviation Void w C : {rms_deviation_Ec:.4f}")
    print(f"Max Deviation Void w C : {max_deviation_Ec:.4f}")


plt.show()


