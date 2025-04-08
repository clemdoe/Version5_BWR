# Données brutes sous forme de chaîne de caractères
raw_data = """
PCOOL:   11137852.0       11112698.0       11087173.0       11061236.0       11034838.0       11007925.0       10980831.0       10945200.0       10881401.0       10800000.0
VCOOL:   9.24568653       9.40701008       9.58030033       9.76725292       9.96992397       10.1908073       10.3912401       11.5481129       12.6550112       13.7449045
DCOOL:   757.130737       744.146790       730.685547       716.700134       702.128662       686.912476       673.661499       606.174988       553.154846       509.292847
TCOOL:   553.513855       560.375366       567.052185       573.524292       579.770935       585.761963       590.559814       590.856445       590.418213       589.856201
EPS:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000      0.108962595      0.197041392      0.270050794
XFL:   0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       0.00000000       1.29405092E-02   2.58223973E-02   3.86294909E-02
TSAT:   592.167847       591.997681       591.824646       591.648499       591.468872       591.285400       591.100342       590.856506       590.418213       589.856201
MUT:   9.49746754E-05   9.21317187E-05   8.93800243E-05   8.67009949E-05   8.40756766E-05   8.14880623E-05   7.93356085E-05   7.92272695E-05   7.93935469E-05   7.96065724E-05
"""

# Séparation des données en lignes et traitement de chaque ligne
data_lines = raw_data.strip().split("\n")


# Traitement de chaque ligne pour extraire les valeurs numériques
def extraire(var_name, line_nb):
    values = data_lines[line_nb].split(":")[1]
    # Séparation des valeurs en une liste de flottants
    for value in values.split():
        var_name.append(float(value))
    return


# Affichage du tableau des valeurs triées
PCOOL=[]
extraire(PCOOL, 0)
extraire(VCOOL, 1)
extraire(DCOOL, 2)
extraire(TCOOL, 3)
extraire(EPS, 4)
extraire(XFL, 5)
extraire(TSAT, 6)
extraire(MUT, 7)




