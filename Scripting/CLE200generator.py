def generate_cle2000_script(num_volumes, total_height):
    # Calcul de la hauteur de chaque volume
    segment_height = total_height / num_volumes
    
    # Génération des positions des limites des volumes
    z_positions = []
    current_height = 0.0
    for i in range(num_volumes):
        current_height += segment_height
        z_positions.append(current_height)
    
    # Génération du script CLE-2000

    script = """***********************************************************
* Input file :  pincell_mphy_thm.x2m                      *
*                                                         *
* Test of a PWR single pincell                            *
*                                                         *
* The aim is to test the THM module of DONJON             *
***********************************************************
 
LINKED_LIST Geom Track Flux Cpo Matex Lib MicroF MacroF Macro2 Fmap
            System Power Burnup Thm DONNEES ;
MODULE      GEO: RESINI: TRIVAT: TRIVAA: GREP: FIND0: NCR: FLUD: FLPOW:
            MACINI: USPLIT: TINST: UTL: DELETE: ABORT: THM: END: ;

PROCEDURE   assertS ;
SEQ_ASCII   _ACompo ;
INTEGER     maxstep := 67 ;
REAL        Fuelpwr := 30.0796 ; ! W/g ou kW/kg ou MW/t
REAL        Mass := 6.464E-3 ; ! kg

*----
*  Calculation options
*----
STRING Dir := "EDI2B" ;
REAL   Tfuel := 900.0 ; (*fuel temperature (K)*)
REAL   Tcool := 550.0 ; (*coolant temperature (K)*)
REAL   Dfuel := 9.7413951 ; (*fuel density (g/cc)*)
REAL   dens_mod_0 := 0.65 ;
REAL   powi := 0.0591122 ;
ECHO "total reactor power=" powi "MW" ;
*----
*  Recover the Multicompo
*----
Cpo := _ACompo ;
UTL: Cpo :: DIR ;

INTEGER MaxR := 10000 ;
INTEGER Iter := 1 ;
REAL keff11 keff12 keff1n ;
REAL Y1 Y2 ROOT YNEW ;
LOGICAL CONV ;"""
    script += f"""
REAL maxh := {total_height} ;
"""
    for i, z in enumerate(z_positions, start=1):
        script += f"REAL z{i} := {z} ;\n"


    script += f"""
REAL Cote      := 1.404  ;
Geom := GEO: :: CAR3D 1 1 {num_volumes}
   X- REFL X+ REFL    Y- REFL Y+ REFL    Z- REFL Z+ REFL
   MESHX 0.0 <<Cote>>
   MESHY 0.0 <<Cote>>
   MESHZ 0.0 """
    
    for i in range(1, num_volumes + 1):
        if i%9 == 0:
            script += f"\n            "
        if i == num_volumes:
            script += f"<<z{i}>>"
        else:
            script += f"<<z{i}>> "
    
    script += """
   MIX
   PLANE 1
      1"""
    for i in range(2, num_volumes + 1):
        script += f"\n   PLANE {i}  SAME 1"

    script += f"""
;

Geom Matex := USPLIT: Geom :: NGRP 2 MAXR <<MaxR>>
               NFUEL 1  FMIX  1
;

Track := TRIVAT: Geom ::
   EDIT 1 MAXR <<MaxR>> MCFD 1 ;

*--
* Fuel map definition
*--
Fmap Matex := RESINI: Matex ::
      ::: GEO: CAR3D 1 1 {num_volumes}
                EDIT  0
                X- REFL X+ REFL    Y- REFL Y+ REFL    Z- REFL Z+ REFL
   MESHX 0.0 <<Cote>>
   MESHY 0.0 <<Cote>>
   MESHZ 0.0 """
    
    for i in range(1, num_volumes + 1):
        if i%9 == 0:
            script += f"\n            "
        if i == num_volumes:
            script += f"<<z{i}>>"
        else:
            script += f"<<z{i}>> "
    
    script += """
   MIX
   PLANE 1
      1"""
    for i in range(2, num_volumes + 1):
        script += f"\n   PLANE {i}  SAME 1"

    script += f"""
;
!
NXNAME '01' NYNAME  'A'
NCOMB 1
B-ZONE 1

ADD-PARAM PNAME 'T-FUEL' PARKEY 'TFA' GLOBAL
ADD-PARAM PNAME 'T-COOL' PARKEY 'TCA' GLOBAL
ADD-PARAM PNAME 'D-FUEL' PARKEY 'DFA' GLOBAL
ADD-PARAM PNAME 'D-COOL' PARKEY 'DCA' GLOBAL
BTYPE INST-BURN
INST-BVAL CHAN 0.0
"""

    chaine = ""
    for i in range(num_volumes):
        if i % 10 == 0:
            chaine += "\n"
        chaine += "1.0 "
    chaine = chaine.strip()  # Supprime l'espace final

    script += f"REACTOR-POW <<powi>> AXIAL-PFORM {chaine} \n"

    script += """SET-PARAM 'T-FUEL' <<Tfuel>>
SET-PARAM 'T-COOL' <<Tcool>>
SET-PARAM 'D-FUEL' <<Dfuel>>
SET-PARAM 'D-COOL' <<dens_mod_0>>
FUEL WEIGHT <<Mass>>
;

UTL: Fmap :: STEP UP PARAM STEP AT 4 DIR IMPR P-NAME * ;
*--
* THM single-stage calculation
*--
Thm Fmap := THM: Fmap ::
    EDIT 100
    FLUID H2O
    FPUISS 1.0
    INLET 10.8E6 (*Pa*) 550.0 (*K*)
    INLET-Q 2.1268E-5 (*m2*) 0.148880 (*inlet mass flow rate kg/s*)
    ASSMB 1 0
    RADIUS 5.6E-3 6.14E-3 6.52E-3 7.02E-3 (* m *)
    RODMESH 8 11
    HGAP 10000.0
    CONDC 0 10.0 KELVIN
    CONDF 0 5.0 KELVIN
    BOWR
    PDROP 1
;

*--
* Dump THM object
*--
UTL: Thm :: DIR DUMP ;

*--
* Cross-section database interpolation
*--
MicroF := NCR: Cpo Fmap ::
       EDIT 2
       MICRO LINEAR
       TABLE Cpo <<Dir>> 'burnup'
         MIX 1 INST-BURN
               SET LINEAR 'burnup' MAP
               SET CUBIC 'DCA' <<dens_mod_0>>
               SET CUBIC 'DCAH' <<dens_mod_0>>
               ADD 'DCA' <<dens_mod_0>> MAP
                        REF 'burnup' SAMEASREF
                        ENDREF
         ENDMIX
  ;
MacroF := MicroF :: STEP UP 'MACROLIB' ;
  
Macro2 Matex := MACINI: Matex MacroF :: FUEL ;
  
*--
* Steady-state diffusion calculation
*--
System := TRIVAA: Macro2 Track ;

Flux := FLUD: System Track :: EDIT 1 ADI 4 ACCE 5 3 ;
System MacroF Macro2 := DELETE: System MacroF Macro2 ;

GREP: Flux :: GETVAL 'K-EFFECTIVE' 1 >>keff11<< ;
ECHO "+++ Burnup= 0.0 Keff=" keff11 ;

*assertS Flux :: 'K-EFFECTIVE' 1 1.354165 ;

ECHO "test pincell_mphy_thm.x2m completed" ;
END: ;
"""
    return script

def access(filename):
    script=f"""#!/bin/sh
if [ $# = 0 ]
   then
   echo "usage: {filename}.access directory" 1>&2
   exit 1
fi
ln -s "$1"/data/pincell_mphy_thm_proc/_ACompo _ACompo
ls -l
echo "pincell_mphy_thm access script terminated"
"""
    # Écrire le script dans un fichier
    with open(f"{filename}.access", 'w') as file:
        file.write(script)

    return script
# Exemple d'utilisation
num_volumes = int(input("Combien de volume de contrôle il faut ?"))  # Nombre de volumes de contrôle
total_height = 155  # Hauteur totale en cm

# Générer le script CLE-2000
cle2000_script = generate_cle2000_script(num_volumes, total_height)
access_script = access(f"THM_{num_volumes}_{int(total_height)}")

# Nom du fichier de sortie
output_filename = f"THM_{num_volumes}_{int(total_height)}"


# Écrire le script dans un fichier
with open(f"{output_filename}.x2m", 'w') as file:
    file.write(f'{cle2000_script}')

print(f"Le script CLE-2000 a été écrit dans le fichier : {output_filename}")

