*DECK THMDFM
      SUBROUTINE THMDFM(ITIME,I,J,K,K0,PINLET,MFLOW,HMAVG,ENT,HD,IFLUID,
     > IHCONV,KHCONV,ISUBM,RADCL,ZF,PHI,XFL,EPS,SLIP,ACOOL,PCH,DZ,TCALO,
     > RHO,RHOLAV,TSCLAD,KWA)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Drift-flux Model for the computation of thermohydraulics parameters in two-phase flow
*
*Copyright:
* Copyright (C) 2025 Ecole Polytechnique de Montreal.
*
*Author(s): 
* M.Bellier
*
*Parameters: input
* ITIME   type of calculation  (0=steady-state; 1=transient).
* I       position of channel alon X-axis
* J       position of channel alon Y-axis
* K       position of channel alon Z-axis
* K0      onser of nuclear boiling point
* PINLET  pressure in Pascal
* MFLOW   massic coolant flow rate in Kg/m^2/s
* HMAVG   averaged enthalpy
* ENT     four values of enthalpy in J/Kg to be used in Gaussian
*         integration
* HD      hydraulic diameter in m
* IFLUID  type of fluid (0=H2O; 1=D2O).
* IHCONV  flag indicating HCONV chosen (0=default/1=user-provided).
* KHCONV  fixed user-provided HCONV value in W/m^2/K.
* ISUBM   subcooling model (0: one-phase; 1: Jens-Lottes model;
*         2: Saha-Zuber model).
* RADCL   outer clad radius in m
* ZF      parameters used to compute heat flux on clad surface in
*         transient cases.
* PHI     heat flow exchanged between clad and fluid in W/m^2.
*         Given in steady-state cases.
* XFL     input coolant flow quality
* EPS     input coolant void fraction
* SLIP    input slip ratio of vapor phase speed to liquid phase speed.
* ACOOL   coolant cross section area in m^2.
* PCH     heating perimeter in m.
* DZ      axial mesh width in m.
*
*Parameters: output
* PHI     heat flow exchanged between clad and fluid in W/m^2.
*         Computed in transient cases.
* XFL     output coolant flow quality
* EPS     output coolant void fraction
* SLIP    output slip ratio of vapor phase speed to liquid phase speed.
* TCALO   coolant temperature in K
* RHO     coolant density in Kg/m^3
* RHOLAV  liquid density in kg/m^3
* TSCLAD  clad temperature in K
* KWA     flow regime (=0: single-phase; =1: subcooled; =2: nucleate
*         boiling; =3 superheated steam)
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER I,J,K,K0,IFLUID,IHCONV,ISUBM,KWA
      REAL PINLET,MFLOW,HMAVG,ENT(4),HD,KHCONV,RADCL,ZF(2),PHI,TCALO,
     > RHO,RHOLAV,TSCLAD,XFL,EPS,SLIP,ACOOL,PCH,DZ
*----
*  LOCAL VARIABLES
*----
      REAL W(4),HL(4),JL,JG
      CHARACTER HSMG*131
      LOGICAL LFIRST
      REAL EPSold
*----
*  SAVE VARIABLES
*----
      SAVE DHSUB,DSAT,W
      DATA W /0.347855,0.652145,0.652145,0.347855/

*réécrit à partir d'ici
*----
* INITIALIZE VARIABLES
*----
      VGJ = 0
      C0 = 1
      VGJprime = 0

*----
*  MAIN LOOP
*----
     I=0
     ERREPS=1

    10 CONTINUE

*----
*  SAVE THE OLD EPSILON VALUE
*----
     EPSold=EPS 
     I = I+1

*----
* TEST ON ERR EPS
*----
    IF (I .GT. 1000) GOTO 20
    IF (ERREPS < 1E-3) GOTO 20


*----
*  COMPUTE PHASES VELOCITIES AND REYNOLDS
*----
    VLIQ = VCOOL - (1/(1- EPS) - RHOLAV/RHO) *VGJprime
    VVAP = VCOOL + RHOLAV/RHO * VGJprime

    REY = RHO * ABS(VCOOL) * HD / (ZMUL*ZMUG/ (ZMUL*(1-EPS) + ZMUG*EPS))
  
*----
*  COMPUTE FLOW QUALITY
*----
#Get the enthalpy values for the liquid and vapor phases at a given pressure
    def getPhasesEnthalpy(self, i):
        P = self.P[i]
        vapor = IAPWS97(P = P*(10**(-6)), x = 1)
        liquid = IAPWS97(P = P*(10**(-6)), x = 0)
        return liquid.h, vapor.h


    IF CORREL.EQ.'EPRI' THEN
      Xs = 0.05
      Xh = 0.025
      Xeq = (HLAVG - HMAVG*0.001) / (HLAVG - hg) !hg à récupérer 
            
      IF (Xeq.GE.0.05) THEN
        XFL = Xeq
      ELSE 
      
        else:
          rhol = self.rholTEMP[i]
          rhog = self.rhogTEMP[i]
          u = self.U[i]
          muf = IAPWS97(P = p*(10**(-6)), x = 1).mu
          Re = rhol * abs(u) * self.D_h[i] / muf

          Cpf = IAPWS97(P = p*(10**(-6)), x = 1).cp
          k_f = IAPWS97(P = p*(10**(-6)), x = 1).k
          Pr = Cpf * muf / k_f

          qdp = self.q__[i] * self.DV / (2 * np.pi * self.rw * self.height)

          ! Calculate heat transfer coefficients
          hb = np.exp(p / 4.35e6) / (22.7)**2 * 1000.0  # W/(m^2·K)
          Chn = 0.2 / 4.0 * self.D_h[i] / self.rf
          hhn = Chn * Re**0.662 * Pr * k_f / (self.D_h[i])
          Cdb = (0.033 *self.areaMatrix[i] / (self.areaMatrix[i] +  np.pi * self.rf**2 + np.pi *self.rw**2) + 0.013)
          hdb = Cdb * Re**0.8 * Pr**0.4 * k_f / (self.D_h[i])

          !Intermediate calculations
          tmp1 = 4.0 * hb * (hdb + hhn)**2
          tmp2 = 2.0 * hdb**2 * (hhn + hdb / 2.0) + 8.0 * qdp * hb * (hdb + hhn)
          tmp3 = qdp * (4.0 * hb * qdp + hdb**2)

          !Calculate characteristic quality xd  ///////////////PROBLEM HERE 
          !Toujours problème ? 
          delta_h = hg - hl
          xd = -Cpf / (delta_h) * ((-tmp2 + np.sqrt(tmp2**2 - 4.0 * tmp1 * tmp3)) / (2.0 * tmp1))

          !Determine quality based on xeq and xd
                    if xeq <= 0.0:
                        if xeq <= -xd:
                            QUALITY = 0.0       
                            print(f'IF 2')
                        else:
                            tmp1 = 1 + xeq / xd
                            tmp2 = tmp1**2
                            QUALITY = xd * tmp2 * (0.1 + 0.087 * tmp1 + 0.05 * tmp2)
                            print(f'IF 3')
                    elif Xh > xd:
                        if xeq >= 2 * xd:
                            QUALITY = xeq
                            print(f'IF 4')
                        else:
                            tmp1 = xeq / xd
                            QUALITY = xd * (0.237 + tmp1 * (0.661 + tmp1 * (0.153 + tmp1 * (-0.01725 - tmp1 * 0.0020625))))
                            print(f'IF 5')
                    else:
                        tmp1 = xeq / Xh
                        tmp2 = xd / Xh
                        tmp3 = 0.237 * tmp2
                        QUALITY = Xh * (tmp3 + tmp1 * (0.661 + tmp1 * (0.5085 - 0.3555 * tmp2 + tmp1 * (tmp3 - 0.25425 + tmp1 * (0.042375 - 0.0444375 * tmp2)))))
                        print(f'IF 6')

                print(f'Quality: {QUALITY}')
                if QUALITY >= 1.0:
                    return 0.99
                else:
                    return QUALITY

    ELSE 
    !hl, hg = self.getPhasesEnthalpy(i) trouver un moyen de les récupérer, solution provisoire = utiliser HGSAT et HLSAT
      H = self.H[i]
      hfg = IAPWS97(P = self.P[i]*(10**(-6)), x = 0).h
      if H*0.001 < hl: !demander à clément: HMAVG plutôt que H, qu'est ce que HFG ? 
        x = (H*0.001 - hl)/(hg - hl)
        return 0
      elif H*0.001 > hg:
        return 1
            
      elif H*0.001 <= hg and H*0.001 >= hl:
        return (H*0.001 - hl)/(hg - hl)
*----
*  COMPUTE DENSITIES
*----
def getDensity(self, i):
        vapor = IAPWS97(P = self.P[i]*(10**(-6)), x = 1)
        liquid = IAPWS97(P = self.P[i]*(10**(-6)), x = 0)
        rho_g = vapor.rho
        rho_l = liquid.rho
        rho = rho_l * (1 - self.voidFractionTEMP[i]) + rho_g * self.voidFractionTEMP[i]
        return rho_l, rho_g, rho

*----
*  COMPUTE VGJ, VGJprime and C0 AFTER CHOSEN CORRELATION
*----
CALL MARIE(VCOOL, DCOOL, PCOOL, MUT, XFL, HD, RHOG, RHOL, EPS, CORREL, VGJ, C0)
#Get the drift velocity for a given cell based on the selected void fraction correlation
    def getVgj_prime(self, i):
        U = self.U[i]
        C0 = self.C0TEMP[i]
        Vgj = self.VgjTEMP[i]
        Vgj_prime = Vgj + (C0-1) * U
        return Vgj_prime

*----
*  COMPUTE HD (hydraulic diameter)
*----
!pourquoi cecalculer le Dh ? Qui est fixe ? Est-ce que c'était un autre paramètre à recalculer ? En fait DHGL ? 
*----
*  COMPUTE NEW EPS VALUE
*----
#Get the void fraction for a given cell using different correlations for drift flux model
    def getVoidFraction(self, i):
        correl = 'paths'
        if correl == 'simple':
            x_th = self.xThTEMP[i]
            rho_l = self.rholTEMP[i]
            rho_g = self.rhogTEMP[i]
            if x_th == 0:
                return 0
            elif x_th == 1:
                return 0.99
            else:
                return (x_th * rho_l)/(x_th * rho_l + (1 - x_th) * rho_g)
        elif correl == 'paths':
            x_th = self.xThTEMP[i]
            rho_l = self.rholTEMP[i]
            rho_g = self.rhogTEMP[i]
            u = self.U[i]
            V_gj = self.VgjTEMP[i]
            C0 = self.C0TEMP[i]
            if x_th == 0:
                return 0
            elif x_th == 1:
                return 0.99
            else:
                return x_th / (C0 * (x_th + (rho_g / rho_l) * (1 - x_th)) + (rho_g * V_gj) / (rho_l * u))
    
*----
*  COMPUTE DELTA BETWEEN EPSold AND EPS
*----
    ERREPS = ABS(ERRold - EPS)
    GOTO 10


*----
* EXIT LOOP
*----
    20 CONTINUE

      IF (I == 1000) THEN
        PRINT *, 'Nombre maximum d''itérations max atteint'
      ELSE
        PRINT *, 'Convergence atteinte à I = ', I
      ENDIF



*----
*  COMPUTE THE STEAM FLOW QUALITY AND LIQUID ENTHALPY
*  Reference: R. T. Lahey Jr. and F. J. Moody, "The thermal hydraulics
*  of a Boiling water nuclear reactor," American Nuclear Society, 1977.
*  Equation (5.177), page 224
*  F2: Thermodynamic quality
*----
      TSCLAD=600.0
      IF(K0.GT.0) TSCLAD=TSAT+DSAT
      XFL0=XFL
      EPS0=EPS
      SLIP0=SLIP
      LFIRST=.TRUE.
      HLAVG=HMAVG
      F2=0.0
      F3=0.0
      IF(K0.GT.0) THEN
        HLV=HGSAT-HLSAT
        IF((HLV.GT.0.0).AND.(DHSUB.GT.0.0)) THEN
          F2=(HMAVG-HLSAT)/HLV
          F3=(DHSUB/HLV)*EXP(-(HMAVG-HLSAT)/DHSUB-1.0)
        ENDIF
        IF(HMAVG.GE.HGSAT) THEN
          XFL=1.0
          EPS=1.0
          SLIP=1.0
          HLAVG=0.0
        ELSE
          IF(ISUBM.EQ.1) THEN
*           Use the Paul Gallet thesis model.
            PI=RHOLAV*CPLAV*(TSCLAD-TL)/(RHOG*HLV)
            XFL=XFL0+PCH*PHI*DZ/(MFLOW*ACOOL*HLV)/(1.0+PI)
          ELSE IF(ISUBM.EQ.2) THEN
*           Use a profile fit model.
            XFL=MAX(XFL0,(F2+F3)/(1.0+F3))
          ENDIF
          HLAVG=MIN(HLSAT,(HMAVG-XFL*HGSAT)/(1.0-XFL))
        ENDIF
*----
*  RECOMPUTE THE LIQUID PROPERTIES
*----
        IF(HLAVG.GT.0.0) THEN
          CALL THMPH(IFLUID,PINLET,HLAVG,RHOL,TL)
          IF(IFLUID.EQ.0) THEN
            CALL THMPT(PINLET,TL,R1,R2,ZKL,ZMUL,CPL)
          ELSE IF(IFLUID.EQ.1) THEN
            CALL THMHPT(PINLET,TL,R1,R2,ZKL,ZMUL,CPL)
          ENDIF
*----
*  COMPUTE THE COOLANT VOID FRACTION AND SLIP RATIO
*  A drift-flux model is proposed by means of the concentration
*  parameter CO and the drift velocity VGJ (their correspondent 
*  correlations are supposed to work fine under different flow regimes).
*----
! Condenser cette étape à la suivante en distinguant trois cas 
! One pahse liquid (garder l'existant)
! Two-phase flow
! Superheated steam 
! OUUUUUU
! Garder cette cellule pour obtenir une première valeur de epsilon 
! et ensuite dans le two-phase flow, itérer sur ce epsilon dans le recalcul de tous les paramètres hydrauliques
          IF(HGSAT.GT.HLSAT) THEN
            CO=1.13
            PR=PINLET/10**6
            SIGM=-7.2391E-6*PR**3+2.8345E-4*PR**2-5.1566E-3*PR+4.2324E-2
            VGJ=1.18*((SIGM*9.81*(RHOL-RHOG))/RHOL**2)**0.25
            F4=CO*((XFL*RHOL)+((1.0-XFL)*RHOG))+(RHOL*RHOG*VGJ/MFLOW)
            EPS=(XFL*RHOL)/F4
            JL=(1.0-XFL)*MFLOW/RHOL
            JG=XFL*MFLOW/RHOG
            IF(EPS.NE.0) SLIP=JG*(1.0-EPS)/(JL*EPS)
          ENDIF
        ELSE
*         superheated steam
          CALL THMPH(IFLUID,PINLET,HMAVG,RHOG,TCALO)
          IF(IFLUID.EQ.0) THEN
            CALL THMPT(PINLET,TCALO,R1,R2,ZKG,ZMUG,CPG)
          ELSE IF(IFLUID.EQ.1) THEN
            CALL THMHPT(PINLET,TCALO,R1,R2,ZKG,ZMUG,CPG)
          ENDIF
        ENDIF
      ENDIF
*----
*  COMPUTE THE FLUID PROPERTIES
*  RHO: fluid density
*  REL: Reynolds number of liquid phase
*  PRL: Prandtl number of liquid phase
*----
      IF(XFL.EQ.0.0) THEN
*       One phase liquid
        TB=TSAT-0.1
        IF(TL.LT.TB) THEN
          TCALO=TL
        ELSE
          TCALO=TB
        ENDIF
        IF(IFLUID.EQ.0) THEN
          CALL THMPT(PINLET,TCALO,R1,R2,ZKONE,ZMUONE,CPONE)
        ELSE IF(IFLUID.EQ.1) THEN
          CALL THMHPT(PINLET,TCALO,R1,R2,ZKONE,ZMUONE,CPONE)
        ENDIF
        RHO=RHOLAV
        REL=MFLOW*HD/ZMUONE
        PRL=ZMUONE*CPONE/ZKONE
      ELSE IF(HMAVG.LT.HGSAT) THEN
*       Two-phase flow
! remplacer ça par une référence à mon tout joli algo DFM ? 
        TCALO=EPS*TSAT+(1.0-EPS)*TL
        ZKONE=ZKL
        CPONE=CPL
        RHO=EPS*RHOG+(1.0-EPS)*RHOL
        REL=MFLOW*(1.0-XFL)*HD/ZMUL
        PRL=ZMUL*CPL/ZKL
      ELSE
*       superheated steam
        RHO=RHOG
        REL=MFLOW*HD/ZMUG
        PRL=ZMUG*CPG/ZKG
      ENDIF
*----
*  THERMAL EXCHANGE BETWEEN CLAD AND FLUID USING THE DITTUS AND BOELTER
*  CORRELATION (SINGLE PHASE) OR CHEN CORRELATION (SATURATED BOILING)
*----
      IF(IHCONV.EQ.0) THEN
        ITER=0
        KWA=99
        DO
          ITER=ITER+1
          IF(ITER.GT.50) THEN
            WRITE(HSMG,'(30HTHMH2O: HCONV FAILURE IN SLICE,I5,1H.)') K
            CALL XABORT(HSMG)
          ENDIF
          HA=0.023*(ZKONE/HD)*REL**0.8*PRL**0.4
          F=1.0
          S=1.0
          IF((XFL.EQ.XFL0).OR.(TSCLAD.LE.TSAT).OR.(KWA.EQ.0)) THEN
*           Single-phase convection. Use Dittus-Boelter correlation
            KWA=0
            HB=0.0
            K0=0
            XFL=XFL0
            EPS=EPS0
            SLIP=SLIP0
          ELSE IF(HMAVG.LT.HGSAT) THEN
*           Subcooled convection. Use Dittus-Boelter and Forster-Zuber
*           correlations
*           XM: Martinelli parameter
*           F: Reynolds number factor
*           S: nucleate boiling suppression factor
*           SIGM: surface tension in N/m
*           HA: Dittus-Boelter coefficient
*           HB: Forster-Zuber coefficient
*
            IF(HMAVG.LT.HLSAT) THEN
              KWA=1
            ELSE
              KWA=2
            ENDIF
            XM=(XFL/(1.0-XFL))**0.9*(RHOL/RHOG)**0.5*(ZMUG/ZMUL)**0.1
            IF(XM.LE.0.100207) THEN
              F=1.0
            ELSE
              F=2.35*(0.213+XM)**0.736
            ENDIF
            RE=REL*F**1.25
            S=1.0/(1.0+2.53E-6*RE**1.17)
            PR=PINLET/10**6
            SIGM=-7.2391E-6*PR**3+2.8345E-4*PR**2-5.1566E-3*PR+4.2324E-2
            HA=0.023*(ZKL/HD)*REL**0.8*PRL**0.4
            DTSAT=TSCLAD-TSAT
            IF(IFLUID.EQ.0) THEN
              CALL THMSAP(PW, TSCLAD)
            ELSE
              CALL THMHSP(PW, TSCLAD)
            ENDIF
            DP=PW-PINLET
*           Forster-Zuber equation
            HLV=HGSAT-HLSAT
            HB=0.00122*ZKL**0.79*CPL**0.45*RHOL**0.49/(ZMUL**0.29*
     >      SIGM**0.5*HLV**0.24*RHOG**0.24)*DTSAT**0.24*DP**0.75
          ELSE
*           Superheated steam. Use Mokry correlation
            KWA=3
            IF(IFLUID.EQ.0) THEN
              CALL THMPT(PINLET,TSCLAD,RHOW,R2,R3,R4,R5)
            ELSE IF(IFLUID.EQ.1) THEN
              CALL THMHPT(PINLET,TSCLAD,RHOW,R2,R3,R4,R5)
            ENDIF
            HA=0.0061*(ZKG/HD)*REL**0.904*PRL**0.684*(RHOW/RHO)**0.564
            HB=0.0
          ENDIF
*         Chen correlation
          HCONV=F*HA+S*HB
          IF(HCONV.LE.0.0) THEN
            WRITE(HSMG,'(34HTHMH2O: DRY OUT REACHED IN CHANNEL,3I5)')
     >      I,J,K
            CALL XABORT(HSMG)
          ENDIF
          IF(ITIME.EQ.0) THEN
            TWAL=(PHI+S*HB*TSAT+F*HA*TCALO)/(S*HB+F*HA)
          ELSE
            ZNUM=ZF(1)+RADCL*S*HB*TSAT+RADCL*F*HA*TCALO
            ZDEN=ZF(2)+RADCL*S*HB+RADCL*F*HA
            TWAL=MAX(273.15,ZNUM/ZDEN)
            PHI=MAX(0.0,(ZF(1)-TWAL*ZF(2))/RADCL)
          ENDIF
          IF(ABS(TSCLAD-TWAL).GT.1.0E-5*TSCLAD) THEN
            TSCLAD=TWAL
          ELSE
            EXIT
          ENDIF
        ENDDO
      ELSE IF(IHCONV.EQ.1) THEN
        IF(ITIME.EQ.0) THEN
          TSCLAD=TCALO+PHI/KHCONV
        ELSE
          RCHC=RADCL*KHCONV
          TSCLAD=MAX(273.15,(ZF(1)+RCHC*TCALO)/(ZF(2)+RCHC))
          PHI=(ZF(1)-TSCLAD*ZF(2))/RADCL
        ENDIF
      ENDIF
*----
*  COMPUTE INITIAL BULK LIQUID ENTHALPY SUBCOOLING DHSUB
*----
      IF((ISUBM.GT.0).AND.(K0.EQ.0).AND.LFIRST) THEN
        DTSUB=0.0
        IF(ISUBM.EQ.1) THEN
*         Bowring correlation
*         Reference: R. W. Bowring, "Physical Model, Based on Bubble
*         Detachment, and Calculation of Steam Voidage in the Subcooled
*         Region of a Heated Channel," OECD Report HPR-10, 1962.
*         Equation3 (3) and (17)
          VC=MFLOW/RHOL
          ETA=14.0+0.1*PINLET/1.01325E+05
          DTSUB=ETA*PHI/VC*1.0E-6
         ELSE IF(ISUBM.EQ.2) THEN
*         Saha-Zuber subcooling model
*         PE: Peclet number
          PE=MFLOW*CPL*HD/ZKL
          IF(PE.LE.70000.0) THEN
            DTSUB=PHI*HD/(455.0*ZKL)
          ELSE
*           reactor conditions
            DTSUB=154.0*PHI/(MFLOW*CPL)
          ENDIF
        ENDIF
        IF(TCALO.GE.TSAT-DTSUB) K0=K
        DSAT=TSCLAD-TCALO-DTSUB
        DHSUB=CPL*DTSUB
        LFIRST=.FALSE.
      ENDIF
      RETURN
      END
