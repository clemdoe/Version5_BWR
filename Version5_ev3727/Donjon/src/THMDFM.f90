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
*  COMPUTE DENSITIES
*----
      CALL THMPX(PINLET,0.0,RHOL,HLSAT,ZKL,ZMUL,CPL)
      CALL THMPX(PINLET,1.0,RHOG,HGSAT,ZKG,ZMUG,CPG)

      RHO = RHOL*(1 - EPS)+ EPS*RHOG
    
*----
*  COMPUTE PHASES VELOCITIES AND REYNOLDS
*----
      VLIQ = VCOOL - (1/(1- EPS) - RHOLAV/RHO) *VGJprime
      VVAP = VCOOL + RHOLAV/RHO * VGJprime
      ZMU = (ZMUL*ZMUG/ (ZMUL*(1-EPS) + ZMUG*EPS))
      REY = RHO * ABS(VCOOL) * HD / ZMU

*----
*  COMPUTE FLOW QUALITY
*----

      IF (HLSAT. GT. HMAVG) THEN 
        XFL = 0
      ELSE IF (HMAVG. GT. HGSAT) THEN
        XFL = 1
      ELSE     
        XFL = (HMAVG - HLSAT)/(HGSAT - HLSAT)
      ENDIF

*----
*  COMPUTE VGJ, VGJprime and C0 AFTER CHOSEN CORRELATION
*----
      CORREL = 'EPRI'
      CALL THMVGJ(VCOOL, DCOOL, PINLET, ZMU, XFL, HD, RHOG, RHOL, EPS, CORREL, VGJ, C0)

      VGJprime = VGJ + (C0-1)*VCOOL

*----
*  COMPUTE HLV
*----
      HLV=HGSAT-HLSAT
*----
*  COMPUTE NEW EPS VALUE
*----
      IF (XFL.EQ.0) THEN 
        EPS = 0
      ELSE IF (XFL.EQ.1) THEN
        EPS = 1
      ELSE 
        EPS = XFL / (C0 * (XFL + (RHOG/RHOL) * (1 - XFL)) + (RHOG * VGJ) / (RHOL * VCOOL))
      ENDIF
    
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
END
