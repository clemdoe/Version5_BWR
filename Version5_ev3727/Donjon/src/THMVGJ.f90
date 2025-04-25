SUBROUTINE THMVGJ(VCOOL, DCOOL, PCOOL, MUT, XFL, HD, RHOG, RHOL, EPS, CORREL, VGJ, C0)
!
!-----------------------------------------------------------------------
!
!Purpose:
! Update the concentration parameter CO and the drift velocity VGJ 
! in the THM model after several correlations to implement the drift flux model 
! in the THM code
!
!Copyright:
! Copyright (C) 2025 Ecole Polytechnique de Montreal
!
!Author(s): M. Bellier
! 04/2025: M.Bellier - Creation
!
!Parameters: input
! XFL     quality of the fluid in the channel
! DCOOL   density of the fluid in the channel
! VCOOL   velocity of the fluid in the channel
! PCOOL   pressure of the fluid in the channel
! MUT     dynamic viscosity of the fluid in the channel
! HD      hydraulic diameter of the channel
! RHOG    density of the vapour in given thermohydraulic conditions
! RHOL    density of the liquid in given thermohydraulic conditions
! ESP     void fraction of the fluid
! CORREL  correlation used to compute VGJ and C0
!
!Parameters: output
! VGJ     drift velocity
! C0      concentration parameter
!
!-----------------------------------------------------------------------
!
    USE GANLIB
    IMPLICIT NONE
!----
!   SUBROUTINE ARGUMENTS
!----
    REAL VCOOL, DCOOL, PCOOL, MUT, XFL, HD 
    REAL EPS, RHOG, RHOL
    CHARACTER CORREL*10
!----
!   LOCAL VARIABLES
!----
    REAL g
    REAL C1, k1, k0, r, PR, SIGM
    REAL VGJ, C0
    REAL REY 
    INTEGER PC

    REY = ABS(VCOOL*DCOOL) * (1.0 - XFL) * HD / MUT !Reynolds
    g =  9.81 !gravity
    PR=PCOOL/10**6 ! PCOOL ou autre valeur de P ? initialement Pinlet
    SIGM=-7.2391E-6*PR**3+2.8345E-4*PR**2-5.1566E-3*PR+4.2324E-2
    !sigma = IAPWS97(P = self.P[i]*(10**(-6)), x = 0).sigma
!----
!   VGJ AND C0 CALCULATION
!----
IF (RHOG.EQ.0) THEN
    C0=0
    VGJ=0
ELSE IF (RHOL.EQ.0) THEN
    C0=0
    VGJ=0
ELSE IF (CORREL.EQ.'HEM1') THEN
    VGJ = 0
    C0 = 1
    
ELSE IF (CORREL.EQ.'CHEXAL') THEN
! Correlation used in previous codes, after the work of Sarra Zoghlami for CANDU reactors
    C0=1.13
    VGJ=1.18*((SIGM*9.81*(RHOL-RHOG))/(RHOL**2))**0.25
    

ELSE IF (CORREL.EQ.'GERAMP') THEN
    PRINT *, 'THMVGJ : GERAMP CORREL USED'
    IF (SIGM.EQ.0) THEN
        VGJ = 0
    ELSE 
        VGJ = (g*SIGM*(RHOL-RHOG)/(RHOG**2))**0.25  
        IF (EPS.GT.0.65) THEN
            VGJ= VGJ*(2.9/0.35)*(1-EPS) 
            C0= 1 + (0.1/0.35)*(1-EPS)
        ELSE
            VGJ = 2.9*VGJ
            C0= 1.1 
        ENDIF
    ENDIF

ELSE IF (CORREL.EQ.'EPRI') THEN
    PRINT *, 'THMVGJ : EPRI CORREL USED'
    VGJ= ((2**0.5)*g*SIGM*(RHOL-RHOG)/(RHOL**2))**0.25 * ((1+EPS)**1.5)
    PC = 22060000
    C1 = (4 * (PC**2))/(PCOOL*(PC - PCOOL))
    k1 = MIN(0.8, 1/(1 + exp(-REY /60000)))
    k0 = k1 + (1-k1) * (RHOG / RHOL)**2
    r = (1+1.57*(RHOG/RHOL))/(1-k1)
    IF (EPS.GT.0) THEN
        C0 = (k0 + (1 - k0)*(EPS**r)*(1 - exp((-1)*C1))/(1 - exp((-1)*C1 * EPS)))**(-1)
    ENDIF

ELSE IF (CORREL.EQ.'MODBESTION') THEN
    PRINT *, 'THMVGJ : MODBESTION CORREL USED'
    VGJ =  0.188 * (((RHOL - RHOG) * g * HD ) / RHOG )*0.5
    C0 = 1.2 - 0.2*(RHOG/RHOL)**0.5

ELSE 
    PRINT *, 'Unknow correlation model, HEM1 used by default'
    VGJ = 0
    C0 = 1
ENDIF
RETURN
END