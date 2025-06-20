*DECK THM
      SUBROUTINE THM(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Simplified thermal-hydraulics module.
*
*Copyright:
* Copyright (C) 2013 Ecole Polytechnique de Montreal.
*
*Author(s): 
* A. Hebert, P. Gallet and V. Salino
* 02/2025: C. HUET - Modifications to include pressure drop calculation
*
*Parameters: input
* NENTRY  number of data structures transfered to this module.
* HENTRY  name of the data structures.
* IENTRY  data structure type where:
*         IENTRY=1 for LCM memory object;
*         IENTRY=2 for XSM file;
*         IENTRY=3 for sequential binary file;
*         IENTRY=4 for sequential ASCII file.
* JENTRY  access permission for the data structure where:
*         JENTRY=0 for a data structure in creation mode;
*         JENTRY=1 for a data structure in modifications mode;
*         JENTRY=2 for a data structure in read-only mode.
* KENTRY  data structure pointer.
*
*Comments:
* The THM: module specification is:
* THERMO MAPFL := THM: [ THERMO ] MAPFL :: (descthm) ;
* where
*   THERMO : name of the \emph{thermo) object that will be created or updated 
*     by the THM: module. Object \emph{thermo} contains thermal-hydraulics 
*     information set or computed by THM: in transient or in permanent 
*     conditions such as the distribution of the enthalpy, the pressure, the 
*     velocity, the density and the temperatures of the coolant for all the 
*     channels in the geometry. It also contains all the values of the fuel 
*     temperatures in transient or in permanent conditions according to the 
*     discretisation chosen for the fuel rods.
*   MAPFL : name of the \emph{map} object containing fuel regions description 
*     and local parameter informations.
*   (descthm) : structure describing the input data to the THM: module. 
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40,IOUT=6,PI=3.141592654,ZKILO=1.0E3)
      PARAMETER(IMAXO=1000,JMAXO=100,KMAXO=100,NMAXO=40,MAXRAD=10)
      PARAMETER(DTEMPR=5.0,DTEMPT=40.0,DPRESS=4.0)
      CHARACTER TEXT*40,TEXT12*12,HSIGN*12,PNAME*12,TXTDIR*12,HSMG*131,
     > UCONDF*12,UCONDC*12,SNAME*32,SCOMP*32,FNAME*32,FCOMP*32
      INTEGER ISTATE(NSTATE),TIMEIT,ITIME
      REAL STATE(NSTATE),DTIME,KHGAP,KHCONV,WTEFF
      REAL POULET,HX(IMAXO),HY(JMAXO),HZ(KMAXO)
      REAL RPRAD(MAXRAD),FPRAD(MAXRAD),TERP(MAXRAD)
      DOUBLE PRECISION DFLOT,DSUM
      LOGICAL LPRAD
      TYPE(C_PTR) IPTHM,IPMAP,JPMAP,KPMAP,JPTHM,KPTHM,LPTHM,MPTHM,
     > KPTHMI,LPTHMI,MPTHMI
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NUM,IREFSC
      REAL, ALLOCATABLE, DIMENSION(:) :: XX,YY,ZZ,BURN,BURN2,PW,FRO,
     1 FNFU,FNTG,FRACPU
      REAL, ALLOCATABLE, DIMENSION(:) :: FPOWER,KCONDF,KCONDC,TIMESR,
     1 TPOWER,PFORM,DTERP
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: XBURN,POW,TCOMB,DCOOL,
     1 TCOOL,TSURF,PCOOL,HCOOL
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: RAD
      DOUBLE PRECISION ARF,ARCI,ARCE,DARF,DARC
      REAL, ALLOCATABLE, DIMENSION(:,:) :: VAL,RVAL
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.NE.2)CALL XABORT('@THM: 2 PARAMETERS EXPECTED.')
      IPTHM=KENTRY(1)
      IPMAP=KENTRY(2)
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2))CALL XABORT('@THM:'
     1 //' LCM OBJECT EXPECTED AT FIRST LHS.')
      IF((IENTRY(2).NE.1).AND.(IENTRY(2).NE.2))CALL XABORT('@THM:'
     1 //' LCM OBJECT EXPECTED AT SECOND LHS.')
      CALL LCMGTC(IPMAP,'SIGNATURE',12,HSIGN)
      IF(HSIGN.NE.'L_MAP')THEN
        TEXT=HENTRY(2)
        CALL XABORT('@THM: SIGNATURE OF '//TEXT//' IS '//HSIGN//
     1  '. L_MAP EXPECTED.')
      ENDIF
*----
*  RECOVER L_MAP STATE-VECTOR
*----
      CALL LCMGET(IPMAP,'STATE-VECTOR',ISTATE)
      NB=ISTATE(1)
      NCH=ISTATE(2)
      NPARM=ISTATE(8)
      NSIMS=ISTATE(13)
      ALLOCATE(FNFU(NCH),FNTG(NCH),FRACPU(NCH))
*----
*  READ DATA
*----
      IMPX=1
      ITIME=0
      DTIME=0.0
      FPUISS=0.974
      CFLUX=2.0E+6
      SPEED=0.0
      TINLET=0.0
      POULET=0.0
      POROS=0.05
      ICONDF=0
      ICONDC=0
      IHGAP=0
      IHCONV=0
      IFRCDI=0
      ISUBM=1
      RC=0.0
      RIG=0.0
      RGG=0.0
      RTG=0.0
      PITCH=0.0
      MAXIT1=50
      MAXIT2=50
      MAXIT3=50
      IFLUID=0
      IFUEL=0
      IGAP=0
      IPRES=0
      IDFM=0
      ERMAXT=1.0
      ERMAXC=1.0E-3
      NFD=5
      NDTOT=8
      NPRAD=0
      TIMEIT=0
      NPOWER=0
      RELAX=1.0
      RTIME=0.0
      WTEFF=5.0/9.0 ! Rowlands weighting factor
      EPSR=0.0
      THETA=0.0
      UCONDF='CELSIUS'
      UCONDC='CELSIUS'
      TPOW=0.0
      FNFU(:NCH)=1.0
      FNTG(:NCH)=0.0
      FRACPU(:NCH)=0.0
      IF(JENTRY(1).EQ.1) THEN
        CALL LCMGET(IPTHM,'STATE-VECTOR',ISTATE)
        IF(ISTATE(1).NE.NCH) CALL XABORT('THM: INVALID STATE VECTOR FO'
     >  //'R IPTHM OPJECT.')
        MAXIT1=ISTATE(3)
        MAXIT2=ISTATE(4)
        MAXIT3=ISTATE(5)
        NFD=ISTATE(6)
        NDTOT=ISTATE(7)
        ITIME=ISTATE(8)
        TIMEIT=ISTATE(9)
        IHGAP=ISTATE(10)
        IHCONV=ISTATE(11)
        ICONDF=ISTATE(12)
        ICONDC=ISTATE(13)
        IFRCDI=ISTATE(14)
        ISUBM=ISTATE(15)
        IF(ICONDF.EQ.1) NCONDF=ISTATE(16)
        IF(ICONDC.EQ.1) NCONDC=ISTATE(17)
        NPRAD=ISTATE(18)
        IFLUID=ISTATE(20)
        IGAP=ISTATE(21)
        IPRES=ISTATE(22)
        IDFM=ISTATE(23)
        CALL LCMGET(IPTHM,'REAL-PARAM',STATE)
        DTIME=STATE(1)
        FPUISS=STATE(2)
        CFLUX=STATE(3)
        SPEED=STATE(4)
        POULET=STATE(5)
        TINLET=STATE(6)
        POROS=STATE(7)
        RC=STATE(8)
        RIG=STATE(9)
        RGG=STATE(10)
        RTG=STATE(11)
        PITCH=STATE(12)
        ERMAXT=STATE(13)
        ERMAXC=STATE(14)
        RELAX=STATE(15)
        RTIME=STATE(16)
        IF(IHGAP.EQ.1) KHGAP=STATE(17)
        IF(IHCONV.EQ.1) KHCONV=STATE(18)
        WTEFF=STATE(19)
        TPOW=STATE(20)
        EPSR=STATE(22)
        THETA=STATE(23)
*----
*  RECOVER CELL-DEPENDENT DATA
*----
        CALL LCMGET(IPTHM,'NB-FUEL',FNFU)
        CALL LCMGET(IPTHM,'NB-TUBE',FNTG)
        CALL LCMGET(IPTHM,'FRACT-PU',FRACPU)
*----
*  RECOVER CONDUCTIVITY INFORMATION ON LCM OBJECT THM
*----
        IF(ICONDF.EQ.1) THEN
          ALLOCATE(KCONDF(NCONDF+3))
          CALL LCMGET(IPTHM,'KCONDF',KCONDF)
          CALL LCMGTC(IPTHM,'UCONDF',12,UCONDF)
        ENDIF
        IF(ICONDC.EQ.1) THEN
          ALLOCATE(KCONDC(NCONDC+1))
          CALL LCMGET(IPTHM,'KCONDC',KCONDC)
          CALL LCMGTC(IPTHM,'UCONDC',12,UCONDC)
        ENDIF
      ENDIF
*----
*  READ INPUT DATA
*----
      IPICK=0
      LPRAD=.FALSE.
   10 CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.EQ.10)GO TO 60
   20 IF(ITYP.NE.3)CALL XABORT('@THM: CHARACTER DATA EXPECTED.')
      IF(TEXT.EQ.'EDIT') THEN
*       Read printing index
        CALL REDGET(ITYP,IMPX,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@THM: INTEGER FOR EDIT EXPECTED.')
      ELSE IF(TEXT.EQ.'TIME') THEN
*       Time at beginning of time-step (s).
        CALL REDGET(ITYP,NITMA,RTIME,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR RTIME EXPECTED.')
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.EQ.2) THEN
          DTIME=FLOT
        ELSE IF(ITYP.EQ.3) THEN
          GO TO 20
        ELSE
          CALL XABORT('@THM: REAL FOR DTIME EXPECTED.')
        ENDIF
        ITIME=1
      ELSE IF(TEXT.EQ.'FLUID') THEN
*       Read fluid type
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.3)CALL XABORT('@THM: CHARACTER FOR FLUID EXPECTED.')
        IF(TEXT.EQ.'H2O') THEN
           IFLUID=0
        ELSE IF(TEXT.EQ.'D2O') THEN
           IFLUID=1
        ELSE IF(TEXT.EQ.'SALT') THEN
           IFLUID=2
           CALL REDGET(ITYP,NITMA,FLOT,SNAME,DFLOT)
           IF(ITYP.NE.3) THEN
             CALL XABORT('@THM: CHARACTER FOR FLUID SALT NAME EXPECTED'
     >       //'.')
           ENDIF
           CALL REDGET(ITYP,NITMA,FLOT,SCOMP,DFLOT)
           IF(ITYP.NE.3) THEN
             CALL XABORT('@THM: CHARACTER FOR FLUID SALT COMPOSITION'
     >       //'EXPECTED.')
           ENDIF
        ELSE
           CALL XABORT('@THM: INVALID FLUID TYPE.')
        ENDIF
      ELSE IF(TEXT.EQ.'FUEL') THEN
*       Read fuel type
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.3)CALL XABORT('@THM: CHARACTER FOR FLUID EXPECTED.')
        IF(TEXT.EQ.'UO2') THEN
           IFUEL=0
        ELSE IF(TEXT.EQ.'SALT') THEN
           IFUEL=1
           CALL REDGET(ITYP,NITMA,FLOT,FNAME,DFLOT)
           IF(ITYP.NE.3) THEN
             CALL XABORT('@THM: CHARACTER FOR FLUID EXPECTED.')
           ENDIF
           CALL REDGET(ITYP,NITMA,FLOT,FCOMP,DFLOT)
           IF(ITYP.NE.3) THEN
             CALL XABORT('@THM: CHARACTER FOR FLUID EXPECTED.')
           ENDIF
        ELSE
           CALL XABORT('@THM: INVALID FUEL TYPE.')
        ENDIF
      ELSE IF(TEXT.EQ.'FPUISS') THEN
*       Coolant power factor
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.EQ.2) THEN
           FPUISS=FLOT
        ELSE
           CALL XABORT('@THM: REAL FOR FPUISS EXPECTED.')
        ENDIF
      ELSE IF(TEXT.EQ.'CRITFL') THEN
*       Critical heat flux (W/m^2)
        CALL REDGET(ITYP,NITMA,CFLUX,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR CFLUX EXPECTED.')
      ELSE IF(TEXT.EQ.'CWSECT') THEN
*       Core coolant section (m^2)
        CALL REDGET(ITYP,NITMA,CWSECT,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR CWSECT EXPECTED.')
*       Coolant flow (m^3/h)
        CALL REDGET(ITYP,NITMA,FLOW,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR FLOW EXPECTED.')
        SPEED=FLOW/(3600.0*CWSECT)
      ELSE IF(TEXT.EQ.'INLET-Q') THEN
*       Core coolant section (m^2)
        IF((POULET.EQ.0.0).OR.(TINLET.EQ.0.0)) CALL XABORT('@THM: INLE'
     >  //'T INFORMATION NOT SET BEFORE USING INLET-Q.')
        CALL REDGET(ITYP,NITMA,CWSECT,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR CWSECT EXPECTED.')
*       Inlet mass flow rate (kg/s)
        CALL REDGET(ITYP,NITMA,QFLUID,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR QFLUID EXPECTED.')
        IF(IFLUID.EQ.0) THEN
          CALL THMPT(POULET,TINLET,RHOL,R2,R3,R4,R5)
        ELSE IF(IFLUID.EQ.1) THEN
          CALL THMHPT(POULET,TINLET,RHOL,R2,R3,R4,R5)
        ELSE IF(IFLUID.EQ.2) THEN
          CALL THMSPT(SNAME,SCOMP,TINLET,RHOL,R2,R3,R4,R5,IMPX)
        ENDIF
        SPEED=QFLUID/(CWSECT*RHOL)
      ELSE IF(TEXT.EQ.'SPEED') THEN
*       Coolant velocity (m/s)
        CALL REDGET(ITYP,NITMA,SPEED,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR SPEED EXPECTED.')
      ELSE IF(TEXT.EQ.'INLET') THEN
*       The POULET and TINLET informations are used to compute initial
*       enthalpy and water density.
*       Outlet pressure (Pa)
        CALL REDGET(ITYP,NITMA,POULET,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR POULET EXPECTED.')
*       Inlet temperature (K)
        CALL REDGET(ITYP,NITMA,TINLET,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR TINLET EXPECTED.')
      ELSE IF(TEXT.EQ.'PUFR') THEN
        ICONDF=0
*       Plutonium mass enrichment
        CALL THMINP('PUFR',NCH,FRACPU)
      ELSE IF(TEXT.EQ.'POROS') THEN
        ICONDF=0
*       Oxyde porosity
        CALL REDGET(ITYP,NITMA,POROS,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR POROS EXPECTED.')
      ELSE IF(TEXT.EQ.'CONDF') THEN
        IF(ICONDF.EQ.1)DEALLOCATE(KCONDF)
        ICONDF=1
*       Fuel conductivity expressed as a function of fuel temperature
*       (function = polynomial + inverse term)
        CALL REDGET(ITYP,NCONDF,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@THM: INTEGER FOR CONDF EXPECTED.')
        IF(NCONDF.LT.0)CALL XABORT('@THM: NCONDF MUST BE LARGER OR '
     >                  //'EQUAL TO 0.')
        ALLOCATE(KCONDF(NCONDF+3))
        DO I=1,NCONDF+1
          CALL REDGET(ITYP,NITMA,KCONDF(I),TEXT,DFLOT)
          IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR KCONDF EXPECTED.')
        ENDDO
        CALL REDGET(ITYP,NITMA,FLOT,TEXT12,DFLOT)
        IF(ITYP.NE.3)CALL XABORT('@THM: CHARACTER DATA EXPECTED (INV, '
     >     //'CELSIUS OR KELVIN) IN CONDF STATEMENT.')
        IF(TEXT12.EQ.'INV') THEN
          CALL REDGET(ITYP,NITMA,KCONDF(NCONDF+2),TEXT,DFLOT)
          IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR INV EXPECTED.')
          CALL REDGET(ITYP,NITMA,KCONDF(NCONDF+3),TEXT,DFLOT)
          IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR REF EXPECTED.')
          CALL REDGET(ITYP,NITMA,FLOT,TEXT12,DFLOT)
        ELSE
          KCONDF(NCONDF+2)=0.0 ! Coefficient for the inverse term
          KCONDF(NCONDF+3)=-273.15 ! Reference for the inverse term
        ENDIF
        IF((TEXT12.NE.'CELSIUS').AND.(TEXT12.NE.'KELVIN')) THEN
          CALL XABORT('@THM: UNIT KEYWORD EXPECTED (CELSIUS OR '
     >                //'KELVIN) IN CONDF STATEMENT.')
        ENDIF
        UCONDF=TEXT12
      ELSE IF(TEXT.EQ.'CONDC') THEN
        IF(ICONDC.EQ.1)DEALLOCATE(KCONDC)
        ICONDC=1
*       Clad conductivity expressed as a polynomial of clad temperature
        CALL REDGET(ITYP,NCONDC,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@THM: INTEGER FOR CONDC EXPECTED.')
        IF(NCONDC.LT.0)CALL XABORT('@THM: NCONDC MUST BE LARGER OR '
     >                  //'EQUAL TO 0.')
        ALLOCATE(KCONDC(NCONDC+1))
        DO I=1,NCONDC+1
          CALL REDGET(ITYP,NITMA,KCONDC(I),TEXT,DFLOT)
          IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR KCONDC EXPECTED.')
        ENDDO
        CALL REDGET(ITYP,NITMA,FLOT,UCONDC,DFLOT)
        IF((ITYP.NE.3).OR.((UCONDC.NE.'CELSIUS').AND.
     >                     (UCONDC.NE.'KELVIN'))) THEN
          CALL XABORT('@THM: UNIT KEYWORD EXPECTED (CELSIUS OR '
     >                //'KELVIN) IN CONDC STATEMENT.')
        ENDIF
      ELSE IF(TEXT.EQ.'HGAP') THEN
        IHGAP=1
*       Fixed, user-chosen value of the HGAP (heat exchange coefficient
*       of the gap) (W/m^2/K)
        CALL REDGET(ITYP,NITMA,KHGAP,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR HGAP EXPECTED.')
      ELSE IF(TEXT.EQ.'HCONV') THEN
        IHCONV=1
*       Fixed, user-chosen value of the HCONV (heat transfer coefficient
*       between clad and fluid) (W/m^2/K)
        CALL REDGET(ITYP,NITMA,KHCONV,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR HCONV EXPECTED.')
      ELSE IF(TEXT.EQ.'TEFF') THEN
*       Surface temperature's weighting factor in effective fuel
*       temperature
        CALL REDGET(ITYP,NITMA,WTEFF,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR TEFF EXPECTED.')
      ELSE IF(TEXT.EQ.'FORCEAVE') THEN
*       Force the use of the average value approximation for fuel
*       conductivity
        IFRCDI=1
      ELSE IF(TEXT.EQ.'MONO') THEN
*       one-phase flow model
        ISUBM=0
      ELSE IF(TEXT.EQ.'BOWR') THEN
*       Bowring's correlation
        ISUBM=1
      ELSE IF(TEXT.EQ.'SAHA') THEN
*       Saha-Zuber correlation
        ISUBM=2
      ELSE IF(TEXT.EQ.'RADIUS') THEN
*       Fuel pellet radius (m)
        CALL REDGET(ITYP,NITMA,RC,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR RC EXPECTED.')
*       Internal clad rod radius (m)
        CALL REDGET(ITYP,NITMA,RIG,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR RIG EXPECTED.')
*       External clad rod radius (m)
        CALL REDGET(ITYP,NITMA,RGG,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR RGG EXPECTED.')
*       Guide tube radius (m)
        CALL REDGET(ITYP,NITMA,RTG,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR RTG EXPECTED.')
      ELSE IF(TEXT.EQ.'ASSMB') THEN
*       Number of active fuel rods
        CALL THMINP('NB-FUEL',NCH,FNFU)
*       Number of guide tubes
        CALL THMINP('NB-TUBE',NCH,FNTG)
      ELSE IF(TEXT.EQ.'CLUSTER') THEN
*       Hexagonal pitch (m)
        CALL REDGET(ITYP,NITMA,PITCH,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR PITCH EXPECTED.')
*       Number of active fuel pins in cluster
        CALL THMINP('NB-FUEL',NCH,FNFU)
      ELSE IF(TEXT.EQ.'CONV') THEN
*       Number of conduction iterations
        CALL REDGET(ITYP,MAXIT1,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@THM: INTEGER FOR MAXIT1 EXPECTED.')
*       Number of center-pellet iterations
        CALL REDGET(ITYP,MAXIT2,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@THM: INTEGER FOR MAXIT2 EXPECTED.')
*       Number of flow iterations
        CALL REDGET(ITYP,MAXIT3,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@THM: INTEGER FOR MAXIT3 EXPECTED.')
*       Temperature maximum error (K)
        CALL REDGET(ITYP,NITMA,ERMAXT,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR ERMAXT EXPECTED.')
*       maximum relative error for the calculation of the properties
*       in the coolant (pressure, enthalpy, density, velocity,...)
        CALL REDGET(ITYP,NITMA,ERMAXC,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR ERMAXC EXPECTED.')
      ELSE IF(TEXT.EQ.'RODMESH') THEN
*       Number of discretisation points in fuel
        CALL REDGET(ITYP,NFD,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@THM: INTEGER FOR NFD EXPECTED.')
*       Number of discretisation points in fuel rod (fuel+cladding)
        CALL REDGET(ITYP,NDTOT,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@THM: INTEGER FOR NDTOT EXPECTED.')
      ELSE IF(TEXT.EQ.'RELAX') THEN
*       Relaxation parameter
        CALL REDGET(ITYP,NITMA,RELAX,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR RELAX EXPECTED.')
      ELSE IF(TEXT.EQ.'RAD-PROF') THEN
*       Set radial power profile
        NPRAD=0
        LPRAD=.TRUE.
   30   CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.EQ.3)GO TO 20
        NPRAD=NPRAD+1
        RPRAD(NPRAD)=FLOT
        IF(NPRAD.GT.MAXRAD) CALL XABORT('@THM: MAXRAD OVERFLOW.')
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR RAD-PROF-X EXPECTED.')
        IF(RPRAD(NPRAD).LT.0.0)CALL XABORT('@THM: R TOO SMALL.')
        IF(RPRAD(NPRAD).GT.RC)CALL XABORT('@THM: R TOO LARGE.')
        CALL REDGET(ITYP,NITMA,FPRAD(NPRAD),TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR RAD-PROF-F EXPECTED.')
        GO TO 30
      ELSE IF(TEXT.EQ.'POWER-LAW') THEN
*       The total power in W generated in the fuel is defined as
*       T-POWER*TIME-LAW(t).
        CALL REDGET(ITYP,NITMA,TPOW,TEXT12,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR T-POWER EXPECTED.')
        CALL REDGET(ITYP,NPOWER,FLOT,TEXT12,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@THM: INTEGER VALUE EXPECTED.')
        ALLOCATE(TIMESR(NPOWER),TPOWER(NPOWER))
        DO I=1,NPOWER
          CALL REDGET(ITYP,NITMA,TIMESR(I),TEXT12,DFLOT)
          IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR TIME EXPECTED.')
          CALL REDGET(ITYP,NITMA,TPOWER(I),TEXT12,DFLOT)
          IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR POWER EXPECTED.')
        ENDDO
        CALL LCMPUT(IPTHM,'TIME-SR1',NPOWER,2,TIMESR)
        CALL LCMPUT(IPTHM,'POWER-SR1',NPOWER,2,TPOWER)
        DEALLOCATE(TPOWER,TIMESR)
      ELSE IF(TEXT.EQ.'F-RUG') THEN
*       Rugosity of the fuel rod
        CALL REDGET(ITYP,NITMA,EPSR,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR F-RUG EXPECTED.')
      ELSE IF(TEXT.EQ.'THETA') THEN
*       Angle of the fuel channel
        CALL REDGET(ITYP,NITMA,THETA,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR THETA EXPECTED.')
      ELSE IF(TEXT.EQ.'PDROP') THEN
*       Pressure drop identification
        CALL REDGET(ITYP,IPRES,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@THM: INTEGER FOR IPRES EXPECTED.')
      ELSE IF(TEXT.EQ.'DFM') THEN
*       Drift Flux Model identification
        CALL REDGET(ITYP,IDFM,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@THM: INTEGER FOR IDFM EXPECTED.')
      ELSE IF(TEXT.EQ.'SET-LOCAL') THEN
*       Reset a global parameter
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.3)CALL XABORT('@THM: CHARACTER NAME EXPECTED.')
        CALL REDGET(ITYP,NITMA,VALUE,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@THM: REAL FOR VALUE EXPECTED.')
        JPMAP=LCMGID(IPMAP,'PARAM')
        DO 40 IPAR=1,NPARM
        KPMAP=LCMGIL(JPMAP,IPAR)
        CALL LCMGTC(KPMAP,'P-NAME',12,PNAME)
        CALL LCMGET(KPMAP,'P-TYPE',ITYPE)
        IF(ITYPE.EQ.1) THEN
          IF(PNAME.EQ.TEXT) THEN
            CALL LCMPUT(KPMAP,'P-VALUE',1,2,VALUE)
            IF(IMPX.GT.0) WRITE(6,500) PNAME,VALUE
            GO TO 10
          ELSE
            GO TO 40
          ENDIF
        ELSE IF(ITYPE.EQ.2) THEN
          CALL XABORT('@THM: CANNOT RESET LOCAL PARAMETER: '//TEXT)
        ENDIF
   40   CONTINUE
        CALL XABORT('@THM: GLOBAL PARAMETER NAME NOT FOUND: '//TEXT)
      ELSE IF(TEXT.EQ.';') THEN
        GO TO 60
      ELSE IF(TEXT.EQ.'PICK') THEN
        IPICK=1
        GO TO 60
      ELSE
        CALL XABORT('@THM: INVALID KEYWORD: '//TEXT//'.')
      ENDIF
      GO TO 10
*----
*  TEST DATA INPUT
*----
   60 IF(TINLET.LE.273.15) CALL XABORT('@THM: INLET TEMPERATURE MUST BE'
     > //' HIGHER THAN 273.15K.')
      IF(SPEED.EQ.0.0) CALL XABORT('@THM: ZERO COOLANT SPEED.')
      IF(POULET.EQ.0.0) CALL XABORT('@THM: ZERO OUTLET PRESSURE.')
      IF(RC.EQ.0.0) CALL XABORT('@THM: ZERO FUEL PELLET RADIUS.')
      IF(RIG.EQ.0.0) CALL XABORT('@THM: ZERO INTERNAL CLAD ROD RADIUS.')
      IF(RGG.EQ.0.0) CALL XABORT('@THM: ZERO EXTERNAL CLAD ROD RADIUS.')
      IF(NDTOT.GT.NMAXO) CALL XABORT('@THM: NFD OVERFLOW, TOO MANY FUE'
     > //'L DOMAINS')
      IF(NDTOT.LT.8) CALL XABORT('@THM: NDTOT MUST AT LEAST BE EQUAL T'
     > //'O 8')
      IF(NFD.LT.4) CALL XABORT('@THM: NFD MUST AT LEAST BE EQUAL TO 4')
      IF(NFD.GE.NDTOT) CALL XABORT('@THM: NFD MUST BE LOWER THAN NDTO'
     > //'T.')
      IF((RELAX.LE.0.0).OR.(RELAX.GT.1.0)) CALL XABORT('@THM: RELAX '
     >  //'PARAMETER EXPECTED BETWEEN 0<RELAX<=1.')
      IF((WTEFF.LT.0.0).OR.(WTEFF.GT.1.0)) CALL XABORT('@THM: WTEFF '
     >  //'PARAMETER EXPECTED BETWEEN 0<=WTEFF<=1.')
      IF(ITIME.EQ.1) RELAX=1.0
      IF((RC.NE.RIG).AND.(IFUEL.EQ.1)) CALL XABORT('@THM: WITH MOLTEN'
     > //' SALT FUEL INNER CLAD RADIUS MUST BE EQUAL TO FUEL RADIUS')
*----
*  PRINT CHANNEL-DEPENDENT DATA
*----
      IF(IMPX.GT.1) THEN
        WRITE(6,'(/28H THM: CHANNEL-DEPENDENT DATA)')
        I1=1
        DO I=1,(NCH-1)/8+1
          I2=I1+7
          IF(I2.GT.NCH) I2=NCH
          WRITE(6,'(//8H CHANNEL,8(I8,6X,1H|))') (J,J=I1,I2)
          WRITE(6,'(8H NB-FUEL,8(F10.2,4X,1H|))') (FNFU(J),J=I1,I2)
          WRITE(6,'(8H NB-TUBE,8(F10.2,4X,1H|))') (FNTG(J),J=I1,I2)
          WRITE(6,'(8H PUFR   ,8(1P,E13.4,2H |))') (FRACPU(J),J=I1,I2)
          I1=I1+8
        ENDDO
      ENDIF
*----
*  SET POWER DISTRIBUTION
*----
      ALLOCATE(FRO(NFD-1))
      IF(NPRAD.EQ.0) THEN
        FRO(:NFD-1)=1.0
      ELSE
        IF(.NOT.LPRAD) THEN
          CALL LCMGET(IPTHM,'RAD-PROF_R',RPRAD)
          CALL LCMGET(IPTHM,'RAD-PROF_F',FPRAD)
        ELSE
          CALL LCMPUT(IPTHM,'RAD-PROF_R',NPRAD,2,RPRAD)
          CALL LCMPUT(IPTHM,'RAD-PROF_F',NPRAD,2,FPRAD)
        ENDIF
        DAR1=0.0
        DELT=0.5*RC**2/REAL(NFD-1)
        DO IM=1,NFD-1
           DAR2=DAR1+DELT
           RADM=SQRT(DAR1+DAR2)
           CALL ALTERP(.FALSE.,NPRAD,RPRAD(1),RADM,.FALSE.,TERP(1))
           DSUM=0.0D0
           DO J=1,NPRAD
             DSUM=DSUM+TERP(J)*FPRAD(J)
           ENDDO
           FRO(IM)=REAL(DSUM)
           DAR1=DAR2
        ENDDO
      ENDIF
      IF(IMPX.GT.1) WRITE(6,480) (FRO(IM),IM=1,NFD-1)
*----
*  RECOVER GEOMAP STATE-VECTOR
*  ISTATE(1): 7 = XYZ, 9 = HEXZ
*  In 3d hexagonal, NY=0, but THM: expects a 3D geometry, so we set
*  NY=1 and  continue.
*---- 
      JPMAP=LCMGID(IPMAP,'GEOMAP')
      CALL LCMGET(JPMAP,'STATE-VECTOR',ISTATE)
      NX=ISTATE(3)
      NY=ISTATE(4)
      NZ=ISTATE(5)
      NEL=ISTATE(6)
      IF((ISTATE(1).EQ.9).AND.(NY.EQ.0)) NY=1
      IF(NX.GT.IMAXO) CALL XABORT('@THM: NX OVERFLOW.')
      IF(NY.GT.JMAXO) CALL XABORT('@THM: NY OVERFLOW.')
      IF(NZ.GT.KMAXO) CALL XABORT('@THM: NZ OVERFLOW.')
*----
*  RECOVER REACTOR MESH IN METER
*  The arrays HX, HY, and HZ contain the mesh size in X-, Y-, and
*  Z-direction and are used to determine the volume of a mesh, i.e.
*  V(I,J,K)=HX(I)*HY(J)*HZ(K)
*  For 3D hexagonal, set HX and HY to the square root of the SA surface
*  SASS
*----
      ALLOCATE(XX(NX+1),YY(NY+1),ZZ(NZ+1))
      IF(ISTATE(1).EQ.7) THEN
        CALL LCMGET(JPMAP,'MESHX',XX)
        CALL LCMGET(JPMAP,'MESHY',YY)
      ENDIF
      CALL LCMGET(JPMAP,'MESHZ',ZZ)
      IF(ISTATE(1).EQ.9) THEN
        CALL LCMGET(JPMAP,'SIDE',SIDE)
        SASS=1.5*SQRT(3.0)*SIDE*SIDE/1.0E4
        DO 70 I=1,NX
        HX(I) = SQRT(SASS)
   70   CONTINUE
        DO 80 I=1,NY
        HY(I) = SQRT(SASS)
   80   CONTINUE
      ELSE
        DO 90 I=1,NX
        HX(I)=(XX(I+1)-XX(I))/100.0
   90   CONTINUE
        DO 100 I=1,NY
        HY(I)=(YY(I+1)-YY(I))/100.0
  100   CONTINUE
      ENDIF
      DO 110 I=1,NZ
      HZ(I)=(ZZ(I+1)-ZZ(I))/100.0
  110 CONTINUE
      DO 120 I=1,NZ+1
      ZZ(I)=ZZ(I)/100.0
  120 CONTINUE
      CALL LCMPUT(IPTHM,'MESHZ',NZ+1,2,ZZ)
      DEALLOCATE(ZZ,YY,XX)
*----
*  RECOVER LOCAL PARAMETER INFORMATION FROM L_MAP OBJECT
*----
      ALLOCATE(NUM(NEL),BURN(NCH*NB),PW(NCH*NB))
      CALL LCMGET(IPMAP,'BMIX',NUM)
      CALL LCMLEN(IPMAP,'BURN-INST',ILONG,ITYLCM)
      IF(ILONG.EQ.NCH*NB) THEN
        CALL LCMGET(IPMAP,'BURN-INST',BURN)
      ELSE
        CALL LCMLEN(IPMAP,'BURN-BEG',ILONG,ITYLCM)
        IF(ILONG.NE.NCH*NB) CALL XABORT('@THM: MISSING BURNUP INFO ON '
     >  //'FUELMAP.')
        ALLOCATE(BURN2(NCH*NB))
        CALL LCMGET(IPMAP,'BURN-BEG',BURN)
        CALL LCMGET(IPMAP,'BURN-END',BURN2)
        DO I=1,NCH*NB
          BURN(I)=(BURN(I)+BURN2(I))/2.0
        ENDDO
        DEALLOCATE(BURN2)
      ENDIF
      CALL LCMLEN(IPTHM,'POWER-SR1',NPOWER,ITYLCM)
      IF(NPOWER.NE.0) THEN
*       USE POWER TIME LAW
        IF(IMPX.GT.0) WRITE(6,*) 'THM: T-POWER = ',TPOW,' W'
        IF(TPOW.EQ.0.0) CALL XABORT('@THM: T-POWER NOT DEFINED.')
        IF(NCH.NE.1) CALL XABORT('@THM: NCH=1 EXPECTED.')
        ALLOCATE(TIMESR(NPOWER),TPOWER(NPOWER),DTERP(NPOWER))
        CALL LCMGET(IPTHM,'TIME-SR1',TIMESR)
        CALL LCMGET(IPTHM,'POWER-SR1',TPOWER)
        IF(ITIME.EQ.0) THEN
          CALL ALTERP(.FALSE.,NPOWER,TIMESR(1),RTIME,.FALSE.,DTERP(1))
        ELSE
          IF(DTIME.EQ.0.0) CALL XABORT('@THM: DTIME NOT DEFINED.')
          CALL ALTERI(.FALSE.,NPOWER,TIMESR(1),RTIME,RTIME+DTIME,
     >    DTERP(1))
          DO J=1,NPOWER
            DTERP(J)=DTERP(J)/DTIME
          ENDDO
        ENDIF
        DPOW=0.0D0
        DO J=1,NPOWER
          DPOW=DPOW+DTERP(J)*TPOWER(J)
        ENDDO
        DPOW=DPOW*TPOW
        DEALLOCATE(DTERP,TPOWER,TIMESR)
        CALL LCMLEN(IPMAP,'AXIAL-FPW',ILONG,ITYLCM)
        IF(ILONG.NE.NB) CALL XABORT('THM: NO AXIAL-FPW ON THE FUELMAP')
        ALLOCATE(PFORM(NB))
        CALL LCMGET(IPMAP,'AXIAL-FPW',PFORM)
        DO I=1,NB
          PW(I)=DPOW*PFORM(I)*1.0E-3
        ENDDO
        DEALLOCATE(PFORM)
      ELSE
*       RECOVER POWER FROM FUELMAP
        CALL LCMGET(IPMAP,'BUND-PW',PW)
      ENDIF
      IF(IMPX.GT.2) THEN
        PTOT=0.0
        DO I=1,NCH*NB
          PTOT=PTOT+PW(I)
        ENDDO
      ENDIF
*----
*  REBUILD LOCAL PARAMETER INFORMATION FOR THM
*----
      ALLOCATE(IREFSC(NCH))
      CALL LCMLEN(IPMAP,'REF-SCHEME',ILONG,ITYLCM)
      IF(ILONG.EQ.NCH) THEN
        CALL LCMGET(IPMAP,'REF-SCHEME',IREFSC)
      ELSE
        IREFSC(:NCH)=1
      ENDIF
      ALLOCATE(XBURN(NZ,NX,NY),POW(NZ,NX,NY))
      XBURN(:NZ,:NX,:NY)=0.0
      POW(:NZ,:NX,:NY)=0.0
      ICH=0
      DO 165 IY=1,NY
      DO 160 IX=1,NX
      IEL=(IY-1)*NX+IX
      DO 130 IZ=1,NZ
      IF(NUM((IZ-1)*NX*NY+IEL).NE.0) GO TO 140
  130 CONTINUE
      GO TO 160
  140 ICH=ICH+1
      IB=0
      DO 150 IZ=1,NZ
      IF(NUM((IZ-1)*NX*NY+IEL).EQ.0) GO TO 150
      IB=IB+1
      IMA=(IB-1)*NCH+ICH
      IF(IREFSC(ICH).GT.0) THEN
        XBURN(IZ,IX,IY)=BURN(IMA)
        POW(IZ,IX,IY)=PW(IMA)*1.0E3
      ELSE
        XBURN(NZ-IZ+1,IX,IY)=BURN(IMA)
        POW(NZ-IZ+1,IX,IY)=PW(IMA)*1.0E3
      ENDIF
  150 CONTINUE
      IF(IB.NE.NB) CALL XABORT('@THM: INVALID NUMBER OF BUNDLES.')
  160 CONTINUE
  165 CONTINUE
      IF(ICH.NE.NCH) CALL XABORT('@THM: INVALID NUMBER OF CHANNELS.')
      DEALLOCATE(PW,BURN)
*----
*  RECOVER AVERAGE FUEL TEMPERATURE FIELD FROM THM OBJECT
*----
      ALLOCATE(TCOMB(NZ,NX,NY))
      TCOMB(:NZ,:NX,:NY)=0.0
      JPMAP=LCMGID(IPMAP,'PARAM')
      DO 220 IPAR=1,NPARM
      KPMAP=LCMGIL(JPMAP,IPAR)
      CALL LCMGTC(KPMAP,'P-NAME',12,PNAME)
      IF(PNAME.EQ.'T-FUEL') THEN
        CALL LCMGET(KPMAP,'P-TYPE',ITYPE)
        ALLOCATE(VAL(NCH,NB))
        IF(ITYPE.EQ.1) THEN
          IF(IMPX.GT.0) WRITE(6,510) 'GLOBAL',PNAME
          CALL LCMGET(KPMAP,'P-VALUE',FLOT)
          DO 175 ICH=1,NCH
          DO 170 IB=1,NB
          VAL(ICH,IB)=FLOT
  170     CONTINUE
  175     CONTINUE
        ELSE IF(ITYPE.EQ.2) THEN
          IF(IMPX.GT.0) WRITE(6,510) 'LOCAL',PNAME
          CALL LCMGET(KPMAP,'P-VALUE',VAL)
        ENDIF
        ICH=0
        DO 215 IY=1,NY
        DO 210 IX=1,NX
        IEL=(IY-1)*NX+IX
        DO 180 IZ=1,NZ
        IF(NUM((IZ-1)*NX*NY+IEL).NE.0) GO TO 190
  180   CONTINUE
        GO TO 210
  190   ICH=ICH+1
        IB=0
        DO 200 IZ=1,NZ
        IF(NUM((IZ-1)*NX*NY+IEL).EQ.0) GO TO 200
        IB=IB+1
        IF(IREFSC(ICH).GT.0) THEN
          TCOMB(IZ,IX,IY)=VAL(ICH,IB)
        ELSE
          TCOMB(NZ-IZ+1,IX,IY)=VAL(ICH,IB)
        ENDIF
  200   CONTINUE
  210   CONTINUE
  215   CONTINUE
        DEALLOCATE(VAL)
      ENDIF
  220 CONTINUE
      DEALLOCATE(IREFSC)
*---- 
*  TEST TO COMPUTE STEADY-STATE OR TRANSIENT CALCULATION
*----
      IF(ITIME.EQ.0) THEN
         GO TO 230
      ELSE IF(ITIME.EQ.1) THEN
         GO TO 310
      ELSE
         CALL XABORT('@THM: UNEXPECTED VALUE FOR ITIME.')
      ENDIF
*----
*  CALL DRIVER FOR STEADY-STATE CALCULATION
*----
*     memory allocation for the steady-state calculation
  230 ALLOCATE(DCOOL(NZ,NX,NY),TCOOL(NZ,NX,NY),TSURF(NZ,NX,NY),
     > PCOOL(NZ,NX,NY),HCOOL(NZ,NX,NY),RAD((NDTOT-1),NZ,NX,NY))
      DCOOL(:NZ,:NX,:NY)=0.0
      TCOOL(:NZ,:NX,:NY)=0.0
      TSURF(:NZ,:NX,:NY)=0.0
      PCOOL(:NZ,:NX,:NY)=0.0
      HCOOL(:NZ,:NX,:NY)=0.0
      RAD(:NDTOT-1,:NZ,:NX,:NY)=0.0
*----
*  COMPUTE FUEL RADII
*----
      ALLOCATE(RVAL((NDTOT-1),NZ))
      IF(JENTRY(1).EQ.0) THEN
        WRITE(6,*)'RC,RIG=',RC,RIG
*CGT THERE IS GAP
        IF(RC.NE.RIG) THEN
          IGAP=0
          ARF=0.5*RC**2  ! at fuel radius
          ARCI=0.5*RIG**2 ! at internal clad radius
          ARCE=0.5*RGG**2 ! at external clad radius
          DARF=ARF/REAL(NFD-1)
          DARC=(ARCE-ARCI)/REAL(NDTOT-NFD-2)
          DO IEL=1,NZ
            RVAL(1,IEL)=0.0
            DO I=1,NFD-1
              RVAL(I+1,IEL)=REAL(SQRT(2.0D0*REAL(I)*DARF))
            ENDDO
            DO I=NFD+1,NDTOT-1
              RVAL(I,IEL)=REAL(SQRT(2.0D0*(ARCI+REAL(I-NFD-1)*DARC)))
            ENDDO
          ENDDO
        ELSE
*CGT NO GAP
          IGAP=1
          ARF=0.5*RC**2  ! at fuel radius
          ARCE=0.5*RGG**2 ! at external clad radius
          DARF=ARF/REAL(NFD-1)
          DARC=(ARCE-ARF)/REAL(NDTOT-NFD-1)
          DO IEL=1,NZ
            RVAL(1,IEL)=0.0
            DO I=1,NFD
              RVAL(I+1,IEL)=REAL(SQRT(2.0D0*REAL(I)*DARF))
            ENDDO
            DO I=NFD+1,NDTOT-1
              RVAL(I,IEL)=REAL(SQRT(2.0D0*(ARF+REAL(I-NFD)*DARC)))
            ENDDO
          ENDDO
        ENDIF
        CALL LCMPUT(IPTHM,'REF-RAD',(NDTOT-1)*NZ,2,RVAL)
        JPTHM=LCMDID(IPTHM,'HISTORY-DATA')
        KPTHM=LCMDID(JPTHM,'TIMESTEP0000')
        LPTHM=LCMLID(KPTHM,'CHANNEL',NCH)
      ELSE
        JPTHM=LCMGID(IPTHM,'HISTORY-DATA')
        KPTHM=LCMGID(JPTHM,'TIMESTEP0000')
        LPTHM=LCMGID(KPTHM,'CHANNEL')
      ENDIF
*----
*  LOOP OVER REACTOR CHANNELS
*----
      ICH=0
      SUMSEC=0.0
      DO 265 IY=1,NY
      DO 260 IX=1,NX
        IEL=IX+(IY-1)*NX
        DO 240 IZ=1,NZ
          IF(NUM((IZ-1)*NX*NY+IEL).NE.0) GO TO 250
  240   CONTINUE
        GO TO 260
  250   ICH=ICH+1
*----
*  COMPUTE HYDRAULICS CONSTANTS
*  SASS:    assembly cross section in m^2
*  RC:      fuel pellet radius in m
*  RTG:     guide tube radius in m
*  ACOOL:   coolant cross section per assembly in m^2
*  RAPCOOL: assembly over coolant volumic ratio
*  RAPFUEL: assembly over fuel volumic ratio
*  FCOOL:   power density fraction in coolant.
*  FFUEL:   power density fraction in fuel.
*  PCH:     heating perimeter in m
*  PM:      perimeter in contact with flow in m
*  HD:      hydraulic diameter of one assembly in m
*  SPEED:   inlet flow velocity in m/s
*----
        IF(FNFU(ICH).EQ.0.0) GO TO 260
        SASS=HX(IX)*HY(IY)
        IF(PITCH.EQ.0.0) THEN
*         PWR ASSEMBLY
          ACOOL=SASS-FNFU(ICH)*PI*RGG*RGG-FNTG(ICH)*PI*RTG*RTG
          RAPCOOL=SASS/ACOOL
          PCH=FNFU(ICH)*2.0*PI*RGG
          PM=PCH+FNTG(ICH)*2.0*PI*RTG
          SUMSEC=SUMSEC+ACOOL
        ELSE
*         CANDU CLUSTER
          ATOTHEX=3.0*PITCH**2.0*(3.0)**0.5/2.0
          ATIGEHEX=3.0*PI*RGG*RGG
          ACOOL=ATOTHEX-ATIGEHEX
          PM=6.0*PI*RGG
          PCH=PM
          RAPCOOL=3.0*SASS/(FNFU(ICH)*ACOOL)
          SUMSEC=SUMSEC+FNFU(ICH)*ACOOL/3.0
        ENDIF
        RAPFUEL=SASS/(FNFU(ICH)*PI*RC*RC)
        FCOOL=(1.0-FPUISS)*RAPCOOL
        FFUEL=FPUISS*RAPFUEL
        HD=4.0*ACOOL/PM
        IF(HD.LE.0.) CALL XABORT('THM: NEGATIVE HYDRAULIC DIAMETER(1).')
*----
*  RECOVER STEADY-STATE RADII
*----
        IF(JENTRY(1).EQ.0) THEN
          RAD(:,:,IX,IY)=RVAL(:,:)
        ELSE IF(JENTRY(1).EQ.1) THEN
          MPTHM=LCMGIL(LPTHM,ICH)
          CALL LCMGET(MPTHM,'RADII',RAD(1,1,IX,IY))
        ENDIF
*----
*  EXECUTION OF THE STEADY-STATE DRIVER PROGRAM
*----
        MPTHM=LCMDIL(LPTHM,ICH)
        CALL THMDRV(MPTHM,IMPX,IX,IY,NZ,XBURN(1,IX,IY),SASS,HZ,CFLUX,
     >  POROS,FNFU(ICH),NFD,NDTOT,IFLUID,SNAME,SCOMP,IGAP,IFUEL,FNAME,
     >  FCOMP,FCOOL,FFUEL,ACOOL,
     >  HD,PCH,RAD(1,1,IX,IY),MAXIT1,MAXIT2,ERMAXT,SPEED,TINLET,POULET,
     >  FRACPU(ICH),ICONDF,NCONDF,KCONDF,UCONDF,ICONDC,NCONDC,KCONDC,
     >  UCONDC,IHGAP,KHGAP,IHCONV,KHCONV,WTEFF,IFRCDI,ISUBM,FRO,
     >  POW(1,IX,IY),IPRES,IDFM,TCOMB(1,IX,IY),DCOOL(1,IX,IY),
     >  TCOOL(1,IX,IY),TSURF(1,IX,IY),HCOOL(1,IX,IY),PCOOL(1,IX,IY))
  260 CONTINUE
  265 CONTINUE
      IF(IMPX.GT.1) WRITE(6,610) SUMSEC,CWSECT
      DEALLOCATE(RVAL)
      CALL LCMPUT(KPTHM,'TIME',1,2,RTIME)
      IF(IMPX.GT.1) WRITE(6,470) 'TIMESTEP0000',RTIME
*----
*  PRINT AVERAGED THERMALHYDRAULICS PROPERTIES OVER THE CORE MAP
*----
      IF(IMPX.GT.1) THEN
        CALL THMAVG(IPMAP,IMPX,NX,NY,NZ,NCH,TCOMB,TSURF,DCOOL,TCOOL,
     >  PCOOL,HCOOL,POW,NSIMS)
      ENDIF
      DEALLOCATE(RAD,HCOOL)
      GO TO 400
*----
* CALL DRIVER FOR TRANSIENT CALCULATION
*----
*     memory allocation for the transient calculation
  310 ALLOCATE(TSURF(NZ,NX,NY),TCOOL(NZ,NX,NY),DCOOL(NZ,NX,NY),
     > PCOOL(NZ,NX,NY))
      TSURF(:NZ,:NX,:NY)=0.0
      TCOOL(:NZ,:NX,:NY)=0.0
      DCOOL(:NZ,:NX,:NY)=0.0
      PCOOL(:NZ,:NX,:NY)=0.0
*----
*  RECOVER TIME INDEX AT INITIAL CONDITIONS
*----
      JPTHM=LCMDID(IPTHM,'HISTORY-DATA')
      KPTHMI=LCMGID(JPTHM,'TIMESTEP0000')
      CALL LCMGET(KPTHMI,'TIME',TIMEPR)
      IF(ABS(RTIME-TIMEPR).LE.1.0E-3*DTIME) THEN
        TIMEIT=1
      ELSE
        DO I=1,TIMEIT
          WRITE(TXTDIR,'(8HTIMESTEP,I4.4)') I
          KPTHMI=LCMGID(JPTHM,TXTDIR)
          CALL LCMGET(KPTHMI,'TIME',TIMEPR)
          IF(ABS(RTIME-TIMEPR).LE.1.0E-3*DTIME) THEN
            TIMEIT=I+1
            GO TO 315
          ENDIF
        ENDDO
        WRITE(HSMG,'(45H@THM: UNABLE TO FIND INITIAL CONDITIONS AT T=,
     >  1P,E14.4,3H S.)') RTIME
        CALL XABORT(HSMG)
      ENDIF
  315 LPTHMI=LCMGID(KPTHMI,'CHANNEL')
      WRITE(TXTDIR,'(8HTIMESTEP,I4.4)') TIMEIT
      KPTHM=LCMDID(JPTHM,TXTDIR)
      LPTHM=LCMLID(KPTHM,'CHANNEL',NCH)
      IF(IMPX.GT.1) WRITE(6,530) TIMEIT,RTIME,RTIME+DTIME
*----
*  LOOP OVER REACTOR CHANNELS
*----
      ICH=0
      DO 355 IY=1,NY
      DO 350 IX=1,NX
        IEL=IX+(IY-1)*NX
        DO 320 IZ=1,NZ
          IF(NUM((IZ-1)*NX*NY+IEL).NE.0) GO TO 330
  320   CONTINUE
        GO TO 350
  330   ICH=ICH+1
*----
*  COMPUTE HYDRAULICS CONSTANTS
*----
        IF(FNFU(ICH).EQ.0.0) GO TO 350
        SASS=HX(IX)*HY(IY)
        IF(PITCH.EQ.0.0) THEN
*         PWR ASSEMBLY
          ACOOL=SASS-FNFU(ICH)*PI*RGG*RGG-FNTG(ICH)*PI*RTG*RTG
          RAPCOOL=SASS/ACOOL
          PCH=FNFU(ICH)*2.0*PI*RGG
          PM=PCH+FNTG(ICH)*2.0*PI*RTG
        ELSE
*         CANDU CLUSTER
          ATOTHEX=3.0*PITCH**2.0*(3.0)**0.5/2.0
          ATIGEHEX=3.0*PI*RGG*RGG
          ACOOL=ATOTHEX-ATIGEHEX
          PM=6.0*PI*RGG
          PCH=PM
          RAPCOOL=3.0*SASS/(FNFU(ICH)*ACOOL)
        ENDIF
        RAPFUEL=SASS/(FNFU(ICH)*PI*RC*RC)
        FCOOL=(1.0-FPUISS)*RAPCOOL
        FFUEL=FPUISS*RAPFUEL
        HD=4.0*ACOOL/PM
        IF(HD.LE.0.) CALL XABORT('THM: NEGATIVE HYDRAULIC DIAMETER(2).')
*----
*  EXECUTION OF THE TRANSIENT DRIVER PROGRAM
*----
        MPTHMI=LCMGIL(LPTHMI,ICH)
        MPTHM=LCMDIL(LPTHM,ICH)
        CALL THMTRS(MPTHMI,MPTHM,IMPX,IX,IY,NZ,XBURN(1,IX,IY),SASS,HZ,
     >  DTIME,CFLUX,POROS,FNFU(ICH),NFD,NDTOT,IFLUID,SNAME,SCOMP,
     >  IGAP,IFUEL,FNAME,FCOMP,   
     >  FCOOL,FFUEL,ACOOL,HD,PCH,MAXIT3,MAXIT1,MAXIT2,ERMAXT,ERMAXC,
     >  SPEED,TINLET,POULET,FRACPU(ICH),ICONDF,NCONDF,KCONDF,UCONDF,
     >  ICONDC,NCONDC,KCONDC,UCONDC,IHGAP,KHGAP,IHCONV,KHCONV,WTEFF,
     >  IFRCDI,ISUBM,FRO,POW(1,IX,IY),TCOMB(1,IX,IY),DCOOL(1,IX,IY),
     >  TCOOL(1,IX,IY),TSURF(1,IX,IY))
  350 CONTINUE
  355 CONTINUE
      CALL LCMPUT(KPTHM,'TIME',1,2,RTIME+DTIME)
      IF(IMPX.GT.1) WRITE(6,470) TXTDIR,RTIME+DTIME
*----
*  RECOVER LOCAL PARAMETER INFORMATION COMPUTED BY THMDRV OR THMTRS
*----
  400 ERRA1=0.0
      ERRA2=0.0
      ERRA3=0.0
      ERRBB=0.0
      ZMINA1=1.0E10
      ZMINA2=1.0E10
      ZMINA3=1.0E10
      ZMINBB=1.0E10
      ZMAXA1=0.0
      ZMAXA2=0.0
      ZMAXA3=0.0
      ZMAXBB=0.0
      RATIOX=0.0
      ALLOCATE(IREFSC(NCH))
      CALL LCMLEN(IPMAP,'REF-SCHEME',ILONG,ITYLCM)
      IF(ILONG.EQ.NCH) THEN
        CALL LCMGET(IPMAP,'REF-SCHEME',IREFSC)
      ELSE
        IREFSC(:NCH)=1
      ENDIF
      JPMAP=LCMGID(IPMAP,'PARAM')
      DO 460 IPAR=1,NPARM
      KPMAP=LCMGIL(JPMAP,IPAR)
      CALL LCMGTC(KPMAP,'P-NAME',12,PNAME)
      IF((PNAME.EQ.'T-FUEL').OR.(PNAME.EQ.'D-COOL').OR.
     1   (PNAME.EQ.'T-COOL').OR.(PNAME.EQ.'T-SURF').OR.
     2   (PNAME.EQ.'P-COOL')) THEN
        CALL LCMGET(KPMAP,'P-TYPE',ITYPE)
        ALLOCATE(VAL(NCH,NB))
        RELAX0=1.0
        IF(ITYPE.EQ.1) THEN
          IF(IMPX.GT.0) WRITE(6,510) 'GLOBAL',PNAME
          CALL LCMGET(KPMAP,'P-VALUE',FLOT)
          DO 415 ICH=1,NCH
          DO 410 IB=1,NB
          VAL(ICH,IB)=FLOT
  410     CONTINUE
  415     CONTINUE
        ELSE IF(ITYPE.EQ.2) THEN
          RELAX0=RELAX
          IF(IMPX.GT.0) WRITE(6,510) 'LOCAL',PNAME
          CALL LCMGET(KPMAP,'P-VALUE',VAL)
        ENDIF
        ICH=0
        DO 455 IY=1,NY
        DO 450 IX=1,NX
        IEL=(IY-1)*NX+IX
        DO 420 IZ=1,NZ
        IF(NUM((IZ-1)*NX*NY+IEL).NE.0) GO TO 430
  420   CONTINUE
        GO TO 450
  430   ICH=ICH+1
        IB=0
        DO 440 IZ=1,NZ
        IF(NUM((IZ-1)*NX*NY+IEL).EQ.0) GO TO 440
        IB=IB+1
        FLOT=0.0
        IF(PNAME.EQ.'T-FUEL') THEN
          IF(IREFSC(ICH).GT.0) THEN
            FLOT=RELAX0*TCOMB(IZ,IX,IY)+(1.0-RELAX0)*VAL(ICH,IB)
          ELSE
            FLOT=RELAX0*TCOMB(NZ-IZ+1,IX,IY)+(1.0-RELAX0)*VAL(ICH,IB)
          ENDIF
          IF(ITIME.EQ.0) ERRA1=MAX(ERRA1,ABS(VAL(ICH,IB)-FLOT))
          ZMINA1=MIN(ZMINA1,FLOT)
          ZMAXA1=MAX(ZMAXA1,FLOT)
        ELSE IF(PNAME.EQ.'D-COOL') THEN
          IF(IREFSC(ICH).GT.0) THEN
            FLOT=RELAX0*DCOOL(IZ,IX,IY)/ZKILO+(1.0-RELAX0)*VAL(ICH,IB)
          ELSE
            FLOT=RELAX0*DCOOL(NZ-IZ+1,IX,IY)/ZKILO+(1.0-RELAX0)
     >      *VAL(ICH,IB)
          ENDIF
          IF(ITIME.EQ.0) ERRA2=MAX(ERRA2,ABS(VAL(ICH,IB)-FLOT))
          ZMINA2=MIN(ZMINA2,FLOT)
          ZMAXA2=MAX(ZMAXA2,FLOT)
        ELSE IF(PNAME.EQ.'T-COOL') THEN
          IF(IREFSC(ICH).GT.0) THEN
            FLOT=RELAX0*TCOOL(IZ,IX,IY)+(1.0-RELAX0)*VAL(ICH,IB)
          ELSE
            FLOT=RELAX0*TCOOL(NZ-IZ+1,IX,IY)+(1.0-RELAX0)*VAL(ICH,IB)
          ENDIF
          IF(ITIME.EQ.0) ERRA3=MAX(ERRA3,ABS(VAL(ICH,IB)-FLOT))
          ZMINA3=MIN(ZMINA3,FLOT)
          ZMAXA3=MAX(ZMAXA3,FLOT)
        ELSE IF(PNAME.EQ.'T-SURF') THEN
          IF(IREFSC(ICH).GT.0) THEN
            FLOT=RELAX0*TSURF(IZ,IX,IY)+(1.0-RELAX0)*VAL(ICH,IB)
          ELSE
            FLOT=RELAX0*TSURF(NZ-IZ+1,IX,IY)+(1.0-RELAX0)*VAL(ICH,IB)
          ENDIF
          IF(ITIME.EQ.0) ERRBB=MAX(ERRBB,ABS(VAL(ICH,IB)-FLOT))
          ZMINBB=MIN(ZMINBB,FLOT)
          ZMAXBB=MAX(ZMAXBB,FLOT)
        ELSE IF(PNAME.EQ.'P-COOL') THEN
          IF(IREFSC(ICH).GT.0) THEN
            FLOT=RELAX0*PCOOL(IZ,IX,IY)+(1.0-RELAX0)*VAL(ICH,IB)
          ELSE
            FLOT=RELAX0*PCOOL(NZ-IZ+1,IX,IY)+(1.0-RELAX0)*VAL(ICH,IB)
          ENDIF
          IF(ITIME.EQ.0) ERRA4=MAX(ERRBB,ABS(VAL(ICH,IB)-FLOT))
          ZMINA4=MIN(ZMINBB,FLOT)
          ZMAXA4=MAX(ZMAXBB,FLOT)
        ELSE
          CALL XABORT('@THM: INVALID PARAMETER TYPE: '// PNAME//'.')
        ENDIF
        VAL(ICH,IB)=FLOT
  440   CONTINUE
  450   CONTINUE
  455   CONTINUE
        ITYPE=2
        CALL LCMPUT(KPMAP,'P-TYPE',1,1,ITYPE)
        CALL LCMPUT(KPMAP,'P-VALUE',NCH*NB,2,VAL)
        CALL LCMLEN(IPMAP,'AXIAL-FPW',JLONG,ITYLCM)
        DD1=0.0
        DD2=0.0
        IF(JLONG.NE.0) THEN
          ALLOCATE(FPOWER(NB))
          IF(JLONG.NE.NB) CALL XABORT('THM: UNABLE TO FIND RECORD AXIA'
     1    //'L-FPW IN THE FUELMAP.')
          CALL LCMGET(IPMAP,'AXIAL-FPW',FPOWER)
          DO ICH=1,NCH
            DO IB=1,NB
              DD1=DD1+VAL(ICH,IB)*FPOWER(IB)**2
              DD2=DD2+FPOWER(IB)**2
            ENDDO
          ENDDO
          DEALLOCATE(FPOWER)
        ELSE
          ALLOCATE(PW(NCH*NB))
          CALL LCMGET(IPMAP,'BUND-PW',PW)
          ITOT=0
          DO IB=1,NB
            DO ICH=1,NCH
              ITOT=ITOT+1
              DD1=DD1+VAL(ICH,IB)*PW(ITOT)**2
              DD2=DD2+PW(ITOT)**2
            ENDDO
          ENDDO
          DEALLOCATE(PW)
        ENDIF
        TMOY0=DD1/DD2
        TEXT12='AVG-'//PNAME(:8)
        CALL LCMLEN(IPTHM,TEXT12,KLONG,ITYLCM)
        IF(((PNAME.EQ.'T-FUEL').OR.(PNAME.EQ.'T-COOL').OR.
     1  (PNAME.EQ.'P-COOL')).AND.(KLONG.GT.0)) THEN
          CALL LCMGET(IPTHM,TEXT12,TMOY0I)
          IF(PNAME.EQ.'T-FUEL') THEN
            RATIO=ABS(TMOY0/DTEMPR-TMOY0I/DTEMPR)
            IF(IMPX.GT.0) WRITE(6,490) TEXT12,TMOY0I,TMOY0,RATIO
            RATIOX=MAX(RATIOX,RATIO)
          ELSE IF(PNAME.EQ.'T-COOL') THEN
            RATIO=ABS(TMOY0/DTEMPT-TMOY0I/DTEMPT)
            IF(IMPX.GT.0) WRITE(6,490) TEXT12,TMOY0I,TMOY0,RATIO
            RATIOX=MAX(RATIOX,RATIO)
          ELSE IF(PNAME.EQ.'P-COOL') THEN
            RATIO=ABS(TMOY0/DPRESS-TMOY0I/DPRESS)
            IF(IMPX.GT.0) WRITE(6,490) TEXT12,TMOY0I,TMOY0,RATIO
            RATIOX=MAX(RATIOX,RATIO)
          ENDIF
        ENDIF
        CALL LCMPUT(IPTHM,TEXT12,1,2,TMOY0)
        DEALLOCATE(VAL)
      ENDIF
      IF(PNAME.EQ.'T-FUEL') THEN
        IF(ITIME.EQ.0) CALL LCMPUT(IPTHM,'ERROR-T-FUEL',1,2,ERRA1)
        CALL LCMPUT(IPTHM,'MIN-T-FUEL',1,2,ZMINA1)
        CALL LCMPUT(IPTHM,'MAX-T-FUEL',1,2,ZMAXA1)
        IF(IMPX.GT.0) WRITE(6,520) 'FUEL TEMPERATURE',ERRA1,'K',
     1  ZMINA1,'K',ZMAXA1,'K'
      ELSE IF(PNAME.EQ.'D-COOL') THEN
        IF(ITIME.EQ.0) CALL LCMPUT(IPTHM,'ERROR-D-COOL',1,2,ERRA2)
        CALL LCMPUT(IPTHM,'MIN-D-COOL',1,2,ZMINA2)
        CALL LCMPUT(IPTHM,'MAX-D-COOL',1,2,ZMAXA2)
        IF(IMPX.GT.0) WRITE(6,520) 'COOLANT DENSITY',ERRA2,'g/cc',
     1  ZMINA2,'g/cc',ZMAXA2,'g/cc'
      ELSE IF(PNAME.EQ.'T-COOL') THEN
        IF(ITIME.EQ.0) CALL LCMPUT(IPTHM,'ERROR-T-COOL',1,2,ERRA3)
        CALL LCMPUT(IPTHM,'MIN-T-COOL',1,2,ZMINA3)
        CALL LCMPUT(IPTHM,'MAX-T-COOL',1,2,ZMAXA3)
        IF(IMPX.GT.0) WRITE(6,520) 'COOLANT TEMPERATURE',ERRA3,'K',
     1  ZMINA3,'K',ZMAXA3,'K'
      ELSE IF(PNAME.EQ.'T-SURF') THEN
        IF(ITIME.EQ.0) CALL LCMPUT(IPTHM,'ERROR-T-SURF',1,2,ERRBB)
        IF(IMPX.GT.0) WRITE(6,520) 'FUEL SURFACE TEMPERATURE',ERRBB,
     1  'K',ZMINBB,'K',ZMAXBB,'K'
      ELSE IF(PNAME.EQ.'P-COOL') THEN
        IF(ITIME.EQ.0) CALL LCMPUT(IPTHM,'ERROR-P-COOL',1,2,ERRA4)
        CALL LCMPUT(IPTHM,'MIN-P-COOL',1,2,ZMINA4)
        CALL LCMPUT(IPTHM,'MAX-P-COOL',1,2,ZMAXA4)
        IF(IMPX.GT.0) WRITE(6,520) 'COOLANT PRESSURE',ERRA4,'Pa',
     1  ZMINA4,'Pa',ZMAXA4,'Pa'
      ENDIF
  460 CONTINUE
      DEALLOCATE(IREFSC)
*----
*  SAVE CONDUCTIVITY INFORMATION ON LCM OBJECT THM
*----
      IF(ICONDF.EQ.1) THEN
        CALL LCMPUT(IPTHM,'KCONDF',NCONDF+3,2,KCONDF)
        CALL LCMPTC(IPTHM,'UCONDF',12,UCONDF)
      ENDIF
      IF(ICONDC.EQ.1) THEN
        CALL LCMPUT(IPTHM,'KCONDC',NCONDC+1,2,KCONDC)
        CALL LCMPTC(IPTHM,'UCONDC',12,UCONDC)
      ENDIF
*----
*  RELEASE MEMORY
*----
      DEALLOCATE(PCOOL,DCOOL,TCOOL,TSURF,TCOMB)
      DEALLOCATE(NUM)
      DEALLOCATE(POW,XBURN)
      DEALLOCATE(FRO)
      IF(ICONDF.EQ.1)DEALLOCATE(KCONDF)
      IF(ICONDC.EQ.1)DEALLOCATE(KCONDC)
*----
*  STATE-VECTOR FOR THM
*----
      HSIGN='L_THM'
      CALL LCMPTC(IPTHM,'SIGNATURE',12,HSIGN)
      ISTATE(:NSTATE)=0
      ISTATE(1)=NCH
      ISTATE(2)=NZ
      ISTATE(3)=MAXIT1
      ISTATE(4)=MAXIT2
      ISTATE(5)=MAXIT3
      ISTATE(6)=NFD
      ISTATE(7)=NDTOT
      ISTATE(8)=ITIME
      ISTATE(9)=TIMEIT
      ISTATE(10)=IHGAP
      ISTATE(11)=IHCONV
      ISTATE(12)=ICONDF
      ISTATE(13)=ICONDC
      ISTATE(14)=IFRCDI
      ISTATE(15)=ISUBM
      IF(ICONDF.EQ.1) ISTATE(16)=NCONDF
      IF(ICONDC.EQ.1) ISTATE(17)=NCONDC
      ISTATE(18)=NPRAD
      ISTATE(19)=NPOWER
      ISTATE(20)=IFLUID
      ISTATE(21)=IGAP
      ISTATE(22)=IPRES
      ISTATE(23)=IDFM
      CALL LCMPUT(IPTHM,'STATE-VECTOR',NSTATE,1,ISTATE)
      STATE(:NSTATE)=0.0
      STATE(1)=DTIME
      STATE(2)=FPUISS
      STATE(3)=CFLUX
      STATE(4)=SPEED
      STATE(5)=POULET
      STATE(6)=TINLET
      STATE(7)=POROS
      STATE(8)=RC
      STATE(9)=RIG
      STATE(10)=RGG
      STATE(11)=RTG
      STATE(12)=PITCH
      STATE(13)=ERMAXT
      STATE(14)=ERMAXC
      STATE(15)=RELAX
      STATE(16)=RTIME
      IF(IHGAP.EQ.1) STATE(17)=KHGAP
      IF(IHCONV.EQ.1) STATE(18)=KHCONV
      STATE(19)=WTEFF
      STATE(20)=TPOW
      STATE(21)=RATIOX
      STATE(22)=EPSR
      STATE(23)=THETA
      CALL LCMPUT(IPTHM,'REAL-PARAM',NSTATE,2,STATE)
      IF(IMPX.GT.0) THEN
         WRITE(6,540) ISTATE(:15),ISTATE(18:22)
         IF(ISTATE(10).EQ.1) WRITE(6,550) (ISTATE(16))
         IF(ISTATE(11).EQ.1) WRITE(6,560) (ISTATE(17))
         WRITE(6,570) STATE(:16),STATE(19),STATE(21:23)
         IF(ISTATE(10).EQ.1) WRITE(6,580) (STATE(17))
         IF(ISTATE(11).EQ.1) WRITE(6,590) (STATE(18))
         IF(ISTATE(19).GT.0) WRITE(6,600) (STATE(20))
      ENDIF
      IF(IMPX.GT.4) CALL LCMLIB(IPTHM)
*----
*  SAVE CELL-DEPENDENT DATA
*----
      CALL LCMPUT(IPTHM,'NB-FUEL',NCH,2,FNFU)
      CALL LCMPUT(IPTHM,'NB-TUBE',NCH,2,FNTG)
      CALL LCMPUT(IPTHM,'FRACT-PU',NCH,2,FRACPU)
      DEALLOCATE(FRACPU,FNTG,FNFU)
*----
*  RECOVER THE VARIATION RATIO AND SAVE IT IN A CLE-2000 VARIABLE
*----
      IF(IPICK.EQ.1) THEN
         CALL REDGET(ITYP,NITMA,FLOT,TEXT12,DFLOT)
         IF(ITYP.NE.-2) CALL XABORT('THM: OUTPUT REAL EXPECTED.')
         ITYP=2
         CALL REDPUT(ITYP,NITMA,RATIOX,TEXT12,DFLOT)
         CALL REDGET(ITYP,NITMA,FLOT,TEXT12,DFLOT)
         IF((ITYP.NE.3).OR.(TEXT12.NE.';')) THEN
           CALL XABORT('THM: ; CHARACTER EXPECTED.')
         ENDIF      
      ENDIF      
      RETURN
*
  470 FORMAT(/11H THM: SAVE ,A,9H AT TIME=,1P,E12.4,3H S.)
  480 FORMAT(/31H THM: RADIAL POWER FORM FACTORS/(1P,10E12.4))
  490 FORMAT(/18H THM: PARAMETER = ,A,1P,E12.4,3H ->,E12.4,7H RATIO=,
     1 E12.4)
  500 FORMAT(/27H THM: SET GLOBAL PARAMETER ,A,2H =,1P E12.4)
  510 FORMAT(/14H THM: RECOVER ,A,13H PARAMETER = ,A,1H.)
  520 FORMAT(/15H THM: ERROR ON ,A,2H =,F12.3,1X,A,13H  MIN VALUE =,
     1 F12.3,1X,A,13H  MAX VALUE =,F12.3,1X,A)
  530 FORMAT(/28H THM: PERFORM TRANSIENT STEP,I5,9H BETWEEN ,1P,E14.4,
     1 4H AND,E14.4,3H S.)
  540 FORMAT(/
     1 14H STATE VECTOR:/
     2 7H NZ    ,I9,27H   (NUMBER OF AXIAL MESHES)/
     3 7H NCH   ,I9,43H   (NUMBER OF CHANNELS IN THE RADIAL PLANE)/
     4 7H MAXIT1,I9,36H   (NUMBER OF CONDUCTION ITERATIONS)/
     5 7H MAXIT2,I9,39H   (NUMBER OF CENTER-PELLET ITERATIONS)/
     6 7H MAXIT3,I9,30H   (NUMBER OF FLOW ITERATIONS)/
     7 7H NFD   ,I9,32H   (NUMBER OF FUEL RADIAL ZONES)/
     8 7H NDTOT ,I9,36H   (NUMBER OF DISCRETISATION POINTS)/
     9 7H ITIME ,I9,21H   (CALCULATION TYPE)/
     1 7H TIMEIT,I9,30H   (TRANSIENT ITERATION INDEX)/
     2 7H IHGAP ,I9,34H   (HGAP FLAG (0=DEFAULT/1=FIXED))/
     3 7H IHCONV,I9,42H   (HCONV FLAG (0=DITTUS-BOELTER/1=FIXED))/
     4 7H ICONDF,I9,46H   (FUEL CONDUCTIVITY FLAG (0=STORA-CHENEBAULT,
     5 54H (UOX), COMETHE (MOX)/1=USER-PROVIDED FUNCTION OF FUEL,
     6 14H TEMPERATURE))/
     7 7H ICONDC,I9,39H   (CLAD CONDUCTIVITY FLAG (0=DEFAULT/1,
     8 47H=USER-PROVIDED POLYNOMIAL OF CLAD TEMPERATURE))/
     9 7H IFRCDI,I9,40H   (FUEL CONDUCTIVITY APPROXIMATION FLAG,
     1 44H (0=DEFAULT/1=AVERAGE APPROXIMATION FORCED))/
     2 7H ISUBM ,I9,47H   (BOILING MODEL FLAG (0=ONE-PHASE/1=BOWRING C,
     3 37HORRELATION/2=SAHA-ZUBER CORRELATION))/
     4 7H NPRAD ,I9,47H   (RADIAL POWER FORM FACTOR (0=FLAT/NUMBER OF ,
     5 8HPOINTS))/
     6 7H NPOWER,I9,36H   (NUMBER OF POINTS IN POWER-TABLE)/
     7 7H IFLUID,I9,32H   (TYPE OF FLUID (0=H2O/1=D2O))/
     8 7H IGAP  ,I9,40H   (GAP IS CONSIDERED (0=GAP/1=NO GAP))/
     9 7H IPRES ,I9,46H   (PRESSURE DROP (0=CONSTANT/1=NON CONSTANT))/
     1 7H IDFM ,I9,35H   (DRIFT FLUX MODEL (0=HEM1/1=EPRI
     2 33H/2=MODEBSTION/3=GERAMP/4=CHEXAL)))
  550 FORMAT(
     1 7H NCONDF,I9,43H   (DEGREE OF FUEL CONDUCTIVITY POLYNOMIAL))
  560 FORMAT(
     1 7H NCONDC,I9,43H   (DEGREE OF CLAD CONDUCTIVITY POLYNOMIAL))
  570 FORMAT(/
     1 12H REAL PARAM:,1P/
     2 7H DTIME ,E12.4,19H   (TIME STEP IN S)/
     3 7H FPUISS,E12.4,25H   (COOLANT POWER FACTOR)/
     4 7H CFLUX ,E12.4,32H   (CRITICAL HEAT FLUX IN W/M^2)/
     5 7H SPEED ,E12.4,28H   (COOLANT VELOCITY IN M/S)/
     6 7H POULET,E12.4,34H   (OUTLET COOLANT PRESSURE IN PA)/
     7 7H TINLET,E12.4,35H   (INLET COOLANT TEMPERATURE IN K)/
     8 7H POROS ,E12.4,19H   (OXYDE POROSITY)/
     9 7H RC    ,E12.4,28H   (FUEL PELLET RADIUS IN M)/
     1 7H RIG   ,E12.4,34H   (INTERNAL CLAD ROD RADIUS IN M)/
     2 7H RGG   ,E12.4,34H   (EXTERNAL CLAD ROD RADIUS IN M)/
     3 7H RTG   ,E12.4,27H   (GUIDE TUBE RADIUS IN M)/
     4 7H PITCH ,E12.4,24H   (HEXAGONAL SIDE IN M)/
     5 7H ERMAXT,E12.4,35H   (TEMPERATURE MAXIMUM ERROR IN K)/
     6 7H ERMAXC,E12.4,32H   (FLOW MAXIMUM RELATIVE ERROR)/
     7 7H RELAX ,E12.4,25H   (RELAXATION PARAMETER)/
     8 7H RTIME ,E12.4,20H   (TIME VALUE IN S)/
     9 7H WTEFF ,E12.4,44H   (SURFACE TEMPERATURE WEIGHTING FACTOR IN ,
     1 27HEFFECTIVE FUEL TEMPERATURE)/
     2 7H RATIOX,E12.4,35H   (MAXIMUM OF VARIABLE VARIATIONS)/
     3 7H EPSR  ,E12.4,34H   (RUGOSITY IN M OF THE FUEL ROD)/
     4 7H THETA ,E12.4,41H   (ANGLE IN RADIANS OF THE FUEL CHANNEL))
  580 FORMAT(7H HGAP  ,1P,E12.4,20H   (HGAP IN W/m^2/K))
  590 FORMAT(7H HCONV ,1P,E12.4,21H   (HCONV IN W/m^2/K))
  600 FORMAT(7H HCONV ,1P,E12.4,22H   (POWER FACTOR IN W))
  610 FORMAT(/37H THM: CORE COOLANT SECTION. COMPUTED=,1P,E9.2,
     1 7H GIVEN=,E9.2,4H m2.)
      END
