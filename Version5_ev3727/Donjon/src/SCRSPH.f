*DECK SCRSPH
      SUBROUTINE SCRSPH(IPMEM,IPMAC,ICAL,IMPX,HEQUI,HMASL,NMIL,NGROUP,
     > ILUPS,SPH,B2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Extract a Macrolib corresponding to an elementary calculation in a
* memory-resident Saphyb.
*
*Copyright:
* Copyright (C) 2011 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IPMEM   pointer to the memory-resident Saphyb object.
* ICAL    index of the elementary calculation being considered.
* IMPX    print parameter (equal to zero for no print).
* HEQUI   keyword of SPH-factor set to be recovered.
* HMASL   keyword of MASL data set to be recovered.
* NMIL    number of mixtures in the elementary calculation.
* NGROUP  number of energy groups in the elementary calculation.
* ILUPS   up-scattering removing flag (=1 to remove up-scattering from
*         output cross-sections).
* B2      imposed buckling.
*
*Parameters: output
* IPMAC   pointer to the Macrolib (L_MACROLIB signature).
* SPH     SPH-factor set extracted from the Saphyb.
* B2      buckling recovered from the Saphyb.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMEM,IPMAC
      INTEGER ICAL,IMPX,NMIL,NGROUP,ILUPS
      REAL SPH(NMIL,NGROUP),B2
      CHARACTER HEQUI*4,HMASL*4
*----
*  LOCAL VARIABLES
*----
      INTEGER, PARAMETER::IOUT=6
      INTEGER, PARAMETER::MAXDIV=3
      INTEGER, PARAMETER::MAXLOC=10
      INTEGER, PARAMETER::MAXREA=25
      INTEGER, PARAMETER::MAXMAC=2
      INTEGER, PARAMETER::NSTATE=40
      REAL DEN,FF,CSCAT
      INTEGER I, J, IR, I0, IAD, IBM, IDF, IFISS, IGMAX, IGMIN, IGR, IL,
     & ILENG, ILOC, ILONG, IMAC, INDX, IPOSDE, IPRC, IREA, IRES, IS2,
     & ISO, ISOKEP, ITRANC, ITYLCM, JGR, NADRX, NCALS, NDATAP, NDATAX,
     & NED,  NISO, NISOTS, NL, NW, NLOC, NMAC, NPARL, NPR, NPRC, NREA,
     & NSURFD, NVDIV, IRENT0, IRENT1
      INTEGER ISTATE(NSTATE),DIMSAP(50)
      REAL VALDIV(MAXDIV)
      LOGICAL LSTRD,LDIFF,LSPH,LMASL
      CHARACTER TEXT12*12,HSMG*131,NOMREA(MAXREA)*12,CM*2,
     1 IDVAL(MAXDIV)*4,LOCTYP(MAXLOC)*4,LOCKEY(MAXLOC)*4,TEXT8*8,
     2 TEXT9*8
      TYPE(C_PTR) JPMAC,KPMAC,JPMEM,KPMEM,LPMEM,MPMEM
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: LOCAD,ISADRX,LENGDX,LENGDP,
     1 IDATA,IHEDI,TOTM,RESM,ISOTS,NOMISO,IPOS,NJJM,IJJM
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: IADRX
      REAL, ALLOCATABLE, DIMENSION(:) :: ENER,XVOLM,FLUXS,RDATA,STR,WRK,
     1 SCAT,GAR,RVALO,CONCES,LAMB,SURF,FMASL
      REAL, ALLOCATABLE, DIMENSION(:,:) :: NWT0,XSB,SIGS0,SIGSB,SURFLX,
     1 WORK,BETAR,INVELS
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: XS,SIGS,SS2DB,CHIRS
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: SS2D
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LXS
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: HADF
      CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:) :: NOMMIL
*----
*  SCRATCH STORAGE ALLOCATION
*   SIGS0    P0 scattering cross sections.
*----
      ALLOCATE(IPOS(NMIL),NJJM(NMIL),IJJM(NMIL),NOMMIL(NMIL))
      ALLOCATE(SIGS0(NMIL,NGROUP),FMASL(NMIL))
      FMASL(:NMIL)=0.0
*----
*  RECOVER SAPHYB CHARACTERISTICS
*----
      CALL LCMLEN(IPMEM,'DIMSAP',ILENG,ITYLCM)
      IF(ILENG.EQ.0) CALL XABORT('SCRSPH: INVALID SAPHYB.')
      CALL LCMGET(IPMEM,'DIMSAP',DIMSAP)
      IF(NMIL.NE.DIMSAP(7)) THEN
         CALL XABORT('SCRSPH: INVALID VALUE OF NMIL.')
      ELSE IF(NGROUP.NE.DIMSAP(20)) THEN
         CALL XABORT('SCRSPH: INVALID VALUE OF NGROUP.')
      ENDIF
      NREA=DIMSAP(4)   ! number of reactions
      NISO=DIMSAP(5)   ! number of particularized isotopes
      NMAC=DIMSAP(6)   ! number of macroscopic sets
      NPARL=DIMSAP(11) ! number of local variables
      NADRX=DIMSAP(18) ! number of address sets
      NCALS=DIMSAP(19) ! number of elementary calculations in the Saphyb
      NPRC=DIMSAP(31)  ! number of delayed neutron precursor groups
      NISOTS=DIMSAP(32) ! number of isotopes in edition tables
      IF(IMPX.GT.1) THEN
        WRITE(IOUT,'(30H SCRSPH: number of reactions =,I3)') NREA
        WRITE(IOUT,'(44H SCRSPH: number of particularized isotopes =,
     1  I4)') NISO
        WRITE(IOUT,'(37H SCRSPH: number of macroscopic sets =,I2)') NMAC
        WRITE(IOUT,'(29H SCRSPH: number of mixtures =,I5)') NMIL
        WRITE(IOUT,'(36H SCRSPH: number of local variables =,I4)') NPARL
        WRITE(IOUT,'(33H SCRSPH: number of address sets =,I4)') NADRX
        WRITE(IOUT,'(33H SCRSPH: number of calculations =,I5)') NCALS
        WRITE(IOUT,'(34H SCRSPH: number of energy groups =,I4)') NGROUP
        WRITE(IOUT,'(37H SCRSPH: number of precursor groups =,I4)') NPRC
        WRITE(IOUT,'(46H SCRSPH: number of isotopes in output tables =,
     1  I4)') NISOTS
      ENDIF
      IF(NREA.GT.MAXREA) CALL XABORT('SCRSPH: MAXREA OVERFLOW')
      IF(NMAC.GT.MAXMAC) CALL XABORT('SCRSPH: MAXMAC OVERFLOW')
      INDX=NISO+NMAC
      IF(INDX.EQ.0) CALL XABORT('SCRSPH: NO CROSS SECTIONS FOUND.')
*----
*  RECOVER INFORMATION FROM constphysiq DIRECTORY.
*----
      ALLOCATE(ENER(NGROUP+1))
      CALL LCMSIX(IPMEM,'constphysiq',1)
      CALL LCMGET(IPMEM,'ENRGS',ENER)
      CALL LCMSIX(IPMEM,' ',2)
      DO IGR=1,NGROUP+1
        ENER(IGR)=ENER(IGR)/1.0E-6
      ENDDO
      CALL LCMPUT(IPMAC,'ENERGY',NGROUP+1,2,ENER)
      DEALLOCATE(ENER)
*----
*  RECOVER INFORMATION FROM contenu DIRECTORY.
*----
      ALLOCATE(TOTM(NMIL),RESM(NMIL))
      CALL LCMSIX(IPMEM,'contenu',1)
      IF(NREA.GT.0) THEN
        CALL LCMGTC(IPMEM,'NOMREA',12,NREA,NOMREA)
        IF(IMPX.GT.1) THEN
          WRITE(IOUT,'(29H SCRSPH: Available reactions:/(1X,10A13))')
     1    (NOMREA(I),I=1,NREA)
        ENDIF
      ENDIF
      CALL LCMGET(IPMEM,'TOTMAC',TOTM)
      CALL LCMGET(IPMEM,'RESMAC',RESM)
      IF(NISO.GT.0) THEN
        ALLOCATE(NOMISO(NISO*2))
        CALL LCMGET(IPMEM,'NOMISO',NOMISO)
      ENDIF
      CALL LCMSIX(IPMEM,' ',2)
*----
*  RECOVER INFORMATION FROM adresses DIRECTORY.
*----
      NL=0
      IF(NADRX.GT.0) THEN
         ALLOCATE(IADRX((NREA+2),(NISO+NMAC),NADRX))
         CALL LCMSIX(IPMEM,'adresses',1)
         CALL LCMGET(IPMEM,'ADRX',IADRX)
         CALL LCMSIX(IPMEM,' ',2)
         DO IAD=1,NADRX
          DO ISO=1,NISO+NMAC
            NL=MAX(NL,IADRX(NREA+1,ISO,IAD))
            NL=MAX(NL,IADRX(NREA+2,ISO,IAD))
          ENDDO
         ENDDO
      ENDIF
      IF(IMPX.GT.1) THEN
        WRITE(IOUT,'(36H SCRSPH: number of legendre orders =,I4)') NL
      ENDIF
*----
*  RECOVER INFORMATION FROM geom DIRECTORY.
*----
      NSURFD=0
      CALL LCMSIX(IPMEM,'geom',1)
      ALLOCATE(XVOLM(NMIL))
      CALL LCMGET(IPMEM,'XVOLMT',XVOLM)
      CALL LCMGTC(IPMEM,'NOMMIL',20,NMIL,NOMMIL)
      CALL LCMLEN(IPMEM,'outgeom',ILONG,ITYLCM)
      IF(ILONG.NE.0) THEN
        CALL LCMSIX(IPMEM,'outgeom',1)
        CALL LCMLEN(IPMEM,'SURF',NSURFD,ITYLCM)
        IF(IMPX.GT.1) THEN
          WRITE(IOUT,'(42H SCRSPH: number of discontinuity factors =,
     1    I4/)') NSURFD
        ENDIF
        CALL LCMSIX(IPMEM,' ',2)
      ENDIF
      ALLOCATE(SURFLX(NSURFD,NGROUP),SURF(NSURFD))
      IF(NSURFD.GT.0) THEN
        CALL LCMSIX(IPMEM,'outgeom',1)
        CALL LCMGET(IPMEM,'SURF',SURF)
        CALL LCMSIX(IPMEM,' ',2)
      ENDIF
      CALL LCMSIX(IPMEM,' ',2)
*----
*  RECOVER INFORMATION FROM caldir DIRECTORY.
*----
      JPMEM=LCMGID(IPMEM,'calc')
      KPMEM=LCMGIL(JPMEM,ICAL)
      CALL LCMSIX(KPMEM,'info',1)
      LSPH=.FALSE.
      LMASL=.FALSE.
      IF(NPARL.GT.0) THEN
        CALL LCMGET(KPMEM,'NLOC',NLOC)
        IF(NLOC.GT.MAXLOC) CALL XABORT('SCRSPH: MAXLOC OVERFLOW')
        CALL LCMGTC(KPMEM,'LOCTYP',4,NLOC,LOCTYP)
        CALL LCMGTC(KPMEM,'LOCKEY',4,NLOC,LOCKEY)
        ALLOCATE(LOCAD(NLOC+1))
        CALL LCMGET(KPMEM,'LOCADR',LOCAD)
        DO ILOC=1,NLOC
          LSPH=LSPH.OR.((LOCTYP(ILOC).EQ.'EQUI').AND.
     1                  (LOCKEY(ILOC).EQ.HEQUI))
          LMASL=LMASL.OR.((LOCTYP(ILOC).EQ.'MASL').AND.
     1                    (LOCKEY(ILOC).EQ.HMASL))
        ENDDO
      ENDIF
      IF((HEQUI.NE.' ').AND.(.NOT.LSPH)) THEN
        WRITE(HSMG,'(46HSCRSPH: UNABLE TO FIND A LOCAL PARAMETER SET O,
     1  25HF TYPE EQUI WITH KEYWORD ,A4,1H.)') HEQUI
        CALL XABORT(HSMG)
      ELSE IF((HMASL.NE.' ').AND.(.NOT.LMASL)) THEN
        WRITE(HSMG,'(46HSCRSPH: UNABLE TO FIND A LOCAL PARAMETER SET O,
     1  25HF TYPE MASL WITH KEYWORD ,A4,1H.)') HMASL
        CALL XABORT(HSMG)
      ENDIF
      ALLOCATE(ISADRX(NMIL),LENGDX(NMIL),LENGDP(NMIL))
      CALL LCMGET(KPMEM,'ISADRX',ISADRX)
      CALL LCMGET(KPMEM,'LENGDX',LENGDX)
      CALL LCMGET(KPMEM,'LENGDP',LENGDP)
      IF(NISOTS.GT.0) THEN
        ALLOCATE(ISOTS(NISOTS*2))
        CALL LCMGET(KPMEM,'ISOTS',ISOTS)
      ENDIF
      CALL LCMSIX(KPMEM,' ',2)
      CALL LCMSIX(KPMEM,'divers',1)
      CALL LCMLEN(KPMEM,'NVDIV',ILONG,ITYLCM)
      IF(ILONG.EQ.0) THEN
        NVDIV=0
      ELSE
        CALL LCMGET(KPMEM,'NVDIV',NVDIV)
      ENDIF
      LSTRD=(B2.EQ.0.0)
      IF(NVDIV.GT.0) THEN
        IF(NVDIV.GT.MAXDIV) CALL XABORT('SCRSPH: MAXDIV OVERFLOW.')
        CALL LCMGTC(KPMEM,'IDVAL',4,NVDIV,IDVAL)
        CALL LCMGET(KPMEM,'VALDIV',VALDIV)
        DO I=1,NVDIV
          IF(IMPX.GT.1) THEN
            WRITE(IOUT,'(9H SCRSPH: ,I3,2X,A,1H=,1P,E13.5)') I,IDVAL(I),
     1      VALDIV(I)
          ENDIF
          IF(IDVAL(I).EQ.'KEFF') THEN
            CALL LCMPUT(IPMAC,'K-EFFECTIVE',1,2,VALDIV(I))
          ELSE IF(IDVAL(I).EQ.'KINF') THEN
            CALL LCMPUT(IPMAC,'K-INFINITY',1,2,VALDIV(I))
          ELSE IF(IDVAL(I).EQ.'B2') THEN
            B2=VALDIV(I)
            LSTRD=(B2.EQ.0.0)
            CALL LCMPUT(IPMAC,'B2  B1HOM',1,2,VALDIV(I))
          ENDIF
        ENDDO
      ENDIF
      CALL LCMSIX(KPMEM,' ',2)
*----
*  ALLOCATE MACROLIB WORKING ARRAYS.
*----
      ALLOCATE(LXS(NREA),NWT0(NMIL,NGROUP),SIGS(NMIL,NGROUP,NL),
     1 SS2D(NMIL,NGROUP,NGROUP,NL),XS(NMIL,NGROUP,NREA))
      NWT0(:NMIL,:NGROUP)=0.0
      SIGS(:NMIL,:NGROUP,:NL)=0.0
      SS2D(:NMIL,:NGROUP,:NGROUP,:NL)=0.0
      XS(:NMIL,:NGROUP,:NREA)=0.0
      LXS(:NREA)=.FALSE.
*----
*  ALLOCATE DELAYED NEUTRON WORKING ARRAYS.
*----
      ALLOCATE(LAMB(NPRC),CHIRS(NGROUP,NPRC,NMIL),BETAR(NPRC,NMIL),
     1 INVELS(NGROUP,NMIL))
      LAMB(:NPRC)=0.0
      CHIRS(:NGROUP,:NPRC,:NMIL)=0.0
      BETAR(:NPRC,:NMIL)=0.0
      INVELS(:NGROUP,:NMIL)=0.0
      CALL LCMSIX(KPMEM,'divers',1)
      CALL LCMLEN(KPMEM,'NPR',ILONG,ITYLCM)
      IF((NPRC.GT.0).AND.(ILONG.EQ.1)) THEN
        CALL LCMGET(KPMEM,'NPR',NPR)
        IF(NPR.NE.NPRC) CALL XABORT('SCRSPH: NPR INCONSISTENCY(1).')
        CALL LCMGET(KPMEM,'LAMBRS',LAMB)
        DO IBM=1,NMIL
          CALL LCMGET(IPMEM,'CHIRS',CHIRS(1,1,IBM))
          CALL LCMGET(IPMEM,'BETARS',BETAR(1,IBM))
          CALL LCMGET(IPMEM,'INVELS',INVELS(1,IBM))
        ENDDO
      ENDIF
      CALL LCMSIX(KPMEM,' ',2)
*----
*  LOOP OVER SAPHYB MIXTURES.
*----
      IF(NADRX.EQ.0) CALL XABORT('SCRSPH: NO ADDRESS SETS AVAILABLE.')
      LPMEM=LCMGID(KPMEM,'mili')
      DO IBM=1,NMIL
        CALL LCMLEL(LPMEM,IBM,ILONG,ITYLCM)
        IF(ILONG.EQ.0) CYCLE
        MPMEM=LCMGIL(LPMEM,IBM)
        IMAC=TOTM(IBM)
        IRES=RESM(IBM)
        IAD=ISADRX(IBM)
        NDATAX=LENGDX(IBM)
        NDATAP=LENGDP(IBM)
        ALLOCATE(FLUXS(NGROUP),RDATA(NDATAX),IDATA(NDATAP))
        CALL LCMGET(MPMEM,'FLUXS',FLUXS)
        CALL LCMGET(MPMEM,'RDATAX',RDATA)
        CALL LCMGET(MPMEM,'IDATAP',IDATA)
        DO I=1,NGROUP
          NWT0(IBM,I)=NWT0(IBM,I)+FLUXS(I)
        ENDDO
        ALLOCATE(SIGSB(NGROUP,NL),SS2DB(NGROUP,NGROUP,NL),
     1  XSB(NGROUP,NREA))
        IF(IMAC.NE.0) THEN
          CALL SPHSXS(NREA,INDX,NADRX,NGROUP,NL,NDATAX,NDATAP,
     1    NISO+IMAC,IAD,IADRX,RDATA,IDATA,NOMREA,SIGSB,SS2DB,XSB,
     2    LXS)
          DO IL=1,NL
            DO I=1,NGROUP
              SIGS(IBM,I,IL)=SIGS(IBM,I,IL)+SIGSB(I,IL)
            ENDDO
          ENDDO
          DO IL=1,NL
            DO J=1,NGROUP
              DO I=1,NGROUP
                SS2D(IBM,I,J,IL)=SS2D(IBM,I,J,IL)+SS2DB(I,J,IL)
              ENDDO
            ENDDO
          ENDDO
          DO IR=1,NREA
            DO I=1,NGROUP
              XS(IBM,I,IR)=XS(IBM,I,IR)+XSB(I,IR)
            ENDDO
          ENDDO
        ELSE IF(NISO.NE.0) THEN
          IF(NISOTS.EQ.0) CALL XABORT('SCRSPH: MISSING CONCES INFO.')
          ALLOCATE(CONCES(NISOTS))
          CALL LCMGET(MPMEM,'CONCES',CONCES)
          DO ISO=1,NISO
            WRITE(TEXT8,'(2A4)') (NOMISO(2*(ISO-1)+I0),I0=1,2)
            ISOKEP=0
            DO IS2=1,NISOTS
              ISOKEP=IS2
              WRITE(TEXT9,'(2A4)') (ISOTS(2*(IS2-1)+I0),I0=1,2)
              IF(TEXT9.EQ.TEXT8) GO TO 10
            ENDDO
            CYCLE
   10       DEN=CONCES(ISOKEP)
            IF(DEN.NE.0.0) THEN
              CALL SPHSXS(NREA,INDX,NADRX,NGROUP,NL,NDATAX,NDATAP,ISO,
     1        IAD,IADRX,RDATA,IDATA,NOMREA,SIGSB,SS2DB,XSB,LXS)
              DO IL=1,NL
                DO I=1,NGROUP
                  SIGS(IBM,I,IL)=SIGS(IBM,I,IL)+DEN*SIGSB(I,IL)
                ENDDO
              ENDDO
              DO IL=1,NL
                DO J=1,NGROUP
                  DO I=1,NGROUP
                    SS2D(IBM,I,J,IL)=SS2D(IBM,I,J,IL)+DEN*SS2DB(I,J,IL)
                  ENDDO
                ENDDO
              ENDDO
              DO IR=1,NREA
                DO I=1,NGROUP
                  XS(IBM,I,IR)=XS(IBM,I,IR)+DEN*XSB(I,IR)
                ENDDO
              ENDDO
            ENDIF
          ENDDO
          DEALLOCATE(CONCES)
          IF(IRES.NE.0) THEN
            CALL SPHSXS(NREA,INDX,NADRX,NGROUP,NL,NDATAX,NDATAP,
     1      NISO+IRES,IAD,IADRX,RDATA,IDATA,NOMREA,SIGSB,SS2DB,XSB,
     2      LXS)
            DO IL=1,NL
              DO I=1,NGROUP
                SIGS(IBM,I,IL)=SIGS(IBM,I,IL)+SIGSB(I,IL)
              ENDDO
            ENDDO
            DO IL=1,NL
              DO J=1,NGROUP
                DO I=1,NGROUP
                  SS2D(IBM,I,J,IL)=SS2D(IBM,I,J,IL)+SS2DB(I,J,IL)
                ENDDO
              ENDDO
            ENDDO
            DO IR=1,NREA
              DO I=1,NGROUP
                XS(IBM,I,IR)=XS(IBM,I,IR)+XSB(I,IR)
              ENDDO
            ENDDO
          ENDIF
        ELSE
          CALL XABORT('SCRSPH: NO MACROSCOPIC SET.')
        ENDIF
        DEALLOCATE(XSB,SS2DB,SIGSB,IDATA,RDATA,FLUXS)
*
*       UP-SCATTERING CORRECTION OF THE MACROLIB.
        IF(ILUPS.EQ.1) THEN
          IRENT0=0
          IRENT1=0
          DO IREA=1,NREA
            IF(NOMREA(IREA).EQ.'TOTALE') IRENT0=IREA
            IF(NOMREA(IREA).EQ.'TOTALE P1') IRENT1=IREA
          ENDDO
          IF(IRENT0.EQ.0) CALL XABORT('SCRSPH: MISSING NTOT0.')
          DO JGR=2,NGROUP
            DO IGR=1,JGR-1 ! IGR < JGR
              FF=NWT0(IBM,JGR)/NWT0(IBM,IGR)
              CSCAT=SS2D(IBM,IGR,JGR,1) ! IGR < JGR
              XS(IBM,IGR,IRENT0)=XS(IBM,IGR,IRENT0)-CSCAT*FF
              XS(IBM,JGR,IRENT0)=XS(IBM,JGR,IRENT0)-CSCAT
              IF((IRENT1.GT.0).AND.(NL.GT.1)) THEN
                CSCAT=SS2D(IBM,IGR,JGR,2)
                XS(IBM,IGR,IRENT1)=XS(IBM,IGR,IRENT1)-CSCAT*FF
                XS(IBM,JGR,IRENT1)=XS(IBM,JGR,IRENT1)-CSCAT
              ENDIF
              DO IL=1,NL
                CSCAT=SS2D(IBM,IGR,JGR,IL)
                SIGS(IBM,IGR,IL)=SIGS(IBM,IGR,IL)-CSCAT*FF
                SIGS(IBM,JGR,IL)=SIGS(IBM,JGR,IL)-CSCAT
                SS2D(IBM,JGR,IGR,IL)=SS2D(IBM,JGR,IGR,IL)-CSCAT*FF
                SS2D(IBM,IGR,JGR,IL)=0.0
              ENDDO
            ENDDO
          ENDDO
        ENDIF
*
        IF(LSPH) THEN
          ALLOCATE(RVALO(LOCAD(NLOC+1)-1))
          CALL LCMGET(MPMEM,'RVALOC',RVALO)
          DO ILOC=1,NLOC
            IF((LOCTYP(ILOC).EQ.'EQUI').AND.(LOCKEY(ILOC).EQ.HEQUI))
     1      THEN
              IF(LOCAD(ILOC+1)-LOCAD(ILOC).NE.NGROUP) THEN
                CALL XABORT('SCRSPH: INVALID NUMBER OF COMPONENTS FOR '
     1          //'SPH FACTORS')
              ENDIF
              DO IGR=1,NGROUP
                SPH(IBM,IGR)=RVALO(LOCAD(ILOC)+IGR-1)
              ENDDO
            ENDIF
          ENDDO
          DEALLOCATE(RVALO)
        ELSE
          SPH(IBM,:NGROUP)=1.0
        ENDIF
        IF(LMASL) THEN
          ALLOCATE(RVALO(LOCAD(NLOC+1)-1))
          CALL LCMGET(MPMEM,'RVALOC',RVALO)
          DO ILOC=1,NLOC
            IF((LOCTYP(ILOC).EQ.'MASL').AND.(LOCKEY(ILOC).EQ.HMASL))
     1      THEN
              IF(LOCAD(ILOC+1)-LOCAD(ILOC).NE.1) THEN
                CALL XABORT('SCRSPH: INVALID NUMBER OF COMPONENTS FOR '
     1          //'MASL')
              ENDIF
              FMASL(IBM)=RVALO(LOCAD(ILOC))
            ENDIF
          ENDDO
          DEALLOCATE(RVALO)
        ENDIF
*
        CALL LCMLEN(MPMEM,'cinetique',ILONG,ITYLCM)
        IF((NPRC.GT.0).AND.(ILONG.NE.0)) THEN
          CALL LCMSIX(MPMEM,'cinetique',1)
          CALL LCMGET(MPMEM,'NPR',NPR)
          IF(NPR.NE.NPRC) CALL XABORT('SCRSPH: NPR INCONSISTENCY(2).')
          CALL LCMGET(MPMEM,'LAMBRS',LAMB)
          CALL LCMGET(MPMEM,'CHIRS',CHIRS(1,1,IBM))
          CALL LCMGET(MPMEM,'BETARS',BETAR(1,IBM))
          CALL LCMGET(MPMEM,'INVELS',INVELS(1,IBM))
          CALL LCMSIX(MPMEM,' ',2)
        ENDIF
*       END OF LOOP OVER SAPHYB MIXTURES
      ENDDO
      IF(NPARL.GT.0) DEALLOCATE(LOCAD)
*----
*  RECOVER DISCONTINUITY FACTOR INFORMATION
*----
      IDF=0
      IF(NSURFD.GT.0) THEN
        IDF=2
        CALL LCMSIX(KPMEM,'outflx',1)
        CALL LCMGET(KPMEM,'SURFLX',SURFLX)
        CALL LCMSIX(KPMEM,' ',2)
        CALL LCMSIX(IPMAC,'ADF',1)
        CALL LCMPUT(IPMAC,'NTYPE',1,1,NSURFD)
        ALLOCATE(HADF(NSURFD),WORK(NMIL,NGROUP))
        DO I=1,NSURFD
          WRITE(HADF(I),'(3HFD_,I5.5)') I
          DO IGR=1,NGROUP
            WORK(:,IGR)=SURFLX(I,IGR)/SURF(I)
          ENDDO
          CALL LCMPUT(IPMAC,HADF(I),NMIL*NGROUP,2,WORK)
        ENDDO
        CALL LCMPTC(IPMAC,'HADF',8,NSURFD,HADF)
        DEALLOCATE(WORK,HADF)
        CALL LCMSIX(IPMAC,' ',2)
      ENDIF
      DEALLOCATE(SURFLX,SURF)
*----
*  IDENTIFY SPECIAL FLUX EDITS
*----
      ALLOCATE(IHEDI(2*NREA))
      NED=0
      DO IREA=1,NREA
        IF(.NOT.LXS(IREA)) CYCLE
        IF(NOMREA(IREA).EQ.'TOTALE') CYCLE
        IF(NOMREA(IREA).EQ.'TOTALE P1') CYCLE
        IF(NOMREA(IREA).EQ.'EXCESS') CYCLE
        IF(NOMREA(IREA).EQ.'FISSION') CYCLE
        IF(NOMREA(IREA).EQ.'SPECTRE') CYCLE
        IF(NOMREA(IREA).EQ.'NU*FISSION') CYCLE
        IF(NOMREA(IREA)(:7).EQ.'ENERGIE') CYCLE
        IF(NOMREA(IREA).EQ.'SELF') CYCLE
        IF(NOMREA(IREA).EQ.'TRANSP-CORR') CYCLE
        IF(NOMREA(IREA).EQ.'FUITES') CYCLE
        IF(NOMREA(IREA).EQ.'DIFFUSION') CYCLE
        IF(NOMREA(IREA).EQ.'TRANSFERT') CYCLE
        NED=NED+1
        READ(NOMREA(IREA),'(2A4)') IHEDI(2*NED-1),IHEDI(2*NED)
      ENDDO
*----
*  STORE MACROLIB.
*----
      CALL LCMPUT(IPMAC,'VOLUME',NMIL,2,XVOLM)
      IF(LMASL) CALL LCMPUT(IPMAC,'MASL',NMIL,2,FMASL)
      IF(NPRC.GT.0) CALL LCMPUT(IPMAC,'LAMBDA-D',NPRC,2,LAMB)
      IFISS=0
      ITRANC=0
      LDIFF=.FALSE.
      NW=0
      ALLOCATE(STR(NMIL),WRK(NMIL))
      SIGS0(:NMIL,:NGROUP)=0.0
      JPMAC=LCMLID(IPMAC,'GROUP',NGROUP)
      DO IGR=1,NGROUP
        STR(:NMIL)=0.0
        KPMAC=LCMDIL(JPMAC,IGR)
        CALL LCMPUT(KPMAC,'FLUX-INTG',NMIL,2,NWT0(1,IGR))
        IF(NPRC.GT.0) THEN
          DO IBM=1,NMIL
            WRK(IBM)=INVELS(IGR,IBM)
          ENDDO
          CALL LCMPUT(KPMAC,'OVERV',NMIL,2,WRK)
        ENDIF
        DO IREA=1,NREA
          IF(.NOT.LXS(IREA)) CYCLE
          IF(NOMREA(IREA).EQ.'TOTALE') THEN
            IF(LSTRD) THEN
              DO IBM=1,NMIL
                STR(IBM)=STR(IBM)+XS(IBM,IGR,IREA)
              ENDDO
            ENDIF
            CALL LCMPUT(KPMAC,'NTOT0',NMIL,2,XS(1,IGR,IREA))
          ELSE IF(NOMREA(IREA).EQ.'TOTALE P1') THEN
            NW=1
            CALL LCMPUT(KPMAC,'NTOT1',NMIL,2,XS(1,IGR,IREA))
          ELSE IF(NOMREA(IREA).EQ.'EXCESS') THEN
*           correct scattering XS with excess XS
            DO IBM=1,NMIL
              SIGS0(IBM,IGR)=SIGS0(IBM,IGR)+XS(IBM,IGR,IREA)
            ENDDO
            CALL LCMPUT(KPMAC,'N2N',NMIL,2,XS(1,IGR,IREA))
          ELSE IF(NOMREA(IREA).EQ.'FISSION') THEN
            CALL LCMPUT(KPMAC,'NFTOT',NMIL,2,XS(1,IGR,IREA))
          ELSE IF(NOMREA(IREA).EQ.'SPECTRE') THEN
            CALL LCMPUT(KPMAC,'CHI',NMIL,2,XS(1,IGR,IREA))
            DO IPRC=1,NPRC
              DO IBM=1,NMIL
                WRK(IBM)=CHIRS(IGR,IPRC,IBM)
              ENDDO
              WRITE(TEXT12,'(A3,I2.2)') 'CHI',IPRC
              CALL LCMPUT(KPMAC,TEXT12,NMIL,2,WRK)
            ENDDO
          ELSE IF(NOMREA(IREA).EQ.'NU*FISSION') THEN
            IFISS=1
            CALL LCMPUT(KPMAC,'NUSIGF',NMIL,2,XS(1,IGR,IREA))
            DO IPRC=1,NPRC
              DO IBM=1,NMIL
                WRK(IBM)=XS(IBM,IGR,IREA)*BETAR(IPRC,IBM)
              ENDDO
              WRITE(TEXT12,'(A6,I2.2)') 'NUSIGF',IPRC
              CALL LCMPUT(KPMAC,TEXT12,NMIL,2,WRK)
            ENDDO
          ELSE IF(NOMREA(IREA).EQ.'ENERGIE') THEN
            CALL LCMPUT(KPMAC,'H-FACTOR',NMIL,2,XS(1,IGR,IREA))
          ELSE IF(NOMREA(IREA).EQ.'SELF') THEN
            CALL LCMPUT(KPMAC,'SIGW00',NMIL,2,XS(1,IGR,IREA))
          ELSE IF(NOMREA(IREA).EQ.'TRANSP-CORR') THEN
            ITRANC=2
            IF(LSTRD) THEN
              DO IBM=1,NMIL
                STR(IBM)=STR(IBM)-XS(IBM,IGR,IREA)
              ENDDO
            ENDIF
            CALL LCMPUT(KPMAC,'TRANC',NMIL,2,XS(1,IGR,IREA))
          ELSE IF(NOMREA(IREA).EQ.'FUITES') THEN
            LDIFF=LSTRD
            IF(.NOT.LSTRD) THEN
              DO IBM=1,NMIL
                LDIFF=LDIFF.OR.(XS(IBM,IGR,IREA).NE.0.0)
                STR(IBM)=XS(IBM,IGR,IREA)/B2
              ENDDO
            ENDIF
          ELSE IF(NOMREA(IREA).EQ.'DIFFUSION') THEN
            DO IL=1,NL
              WRITE(CM,'(I2.2)') IL-1
              IF(IL.EQ.1) THEN
                DO IBM=1,NMIL
                  SIGS0(IBM,IGR)=SIGS0(IBM,IGR)+SIGS(IBM,IGR,IL)
                ENDDO
              ELSE
                CALL LCMPUT(KPMAC,'SIGS'//CM,NMIL,2,SIGS(1,IGR,IL))
              ENDIF
            ENDDO
          ELSE IF(NOMREA(IREA).EQ.'TRANSFERT') THEN
            ALLOCATE(SCAT(NGROUP*NMIL),GAR(NMIL))
            DO IL=1,NL
              WRITE(CM,'(I2.2)') IL-1
              IPOSDE=0
              DO IBM=1,NMIL
                IPOS(IBM)=IPOSDE+1
                IGMIN=IGR
                IGMAX=IGR
                DO JGR=NGROUP,1,-1
                  IF(SS2D(IBM,IGR,JGR,IL).NE.0.0) THEN
                    IGMIN=MIN(IGMIN,JGR)
                    IGMAX=MAX(IGMAX,JGR)
                  ENDIF
                ENDDO
                IJJM(IBM)=IGMAX
                NJJM(IBM)=IGMAX-IGMIN+1
                DO JGR=IGMAX,IGMIN,-1
                  IPOSDE=IPOSDE+1
                  SCAT(IPOSDE)=SS2D(IBM,IGR,JGR,IL)
                ENDDO
                GAR(IBM)=SCAT(IPOS(IBM)+IJJM(IBM)-IGR)
              ENDDO
              CALL LCMPUT(KPMAC,'SCAT'//CM,IPOSDE,2,SCAT)
              CALL LCMPUT(KPMAC,'NJJS'//CM,NMIL,1,NJJM)
              CALL LCMPUT(KPMAC,'IJJS'//CM,NMIL,1,IJJM)
              CALL LCMPUT(KPMAC,'IPOS'//CM,NMIL,1,IPOS)
              CALL LCMPUT(KPMAC,'SIGW'//CM,NMIL,2,GAR)
            ENDDO
            DEALLOCATE(GAR,SCAT)
          ELSE
            CALL LCMPUT(KPMAC,NOMREA(IREA),NMIL,2,XS(1,IGR,IREA))
          ENDIF
        ENDDO
        IF(LSTRD) THEN
          IF((ITRANC.EQ.0).AND.(NL.GT.1)) THEN
*           Apollo-type transport correction
            DO IBM=1,NMIL
              STR(IBM)=STR(IBM)-SIGS(IBM,IGR,2)
            ENDDO
          ENDIF
          DO IBM=1,NMIL
            STR(IBM)=1.0/(3.0*STR(IBM))
          ENDDO
          LDIFF=.TRUE.
        ENDIF
        IF((ITRANC.EQ.0).AND.(NL.GT.1)) THEN
*         Apollo-type transport correction
          IF(IGR.EQ.NGROUP) ITRANC=2
          CALL LCMPUT(KPMAC,'TRANC',NMIL,2,SIGS(1,IGR,2))
        ENDIF
        IF(LDIFF) CALL LCMPUT(KPMAC,'DIFF',NMIL,2,STR)
      ENDDO
      DEALLOCATE(WRK,STR)
*----
*  RELEASE MEMORY
*----
      DEALLOCATE(INVELS,BETAR,CHIRS,LAMB,LXS,XS,SS2D,SIGS,NWT0,LENGDP,
     1 LENGDX,ISADRX,XVOLM)
      IF(NISOTS.GT.0) DEALLOCATE(ISOTS)
      IF(NADRX.GT.0) DEALLOCATE(IADRX)
      IF(NISO.GT.0) DEALLOCATE(NOMISO)
      DEALLOCATE(RESM,TOTM)
*----
*  SAVE SCATTERING P0 INFO
*----
      DO IGR=1,NGROUP
        KPMAC=LCMDIL(JPMAC,IGR)
        CALL LCMPUT(KPMAC,'SIGS00',NMIL,2,SIGS0(1,IGR))
      ENDDO
*----
*  WRITE STATE VECTOR
*----
      TEXT12='L_MACROLIB'
      CALL LCMPTC(IPMAC,'SIGNATURE',12,TEXT12)
      ISTATE(:NSTATE)=0
      ISTATE(1)=NGROUP
      ISTATE(2)=NMIL
      ISTATE(3)=NL ! 1+scattering anisotropy
      ISTATE(4)=IFISS
      ISTATE(5)=NED
      ISTATE(6)=ITRANC
      ISTATE(7)=NPRC
      IF(LDIFF) ISTATE(9)=1
      ISTATE(10)=NW
      ISTATE(12)=IDF
      CALL LCMPUT(IPMAC,'STATE-VECTOR',NSTATE,1,ISTATE)
      IF(NED.GT.0) CALL LCMPUT(IPMAC,'ADDXSNAME-P0',2*NED,3,IHEDI)
      DEALLOCATE(IHEDI)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FMASL,SIGS0)
      DEALLOCATE(NOMMIL,IJJM,NJJM,IPOS)
      RETURN
      END
