*DECK LIBEAR
      SUBROUTINE LIBEAR(CFILNA,MAXR,NEL,NMDEPL,ITNAM,ITZEA,KPAX,BPAX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read depletion data on an APOLIB-2 formatted library.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): A. Hebert
*
*Parameters: input
* CFILNA  APOLIB-2 file name.
* MAXR    number of reaction types.
* NEL     number of isotopes on library.
* NMDEPL  names of reactions:
*           NMDEPL(1)='DECAY'; NMDEPL(2)='NFTOT';
*           NMDEPL(3)='NG'   ; NMDEPL(4)='N2N';
*           etc.
*
*Parameters: output
* ITNAM   reactive isotope names in chain.
* ITZEA   6-digit nuclide identifier:
*         atomic number z*10000 (digits) + mass number a*10 +
*         energy state (0 = ground state, 1 = first state, etc.).
* KPAX    complete reaction type matrix.
* BPAX    complete branching ratio matrix.
*
*-----------------------------------------------------------------------
*
*----
*  INPUT FORMAT
*----
*    LIB: APLIB2 FIL: CFILNA CHAIN
*    [[ hnamson
*    [ FROM  [[ { DECAY | reaction } yield hnampar ]] ]
*    ]]
*    ENDCHAIN
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER CFILNA*(*),NMDEPL(MAXR)*8
      INTEGER MAXR,NEL,ITNAM(3,NEL),ITZEA(NEL),KPAX(NEL+MAXR,NEL)
      REAL BPAX(NEL+MAXR,NEL)
*
      EXTERNAL LIBA21
      INTEGER ISFICH(3),NITCA(5)
      PARAMETER (IOUT=6)
      CHARACTER TEXT20*20,NOMOBJ*20,TEXT8*8,TEXT12*12,TYPOBJ*8,TYPSEG*8,
     > HNISOR*20,HITNAM*20,HSMG*131,TEXT16*16
      LOGICAL LPHEAD,LPCONS,LPNUMF,LTEST,LPFIX
      DOUBLE PRECISION DBLINP
*----
*  SCRATCH STORAGE ALLOCATION
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: VINTE,NOMOB,KDS,LGS,ITCARO,
     1 NOM,IA,IZ,NFG,IZSECT,ISECTT,IKEEP
      REAL, ALLOCATABLE, DIMENSION(:) :: GAMMA
      TYPE(C_PTR) ICHDIM_PTR,ICHTYP_PTR,ICHDKL_PTR,TSEGM_PTR
      INTEGER, POINTER, DIMENSION(:) :: ICHDIM,ICHTYP,ICHDKL,ITSEGM
      REAL, POINTER, DIMENSION(:) :: RTSEGM
*
      INTEGER TKCARO(31)
      SAVE TKCARO
      DATA TKCARO /
     &   0,   1,   2,   3,  4,   5,  6,  30,   7,  -8,
     &   9, -10,  11, -12, 13, -14, 15,  16, -17,  18,
     & -19,  20, -21,  22, 23, -24, 25, -26,  27, -28,
     &  29   /
*----
*  OPEN AND PROBE THE APOLIB-2 FILE.
*----
      CALL AEXTPA(CFILNA,ISFICH)
      IADRES=ISFICH(1)
      NBOBJ=ISFICH(2)
      LBLOC=ISFICH(3)
      IUNIT=KDROPN(CFILNA,2,4,LBLOC)
      IF(IUNIT.LE.0) THEN
         TEXT12=CFILNA
         CALL XABORT('LIBEAR: APOLLO-2 LIBRARY '//TEXT12//' CANNOT B'//
     1   'E OPENED')
      ENDIF
*----
*  INDEX THE APOLIB-2 FILE.
*----
      ALLOCATE(VINTE(2*NBOBJ))
      CALL AEXDIR(IUNIT,LBLOC,VINTE,IADRES,2*NBOBJ)
      IDKNO=1-TKCARO(14)
      IDKTY=1-TKCARO(21)
      IDKDS=1-TKCARO(10)
      IDKTS=1-TKCARO(23)
      IDKNS=TKCARO(2)+1
      IDKLS=TKCARO(8)
*
      NSEGM=0
      NMGY=0
      NISOT=0
      ALLOCATE(NOMOB(5*(NBOBJ-3)),KDS(NBOBJ-3),LGS(NBOBJ-3))
      LPHEAD=.FALSE.
      LPCONS=.FALSE.
      LPNUMF=.FALSE.
      DO 80 IOBJ=3,NBOBJ
      IDKOBJ=VINTE(2*IOBJ-1)
      LGSEG=VINTE(2*IOBJ)+1
      ALLOCATE(ITCARO(LGSEG))
      CALL AEXDIR(IUNIT,LBLOC,ITCARO,IDKOBJ,LGSEG)
      IDK=ITCARO(IDKNO)
      CALL AEXCPC(IDK,20,ITCARO(1),NOMOBJ)
      IDK=ITCARO(IDKTY)
      CALL AEXCPC(IDK,8,ITCARO(1),TYPOBJ)
      JDKDS=ITCARO(IDKDS)
      JDKTS=ITCARO(IDKTS)
      NS=ITCARO(IDKNS)
      IF(TYPOBJ.EQ.'APOLIB') THEN
         DO 70 IS=1,NS
         IDK=JDKTS+8*(IS-1)
         CALL AEXCPC(IDK,8,ITCARO(1),TYPSEG)
         LTESTS=ITCARO(IDKLS+IS)
         IF(LTESTS.LE.0) GO TO 70
         JDKS=ITCARO(JDKDS+IS)
         CALL AEXTRT(LIBA21,TYPSEG,NBRTYP,ICHDIM_PTR,ICHTYP_PTR,
     1   ICHDKL_PTR)
         CALL C_F_POINTER(ICHDIM_PTR,ICHDIM,(/ NBRTYP /))
         CALL C_F_POINTER(ICHTYP_PTR,ICHTYP,(/ NBRTYP /))
         CALL C_F_POINTER(ICHDKL_PTR,ICHDKL,(/ NBRTYP /))
         TSEGM_PTR=LCMARA(LTESTS+1)
         CALL C_F_POINTER(TSEGM_PTR,ITSEGM,(/ LTESTS+1 /))
         CALL C_F_POINTER(TSEGM_PTR,RTSEGM,(/ LTESTS+1 /))
         CALL AEXDIR(IUNIT,LBLOC,ITSEGM,JDKS,LTESTS+1)
         IF(TYPSEG.EQ.'PHEAD') THEN
            LPHEAD=.TRUE.
            CALL AEXGNV(3,ITSEGM,ICHDIM,ICHTYP,
     1      ICHDKL,IDK,NV)
            IF(NV.EQ.0) THEN
               TEXT12=CFILNA
               CALL XABORT('LIBEAR: NO ISOTOPES PRESENT ON APOLIB-2 '//
     1         'FILE NAMED: '//TEXT12)
            ENDIF
            NISOT=NV/20
            IF(NISOT.NE.NEL) CALL XABORT('LIBEAR: INVALID NEL.')
            ALLOCATE(NOM(5*NISOT))
            DO 20 ISO=1,NISOT
            ISO2=(ISO-1)*5+1
            CALL AEXCPC(0,20,ITSEGM(IDK+ISO2-1),HNISOR)
            CALL LCMCAR(HNISOR,.TRUE.,NOM(ISO2))
            READ(HNISOR,'(3A4)') (ITNAM(II,ISO),II=1,3)
   20       CONTINUE
         ELSE IF(TYPSEG.EQ.'PCONST') THEN
            LPCONS=.TRUE.
            CALL AEXGNV(1,ITSEGM,ICHDIM,ICHTYP,ICHDKL,IDK,NV)
            ALLOCATE(IA(NV))
            DO 30 I=1,NV
            IA(I)=ITSEGM(IDK+I-1)
   30       CONTINUE
            CALL AEXGNV(3,ITSEGM,ICHDIM,ICHTYP,ICHDKL,IDK,NV)
            ALLOCATE(IZ(NV))
            DO 40 I=1,NV
            IZ(I)=ITSEGM(IDK+I-1)
   40       CONTINUE
            CALL AEXGNV(5,ITSEGM,ICHDIM,ICHTYP,ICHDKL,IDK,NV)
            ALLOCATE(NFG(NV))
            DO 50 I=1,NV
            NFG(I)=ITSEGM(IDK+I-1)
   50       CONTINUE
         ELSE IF(TYPSEG.EQ.'PNUMF') THEN
            LPNUMF=.TRUE.
            CALL AEXGNV(1,ITSEGM,ICHDIM,ICHTYP,ICHDKL,IDK,NBFISS)
            NBFISS=NBFISS/7
            CALL AEXGNV(4,ITSEGM,ICHDIM,ICHTYP,ICHDKL,IDK,NBPF)
            NBPF=NBPF/7
            CALL AEXGNV(7,ITSEGM,ICHDIM,ICHTYP,ICHDKL,IDK,NV)
            ALLOCATE(GAMMA(NV))
            DO 60 I=1,NV
            GAMMA(I)=RTSEGM(IDK+I-1)
   60       CONTINUE
            NMGY=NV/(NBFISS*NBPF)
         ENDIF
         CALL LCMDRD(TSEGM_PTR)
         CALL LCMDRD(ICHDIM_PTR)
         CALL LCMDRD(ICHTYP_PTR)
         CALL LCMDRD(ICHDKL_PTR)
   70    CONTINUE
      ELSE IF(TYPOBJ.EQ.'APOLIBE') THEN
         NSEGM=NSEGM+1
         ISO2=(NSEGM-1)*5+1
         CALL LCMCAR(NOMOBJ,.TRUE.,NOMOB(ISO2))
         KDS(NSEGM)=IDKOBJ
         LGS(NSEGM)=LGSEG
      ELSE
         CALL XABORT('LIBEAR: WEIRD SEGMENT TYPE: '//TYPOBJ//' (1).')
      ENDIF
      DEALLOCATE(ITCARO)
   80 CONTINUE
      IF(.NOT.LPHEAD) CALL XABORT('LIBEAR: PHEAD SEGMENT NOT FOUND.')
      IF(.NOT.LPCONS) CALL XABORT('LIBEAR: PCONST SEGMENT NOT FOUND.')
      IF(.NOT.LPNUMF) CALL XABORT('LIBEAR: PNUMF SEGMENT NOT FOUND.')
      DEALLOCATE(VINTE)
*----
*  SET THE CORRESPONDANCE BETWEEN THE APOLIB AND THE LIST OF ISOTOPES.
*----
      KISEG2=0
      DO 260 ISO=1,NISOT
      ISO2=(ISO-1)*5+1
      CALL LCMCAR(TEXT16,.FALSE.,NOM(ISO2))
      TEXT20='ISOTOP'//TEXT16(:14)
      CALL LCMCAR(TEXT20,.TRUE.,NITCA(1))
      DO 90 ISEG=1,NSEGM
      ISEG2=(ISEG-1)*5+1
      IF(NITCA(1).EQ.NOMOB(ISEG2)) THEN
         IF(NITCA(2).EQ.NOMOB(ISEG2+1)) THEN
            IF(NITCA(3).EQ.NOMOB(ISEG2+2)) THEN
               IF(NITCA(4).EQ.NOMOB(ISEG2+3)) THEN
                  IF(NITCA(5).EQ.NOMOB(ISEG2+4)) THEN
                     KISEG2=ISEG
                     GO TO 100
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
   90 CONTINUE
      WRITE (HSMG,500) HNISOR,CFILNA
      CALL XABORT(HSMG)
  100 KISEG3=0
      TEXT20='PHYSIQ'//TEXT16(:14)
      CALL LCMCAR(TEXT20,.TRUE.,NITCA(1))
      DO 110 ISEG=1,NSEGM
      ISEG2=(ISEG-1)*5+1
      IF(NITCA(1).EQ.NOMOB(ISEG2)) THEN
         IF(NITCA(2).EQ.NOMOB(ISEG2+1)) THEN
            IF(NITCA(3).EQ.NOMOB(ISEG2+2)) THEN
               IF(NITCA(4).EQ.NOMOB(ISEG2+3)) THEN
                  IF(NITCA(5).EQ.NOMOB(ISEG2+4)) THEN
                     KISEG3=ISEG
                     GO TO 120
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
  110 CONTINUE
*----
*  ACTIVATION OF CORRESPONDING 'ISOTOP'//NAME SEGMENT.
*----
  120 IDKOBJ=KDS(KISEG2)
      LGSEG=LGS(KISEG2)
      ALLOCATE(ITCARO(LGSEG))
      CALL AEXDIR(IUNIT,LBLOC,ITCARO,IDKOBJ,LGSEG)
      IDK=ITCARO(IDKNO)
      CALL AEXCPC(IDK,20,ITCARO(1),NOMOBJ)
      IDK=ITCARO(IDKTY)
      CALL AEXCPC(IDK,8,ITCARO(1),TYPOBJ)
      JDKDS=ITCARO(IDKDS)
      JDKTS=ITCARO(IDKTS)
      NS=ITCARO(IDKNS)
      NSETOT=0
      NPHY=0
*----
*  RECOVER THE INFINITE DILUTION CROSS SECTION NUMEROTATION.
*----
      LPFIX=.FALSE.
      DO 160 IS=1,NS
      IDK=JDKTS+8*(IS-1)
      CALL AEXCPC(IDK,8,ITCARO(1),TYPSEG)
      LTESTS=ITCARO(IDKLS+IS)
      IF(LTESTS.LE.0) GO TO 160
      JDKS=ITCARO(JDKDS+IS)
      CALL AEXTRT(LIBA21,TYPSEG,NBRTYP,ICHDIM_PTR,ICHTYP_PTR,ICHDKL_PTR)
      CALL C_F_POINTER(ICHDIM_PTR,ICHDIM,(/ NBRTYP /))
      CALL C_F_POINTER(ICHTYP_PTR,ICHTYP,(/ NBRTYP /))
      CALL C_F_POINTER(ICHDKL_PTR,ICHDKL,(/ NBRTYP /))
      TSEGM_PTR=LCMARA(LTESTS+1)
      CALL C_F_POINTER(TSEGM_PTR,ITSEGM,(/ LTESTS+1 /))
      CALL C_F_POINTER(TSEGM_PTR,RTSEGM,(/ LTESTS+1 /))
      CALL AEXDIR(IUNIT,LBLOC,ITSEGM,JDKS,LTESTS+1)
      IF(TYPSEG.EQ.'PFIX') THEN
         LPFIX=.TRUE.
*        NG ENERGY.
         CALL AEXGNV(11,ITSEGM,ICHDIM,ICHTYP,ICHDKL,IDK,NV)
         IF(NV.NE.0) THEN
            IF(RTSEGM(IDK).NE.0.0) THEN
               KPAX(NEL+3,ISO)=1
               BPAX(NEL+3,ISO)=RTSEGM(IDK)
            ENDIF
         ENDIF
*        AVAILABLE CROSS SECTION TYPES.
         CALL AEXGNV(12,ITSEGM,ICHDIM,ICHTYP,ICHDKL,IDK,NSECTT)
         ALLOCATE(IZSECT(NSECTT))
         NSETOT=0
         NPHY=MAX(0,NSECTT-5)
         DO 130 I=1,NSECTT
         IZSECT(I)=ITSEGM(IDK+I-1)
         IF((IZSECT(I).NE.0).AND.(I.LE.5)) NSETOT=NSETOT+1
  130    CONTINUE
*        FISSION ENERGIES.
         CALL AEXGNV(20,ITSEGM,ICHDIM,ICHTYP,ICHDKL,IDK,NMGEF)
         IF(NMGEF.NE.0) THEN
            IF(RTSEGM(IDK+NMGEF-1).NE.0.0) THEN
               KPAX(NEL+2,ISO)=1
               BPAX(NEL+2,ISO)=RTSEGM(IDK+NMGEF-1)
            ENDIF
         ENDIF
*        RADIOACTIVE DECAY CONSTANTS.
         CALL AEXGNV(22,ITSEGM,ICHDIM,ICHTYP,ICHDKL,IDK,NCHANN)
         SUM=0.0
         DO 140 I=1,NCHANN
         SUM=SUM+RTSEGM(IDK+I-1)
  140    CONTINUE
         IF(SUM.NE.0.0) BPAX(NEL+1,ISO)=SUM*1.0E8
*        X-S NAMES.
         CALL AEXGNV(26,ITSEGM,ICHDIM,ICHTYP,ICHDKL,IDK,NV)
         IF(NV/8.NE.NSECTT) CALL XABORT('LIBEAR: INVALID TYPSECT.')
         ALLOCATE(ISECTT(2*NSECTT))
         III=0
         DO 150 I=1,NSECTT
         I2=(I-1)*2+1
         IF(IZSECT(I).NE.0) THEN
            III=III+1
            I3=(III-1)*2+1
            CALL AEXCPC(0,8,ITSEGM(IDK+I2-1),TEXT8)
            CALL LCMCAR(TEXT8,.TRUE.,ISECTT(I3))
         ENDIF
  150    CONTINUE
      ENDIF
      CALL LCMDRD(TSEGM_PTR)
      CALL LCMDRD(ICHDIM_PTR)
      CALL LCMDRD(ICHTYP_PTR)
      CALL LCMDRD(ICHDKL_PTR)
  160 CONTINUE
      IF(.NOT.LPFIX) CALL XABORT('LIBEAR: NO PFIX SEGMENT.')
*----
*  TEST THE INFINITE DILUTION CROSS SECTIONS.
*----
      ITSEC=0
      DO 210 IS=1,NS
      LTESTS=ITCARO(IDKLS+IS)
      IF(LTESTS.LE.0) GO TO 210
      IDK=JDKTS+8*(IS-1)
      CALL AEXCPC(IDK,8,ITCARO(1),TYPSEG)
      JDKS=ITCARO(JDKDS+IS)
      CALL AEXTRT(LIBA21,TYPSEG,NBRTYP,ICHDIM_PTR,ICHTYP_PTR,ICHDKL_PTR)
      CALL C_F_POINTER(ICHDIM_PTR,ICHDIM,(/ NBRTYP /))
      CALL C_F_POINTER(ICHTYP_PTR,ICHTYP,(/ NBRTYP /))
      CALL C_F_POINTER(ICHDKL_PTR,ICHDKL,(/ NBRTYP /))
      TSEGM_PTR=LCMARA(LTESTS+1)
      CALL C_F_POINTER(TSEGM_PTR,ITSEGM,(/ LTESTS+1 /))
      CALL C_F_POINTER(TSEGM_PTR,RTSEGM,(/ LTESTS+1 /))
      CALL AEXDIR(IUNIT,LBLOC,ITSEGM,JDKS,LTESTS+1)
      IF(TYPSEG.EQ.'PSECT') THEN
*        RECOVER A VECTOR CROSS SECTION.
         ITSEC=ITSEC+1
         IF(ITSEC.GT.NSETOT) GO TO 200
         CALL AEXGNV(1,ITSEGM,ICHDIM,ICHTYP,
     1   ICHDKL,IDK,NV)
         I3=(ITSEC-1)*2+1
         CALL LCMCAR(TEXT8,.FALSE.,ISECTT(I3))
         IF(TEXT8.EQ.'SIGA') THEN
            LTEST=.FALSE.
            DO 170 IG=1,NV
            LTEST=LTEST.OR.(RTSEGM(IDK+IG-1).NE.0.0)
  170       CONTINUE
            IF(LTEST) KPAX(NEL+3,ISO)=1
         ELSE IF(TEXT8.EQ.'NEXCESS') THEN
            LTEST=.FALSE.
            DO 180 IG=1,NV
            LTEST=LTEST.OR.(RTSEGM(IDK+IG-1).NE.0.0)
  180       CONTINUE
            IF(LTEST) KPAX(NEL+4,ISO)=1
         ELSE IF(TEXT8.EQ.'SIGF') THEN
            LTEST=.FALSE.
            DO 190 IG=1,NV
            LTEST=LTEST.OR.(RTSEGM(IDK+IG-1).NE.0.0)
  190       CONTINUE
            IF(LTEST) KPAX(NEL+2,ISO)=1
         ENDIF
      ENDIF
  200 CALL LCMDRD(TSEGM_PTR)
      CALL LCMDRD(ICHDIM_PTR)
      CALL LCMDRD(ICHTYP_PTR)
      CALL LCMDRD(ICHDKL_PTR)
  210 CONTINUE
*----
*  TEST THE PRODUCTION REACTIONS.
*----
      IF(NPHY.GE.1) THEN
         IF(KISEG3.EQ.0) CALL XABORT('LIBEAR: INVALID PRODUCTION X-S ('
     1   //'1).')
         IDKOBJ=KDS(KISEG3)
         LGSEG=LGS(KISEG3)
         CALL AEXDIR(IUNIT,LBLOC,ITCARO,IDKOBJ,LGSEG)
         LDKDS=ITCARO(IDKDS)
         LDKTS=ITCARO(IDKTS)
         NS=ITCARO(IDKNS)
         IF(NS.NE.NPHY) CALL XABORT('LIBEAR: INVALID PRODUCTION X-S(2)'
     1   //'.')
         DO 240 IPHY=1,NPHY
         IDK=LDKTS+8*(IPHY-1)
         CALL AEXCPC(IDK,8,ITCARO(1),TYPSEG)
         IF(TYPSEG.NE.'PSECT') CALL XABORT('LIBEAR: INVALID PRODUCTION'
     1   //' X-S(3).')
         LNGS=ITCARO(IDKLS+IPHY)
         IF(LNGS.LE.0) CALL XABORT('LIBEAR: INVALID PRODUCTION X-S(4).')
         LDKS=ITCARO(LDKDS+IPHY)
         CALL AEXTRT(LIBA21,TYPSEG,NBRTYP,ICHDIM_PTR,ICHTYP_PTR,
     1   ICHDKL_PTR)
         CALL C_F_POINTER(ICHDIM_PTR,ICHDIM,(/ NBRTYP /))
         CALL C_F_POINTER(ICHTYP_PTR,ICHTYP,(/ NBRTYP /))
         CALL C_F_POINTER(ICHDKL_PTR,ICHDKL,(/ NBRTYP /))
         TSEGM_PTR=LCMARA(LNGS+1)
         CALL C_F_POINTER(TSEGM_PTR,ITSEGM,(/ LNGS+1 /))
         CALL C_F_POINTER(TSEGM_PTR,RTSEGM,(/ LNGS+1 /))
         CALL AEXDIR(IUNIT,LBLOC,ITSEGM,LDKS,LNGS+1)
         CALL AEXGNV(1,ITSEGM,ICHDIM,ICHTYP,ICHDKL,IDK,NV)
         CALL LCMDRD(ICHDIM_PTR)
         CALL LCMDRD(ICHTYP_PTR)
         CALL LCMDRD(ICHDKL_PTR)
         I3=(NSETOT+IPHY-1)*2+1
         CALL LCMCAR(TEXT8,.FALSE.,ISECTT(I3))
         IF(TEXT8.EQ.'CREA-P') THEN
            LTEST=.FALSE.
            DO 220 IG=1,NV
            LTEST=LTEST.OR.(RTSEGM(IDK+IG-1).NE.0.0)
  220       CONTINUE
            IF(LTEST) KPAX(NEL+8,ISO)=1
         ELSE IF(TEXT8.EQ.'CREA-H2') THEN
            LTEST=.FALSE.
            DO 225 IG=1,NV
            LTEST=LTEST.OR.(RTSEGM(IDK+IG-1).NE.0.0)
  225       CONTINUE
            IF(LTEST) KPAX(NEL+11,ISO)=1
         ELSE IF(TEXT8.EQ.'CREA-H3') THEN
            LTEST=.FALSE.
            DO 230 IG=1,NV
            LTEST=LTEST.OR.(RTSEGM(IDK+IG-1).NE.0.0)
  230       CONTINUE
            IF(LTEST) KPAX(NEL+12,ISO)=1
         ENDIF
         CALL LCMDRD(TSEGM_PTR)
  240    CONTINUE
      ENDIF
      DEALLOCATE(ITCARO)
*----
*  SET OTHER INFORMATION.
*----
      ITZEA(ISO)=IZ(ISO)*10000+IA(ISO)*10
      IPF=NFG(ISO)
      IF(IPF.LT.0) THEN
         KPAX(NEL+2,ISO)=-1
         DO 250 JSO=1,NISOT
         IFI=NFG(JSO)
         IF(IFI.GT.0) THEN
            IOFSET=((-IPF-1)*NBFISS+(IFI-1))*NMGY+NMGY
            BPAX(ISO,JSO)=GAMMA(IOFSET)
            IF(BPAX(ISO,JSO).NE.0.0) KPAX(ISO,JSO)=2
         ENDIF
  250    CONTINUE
      ENDIF
      DEALLOCATE(ISECTT,IZSECT)
  260 CONTINUE
*
      DEALLOCATE(LGS,KDS,NOMOB,GAMMA,NFG,IZ,IA,NOM)
      IERR=KDRCLS(IUNIT,1)
      IF(IERR.LT.0) THEN
         TEXT12=CFILNA
         CALL XABORT('LIBEAR: APOLLO-2 LIBRARY '//TEXT12//' CANNOT B'//
     1   'E CLOSED')
      ENDIF
*----
*  RECOVER INFORMATION FROM INPUT DATA STREAM.
*----
      ALLOCATE(IKEEP(NEL))
      IKEEP(:NEL)=0
      TEXT12=' '
      CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
      IF(INDIC.NE.3.OR.TEXT12.NE.'CHAIN')
     >  CALL XABORT('LIBEAR: KEYWORD CHAIN MISSING')
      CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
      DO 340 IEL=1,NEL
      IF(TEXT12.EQ.'ENDCHAIN') GO TO 350
      IF(INDIC.NE.3) CALL XABORT('LIBEAR: ISOTOPE NAME hnamson MISSING')
      I1=INDEX(TEXT12,'_')
      HNISOR=' '
      IF(I1.EQ.0) THEN
         HNISOR(:12)=TEXT12
      ELSE
         HNISOR(:I1-1)=TEXT12(:I1-1)
      ENDIF
      IDEPL=0
      DO 270 JEL=1,NEL
      WRITE(TEXT12,'(3A4)') (ITNAM(II,JEL),II=1,3)
      I1=INDEX(TEXT12,'_')
      HITNAM=' '
      IF(I1.EQ.0) THEN
         HITNAM(:12)=TEXT12
       ELSE
         HITNAM(:I1-1)=TEXT12(:I1-1)
      ENDIF
      IF(HNISOR.EQ.HITNAM) THEN
         IDEPL=JEL
         GO TO 280
      ENDIF
  270 CONTINUE
      WRITE(HSMG,'(25HLIBEAR: MISSING ISOTOPE '',A12,5H''(1).)')
     > HNISOR
      CALL XABORT(HSMG)
  280 IKEEP(IDEPL)=1
      IF(BPAX(NEL+1,IDEPL).NE.0.0) KPAX(NEL+1,IDEPL)=1
      CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
      IF(INDIC.NE.3) CALL XABORT('LIBEAR: REACTION TYPE EXPECTED')
      IF(TEXT12.EQ.'FROM') THEN
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DBLINP)
  290    IF(INDIC.NE.3) CALL XABORT('LIBEAR: REACTION TYPE EXPECTED')
         DO 330 IREAC=1,MAXR
         RRAT=1.0
         IF(TEXT12.EQ.NMDEPL(IREAC)) THEN
            DO 320 JEL=1,NEL
            CALL REDGET(INDIC,ISOT,RRAT,TEXT12,DBLINP)
            IF(INDIC.NE.2) GO TO 290
            CALL REDGET(INDIC,ISOT,FLOTT,TEXT12,DBLINP)
            IF(INDIC.NE.3) CALL XABORT('LIBEAR: ISOTOPE NAME HNAMPAR '
     >      //'MISSING')
            I1=INDEX(TEXT12,'_')
            TEXT20=' '
            IF(I1.EQ.0) THEN
               TEXT20(:12)=TEXT12
            ELSE
               TEXT20(:I1-1)=TEXT12(:I1-1)
            ENDIF
            JDEPL=0
            DO 300 JREL=1,NEL
            WRITE(TEXT12,'(3A4)') (ITNAM(II,JREL),II=1,3)
            I1=INDEX(TEXT12,'_')
            HITNAM=' '
            IF(I1.EQ.0) THEN
               HITNAM(:12)=TEXT12
            ELSE
               HITNAM(:I1-1)=TEXT12(:I1-1)
            ENDIF
            IF(TEXT20.EQ.HITNAM) THEN
               JDEPL=JREL
               GO TO 310
            ENDIF
  300       CONTINUE
            WRITE(HSMG,'(25HLIBEAR: MISSING ISOTOPE '',A12,5H''(2).)')
     >      TEXT20
            CALL XABORT(HSMG)
  310       KPAX(IDEPL,JDEPL)=IREAC
            BPAX(IDEPL,JDEPL)=RRAT
  320       CONTINUE
            CALL XABORT('LIBEAR: TO MANY PARENT ISOTOPES')
         ENDIF
  330    CONTINUE
      ENDIF
  340 CONTINUE
      IF(INDIC.NE.3.OR.TEXT12.NE.'ENDCHAIN')
     >  CALL XABORT('LIBEAR: KEYWORD ENDCHAIN MISSING')
  350 DO 380 JEL=1,NEL
      IF(IKEEP(JEL).EQ.0) THEN
         DO 360 IREAC=1,NEL+MAXR
         KPAX(IREAC,JEL)=0
  360    CONTINUE
         DO 370 IEL=1,NEL
         KPAX(JEL,IEL)=0
  370    CONTINUE
      ENDIF
  380 CONTINUE
      DEALLOCATE(IKEEP)
      RETURN
*
  500 FORMAT(26HLIBEAR: MATERIAL/ISOTOPE ',A20,20H' IS MISSING ON APOL,
     > 15HIB-2 FILE NAME ,A12,1H.)
      END
