*DECK MCT
      SUBROUTINE MCT(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Monte-Carlo method based on NXT geometry analysis.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier and B. Arsenault
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1) creation type(L_MC);
*         HENTRY(2) read-only or modification type(L_TRACK);
*         HENTRY(3) read-only type(L_LIBRARY) or type(L_MACROLIB).
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file.
* JENTRY  access of each LCM object or file:
*         =0 the LCM object or file is created;
*         =1 the LCM object or file is open for modifications;
*         =2 the LCM object or file is open in read-only mode.
* KENTRY  LCM object address or file unit number.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPOUT,IPTRK,IPLIB
      INTEGER    NSTATE,IOUT,ITC,NSRCK,IKZ,KCT,NGRP,NFREG,NBMIX,NMIX,
     1           NFM,NL,NDEL,NED,ISEED,IPRINT,NBSCO,NMERGE,NGCOND
      PARAMETER (NSTATE=40,IOUT=6)
      INTEGER    ISTATE(NSTATE),GSTATE(NSTATE)
      DOUBLE PRECISION XYZL(2,3),KEFF,REKEFF
      CHARACTER  NAMREC*12,HSIGN*12
      LOGICAL    MODIF,LN2N
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INMIX,NAMEAD
      REAL, ALLOCATABLE, DIMENSION(:) :: XST,XSS,XSSNN,XNUFI,XCHI,XSN2N,
     < XSN3N,XSEDI
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.NE.3) CALL XABORT('MCT: THREE PARAMETERS EXPECTED.')
*     output table in creation or modification mode
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2))
     1  CALL XABORT('MCT: LCM OBJECT EXPECTED AT LHS (1).')
      IF((JENTRY(1).NE.0).AND.(JENTRY(1).NE.1)) 
     1  CALL XABORT('MCT: ENTRY IN CREATION '//
     2              'OR MODIFICATION MODE EXPECTED.')
      IPOUT=KENTRY(1)
      MODIF=(JENTRY(1).EQ.1)
*     tracking table in read-only mode
      IF((IENTRY(2).NE.1).AND.(IENTRY(2).NE.2))
     1  CALL XABORT('MCT: LCM OBJECT EXPECTED AT LHS (2).')
      IF((JENTRY(2).NE.1).AND.(JENTRY(2).NE.2)) CALL XABORT('MCT: ENTR'
     1  //'Y IN READ-ONLY OR MODIFICATION MODE EXPECTED.')
      IPTRK=KENTRY(2)
      CALL LCMGTC(IPTRK,'SIGNATURE',12,HSIGN)
      IF(HSIGN.NE.'L_TRACK') THEN
         NAMREC=HENTRY(2)
         CALL XABORT('MCT: INVALID SIGNATURE FOR '//NAMREC)
      ENDIF
*     xs library in read-only mode      
      IF((IENTRY(3).NE.1).AND.(IENTRY(3).NE.2))
     1  CALL XABORT('MCT: LCM OBJECT EXPECTED AT LHS (3).')
      IF(JENTRY(3).NE.2) 
     1  CALL XABORT('MCT: ENTRY IN READ-ONLY MODE EXPECTED (2).')
      IPLIB=KENTRY(3)
      CALL LCMGTC(IPLIB,'SIGNATURE',12,HSIGN)
      IF(HSIGN.EQ.'L_LIBRARY') THEN
         CALL LCMSIX(IPLIB,'MACROLIB',1)
         CALL LCMGTC(IPLIB,'SIGNATURE',12,HSIGN)
         IF(HSIGN.NE.'L_MACROLIB') THEN
            NAMREC=HENTRY(3)
            CALL XABORT('MCT: INVALID SIGNATURE FOR '//NAMREC)
         ENDIF
      ELSE IF(HSIGN.NE.'L_MACROLIB') THEN
         NAMREC=HENTRY(3)
         CALL XABORT('MCT: INVALID SIGNATURE FOR '//NAMREC)
      ENDIF
*----
*  INITIALIZE OUTPUT TABLE 
*----
      IF(MODIF) THEN
         CALL LCMGTC(IPOUT,'SIGNATURE',12,HSIGN)
         IF(HSIGN.NE.'L_MC') THEN
            NAMREC=HENTRY(1)
            CALL XABORT('MCT: INVALID SIGNATURE FOR '//NAMREC)
         ENDIF
         CALL LCMGET(IPOUT,'STATE-VECTOR',ISTATE)
      ELSE
         HSIGN='L_MC'
         CALL LCMPTC(IPOUT,'SIGNATURE',12,HSIGN)
         ISTATE(:NSTATE)=0
      ENDIF
*----
*  READ INPUT PARAMETERS
*----
      GSTATE(:NSTATE)=0
      CALL LCMGET(IPLIB,'STATE-VECTOR',GSTATE)
      NGRP = GSTATE(1)
      GSTATE(:NSTATE)=0
      CALL LCMGET(IPTRK,'STATE-VECTOR',GSTATE)
      NFREG= GSTATE(1)
      NBMIX= GSTATE(4)
      ALLOCATE(INMIX(NFREG))
      CALL LCMGET(IPTRK,'MATCOD',INMIX)
      CALL MCTGET(IPOUT,NGRP,NFREG,NBMIX,INMIX,IPRINT)
      DEALLOCATE(INMIX)
      CALL LCMGET(IPOUT,'STATE-VECTOR',ISTATE)
      ISEED=ISTATE(4)
      LN2N=(ISTATE(5).EQ.1)
*----
*  VALIDATE THE OPTIONS SPECIFIED FOR THE KCODE CARD
*----
      IF ((ISTATE(1).GT.0).AND.
     <    (ISTATE(2).GT.0).AND.
     <    (ISTATE(3).GT.0).AND.
     <    (ISTATE(2).LE.ISTATE(3))) THEN
         NSRCK=ISTATE(1)
         IKZ  =ISTATE(2)
         KCT  =ISTATE(3)
       ELSE
         CALL XABORT('MCT: INVALID PARAMETERS SPECIFIED FOR KCODE CARD')
      ENDIF
*----
*  RECOVER MACROSCOPIC CROSS SECTIONS.
*----
      GSTATE(:NSTATE)=0
      CALL LCMGET(IPLIB,'STATE-VECTOR',GSTATE)
      NGRP = GSTATE(1)
      NMIX = GSTATE(2)
      NL   = GSTATE(3)
      NFM  = GSTATE(4)
      NED  = GSTATE(5)
      NDEL = GSTATE(7)
      ALLOCATE(NAMEAD(2*NED))
      IF(NED.GT.0) CALL LCMGET(IPLIB,'ADDXSNAME-P0',NAMEAD)
      ALLOCATE(XST(NMIX*NGRP),XSS(NMIX*NGRP*NL),
     > XSSNN(NGRP*NGRP*NMIX*NL),XNUFI(NFM*NMIX*NGRP*(1+NDEL)),
     > XCHI(NFM*NMIX*NGRP*(1+NDEL)),XSN2N(NMIX*NGRP),XSN3N(NMIX*NGRP),
     > XSEDI(NMIX*NGRP*NED))
      CALL MCTLIB(IPLIB,NMIX,NGRP,NL,NFM,NDEL,NED,NAMEAD,LN2N,
     <            XST,XSS,XSSNN,XNUFI,XCHI,XSN2N,XSN3N,XSEDI)
*----
*  POWER ITERATION WITH THE MONTE-CARLO METHOD IN 1D/2D/3D CARTESIAN
*  GEOMETRY.
*----
      GSTATE(:NSTATE)=0
      CALL LCMGET(IPOUT,'STATE-VECTOR',GSTATE)
      NMERGE=GSTATE(7)
      NGCOND=GSTATE(8)
      NBSCO=5+NGCOND*NL+2*NFM*(1+NDEL)+NED
      CALL MCTFLX(IPTRK,IPOUT,IPRINT,NMIX,NGRP,NL,NFM,NDEL,NED,
     <            NAMEAD,XST,XSS,XSSNN,XNUFI,XCHI,XSN2N,XSN3N,
     <            XSEDI,NSRCK,IKZ,KCT,ISEED,XYZL,NBSCO,NMERGE,
     <            NGCOND,KEFF,REKEFF)
*
      DEALLOCATE(XSEDI,XSN3N,XSN2N,XCHI,XNUFI,XSSNN,XSS,XST,NAMEAD)
*----
*  RESET LIBRARY ON ROOT LEVEL
*----
      CALL LCMSIX(IPLIB,' ',0)
*----
*  SAVE KEFF INFORMATION ON MC OBJECT
*----
      CALL LCMPUT(IPOUT,'K-EFFECTIVE',1,2,REAL(KEFF))
      CALL LCMPUT(IPOUT,'K-EFFECTI-SD',1,2,REAL(REKEFF))
      IF(IPRINT.GT.0) WRITE(6,100) (ISTATE(ITC),ITC=1,9)
      RETURN
*
  100 FORMAT(/8H OPTIONS/8H -------/
     1 7H NSRCK ,I8,34H   (NUMBER OF HISTORIES PER CYCLE)/
     2 7H IKZ   ,I8,29H   (NUMBER OF CYCLES TO SKIP)/
     3 7H KCT   ,I8,27H   (TOTAL NUMBER OF CYCLES)/
     4 7H ISEED ,I8,45H   (INITIAL SEED FOR RANDOM NUMBER GENERATOR)/
     5 7H IN2N  ,I8,24H   (N2N PROCESSING FLAG)/
     6 7H ITALLY,I8,20H   (TYPE OF TALLIES)/
     7 7H NMERGE,I8,44H   (NUMBER OF HOMOGENIZED MIXTURES IN TALLY)/
     8 7H NGCOND,I8,40H   (NUMBER OF CONDENSED GROUPS IN TALLY)/
     9 7H NREG  ,I8,34H   (NUMBER OF REGIONS IN GEOMETRY))
      END
