*DECK KININI
      SUBROUTINE KININI(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Initialize the space-time kinetics parameters.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal.
*
*Author(s): D. Sekki
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1): create type(L_KINET);
*         HENTRY(2): read-only type(L_MACROLIB);
*         HENTRY(3): read-only type(L_TRACK);
*         HENTRY(4): read-only type(L_SYSTEM);
*         HENTRY(5): read-only type(L_FLUX).
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
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      CHARACTER HENTRY(NENTRY)*12
      TYPE(C_PTR) KENTRY(NENTRY)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      INTEGER ISTATE(NSTATE)
      CHARACTER TEXT12*12,HSIGN*12,CMODUL*12,HSMG*131
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.NE.5)CALL XABORT('@KININI: INVALID NUMBER OF MODULE PA'
     1 //'RAMETERS.')
      DO IEN=2,NENTRY
       IF((IENTRY(IEN).NE.1).AND.(IENTRY(IEN).NE.2))
     1   CALL XABORT('@KININI: LCM OBJECTS EXPECTED AT RHS')
       IF(JENTRY(IEN).NE.2)CALL XABORT('@KININI: LCM OBJEC'
     1  //'TS IN READ-ONLY MODE EXPECTED AT RHS.')
      ENDDO
*     L_KINET
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2))
     1 CALL XABORT('@KININI: LCM OBJECT EXPECTED AT LHS.')
      IF(JENTRY(1).NE.0)CALL XABORT('@KININI: L_KINET IN'
     1 //' CREATE MODE EXPECTED.')
      HSIGN='L_KINET'
      CALL LCMPTC(KENTRY(1),'SIGNATURE',12,HSIGN)
*     L_MACROLIB
      CALL LCMGTC(KENTRY(2),'SIGNATURE',12,HSIGN)
      IF(HSIGN.NE.'L_MACROLIB')THEN
        TEXT12=HENTRY(2)
        CALL XABORT('@KININI: SIGNATURE OF '//TEXT12//' IS '
     1  //HSIGN//'. L_MACROLIB EXPECTED.')
      ENDIF
*     L_TRACK
      CALL LCMGTC(KENTRY(3),'SIGNATURE',12,HSIGN)
      IF(HSIGN.NE.'L_TRACK')THEN
        TEXT12=HENTRY(3)
        CALL XABORT('@KININI: SIGNATURE OF '//TEXT12//' IS '
     1  //HSIGN//'. L_TRACK EXPECTED.')
      ENDIF
*     L_SYSTEM
      CALL LCMGTC(KENTRY(4),'SIGNATURE',12,HSIGN)
      IF(HSIGN.NE.'L_SYSTEM')THEN
        TEXT12=HENTRY(4)
        CALL XABORT('@KININI: SIGNATURE OF '//TEXT12//' IS '
     1  //HSIGN//'. L_SYSTEM EXPECTED.')
      ENDIF
      CALL LCMGTC(KENTRY(4),'LINK.MACRO',12,TEXT12)
      IF(HENTRY(2).NE.TEXT12) THEN
        WRITE(HSMG,'(40H@KININI: INVALID MACROLIB OBJECT NAME ='',A12,
     1  18H'', EXPECTED NAME='',A12,2H''.)') HENTRY(2),TEXT12
        CALL XABORT(HSMG)
      ENDIF
      CALL LCMGTC(KENTRY(4),'LINK.TRACK',12,TEXT12)
      IF(HENTRY(3).NE.TEXT12) THEN
        WRITE(HSMG,'(40H@KININI: INVALID TRACKING OBJECT NAME ='',A12,
     1  18H'', EXPECTED NAME='',A12,2H''.)') HENTRY(3),TEXT12
        CALL XABORT(HSMG)
      ENDIF
*     L_FLUX
      CALL LCMGTC(KENTRY(5),'SIGNATURE',12,HSIGN)
      IF(HSIGN.NE.'L_FLUX')THEN
        TEXT12=HENTRY(5)
        CALL XABORT('@KININI: SIGNATURE OF '//TEXT12//' IS '
     1  //HSIGN//'. L_FLUX EXPECTED.')
      ENDIF
*----
*  OBJECTS VALIDATION
*----
      ISTATE(:NSTATE)=0
      CALL LCMGET(KENTRY(5),'STATE-VECTOR',ISTATE)
      NGR=ISTATE(1)
      NUN=ISTATE(2)
      ISTATE(:NSTATE)=0
      CALL LCMGET(KENTRY(2),'STATE-VECTOR',ISTATE)
      IF(ISTATE(1).NE.NGR)CALL XABORT('@KININI: INVALID NU'
     1 //'MBER OF ENERGY GROUPS IN L_MACROLIB OR IN L_FLUX.')
      NBM=ISTATE(2)
      NBFIS=ISTATE(4)
      NDG=ISTATE(7)
      ISTATE(:NSTATE)=0
      CALL LCMGET(KENTRY(3),'STATE-VECTOR',ISTATE)
      IF(ISTATE(2).NE.NUN)CALL XABORT('@KININI: INVALID TOTAL'
     1 //' NUMBER OF UNKNOWNS IN L_FLUX OR IN L_TRACK.')
      IF(ISTATE(4).GT.NBM) THEN
         WRITE(HSMG,'(46H@KININI: THE NUMBER OF MIXTURES IN THE TRACKIN,
     1   3HG (,I5,50H) IS GREATER THAN THE NUMBER OF MIXTURES IN THE MA,
     2   8HCROLIB (,I5,2H).)') ISTATE(4),NBM
         CALL XABORT(HSMG)
      ENDIF
      ITYPE=ISTATE(6)
      CALL LCMGTC(KENTRY(3),'TRACK-TYPE',12,CMODUL)
*
      IF(CMODUL.EQ.'BIVAC')THEN
        IF((ITYPE.NE.1).AND.(ITYPE.NE.2).AND.(ITYPE.NE.3).AND.
     1     (ITYPE.NE.4).AND.(ITYPE.NE.5).AND.(ITYPE.NE.6).AND.
     2     (ITYPE.NE.8))CALL XABORT('@KININI: TYPE OF GEOMETR'
     3     //'Y NOT COMPATIBLE WITH BIVAC TRACKING-TYPE.')
      ELSEIF(CMODUL.EQ.'TRIVAC')THEN
        IF((ITYPE.NE.1).AND.(ITYPE.NE.2).AND.(ITYPE.NE.3).AND.
     1     (ITYPE.NE.5).AND.(ITYPE.NE.6).AND.(ITYPE.NE.7).AND.
     2     (ITYPE.NE.8).AND.(ITYPE.NE.9))CALL XABORT('@KININI'
     3     //': TYPE OF GEOMETRY NOT COMPATIBLE WITH TRIVAC T'
     4     //'RACKING-TYPE.')
      ENDIF
      NEL=ISTATE(1)
      CALL LCMPTC(KENTRY(1),'TRACK-TYPE',12,CMODUL)
      CALL KINRD1(NENTRY,KENTRY,CMODUL,NGR,NBM,NBFIS,NEL,NUN,NDG)
      RETURN
      END
