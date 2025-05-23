*DECK FMTGIS
      SUBROUTINE FMTGIS(IPRINT,NBISO,NAMISO,MISPRT,NAMRD,
     >                  NOPT,IOPT,ISOPRT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read and process isotopes to print.
*
*Copyright:
* Copyright (C) 2017 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IPRINT  print level.
* NBISO   number of isotopes on BURNUP.
* NAMISO  names of isotopes on BURNUP.
* NOPT    number of options.
* IOPT    processing option.
*
*Parameters: output
* MISPRT  number of isotopes to print.
* NAMRD   isotopes names to process.
* ISOPRT  isotopes print option.
*
*----------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          IPRINT,NBISO
      INTEGER          NAMISO(3,NBISO)
      INTEGER          MISPRT,NOPT,IOPT(NOPT)
      INTEGER          NAMRD(2,NBISO),ISOPRT(NBISO)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='FMTGIS')
*----
*  Variables for input via REDGET
*----
      INTEGER          ITYPLU,INTLIR
      CHARACTER        CARLIR*12
      REAL             REALIR
      DOUBLE PRECISION DBLLIR
*----
*  Local variables
*----
      INTEGER          ISOR,ISOT,II,KISPRT
*----
*  Get data from input file
*----
      IF(IPRINT .GE. 1) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      IF(IOPT(2).LT. 0) THEN
        IOPT(2)=-IOPT(2)
        ISOR=0
 100    CONTINUE
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .EQ. 10) GO TO 105
        IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >  ': Read error -- Character variable expected.')
        IF(CARLIR .EQ. ';') THEN
          GO TO 105
        ELSE
          ISOR=ISOR+1
          READ(CARLIR,'(2A4)') (NAMRD(II,ISOR),II=1,2)
        ENDIF
        GO TO 100
 105    CONTINUE
      ELSE
        ISOR=-1
      ENDIF
*----
*  All isotopes specified.
*  Set print flag
*----
      MISPRT=ISOR
      ISOPRT(:NBISO)=0
      IF(MISPRT .EQ. 0) THEN
        KISPRT=0
        DO ISOT=1,NBISO
          DO ISOR=1,KISPRT
            IF(NAMISO(1,ISOT) .EQ. NAMRD(1,ISOR) .AND.
     >         NAMISO(2,ISOT) .EQ. NAMRD(2,ISOR)) THEN
              ISOPRT(ISOT)=ISOR
              GO TO 115
            ENDIF
          ENDDO
          KISPRT=KISPRT+1
          NAMRD(1,KISPRT)=NAMISO(1,ISOT)
          NAMRD(2,KISPRT)=NAMISO(2,ISOT)
          ISOPRT(ISOT)=KISPRT
 115      CONTINUE
        ENDDO
        MISPRT=KISPRT
      ELSE
        IF(MISPRT.GT.0) THEN
          DO ISOT=1,NBISO
            DO ISOR=1,MISPRT
              IF(NAMISO(1,ISOT) .EQ. NAMRD(1,ISOR) .AND.
     >           NAMISO(2,ISOT) .EQ. NAMRD(2,ISOR)) THEN
                ISOPRT(ISOT)=ISOR
                GO TO 125
              ENDIF
            ENDDO
 125        CONTINUE
          ENDDO
        ENDIF
      ENDIF
      IF(IPRINT .GE. 1) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
*----
*  Processing finished, return
*----
      RETURN
*----
*  FORMATS
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
      END
