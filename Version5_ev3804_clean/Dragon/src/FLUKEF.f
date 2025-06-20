*DECK FLUKEF
      SUBROUTINE FLUKEF(IPRT,IPMACR,NGRP,NREG,NUNKNO,NMAT,NIFIS,NANIS,
     1 MATCOD,VOL,KEYFLX,XSTRC,XSDIA,XSNUF,XSCHI,NMERG,IMERG,DIFHET,
     2 FLUX,B2,ILEAK,LEAKSW,OLDBIL,AKEFF,AFLNOR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the effective multiplication factor.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* IPRT    print flag.
* IPMACR  pointer to the macrolib LCM object.
* NGRP    number of energy groups.
* NREG    number of regions.
* NUNKNO  number of unknowns per energy group including spherical
*         harmonic terms, interface currents and fundamental
*         currents.
* NMAT    number of mixtures.
* NIFIS   number of fissile isotopes.
* NANIS   maximum cross section Legendre order.
* MATCOD  mixture indices.
* VOL     volumes.
* KEYFLX  index of region flux components in unknown vector.
* XSTRC   transport-corrected macroscopic total cross sections.
* XSDIA   transport-corrected macroscopic within-group scattering cross
*         sections.
* XSNUF   nu*macroscopic fission cross sections.
* XSCHI   fission spectrum.
* NMERG   number of leakage zones.
* IMERG   leakage zone index in each material mixture zone.
* DIFHET  heterogeneous leakage coefficients.
* FLUX    neutron flux.
* B2      directionnal bucklings.
* ILEAK   method used to include DB2 effect:
*         <5 uniform DB2 model;
*         =5 Todorova-type isotropic streaming model;
*         =6 Ecco-type isotropic streaming model;
*         >6 Tibere anisotropic streaming model.
* LEAKSW  leakage flag (=.true. if leakage is present on the outer
*         surface).
*
*Parameters: input/output
* OLDBIL  previous norm of the flux on input and 
*         new norm of the flux at output.
*
*Parameters: output
* AKEFF   effective multiplication factor.
* AFLNOR  flux normatization factor.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMACR
      INTEGER IPRT,NGRP,NREG,NUNKNO,NMAT,NIFIS,NANIS,MATCOD(NREG),
     1 KEYFLX(NREG),NMERG,IMERG(NMAT),ILEAK
      REAL VOL(NREG),XSTRC(0:NMAT,NGRP),XSDIA(0:NMAT,0:NANIS,NGRP),
     1 XSNUF(0:NMAT,NIFIS,NGRP),XSCHI(0:NMAT,NIFIS,NGRP),
     2 DIFHET(NMERG,NGRP),FLUX(NUNKNO,NGRP),B2(4)
      DOUBLE PRECISION OLDBIL,AKEFF,AFLNOR
      LOGICAL LEAKSW
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMACR,KPMACR
      DOUBLE PRECISION LOSS,SUMCHI,PROD,FISONE,PHIC,FISOUR,AKINV,
     1 DZERO,DONE
      PARAMETER (DZERO=0.0D0, DONE=1.0D0)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: XSCAT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NJJ(0:NMAT),IJJ(0:NMAT),IPOS(0:NMAT))
      ALLOCATE(XSCAT(0:NMAT*NGRP))
*
      FISOUR=DZERO
      SUMCHI=DZERO
      LOSS=DZERO
      AKINV=DZERO
      JPMACR=LCMGID(IPMACR,'GROUP')
      DO 40 IGRP=1,NGRP
      KPMACR=LCMGIL(JPMACR,IGRP)
      DO 15 IREG=1,NREG
      IND=KEYFLX(IREG)
      IF(IND.EQ.0) GO TO 15
      IBM=MATCOD(IREG)
      IF(IBM.EQ.0) GO TO 15
      PHIC=FLUX(IND,IGRP)*VOL(IREG)
      LOSS=LOSS+XSTRC(IBM,IGRP)*PHIC
      IF((ILEAK.GE.1).AND.(ILEAK.LE.5)) THEN
         INM=IMERG(IBM)
         IF(INM.GT.0) LOSS=LOSS+B2(4)*DIFHET(INM,IGRP)*PHIC
      ELSE IF(ILEAK.EQ.6) THEN
         LOSS=LOSS+B2(4)*FLUX(NUNKNO/2+IND,IGRP)*VOL(IREG)
      ELSE IF(ILEAK.GE.7) THEN
         LOSS=LOSS+B2(1)*FLUX(NUNKNO/4+IND,IGRP)*VOL(IREG)
     1            +B2(2)*FLUX(NUNKNO/2+IND,IGRP)*VOL(IREG)
     2            +B2(3)*FLUX(3*NUNKNO/4+IND,IGRP)*VOL(IREG)
      ENDIF
      DO 10 IFIS=1,NIFIS
      FISOUR=FISOUR+XSNUF(IBM,IFIS,IGRP)*PHIC
      AKINV=AKINV+XSNUF(IBM,IFIS,IGRP)*PHIC
   10 CONTINUE
   15 CONTINUE
*
      CALL LCMGET(KPMACR,'NJJS00',NJJ(1))
      CALL LCMGET(KPMACR,'IJJS00',IJJ(1))
      CALL LCMGET(KPMACR,'IPOS00',IPOS(1))
      CALL LCMGET(KPMACR,'SCAT00',XSCAT(1))
      DO 30 IREG=1,NREG
      IBM=MATCOD(IREG)
      IF(IBM.GT.0) THEN
         IND=KEYFLX(IREG)
         JGRP=IJJ(IBM)
         DO 20 JND=1,NJJ(IBM)
         IF(JGRP.EQ.IGRP) THEN
           LOSS=LOSS-XSDIA(IBM,0,IGRP)*FLUX(IND,IGRP)*VOL(IREG)
         ELSE
           LOSS=LOSS-XSCAT(IPOS(IBM)+JND-1)*FLUX(IND,JGRP)*VOL(IREG)
         ENDIF
         JGRP=JGRP-1
   20    CONTINUE
      ENDIF
   30 CONTINUE
   40 CONTINUE
*
      IF(AKINV.NE.0.0) THEN
         AKINV=DONE/AKINV
         DO 70 IREG=1,NREG
         IND=KEYFLX(IREG)
         IF(IND.EQ.0) GO TO 70
         IBM=MATCOD(IREG)
         DO 65 IFIS=1,NIFIS
         FISONE=DZERO
         DO 50 IGRP=1,NGRP
         FISONE=FISONE+XSNUF(IBM,IFIS,IGRP)*FLUX(IND,IGRP)
   50    CONTINUE
         DO 60 IGRP=1,NGRP
         SUMCHI=SUMCHI+AKINV*XSCHI(IBM,IFIS,IGRP)*FISONE*VOL(IREG)
   60    CONTINUE
   65    CONTINUE
   70    CONTINUE
      ENDIF
*
      PROD=SUMCHI*FISOUR
      IF(PROD.GT.DZERO) THEN
         AFLNOR=DONE/PROD
      ELSE
         AFLNOR=DONE
      ENDIF
      IF(LEAKSW) THEN
         IF(OLDBIL.GT.DZERO) THEN
            AKEFF=PROD/OLDBIL
         ELSE
            AKEFF=PROD
         ENDIF
      ELSE
         IF(LOSS.GT.DZERO) THEN
            AKEFF=PROD/LOSS
            AFLNOR=AFLNOR*LOSS
         ELSE
            AKEFF=PROD
         ENDIF
      ENDIF
      OLDBIL=PROD*AFLNOR
      IF(IPRT.GT.2) THEN
         WRITE(6,*)
         WRITE(6,*) ' ************ OLDBIL=',OLDBIL
         WRITE(6,*) ' ************ PROD  =',PROD
         WRITE(6,*) ' ************ LOSS  =',LOSS
         WRITE(6,*) ' ************ AFLNOR=',AFLNOR
         WRITE(6,*) ' ************ AKEFF =',AKEFF
         WRITE(6,*)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(XSCAT)
      DEALLOCATE(IPOS,IJJ,NJJ)
      RETURN
      END
