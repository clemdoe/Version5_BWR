*DECK LIBFQD
      SUBROUTINE LIBFQD(MAXNOR,LPART,MAXTRA,HNAMIS,IPLIB,NGRO,NL,NED,
     1 NDEL,NDIL,IGRMIN,IGRMAX,LBIN,NFS,IMPX,LSCAT,LSIGF,LADD,DILUT,
     2 FLUX,TOTAL,SIGF,SIGS,SCAT,SADD,ZDEL,EBIN,SIGTF,SIGSF,SIGFF,
     3 AWR,ISMIN,ISMAX,GOLD,IPRECI,NOR,LBSIGF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute NOR-point Calendf-type probability tables;
* compute physical and/or slowing-down correlated probability tables;
* component of the Ribon extended method.
*
*Copyright:
* Copyright (C) 2003 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* MAXNOR  first dimension of matrix PRTSIG. Equal to the maximum order
*         of a probability table.
* LPART   maximum scattering bandwidth for the isotope.
* MAXTRA  maximum number of energy bins of size DELI.
* HNAMIS  name of the isotope.
* IPLIB   pointer to the internal library (L_LIBRARY signature).
* NGRO    number of energy groups.
* NL      number of Legendre orders required in the calculation
*         (NL=1 or higher).
* NED     number of extra vector edits.
* NDEL    number of delayed neutron precursor groups.
* NDIL    number of finite dilutions.
* IGRMIN  first self-shielded group with BIN data.
* IGRMAX  most thermal self-shielded group with BIN data.
* LBIN    number of fine (bin) energy groups.
* NFS     number of fine energy groups in each group.
* IMPX    print flag (equal to zero for no print).
* LSCAT   anisotropy flag (=.true. if a given Legendre order of the
*         scattering cross section exists).
* LSIGF   fission flag (=.true. if the isotope can fission).
* LADD    additional cross section flag (=.true. if a given additional
*         cross section exists).
* DILUT   dilutions.
* FLUX    weighting flux.
* TOTAL   total cross sections.
* SIGF    nu*fission cross sections.
* SIGS    scattering cross sections.
* SCAT    scattering transfer matrices (sec,prim,Legendre,dilution).
* SADD    additional cross sections.
* ZDEL    delayed nu-sigf cross sections.
* EBIN    fine group energy limits in EV.
* SIGTF   microscopic total x-sections in the fine groups.
* SIGSF   microscopic P0 scattering x-sections in the fine groups.
* AWR     mass ratio for current isotope.
* SIGFF   microscopic nu*fission cross sections in the fine groups.
* ISMIN   minimum secondary group corresponding to each primary group.
* ISMAX   maximum secondary group corresponding to each primary group.
* GOLD    method flag: =-998.0 to use the CALENDF approach; =-999.0 to
*         use the Ribon extended approach; =1.0 to use the ST model.
* IPRECI  accuracy index for probability tables in CALENDF.
* LBSIGF  autolib (bin) fission data flag.
*
*Parameters: output
* NOR     number of subgroups in each group.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB
      INTEGER MAXNOR,LPART,MAXTRA,NGRO,NL,NED,NDEL,NDIL,IGRMIN,IGRMAX,
     1 LBIN,NFS(NGRO),IMPX,ISMIN(NL,NGRO),ISMAX(NL,NGRO),IPRECI,
     2 NOR(NGRO)
      REAL DILUT(NDIL+1),FLUX(NGRO,NDIL+1),TOTAL(NGRO,NDIL+1),
     1 SIGF(NGRO,NDIL+1),SIGS(NGRO,NL,NDIL+1),SCAT(NGRO,NGRO,NL,NDIL+1),
     2 SADD(NGRO,NED,NDIL+1),ZDEL(NGRO,NDEL,NDIL+1),EBIN(LBIN+1),
     3 SIGTF(LBIN),SIGSF(LBIN),SIGFF(LBIN),AWR,GOLD(NGRO)
      LOGICAL LSCAT(NL),LSIGF,LADD(NED),LBSIGF
      CHARACTER HNAMIS*12
*----
*  LOCAL VARIABLES
*----
      PARAMETER (IALTER=0,MAXDIL=65)
      TYPE(C_PTR) JPLIB,KPLIB
      CHARACTER HSMG*131,TAG*2
      LOGICAL LNORAJ,LPHYS,LCALEN,LRIBON,LDIL(MAXDIL+1),LPTMC
      INTEGER IPERD(MAXDIL+1)
      REAL XSREF(MAXDIL),TEST(8),DILUT2(MAXDIL+1),TSCAT(20,MAXDIL),
     1     DIFFS(20,MAXDIL)
      DOUBLE PRECISION SIGTI2,SIGAI2,SIGTIN,SIGAIN,DELMAC,T,TF,T0,T1,
     1 T2,ACCUM1,ACCUM2,ACCUM3,ACCUM4,FACT(MAXDIL),BB(MAXDIL+1)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISM
      REAL, ALLOCATABLE, DIMENSION(:) :: WSLD,DELTA,UUU,STIS,SIGAF,
     1 PRTSIW,PRTABS,GAR,SEFR
      REAL, ALLOCATABLE, DIMENSION(:,:) :: PRTSIG,PRI,PRTSIG1,PRTSIG2,
     1 SCAT00,PRTRS
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: PRTPH
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: PHIMT,CC,MOMT,MOMP
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: MATRIX,WORK,RSTAR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: PHI
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(ISM(2,NL))
      ALLOCATE(PRTSIG(MAXNOR,3+NL+NL*LPART+NED+NDEL),WSLD(MAXNOR**2),
     1 DELTA(LBIN),UUU(LBIN),STIS(LBIN),SIGAF(LBIN),
     2 PRI(MAXTRA,NL),PRTPH(MAXNOR,NDIL,NL),PRTSIW(MAXNOR),
     3 PRTABS(MAXNOR),PRTRS(MAXNOR,NDIL+1),
     4 PRTSIG1(MAXNOR,3+NL+NL*LPART+NED+NDEL),SCAT00(LBIN,NGRO),
     5 PRTSIG2(MAXNOR,3+NL+NL*LPART+NED+NDEL),GAR(LBIN))
      ALLOCATE(PHIMT(MAXNOR),MATRIX(MAXNOR,MAXNOR+1),
     1 PHI(LBIN,NDIL,NL),WORK(NDIL+1,MAXNOR),CC(MAXNOR),
     2 RSTAR(LBIN,NDIL+1))
*
      IF(NDIL.GT.MAXDIL) CALL XABORT('LIBFQD: MAXDIL OVERFLOW.')
*----
*  NORMALIZE THE BIN-TYPE DATA AND COMPUTE DELTA AND SIGAF.
*----
      IBIN=0
      DELMIN=1.0E10
      DO 40 IGRP=IGRMIN,IGRMAX
      SIGTIN=0.0D0
      SIGAIN=0.0D0
      SIGSIN=0.0D0
      SIGFIN=0.0D0
      SIGTI2=0.0D0
      SIGAI2=0.0D0
      SIGSI2=0.0D0
      SIGFI2=0.0D0
      DO 20 IGF=1,NFS(IGRP)
      DELM=LOG(EBIN(IBIN+IGF)/EBIN(IBIN+IGF+1))
      DELMIN=MIN(DELMIN,DELM)
      SIGTIN=SIGTIN+SIGTF(IBIN+IGF)*DELM
      SIGAIN=SIGAIN+(SIGTF(IBIN+IGF)-SIGSF(IBIN+IGF))*DELM
      SIGSIN=SIGSIN+SIGSF(IBIN+IGF)*DELM
      IF(LBSIGF) SIGFIN=SIGFIN+SIGFF(IBIN+IGF)*DELM
      SIGTF(IBIN+IGF)=MAX(0.002,SIGTF(IBIN+IGF))
      SIGAF(IBIN+IGF)=SIGTF(IBIN+IGF)-SIGSF(IBIN+IGF)
      SIGTI2=SIGTI2+SIGTF(IBIN+IGF)*DELM
      SIGAI2=SIGAI2+SIGAF(IBIN+IGF)*DELM
      SIGSI2=SIGSI2+SIGSF(IBIN+IGF)*DELM
      IF(LBSIGF) SIGFI2=SIGFI2+SIGFF(IBIN+IGF)*DELM
      UUU(IBIN+IGF)=LOG(EBIN(1)/EBIN(IBIN+IGF+1))
      DELTA(IBIN+IGF)=DELM
   20 CONTINUE
      DO 30 IGF=1,NFS(IGRP)
      SIGTF(IBIN+IGF)=SIGTF(IBIN+IGF)*REAL(SIGTIN/SIGTI2)
      SIGSF(IBIN+IGF)=SIGSF(IBIN+IGF)*REAL(SIGSIN/SIGSI2)
      IF(LBSIGF) SIGFF(IBIN+IGF)=SIGFF(IBIN+IGF)*(SIGFIN/SIGFI2)
      SIGAF(IBIN+IGF)=SIGAF(IBIN+IGF)*REAL(SIGAIN/SIGAI2)
   30 CONTINUE
      IBIN=IBIN+NFS(IGRP)
   40 CONTINUE
*----
*  ASSUME THAT THE ELEMENTARY LETHARGY WIDTH DELI IS A RATIONAL FRACTION
*  OF THE LETHARGY UNIT. CHECK THIS ASSUMPTION.
*----
      CALL LCMLEN(IPLIB,'BIN-DELI',LENGT,ITYLCM)
      IF((LENGT.EQ.1).AND.(ITYLCM.EQ.2)) THEN
        CALL LCMGET(IPLIB,'BIN-DELI',DELI)
      ELSE
        DELI=1.0/REAL(INT(1.00001/DELMIN))
      ENDIF
      IBIN=0
      ERR=0.0
      DO 60 IGRP=IGRMIN,IGRMAX
      DO 50 IGF=1,NFS(IGRP)
      LARGH=INT(DELTA(IBIN+IGF)/DELI+0.1)
      ERR=MAX(ERR,ABS(DELTA(IBIN+IGF)/DELI-REAL(LARGH)))
   50 CONTINUE
      IBIN=IBIN+NFS(IGRP)
   60 CONTINUE
      IF((IMPX.GT.0).OR.(ERR.GT.0.05)) THEN
         WRITE(6,'(/47H LIBFQD: THE ELEMENTARY LETHARGY WIDTH OF ISOTO,
     1   4HPE '',A12,11H'' IS SET TO,1P,E11.4,6H. ERR=,E10.3)') HNAMIS,
     2   DELI,ERR
      ENDIF
      IF(ERR.GT.0.05) THEN
         WRITE(HSMG,'(45HLIBFQD: UNABLE TO SET THE ELEMENTARY LETHARGY,
     1   20H WIDTH FOR ISOTOPE '',A12,2H''.)') HNAMIS
         WRITE(6,'(A)') HSMG
      ENDIF
*----
*  COMPUTE THE PRI ARRAY FOR VARIOUS LEGENDRE ORDERS.
*----
      DO 70 IL=1,NL
      CALL LIBPRI(MAXTRA,DELI,AWR,IALTER,IL-1,NEXT,PRI(1,IL))
   70 CONTINUE
*----
*  COMPUTE AUTOLIB CROSS SECTIONS FOR THE PO SCATTERING MATRIX 
*----
      SCAT00(:LBIN,:NGRO)=0.0
      LLL=0
      DO IGRP=IGRMIN,IGRMAX
         GAR0=0.0
         DO LI=1,NFS(IGRP)
            LLL=LLL+1
            GAR0=GAR0+DELTA(LLL)
            GAR(:LBIN)=0.0
            CALL LIBECT(MAXTRA,LLL,PRI(1,1),UUU,DELI,DELTA,NEXT,1,
     1      MML,STIS)
            DO I=1,MML
               LLJ=LLL-I+1
               GAR(LLJ)=STIS(I)*SIGSF(LLJ)*DELTA(LLJ)/DELTA(LLL)
            ENDDO
            LLJ=0
            DO JGRP=IGRMIN,IGRMAX
               DO LJ=1,NFS(JGRP)
                  LLJ=LLJ+1
                  SCAT00(LLJ,IGRP)=SCAT00(LLJ,IGRP)+GAR(LLJ)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
*----
*  MAIN LOOP OVER THE COARSE ENERGY GROUPS.
*----
      NOR(:NGRO)=0
      CALL LCMSIX(IPLIB,'PT-TABLE',1)
      CALL LCMPUT(IPLIB,'NDEL',1,1,NDEL)
      IBIN=0
      JPLIB=LCMLID(IPLIB,'GROUP-PT',NGRO)
*     ------------------
      DO 810 IGRP=1,NGRO
*----
*  REMOVE BADLY BEHAVED COLLOCATIONS POINTS.
*----
      MDIL=NDIL
      LDIL(:NDIL+1)=.TRUE.
      DO 90 IDIL=NDIL,1,-1
      IF(DILUT(IDIL).LT.1.0) THEN
         MDIL=MDIL-1
         LDIL(IDIL)=.FALSE.
      ELSE IF((DILUT(IDIL).GT.1.0E5).AND.(DILUT(IDIL).LT.1.0E10)) THEN
         MDIL=MDIL-1
         LDIL(IDIL)=.FALSE.
      ELSE IF(TOTAL(IGRP,IDIL).LE.0.0) THEN
         MDIL=MDIL-1
         LDIL(IDIL)=.FALSE.
      ENDIF
   90 CONTINUE
      IDD=0
      DO 100 IDIL=1,NDIL+1
      IF(LDIL(IDIL)) THEN
         IDD=IDD+1
         DILUT2(IDD)=DILUT(IDIL)
         IPERD(IDD)=IDIL
      ENDIF
  100 CONTINUE
      IF(IDD.NE.MDIL+1) CALL XABORT('LIBFQD: INTERNAL ERROR.')
*
      LCALEN=(NFS(IGRP).GT.0).AND.(GOLD(IGRP).EQ.-998.0)
      LRIBON=(NFS(IGRP).GT.0).AND.(GOLD(IGRP).EQ.-999.0)
      LPTMC=(NFS(IGRP).GT.0).AND.(GOLD(IGRP).EQ.-1000.0)
*----
*  ACTIVE SPM IN GROUPS IGRMAX-1 and IGRMAX
*----
      IF(LPTMC.AND.(IGRP.GE.IGRMAX-1)) THEN
         LPTMC=.FALSE.
         LCALEN=.TRUE.
      ENDIF
      LPHYS=(.NOT.LCALEN).AND.(.NOT.LRIBON).AND.(.NOT.LPTMC)
*      
      PRTSIG(:MAXNOR,:3+NL+NL*LPART+NED+NDEL)=0.0
*
      IF(IMPX.GT.1) THEN
         WRITE(6,'(/25H LIBFQD: PROCESSING GROUP,I4,14H FOR ISOTOPE '',
     1   A12,2H''.)') IGRP,HNAMIS
      ENDIF
      DO 110 IDIL=1,MDIL
      JDIL=IPERD(IDIL)
      XSREF(IDIL)=TOTAL(IGRP,JDIL)-SIGS(IGRP,1,JDIL)
  110 CONTINUE
*----
*  COMPUTE THE RESONANT FLUX BY SOLVING A SLOWING-DOWN EQUATION. COMPUTE
*  STIS USING LIBECT AND COMPUTE THE FINE FLUX. NORMALIZE THE RESONANT
*  FLUX TO THE DILUTION-DEPENDENT NJOY COLLISION RATES.
*----
      IF(LRIBON) THEN
         RSTAR(:LBIN,:NDIL+1)=0.0D0
         DO 142 IDIL=1,MDIL
         T0=0.0D0
         DELMAC=0.0D0
         DO 130 IGF=1,NFS(IGRP)
         DELM=DELTA(IBIN+IGF)
         DELMAC=DELMAC+DELM
         CALL LIBECT(MAXTRA,IBIN+IGF,PRI(1,1),UUU,DELI,DELTA,NEXT,1,
     1   MML,STIS)
         PHI(IBIN+IGF,IDIL,1)=DILUT2(IDIL)
         DO 120 J=2,MML
         JGF=IBIN+IGF-J+1
         PHI(IBIN+IGF,IDIL,1)=PHI(IBIN+IGF,IDIL,1)+DBLE(STIS(J)*
     1   (SIGTF(JGF)-SIGAF(JGF))*DELTA(JGF)/DELTA(IBIN+IGF))*
     2   PHI(JGF,IDIL,1)
  120    CONTINUE
         PHI(IBIN+IGF,IDIL,1)=PHI(IBIN+IGF,IDIL,1)/(SIGTF(IBIN+IGF)+
     1   DILUT2(IDIL)-DBLE(STIS(1)*(SIGTF(IBIN+IGF)-SIGAF(IBIN+IGF))))
         T0=T0+PHI(IBIN+IGF,IDIL,1)*SIGTF(IBIN+IGF)*DELM
  130    CONTINUE
*
         JDIL=IPERD(IDIL)
         FACT(IDIL)=FLUX(IGRP,JDIL)*TOTAL(IGRP,JDIL)*DELMAC/T0
         DO 141 IL=2,NL
         DO 140 IGF=1,NFS(IGRP)
         BONDAR=SIGTF(IBIN+IGF)+DILUT2(IDIL)
         PHI(IBIN+IGF,IDIL,IL)=PHI(IBIN+IGF,IDIL,IL-1)/BONDAR
  140    CONTINUE
  141    CONTINUE
  142    CONTINUE
*
*        COMPUTE THE FINE-GROUP SLOWING-DOWN SOURCE.  
         DO 152 IGF=1,NFS(IGRP)
         CALL LIBECT(MAXTRA,IBIN+IGF,PRI(1,1),UUU,DELI,DELTA,NEXT,1,
     1   MML,STIS)
         DO 151 J=1,MML
         JGF=IBIN+IGF-J+1
         ACCUM1=DBLE(STIS(J)*(SIGTF(JGF)-SIGAF(JGF))*DELTA(JGF)/
     1   DELTA(IBIN+IGF))
         DO 150 IDIL=1,MDIL+1
         IF(IDIL.LE.MDIL) THEN
            RSTAR(IBIN+IGF,IDIL)=RSTAR(IBIN+IGF,IDIL)+ACCUM1*
     1      PHI(JGF,IDIL,1)
         ELSE
            RSTAR(IBIN+IGF,IDIL)=RSTAR(IBIN+IGF,IDIL)+ACCUM1
         ENDIF
  150    CONTINUE
  151    CONTINUE
  152    CONTINUE
         DO 162 IDIL=1,MDIL
         DO 161 IGF=1,NFS(IGRP)
         RSTAR(IBIN+IGF,IDIL)=RSTAR(IBIN+IGF,IDIL)*FACT(IDIL)
         DO 160 IL=1,NL
         PHI(IBIN+IGF,IDIL,IL)=PHI(IBIN+IGF,IDIL,IL)*FACT(IDIL)
  160    CONTINUE
  161    CONTINUE
  162    CONTINUE
      ENDIF
*----
*  TEST FINE FLUX.
*----
      IF((IMPX.GT.5).AND.LRIBON) THEN
         WRITE(6,910) IGRP,HNAMIS
         DO 240 IDIL=1,MDIL
         DELMAC=0.0D0
         TF=0.0D0
         T0=0.0D0
         T1=0.0D0
         T2=0.0D0
         DO 230 IGF=1,NFS(IGRP)
         DELM=DELTA(IBIN+IGF)
         DELMAC=DELMAC+DELM
         TF=TF+PHI(IBIN+IGF,IDIL,1)*DELM
         T0=T0+PHI(IBIN+IGF,IDIL,1)*SIGTF(IBIN+IGF)*DELM
         T1=T1+PHI(IBIN+IGF,IDIL,1)*SIGAF(IBIN+IGF)*DELM
         T2=T2+RSTAR(IBIN+IGF,IDIL)*DELM
  230    CONTINUE
         JDIL=IPERD(IDIL)
         BTOT=TOTAL(IGRP,JDIL)*FLUX(IGRP,JDIL)
         BABS=(TOTAL(IGRP,JDIL)-SIGS(IGRP,1,JDIL))*FLUX(IGRP,JDIL)
         WRITE(6,'(1X,I5,1P,8E12.4)') IDIL,T0/DELMAC,BTOT,T1/DELMAC,
     1   BABS,T2/DELMAC,(T0+DILUT2(IDIL)*TF)/DELMAC-DILUT2(IDIL),
     2   TF/DELMAC,((T2/DELMAC)+DILUT2(IDIL))/((T0/TF)+DILUT2(IDIL))
  240    CONTINUE
         DELMAC=0.0D0
         T0=0.0D0
         T1=0.0D0
         DO 250 IGF=1,NFS(IGRP)
         DELM=DELTA(IBIN+IGF)
         DELMAC=DELMAC+DELM
         T0=T0+SIGTF(IBIN+IGF)*DELM
         T1=T1+SIGAF(IBIN+IGF)*DELM
  250    CONTINUE
         BTOT=TOTAL(IGRP,NDIL+1)
         BABS=TOTAL(IGRP,NDIL+1)-SIGS(IGRP,1,NDIL+1)
         WRITE(6,'(3X,3HINF,1P,4E12.4)') T0/DELMAC,BTOT,T1/DELMAC,BABS
      ENDIF
*----
*  PROCESS CLASSICAL PROBABILITY TABLE INFORMATION IN TOTAL XS.
*----
      LNORAJ=.TRUE.
      ERROR1=0.0
  260 NPAR=1
      IF(LPHYS) THEN
         NPART=3+NL+NED+NDEL
         DO 270 IL=1,NL
         NPART=NPART+MAX(ISMAX(IL,IGRP)-ISMIN(IL,IGRP)+1,0)
  270    CONTINUE
         IF(NPART.GT.3+NL+NL*LPART+NED+NDEL) CALL XABORT('LIBFQD: BUG.')
         CALL LIBTAB (IGRP,NGRO,NL,NDIL,NPART,NED,NDEL,HNAMIS,IMPX,
     1   LSCAT,LSIGF,LADD,DILUT,TOTAL,SIGF,SIGS,SCAT,SADD,ZDEL,1.0,
     2   ISMIN,ISMAX,NOR(IGRP),PRTSIG)
         DO 280 JNOR=1,NOR(IGRP)
         PRTABS(JNOR)=PRTSIG(JNOR,2)-PRTSIG(JNOR,4)
  280    CONTINUE
         GO TO 780
      ELSE IF(LCALEN.OR.LRIBON) THEN
         ALLOCATE(MOMT(2*MAXNOR),MOMP(MAXNOR),SEFR((NPAR+2)*MDIL))
         CALL LIBMOM(NFS(IGRP),MDIL,NPAR,DELTA(IBIN+1),SIGTF(IBIN+1),
     1   SIGAF(IBIN+1),SIGTF(IBIN+1),MAXNOR,DILUT2,MOMT,
     2   MOMP,SEFR)
*
         CALL LIBCAT(MAXNOR,NPAR,MDIL,MOMT,MOMP,IPRECI,LNORAJ,DILUT2,
     1   SEFR,NOR(IGRP),PRTSIG,ERRBST)
         ERROR1=ERRBST
*
         DEALLOCATE(SEFR,MOMP,MOMT)
         DO 285 JNOR=1,NOR(IGRP)
         PRTABS(JNOR)=PRTSIG(JNOR,3)! absorption
         PRTSIG(JNOR,3)=0.0
  285    CONTINUE
*---
      ELSE IF(LPTMC) THEN
         IF(LBSIGF) NPAR=2
         ALLOCATE(MOMT(2*MAXNOR),MOMP(NPAR*MAXNOR),SEFR((NPAR+2)*MDIL))
*        CALENDF PT FOR SIGT, SIGS AND NUSIGF
         CALL LIBMOM(NFS(IGRP),MDIL,NPAR,DELTA(IBIN+1),SIGTF(IBIN+1),
     1   SIGSF(IBIN+1),SIGFF(IBIN+1),MAXNOR,DILUT2,MOMT,MOMP,SEFR)
         CALL LIBCAT(MAXNOR,NPAR,MDIL,MOMT,MOMP,IPRECI,LNORAJ,DILUT2,
     1   SEFR,NOR(IGRP),PRTSIG1,ERRBST)
         ERROR1=ERRBST
         DEALLOCATE(SEFR,MOMP,MOMT)
*
         DO INOR=1,NOR(IGRP)
            PRTSIG(INOR,1)=PRTSIG1(INOR,1)!weight
            PRTSIG(INOR,2)=PRTSIG1(INOR,2)!total
            IF(LBSIGF) THEN
               PRTSIG(INOR,3)=PRTSIG1(INOR,4)!fission
            ELSE
               PRTSIG(INOR,3)=0.0
            ENDIF
            PRTSIG(INOR,4)=PRTSIG1(INOR,3)!scattering
            PRTABS(INOR)=PRTSIG(INOR,2)-PRTSIG(INOR,4)! absorption
         ENDDO
*
         IOF2=4+NL
         DO IL=1,NL
            DO JGRP=ISMIN(IL,IGRP),ISMAX(IL,IGRP)
               NPAR2=1
               ALLOCATE(MOMT(2*MAXNOR),MOMP(NPAR*MAXNOR),
     1         SEFR((NPAR+2)*MDIL))
               CALL LIBMOM(NFS(IGRP),MDIL,NPAR2,DELTA(IBIN+1),
     1         SIGTF(IBIN+1),SCAT00(IBIN+1,JGRP),SIGTF(IBIN+1),MAXNOR,
     2         DILUT2,MOMT,MOMP,SEFR)
               LNORAJ=.FALSE.
               CALL LIBCAT(MAXNOR,NPAR2,MDIL,MOMT,MOMP,IPRECI,LNORAJ,
     1         DILUT2,SEFR,NOR(IGRP),PRTSIG2,ERRBST)
               ERROR1=ERRBST
               DEALLOCATE(SEFR,MOMP,MOMT)
               DO INOR=1,NOR(IGRP)
                  PRTSIG(INOR,IOF2)=PRTSIG2(INOR,3)
               ENDDO
               IOF2=IOF2+1
            ENDDO
         ENDDO
         NPAR=MAX(NPAR,NPAR2)
      ENDIF
      IF(NOR(IGRP).EQ.0) THEN
         CALL XABORT('LIBFQD: NO SUBGROUPS.')
      ELSE IF(NOR(IGRP).GT.MDIL) THEN
         LNORAJ=.FALSE.
         NOR(IGRP)=MDIL
         GO TO 260
      ENDIF
*----
*  REMOVING SMALL PROBABILITIES.
*----
      INOR=0
  290 INOR=INOR+1
      IF(INOR.GT.NOR(IGRP)) GO TO 310
      IF(ABS(PRTSIG(INOR,1)).LT.1.0E-10) THEN
         DO 305 JNOR=INOR+1,NOR(IGRP)
         DO 300 J=1,NPAR+2
         PRTSIG(JNOR-1,J)=PRTSIG(JNOR,J)
  300    CONTINUE
  305    CONTINUE
         INOR=INOR-1
         NOR(IGRP)=NOR(IGRP)-1
      ENDIF
      GO TO 290
*
  310 IF(LRIBON.AND.(IMPX.GT.3)) THEN
         WRITE(6,'(/7X,11HPROBABILITY,7X,5HTOTAL,2X,10HABSORPTION)')
         TEST(:3)=0.0
         DO 320 INOR=1,NOR(IGRP)
         TEST(1)=TEST(1)+PRTSIG(INOR,1)
         TEST(2)=TEST(2)+PRTSIG(INOR,1)*PRTSIG(INOR,2)
         TEST(3)=TEST(3)+PRTSIG(INOR,1)*PRTSIG(INOR,3)
         WRITE(6,'(1X,I5,1P,3E12.4)') INOR,(PRTSIG(INOR,J),J=1,3)
  320    CONTINUE
         WRITE(6,'(6H CHECK,1P,3E12.4)') (TEST(J),J=1,3)
         TEST(:3)=0.0
         DO 330 I=1,NFS(IGRP)
         TEST(1)=TEST(1)+DELTA(IBIN+I)
         TEST(2)=TEST(2)+SIGTF(IBIN+I)*DELTA(IBIN+I)
         TEST(3)=TEST(3)+SIGAF(IBIN+I)*DELTA(IBIN+I)
  330    CONTINUE
         DO 340 J=2,3
         TEST(J)=TEST(J)/TEST(1)
  340    CONTINUE
         TEST(1)=1.0
         WRITE(6,'(6H EXACT,1P,3E12.4)') (TEST(J),J=1,3)
      ENDIF
*----
*  COMPUTE THE REFERENCE SELF-SHIELDED CROSS SECTIONS AT SELECTED
*  VALUES OF THE DILUTION FOR AN HOMOGENEOUS MEDIA. SECOL-TYPE
*  APPROXIMATION.
*----
      IF(IBIN+NFS(IGRP).GT.LBIN) CALL XABORT('LIBFQD: PHI OVERFLOW.')
*
      DO 405 IDIL=1,MDIL
      DO 400 IL=1,NL
      IF(LPHYS.OR.LCALEN.OR.LPTMC) THEN
*        USE A BONDARENKO RESONANT FLUX.
         T0=0.0D0
         DO 350 INOR=1,NOR(IGRP)
         BONDAR=(DILUT2(IDIL)+PRTSIG(INOR,2))**IL
         PRTPH(INOR,IDIL,IL)=DILUT2(IDIL)/BONDAR
         T0=T0+PRTSIG(INOR,1)*PRTPH(INOR,IDIL,1)
  350    CONTINUE
         IF(IL.EQ.1) BB(IDIL)=FLUX(IGRP,IPERD(IDIL))/T0
         DO 360 INOR=1,NOR(IGRP)
         PRTPH(INOR,IDIL,IL)=PRTPH(INOR,IDIL,IL)*REAL(BB(IDIL))
  360    CONTINUE
      ELSE
*        COMPUTE THE BASE POINTS OF THE RESONANT FLUX.
         JINI=(1-NOR(IGRP))/2
         PHIMT(:NOR(IGRP))=0.0D0
         DELMAC=0.0D0
         DO 385 IGF=1,NFS(IGRP)
         DELM=DELTA(IBIN+IGF)
         SIGT=SIGTF(IBIN+IGF)
         DELMAC=DELMAC+DELM
         T0=PHI(IBIN+IGF,IDIL,IL)*DELM
         T=T0
         DO 370 INOR=1-JINI,NOR(IGRP)
         PHIMT(INOR)=PHIMT(INOR)+T
         T=T*SIGT
  370    CONTINUE
         T=T0/SIGT
         DO 380 INOR=-JINI,1,-1
         PHIMT(INOR)=PHIMT(INOR)+T
         T=T/SIGT
  380    CONTINUE
  385    CONTINUE
         DO 390 INOR=1,NOR(IGRP)
         PHIMT(INOR)=PHIMT(INOR)/DELMAC
  390    CONTINUE
         CALL LIBMPA(NOR(IGRP),JINI,PRTSIG(1,1),PRTSIG(1,2),PHIMT,
     1   PRTPH(1,IDIL,IL))
      ENDIF
  400 CONTINUE
  405 CONTINUE
*----
*  COMPUTE THE BASE POINTS OF THE SLOWING-DOWN SOURCE.
*----
      IF(LRIBON) THEN
         JINI=-NOR(IGRP)/2
         DO 440 IDIL=1,MDIL+1
         PHIMT(:NOR(IGRP))=0.0D0
         DELMAC=0.0D0
         DO 425 IGF=1,NFS(IGRP)
         DELM=DELTA(IBIN+IGF)
         SIGT=SIGTF(IBIN+IGF)
         DELMAC=DELMAC+DELM
         T0=RSTAR(IBIN+IGF,IDIL)*DELM
         T=T0
         DO 410 INOR=1-JINI,NOR(IGRP)
         PHIMT(INOR)=PHIMT(INOR)+T
         T=T*SIGT
  410    CONTINUE
         T=T0/SIGT
         DO 420 INOR=-JINI,1,-1
         PHIMT(INOR)=PHIMT(INOR)+T
         T=T/SIGT
  420    CONTINUE
  425    CONTINUE
         DO 430 INOR=1,NOR(IGRP)
         PHIMT(INOR)=PHIMT(INOR)/DELMAC
  430    CONTINUE
         CALL LIBMPA(NOR(IGRP),JINI,PRTSIG(1,1),PRTSIG(1,2),PHIMT,
     1   PRTRS(1,IDIL))
  440    CONTINUE
      ENDIF
*----
*  NORMALIZATION OF THE FLUX-RELATED BASE POINTS. THIS NORMALIZATION
*  PERMITS TO RE-OBTAIN THE BASE POINTS IN TOTAL X-SECTION IF THE RMS
*  APPROACH IS APPLIED TO THE PRTPH MATRIX.
*----
      DO 490 IDIL=1,MDIL
      T0=0.0D0
      DO 450 INOR=1,NOR(IGRP)
      T0=T0+PRTSIG(INOR,1)*PRTSIG(INOR,2)*PRTPH(INOR,IDIL,1)
  450 CONTINUE
      JDIL=IPERD(IDIL)
      FACTOR=FLUX(IGRP,JDIL)*TOTAL(IGRP,JDIL)/REAL(T0)
      DO 465 IL=1,NL
      DO 460 INOR=1,NOR(IGRP)
      PRTPH(INOR,IDIL,IL)=PRTPH(INOR,IDIL,IL)*FACTOR
  460 CONTINUE
  465 CONTINUE
*
      IF(IMPX.GT.9) THEN
         WRITE(6,'(/7X,11HPROBABILITY,3X,9HFINE-FLUX,4X,9HDILUTION=,1P,
     1   E8.1,5H BARN)') DILUT2(IDIL)
         TEST(1)=0.0
         TEST(2)=0.0
         TEST(3)=0.0
         DO 470 INOR=1,NOR(IGRP)
         PGAR=PRTPH(INOR,IDIL,1)
         TEST(1)=TEST(1)+PRTSIG(INOR,1)
         TEST(2)=TEST(2)+PRTSIG(INOR,1)*PGAR
         TEST(3)=TEST(3)+PRTSIG(INOR,1)*PRTSIG(INOR,2)*PGAR
         WRITE(6,'(1X,I5,1P,2E12.4)') INOR,PRTSIG(INOR,1),PGAR
  470    CONTINUE
         TEST(3)=TEST(3)/TEST(2)
         TEST(2)=TEST(2)/TEST(1)
         WRITE(6,'(6H CHECK,1P,3E12.4)') (TEST(J),J=1,3)
         IF(LRIBON) THEN
            TEST(1)=0.0
            TEST(2)=0.0
            TEST(3)=0.0
            DO 480 IGF=1,NFS(IGRP)
            DELM=DELTA(IBIN+IGF)
            TEST(1)=TEST(1)+DELM
            TEST(2)=TEST(2)+REAL(PHI(IBIN+IGF,IDIL,1))*DELM
            TEST(3)=TEST(3)+REAL(PHI(IBIN+IGF,IDIL,1))*SIGTF(IBIN+IGF)*
     1      DELM
  480       CONTINUE
            TEST(3)=TEST(3)/TEST(2)
            TEST(2)=TEST(2)/TEST(1)
            TEST(1)=1.0
            WRITE(6,'(6H EXACT,1P,3E12.4)') (TEST(J),J=1,3)
         ENDIF
      ENDIF
  490 CONTINUE
*----
*  USE A ROOT MEAN SQUARE TECHNIQUE TO FIND BASE POINTS OF THE
*  SCATTERING XS VECTOR AND MATRIX CORRELATED TO THE TOTAL XS IN
*  GROUP IGRP. NOTE: PRTPH(INOR,IDIL,1) IS USED INSTEAD OF
*  PRTPH(INOR,IDIL,IL) ON LINE LABELED 500 IN ORDER TO BE CONSISTENT
*  WITH USSIT0 AND USSIT1. THIS MAY CHANGE IN FUTURE.
*----
      IF(LPTMC) GO TO 780
      IOF1=4
      IOF2=NL+4
      DO 560 IL=1,NL
      IF(LSCAT(IL)) THEN
         DO 505 INOR=1,NOR(IGRP)
         WORK(MDIL+1,INOR)=1.0D0
         DO 500 IDIL=1,MDIL
         WORK(IDIL,INOR)=PRTPH(INOR,IDIL,1)
  500    CONTINUE
  505    CONTINUE
         CALL ALST2F(NDIL+1,MDIL+1,NOR(IGRP),WORK,PHIMT)
         DO 510 IDIL=1,MDIL+1
         JDIL=IPERD(IDIL)
         BB(IDIL)=SIGS(IGRP,IL,JDIL)*FLUX(IGRP,JDIL)
  510    CONTINUE
         CALL ALST2S(NDIL+1,MDIL+1,NOR(IGRP),WORK,PHIMT,BB,CC)
         DO 520 INOR=1,NOR(IGRP)
         PRTSIG(INOR,IOF1)=REAL(CC(INOR))/PRTSIG(INOR,1)
  520    CONTINUE
         DO 550 JGRP=ISMIN(IL,IGRP),ISMAX(IL,IGRP)
         DO 530 IDIL=1,MDIL+1
         JDIL=IPERD(IDIL)
         BB(IDIL)=SCAT(JGRP,IGRP,IL,JDIL)*FLUX(IGRP,JDIL)
  530    CONTINUE
         CALL ALST2S(NDIL+1,MDIL+1,NOR(IGRP),WORK,PHIMT,BB,CC)
         DO 540 INOR=1,NOR(IGRP)
         PRTSIG(INOR,IOF2)=REAL(CC(INOR))/PRTSIG(INOR,1)
  540    CONTINUE
         IOF2=IOF2+1
  550    CONTINUE
      ENDIF
      IOF1=IOF1+1
  560 CONTINUE
*----
*  COMPUTE THE ROOT MEAN SQUARE COEFFICIENT MATRIX FOR P0 FLUX.
*----
      DO 575 INOR=1,NOR(IGRP)
      WORK(MDIL+1,INOR)=1.0D0
      DO 570 IDIL=1,MDIL
      WORK(IDIL,INOR)=PRTPH(INOR,IDIL,1)
  570 CONTINUE
  575 CONTINUE
      CALL ALST2F(NDIL+1,MDIL+1,NOR(IGRP),WORK,PHIMT)
*----
*  USE A ROOT MEAN SQUARE TECHNIQUE TO FIND BASE POINTS OF THE
*  ABSORPTION XS CORRELATED TO THE TOTAL XS IN GROUP IGRP.
*----
      DO 580 IDIL=1,MDIL+1
      JDIL=IPERD(IDIL)
      BB(IDIL)=(TOTAL(IGRP,JDIL)-SIGS(IGRP,1,JDIL))*FLUX(IGRP,JDIL)
  580 CONTINUE
      CALL ALST2S(NDIL+1,MDIL+1,NOR(IGRP),WORK,PHIMT,BB,CC)
      DO 590 INOR=1,NOR(IGRP)
      PRTABS(INOR)=REAL(CC(INOR))/PRTSIG(INOR,1)
  590 CONTINUE
*----
*  USE A ROOT MEAN SQUARE TECHNIQUE TO FIND BASE POINTS OF THE NU*SIGF
*  XS CORRELATED TO THE TOTAL XS IN GROUP IGRP.
*----
      IF(LSIGF) THEN
         DO 600 IDIL=1,MDIL+1
         JDIL=IPERD(IDIL)
         BB(IDIL)=SIGF(IGRP,JDIL)*FLUX(IGRP,JDIL)
  600    CONTINUE
         CALL ALST2S(NDIL+1,MDIL+1,NOR(IGRP),WORK,PHIMT,BB,CC)
         DO 610 INOR=1,NOR(IGRP)
         PRTSIG(INOR,3)=REAL(CC(INOR))/PRTSIG(INOR,1)
  610    CONTINUE
      ENDIF
*----
*  USE A ROOT MEAN SQUARE TECHNIQUE TO FIND BASE POINTS OF THE
*  ADDITIONAL XS CORRELATED TO THE TOTAL XS IN GROUP IGRP.
*----
      DO 640 IED=1,NED
      IF(LADD(IED)) THEN
         DO 620 IDIL=1,MDIL+1
         JDIL=IPERD(IDIL)
         BB(IDIL)=SADD(IGRP,IED,JDIL)*FLUX(IGRP,JDIL)
  620    CONTINUE
         CALL ALST2S(NDIL+1,MDIL+1,NOR(IGRP),WORK,PHIMT,BB,CC)
         DO 630 INOR=1,NOR(IGRP)
         PRTSIG(INOR,IOF2)=REAL(CC(INOR))/PRTSIG(INOR,1)
  630    CONTINUE
      ENDIF
      IOF2=IOF2+1
  640 CONTINUE
*----
*  USE A ROOT MEAN SQUARE TECHNIQUE TO FIND BASE POINTS OF THE DELAYED
*  NU*SIGF XS CORRELATED TO THE TOTAL XS IN GROUP IGRP.
*----
      IF(LSIGF) THEN
        DO 670 IDEL=1,NDEL
        DO 650 IDIL=1,MDIL+1
        JDIL=IPERD(IDIL)
        BB(IDIL)=ZDEL(IGRP,IDEL,JDIL)*FLUX(IGRP,JDIL)
  650   CONTINUE
        CALL ALST2S(NDIL+1,MDIL+1,NOR(IGRP),WORK,PHIMT,BB,CC)
        DO 660 INOR=1,NOR(IGRP)
        PRTSIG(INOR,IOF2)=REAL(CC(INOR))/PRTSIG(INOR,1)
  660   CONTINUE
        IOF2=IOF2+1
  670   CONTINUE
      ENDIF
*----
*  USE A ROOT MEAN SQUARE TECHNIQUE TO FIND THE ELEMENTS OF THE
*  SLOWING-DOWN RELATED CORRELATED WEIGHT MATRIX AND SECONDARY
*  SCATTERING XS IN GROUP IGRP.
*----
      IF(LPHYS.OR.LCALEN.OR.LPTMC) THEN
        DO 685 INOR=1,NOR(IGRP)
        PRTSIW(INOR)=PRTSIG(INOR,4)
        DO 680 JNOR=1,NOR(IGRP)
        WSLD((INOR-1)*NOR(IGRP)+JNOR)=PRTSIG(INOR,1)*PRTSIG(JNOR,1)
  680   CONTINUE
  685   CONTINUE
      ELSE
        DO 705 INOR=1,NOR(IGRP)
        DO 690 IDIL=1,MDIL+1
        BB(IDIL)=PRTRS(INOR,IDIL)
  690   CONTINUE
        CALL ALST2S(NDIL+1,MDIL+1,NOR(IGRP),WORK,PHIMT,BB,CC)
        DO 700 I=1,NOR(IGRP)
        WSLD((I-1)*NOR(IGRP)+INOR)=REAL(CC(I))*PRTSIG(INOR,1)
  700   CONTINUE
  705   CONTINUE
*
        DO 730 J=1,NOR(IGRP)
        T0=0.0D0
        DO 710 I=1,NOR(IGRP)
        T0=T0+WSLD((J-1)*NOR(IGRP)+I)
  710   CONTINUE
        DO 720 I=1,NOR(IGRP)
        WSLD((J-1)*NOR(IGRP)+I)=
     1                   REAL(WSLD((J-1)*NOR(IGRP)+I)*(PRTSIG(J,1)/T0))
  720   CONTINUE
        PRTSIW(J)=REAL(T0)/PRTSIG(J,1)
  730   CONTINUE
      ENDIF
      IDOMAX=0
      EROLD1=1.0E10
      EROLD2=1.0E10
      IF(LCALEN.OR.LPTMC) GO TO 780
*----
*  SOLVE SUBGROUP FORM OF THE SLOWING-DOWN EQUATION FOR AN HOMOGENEOUS
*  MIXTURE AT SELECTED DILUTIONS.
*----
      ERROR1=-9999.0
      ERROR2=-9999.0
      IDMAX=0
      DO 770 IDIL=1,MDIL
      DO 750 I=1,NOR(IGRP)
      MATRIX(I,NOR(IGRP)+1)=PRTSIG(I,1)*DILUT2(IDIL)
      DO 740 J=1,NOR(IGRP)
      MATRIX(I,J)=-WSLD((J-1)*NOR(IGRP)+I)*PRTSIW(J)
  740 CONTINUE
      MATRIX(I,I)=MATRIX(I,I)+(PRTSIG(I,2)+DILUT2(IDIL))*PRTSIG(I,1)
  750 CONTINUE
      CALL ALSBD(NOR(IGRP),1,MATRIX,IER,MAXNOR)
      IF(IER.NE.0) CALL XABORT('LIBFQD: SINGULAR MATRIX(2).')
*----
*  TEST THE ACCURACY OF THE PROBABILITY TABLES FOR THIS ENERGY GROUP.
*----
      ACCUM1=0.0D0
      ACCUM2=0.0D0
      ACCUM3=0.0D0
      ACCUM4=0.0D0
      DO 760 I=1,NOR(IGRP)
      ACCUM1=ACCUM1+PRTSIG(I,1)*PRTABS(I)*MATRIX(I,NOR(IGRP)+1)
      ACCUM2=ACCUM2+PRTSIG(I,1)*MATRIX(I,NOR(IGRP)+1)
      ACCUM3=ACCUM3+PRTSIG(I,1)*PRTABS(I)*PRTPH(I,IDIL,1)
      ACCUM4=ACCUM4+PRTSIG(I,1)*PRTPH(I,IDIL,1)
  760 CONTINUE
      ACCUM1=ACCUM1/ACCUM2
      ACCUM3=ACCUM3/ACCUM4
      IF(ABS(ACCUM1-XSREF(IDIL))/ABS(XSREF(IDIL)).GT.ERROR1) THEN
         EROLD1=ERROR1
         EROLD2=ERROR2
         IDOMAX=IDMAX
         ERROR1=ABS(REAL(ACCUM1)-XSREF(IDIL))/ABS(XSREF(IDIL))
         ERROR2=ABS(REAL(ACCUM3)-XSREF(IDIL))/ABS(XSREF(IDIL))
         IDMAX=IDIL
      ELSE IF(ABS(REAL(ACCUM1)-XSREF(IDIL))/ABS(XSREF(IDIL)).GT.EROLD1)
     1  THEN
         EROLD1=ABS(REAL(ACCUM1)-XSREF(IDIL))/ABS(XSREF(IDIL))
         EROLD2=ABS(REAL(ACCUM3)-XSREF(IDIL))/ABS(XSREF(IDIL))
         IDOMAX=IDIL
      ENDIF
  770 CONTINUE
      IF(IMPX.GT.1) THEN
         TAG='=>'
         IF(LPHYS) TAG='--'
         IF(LCALEN) TAG='=='
         IF(LPTMC) TAG='>>'
         WRITE(6,900) TAG,IGRP,NOR(IGRP),ERROR1*100.0,ERROR2*100.0,
     1   DILUT2(IDMAX),EROLD1*100.0,EROLD2*100.0,DILUT2(IDOMAX)
      ENDIF
      IF(ERROR1.GT.0.01) THEN
         WRITE(HSMG,'(42HLIBFQD: UNABLE TO COMPUTE THE PROBABILITY ,
     1   15HTABLES IN GROUP,I4,17H. TABLE ACCURACY=,1P,E9.2,2H %,
     2   10H ISOTOPE='',A12,2H''.)') IGRP,ERROR1*100.0,HNAMIS
         WRITE(6,'(1X,A)') HSMG
      ENDIF
*
  780 IF((IMPX.GT.2).AND.(NOR(IGRP).GT.1)) THEN
         WRITE(6,'(/7H GROUP=,I4,16H TABLE ACCURACY=,1P,E9.2,2H %)')
     1   IGRP,ERROR1*100.0
         WRITE(6,'(/7X,11HPROBABILITY,7X,5HTOTAL,2X,10HABSORPTION,2X,
     1   10HNU-FISSION,2X,10HSCATTERING,12(1H.))')
         TEST(:8)=0.0
         IOF=NL+IGRP-ISMIN(1,IGRP)
         JMIN=5
         JMAX=MIN(JMIN+ISMAX(1,IGRP)-ISMIN(1,IGRP),8)
         DO 790 JNOR=1,NOR(IGRP)
         TEST(1)=TEST(1)+PRTSIG(JNOR,1)
         TEST(2)=TEST(2)+PRTSIG(JNOR,1)*PRTSIG(JNOR,2)
         TEST(3)=TEST(3)+PRTSIG(JNOR,1)*PRTABS(JNOR)
         TEST(4)=TEST(4)+PRTSIG(JNOR,1)*PRTSIG(JNOR,3)
         DO J=JMIN,JMAX
            TEST(J)=TEST(J)+PRTSIG(JNOR,1)*PRTSIG(JNOR,IOF+J-1)
         ENDDO
         WRITE(6,'(1X,I5,1P,8E12.4)') JNOR,(PRTSIG(JNOR,J),J=1,2),
     1   PRTABS(JNOR),PRTSIG(JNOR,3),(PRTSIG(JNOR,IOF+J-1),J=JMIN,JMAX)
  790    CONTINUE
         WRITE(6,'(6H CHECK,1P,8E12.4)') (TEST(J),J=1,JMAX)
         TEST(:8)=0.0
         TEST(1)=1.0
         TEST(2)=TOTAL(IGRP,NDIL+1)
         TEST(3)=TOTAL(IGRP,NDIL+1)-SIGS(IGRP,1,NDIL+1)
         TEST(4)=SIGF(IGRP,NDIL+1)
         DO J=JMIN,JMAX
            TEST(J)=SCAT(IGRP+J-5,IGRP,1,NDIL+1)
         ENDDO
         WRITE(6,'(6H EXACT,1P,8E12.4)') (TEST(I),I=1,JMAX)
         TEST(:8)=0.0
*---
* CHECK POINT BASES OF THE SCATTERING MATRIX
*---
         IF(IGRP.GE.IGRMIN.AND.IGRP.LT.IGRMAX) THEN
            DIFFS(:20,:MAXDIL)=0.0
            DO IPART=4+NL,4+NL+ISMAX(1,IGRP)-ISMIN(1,IGRP)
            DO IDIL=1,NDIL+1
            TEST(:8)=0.0
            DO INOR=1,NOR(IGRP)
            TEST(1)=TEST(1)+PRTSIG(INOR,1)*PRTSIG(INOR,IPART)/
     1          (PRTSIG(INOR,2)+DILUT2(IDIL))
            TEST(2)=TEST(2)+PRTSIG(INOR,1)/(PRTSIG(INOR,2)+DILUT2(IDIL))
            ENDDO
            TSCAT(IPART,IDIL)=TEST(1)/TEST(2)
            DIFFS(IPART,IDIL)=(TSCAT(IPART,IDIL)-
     1      SCAT(IGRP+IPART-4-NL,IGRP,1,IDIL))/
     2      SCAT(IGRP+IPART-4-NL,IGRP,1,IDIL)
            ENDDO
            WRITE(6,*)'SCATTERING MATRIX COEFFICIENTS FOR IDIL=1,NDIL+1'
            WRITE(6,'(11H SECONDARY ,1P,I3,9H PRIMARY ,1P,I3)')
     1      IGRP+IPART-4-NL,IGRP
            WRITE(6,*)'CALENDF'
            WRITE(6,'(5(4X,F12.4))') (TSCAT(IPART,IDIL),IDIL=1,NDIL+1)
            WRITE(6,*)'NJOY'
            WRITE(6,'(5(4X,F12.4))') (SCAT(IGRP+IPART-4-NL,
     1          IGRP,1,IDIL),IDIL=1,NDIL+1)
            WRITE(*,*)'RELATIVE DIFFERENCE (%)'
            WRITE(6,'(5(4X,F12.4))')
     1          (1.E2*DIFFS(IPART,IDIL),IDIL=1,NDIL+1)
            ENDDO
         ENDIF
      ENDIF
      IF(NOR(IGRP).GT.1) THEN
*        SAVE THE PROBABILITY TABLE INTO IPLIB.
         KPLIB=LCMDIL(JPLIB,IGRP)
         NPART=3+NL+NED+NDEL
         DO 800 IL=1,NL
         ISM(1,IL)=ISMIN(IL,IGRP)
         ISM(2,IL)=ISMAX(IL,IGRP)
         NPART=NPART+MAX(0,(ISMAX(IL,IGRP)-ISMIN(IL,IGRP)+1))
  800    CONTINUE
         CALL LCMPUT(KPLIB,'PROB-TABLE',NPART*MAXNOR,2,PRTSIG)
         IF(LRIBON) THEN
            CALL LCMPUT(KPLIB,'SIGQT-SIGS',NOR(IGRP),2,PRTSIW)
            CALL LCMPUT(KPLIB,'SIGQT-SLOW',NOR(IGRP)**2,2,WSLD)
         ELSE IF(LCALEN.OR.LPTMC) THEN
            IOF=NL+IGRP-ISMIN(1,IGRP)
            CALL LCMPUT(KPLIB,'SIGQT-SIGS',NOR(IGRP),2,PRTSIG(1,IOF+4))
         ENDIF
         CALL LCMPUT(KPLIB,'ISM-LIMITS',2*NL,1,ISM)
      ENDIF
      IBIN=IBIN+NFS(IGRP)
  810 CONTINUE
      CALL LCMPUT(IPLIB,'NOR',NGRO,1,NOR)
      CALL LCMSIX(IPLIB,' ',2)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(RSTAR,CC,WORK,PHI,MATRIX,PHIMT)
      DEALLOCATE(GAR,PRTSIG2,SCAT00,PRTSIG1,PRTRS,PRTABS,PRTSIW,PRTPH,
     1 PRI,SIGAF,STIS,UUU,DELTA,WSLD,PRTSIG)
      DEALLOCATE(ISM)
      RETURN
*
  900 FORMAT(/9H LIBFQD: ,A2,6HGROUP=,I3,7H ORDER=,I2,7H ERROR=,1P,
     1 E9.2,4H % (,E9.2,15H %) AT DILUTION,E10.3,5H BARN/29X,
     2 7H ERROR=,1P,E9.2,4H % (,E9.2,15H %) AT DILUTION,E10.3,6H BARN.)
  910 FORMAT(/32H LIBFQD: TEST FINE FLUX IN GROUP,I5,14H FOR ISOTOPE ',
     1 A12,2H':/9H DILUTION,16X,5HTOTAL,14X,10HABSORPTION,12X,
     2 12HSLOWING-DOWN,20X,4HFLUX/11X,7HAUTOLIB,8X,4HNJOY,5X,7HAUTOLIB,
     3 8X,4HNJOY)
      END
