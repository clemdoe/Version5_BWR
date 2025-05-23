*DECK MCGFCF
      SUBROUTINE MCGFCF(SUBFFI,SUBFFA,SUBLDC,SUBSCH,IFTRAK,NBTR,NMAX,
     1                  NDIM,KPN,K,NREG,M,NGEFF,NANGL,NMU,NLF,NFUNL,
     2                  NMOD,NLFX,NLIN,NFUNLX,KEYFLX,KEYCUR,NZON,NCONV,
     3                  CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,S,SIGAL,ISGNR,IDIR,
     4                  NSOUT,NBATCH,XSI,PHI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Flux integration upon the non-cyclic tracking.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier
*
*Parameters: input
* SUBFFI  flux integration subroutine with isotropic source.
* SUBFFA  flux integration subroutine with anisotropic source.
* SUBLDC  flux integration subroutine with linear-discontinuous source.
* SUBSCH  track coefficients calculation subroutine.
* IFTRAK  tracking file unit number.
* NBTR    total number of tracking lines.
* NMAX    maximum number of elements in a track.
* NDIM    number of dimensions for the geometry.
* KPN     total number of unknowns in vectors PHI.
* K       total number of volumes for which specific values
*         of the neutron flux and reactions rates are required.
* NREG    number of volumes.
* M       number of material mixtures.
* NGEFF   number of groups to process.
* NANGL   number of tracking angles in the tracking file.
* NMU     order of the polar quadrature in 2D / 1 in 3D.
* NLF     number of Legendre orders for the flux.
* NFUNL   number of moments of the flux (in 2D: NFUNL=NLF*(NLF+1)/2).
* NMOD    first dimension of ISGNR.
* NLFX    scattering anisotropy used to compute spherical harmonics.
* NLIN    linear discontinuous flag (=1 SC/DD0; =3 LDC/DD1).
* NFUNLX  number of spherical harmonics components.
* KEYFLX  position of flux elements in PHI vector.
* KEYCUR  position of current elements in PHI vector.
* NZON    index-number of the mixture type assigned to each volume.
* NCONV   logical array of convergence status for each group (.TRUE.
*         not converged).
* CAZ0    cosines of the tracking polar angles in 3D.
* CAZ1    first cosines of the different tracking azimuthal angles.
* CAZ2    second cosines of the different tracking azimuthal angles.
* CPO     cosines of the different tracking polar angles in 2D.
* ZMU     polar quadrature set in 2D.
* WZMU    polar quadrature set in 2D.
* S       total source vector components.
* SIGAL   total cross-section and albedo array.
* ISGNR   sign of correction.
* IDIR    direction of fundamental current for TIBERE with MoC 
*         (=0,1,2,3). 
* NSOUT   number of outer surfaces.
* NBATCH  number of tracks processed in each OpenMP core (default: =1).
* XSI     x,y and z component of the shape parameter for TIBERE. 
*
*Parameters: input/output
* PHI     vector containing the zonal scalar flux.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NGEFF,K,KPN,M,NMAX,NDIM,NMU,NZON(K),NLF,NFUNL,NMOD,
     1 NLFX,NLIN,NFUNLX,NREG,KEYFLX(NREG,NLIN,NFUNL),KEYCUR(K-NREG),
     2 IFTRAK,NBTR,NANGL,ISGNR(NMOD,NFUNLX),IDIR,NSOUT,NBATCH
      REAL CPO(NMU),ZMU(NMU),WZMU(NMU),SIGAL(-6:M,NGEFF)
      DOUBLE PRECISION CAZ0(NANGL),CAZ1(NANGL),CAZ2(NANGL),
     1 PHI(KPN,NGEFF),S(KPN,NGEFF),XSI(NSOUT)
      LOGICAL NCONV(NGEFF)
      EXTERNAL SUBFFI,SUBFFA,SUBLDC,SUBSCH
*----
*  LOCAL VARIABLES
*----
      INTEGER I,II,ILINE,IMU,IANG0,NOMP,INDP,NOMM,INDM,NOMI,JF,IND,
     1 NSUB,INDX,INDY,IREG,I0,NDFUNLX,IBATCH,IL1
      REAL XMUANG(1)
      DOUBLE PRECISION WEIGHT,Q0,Q1,Q0X,Q1X,Q0Y,Q1Y,ZMUI,OMEGA2(3),ZZZ
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NSEG,IANG
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NOM
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: RHARM,TRHAR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: T2D,WEITF,B,FLUX
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: HTF,COEFI,FLUV,
     1 DFLUV
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: PHIV,DPHIV
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: STOT,DSTOT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NSEG(NBATCH),WEITF(NBATCH),IANG(NBATCH),NOM(NMAX,NBATCH),
     1 HTF(NMAX,NBATCH),FLUX(KPN))
*---
* Compute flux and currents for this tracking line
*---
      PHI(:KPN,:NGEFF)=0.0D0
      IF((NLF.EQ.1).AND.(NLIN.EQ.1)) THEN
*     --------------------
*     Isotropic Scattering
*     --------------------
         ALLOCATE(B(2*NMAX))
         IF(NDIM.EQ.3) THEN
*        ---
*        3D calculation -> no loop over the polar angle
*        ---
         DO IBATCH=1,(NBTR-1)/NBATCH+1
         DO ILINE=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NBTR)
           IL1=ILINE-(IBATCH-1)*NBATCH
           READ(IFTRAK) NSUB,NSEG(IL1),WEITF(IL1),IANG(IL1),
     1     (NOM(I,IL1),I=1,NSEG(IL1)),(HTF(I,IL1),I=1,NSEG(IL1))
           IF(NSUB.NE.1) CALL XABORT('MCGFCF: NSUB.NE.1.')
         ENDDO
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(IL1,OMEGA2,ZZZ,FLUX,ILINE,B)
         DO II=1,NGEFF
           IF(NCONV(II)) THEN
             FLUX(:KPN)=0.0D0
             DO ILINE=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NBTR)
               IL1=ILINE-(IBATCH-1)*NBATCH
*                  MCGFFIR: 'Source Term Isolation' Strategy turned on
*                  MCGFFIS: 'Source Term Isolation' Strategy turned off
*                  MCGFFIT: 'MOCC/MCI' Iterative Strategy
               OMEGA2(3)=CAZ0(IANG(IL1))*CAZ0(IANG(IL1))
               ZZZ=1.0D0/SQRT(1.0D0-OMEGA2(3))
               OMEGA2(1)=3.0D0*(CAZ1(IANG(IL1))/ZZZ)**2
               OMEGA2(2)=3.0D0*(CAZ2(IANG(IL1))/ZZZ)**2
               OMEGA2(3)=3.0D0*OMEGA2(3)
               CALL SUBFFI(SUBSCH,K,KPN,M,NSEG(IL1),HTF(1,IL1),
     1              NOM(1,IL1),NZON,SIGAL(0,II),S(1,II),NREG,KEYFLX,
     2              KEYCUR,FLUX,B,WEITF(IL1),OMEGA2,IDIR,NSOUT,XSI)
             ENDDO ! ILINE
             PHI(:KPN,II)=PHI(:KPN,II)+FLUX(:KPN)
           ENDIF
         ENDDO ! II
*$OMP END PARALLEL DO
         ENDDO ! IBATCH
         ELSE
*        ---
*        2D calculation -> loop over the polar angle
*        ---
         ALLOCATE(T2D(NMAX))
         DO IBATCH=1,(NBTR-1)/NBATCH+1
         DO ILINE=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NBTR)
           IL1=ILINE-(IBATCH-1)*NBATCH
           READ(IFTRAK) NSUB,NSEG(IL1),WEITF(IL1),IANG(IL1),
     1     (NOM(I,IL1),I=1,NSEG(IL1)),(HTF(I,IL1),I=1,NSEG(IL1))
           IF(NSUB.NE.1) CALL XABORT('MCGFCF: NSUB.NE.1.')
         ENDDO
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(IL1,OMEGA2,ZMUI,WEIGHT,I0,T2D,FLUX,IMU,ILINE,B)
         DO II=1,NGEFF            
           IF(NCONV(II)) THEN
             FLUX(:KPN)=0.0D0
             DO ILINE=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NBTR)
               IL1=ILINE-(IBATCH-1)*NBATCH
               DO IMU=1,NMU
                 ZMUI=ZMU(IMU)
                 OMEGA2(1)=3.0D0*(CAZ1(IANG(IL1))/ZMUI)**2
                 OMEGA2(2)=3.0D0*(CAZ2(IANG(IL1))/ZMUI)**2
                 OMEGA2(3)=3.0D0*(1.0-1.0/ZMUI**2)
                 WEIGHT=WEITF(IL1)*DBLE(WZMU(IMU))
                 DO I0=2,NSEG(IL1)-1
                   T2D(I0)=HTF(I0,IL1)*ZMUI
                 ENDDO
                 CALL SUBFFI(SUBSCH,K,KPN,M,NSEG(IL1),T2D,NOM(1,IL1),
     1                 NZON,SIGAL(0,II),S(1,II),NREG,KEYFLX,KEYCUR,
     2                 FLUX,B,WEIGHT,OMEGA2,IDIR,NSOUT,XSI)
               ENDDO
             ENDDO ! ILINE
             PHI(:KPN,II)=PHI(:KPN,II)+FLUX(:KPN)
           ENDIF
         ENDDO ! II
*$OMP END PARALLEL DO
         ENDDO ! IBATCH
         DEALLOCATE(T2D)
         ENDIF
         DEALLOCATE(B)
      ELSE IF(NLIN.EQ.1) THEN
*     ----------------------
*     Anisotropic Scattering
*     ----------------------
         ALLOCATE(STOT(NMAX,NMU,2),B(2*NMAX))
         ALLOCATE(RHARM(NMU,NFUNL,NBATCH),TRHAR(NMU,NFUNL,2))
         IANG0=0
         IF(NDIM.EQ.3) THEN
*        ---
*        3D calculation -> no loop over the polar angle
*        ---
         DO IBATCH=1,(NBTR-1)/NBATCH+1
         DO ILINE=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NBTR)
           IL1=ILINE-(IBATCH-1)*NBATCH
           READ(IFTRAK) NSUB,NSEG(IL1),WEITF(IL1),IANG(IL1),
     1     (NOM(I,IL1),I=1,NSEG(IL1)),(HTF(I,IL1),I=1,NSEG(IL1))
           IF(NSUB.NE.1) CALL XABORT('MCGFCF: NSUB.NE.1.')
           IF(IANG(IL1).NE.IANG0) THEN
             IANG0=IANG(IL1)
             XMUANG(1)=REAL(CAZ0(IANG(IL1)))
             CALL MOCCHR(3,NLF-1,NFUNL,1,XMUANG(1),CAZ1(IANG(IL1)),
     1                   CAZ2(IANG(IL1)),RHARM(1,1,IL1))
           ELSE IF(IL1.EQ.1) THEN
             RHARM(:NMU,:NFUNL,IL1)=RHARM(:NMU,:NFUNL,NBATCH)
           ELSE
             RHARM(:NMU,:NFUNL,IL1)=RHARM(:NMU,:NFUNL,IL1-1)
           ENDIF
         ENDDO
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(IL1,TRHAR,NOMP,INDP,NOMM,INDM,STOT,NOMI,Q0,Q1,JF,IND)
*$OMP2 PRIVATE(FLUX,I0,ILINE,B)
         DO II=1,NGEFF            
           IF(NCONV(II)) THEN
             FLUX(:KPN)=0.0D0
             DO ILINE=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NBTR)
               IL1=ILINE-(IBATCH-1)*NBATCH
               DO 10 JF=1,NFUNL
                 TRHAR(1,JF,1)=ISGNR(1,JF)*RHARM(1,JF,IL1)
                 TRHAR(1,JF,2)=ISGNR(NMOD,JF)*RHARM(1,JF,IL1)
 10            CONTINUE
               STOT(:NMAX,:NMU,:2)=0.0D0
*              incoming flux in + direction
               NOMP=NOM(1,IL1)
               INDP=KEYCUR(-NOMP)
*              incoming flux in - direction
               NOMM=NOM(NSEG(IL1),IL1)
               INDM=KEYCUR(-NOMM)
               STOT(1,1,1)=WEITF(IL1)*S(INDP,II)
               STOT(NSEG(IL1),1,2)=WEITF(IL1)*S(INDM,II)
*              regional sources   
               DO I0=2,NSEG(IL1)-1
                 NOMI=NOM(I0,IL1)
                 Q0=0.0D0
                 Q1=0.0D0
                 DO JF=1,NFUNL
                   IND=KEYFLX(NOMI,1,JF)         
                   Q0=Q0+S(IND,II)*TRHAR(1,JF,1)
                   Q1=Q1+S(IND,II)*TRHAR(1,JF,2)
                 ENDDO                       
                 STOT(I0,1,1)=WEITF(IL1)*Q0
                 STOT(I0,1,2)=WEITF(IL1)*Q1
               ENDDO
*                   MCGFFAR: 'Source Term Isolation' Strategy turned on
*                   MCGFFAS: 'Source Term Isolation' Strategy turned off
*                   MCGFFAT: 'MOCC/MCI' Iterative Strategy
               CALL SUBFFA(SUBSCH,K,KPN,M,NSEG(IL1),HTF(1,IL1),
     1              NOM(1,IL1),NZON,SIGAL(0,II),STOT(1,1,1),
     2              STOT(1,1,2),NREG,1,NLF,NFUNL,TRHAR,KEYFLX,
     3              KEYCUR,1,FLUX,B)
             ENDDO ! ILINE
             PHI(:KPN,II)=PHI(:KPN,II)+FLUX(:KPN)
           ENDIF
         ENDDO ! II
*$OMP END PARALLEL DO
         ENDDO ! IBATCH
*        ---
         ELSE
*        ---
*        2D calculation -> loop over the polar angle
*        ---
         ALLOCATE(T2D(NMAX))
         IANG0=0
         DO IBATCH=1,(NBTR-1)/NBATCH+1
         DO ILINE=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NBTR)
           IL1=ILINE-(IBATCH-1)*NBATCH
           READ(IFTRAK) NSUB,NSEG(IL1),WEITF(IL1),IANG(IL1),
     1     (NOM(I,IL1),I=1,NSEG(IL1)),(HTF(I,IL1),I=1,NSEG(IL1))
           IF(NSUB.NE.1) CALL XABORT('MCGFCF: NSUB.NE.1.')
           IF(IANG(IL1).NE.IANG0) THEN
             IANG0=IANG(IL1)
             CALL MOCCHR(2,NLF-1,NFUNL,NMU,CPO(1),CAZ1(IANG(IL1)),
     1                   CAZ2(IANG(IL1)),RHARM(1,1,IL1))
           ELSE IF(IL1.EQ.1) THEN
             RHARM(:NMU,:NFUNL,IL1)=RHARM(:NMU,:NFUNL,NBATCH)
           ELSE
             RHARM(:NMU,:NFUNL,IL1)=RHARM(:NMU,:NFUNL,IL1-1)
           ENDIF
         ENDDO
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(IL1,TRHAR,NOMP,INDP,NOMM,INDM,NOMI,IMU,WEIGHT,STOT,I0)
*$OMP2 PRIVATE(JF,Q0,Q1,ZMUI,T2D,FLUX,ILINE,B)
         DO II=1,NGEFF            
           IF(NCONV(II)) THEN
             FLUX(:KPN)=0.0D0
             DO ILINE=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NBTR)
               IL1=ILINE-(IBATCH-1)*NBATCH
               DO 25 JF=1,NFUNL
               DO 20 IMU=1,NMU
                 TRHAR(IMU,JF,1)=ISGNR(1,JF)*RHARM(IMU,JF,IL1)
                 TRHAR(IMU,JF,2)=ISGNR(NMOD,JF)*RHARM(IMU,JF,IL1)
 20            CONTINUE
 25            CONTINUE
               STOT(:NMAX,:NMU,:2)=0.0D0
*              incoming flux in + direction
               NOMP=NOM(1,IL1)
               INDP=KEYCUR(-NOMP)
*              incoming flux in - direction
               NOMM=NOM(NSEG(IL1),IL1)
               INDM=KEYCUR(-NOMM)
               DO IMU=1,NMU
                 WEIGHT=WEITF(IL1)*DBLE(WZMU(IMU))
                 STOT(1,IMU,1)=WEIGHT*S(INDP,II)
                 STOT(NSEG(IL1),IMU,2)=WEIGHT*S(INDM,II)
               ENDDO
*              regional sources               
               DO I0=2,NSEG(IL1)-1
                 NOMI=NOM(I0,IL1)
                 DO IMU=1,NMU
                   Q0=0.0D0
                   Q1=0.0D0
                   WEIGHT=WEITF(IL1)*DBLE(WZMU(IMU))
                   DO JF=1,NFUNL
                     IND=KEYFLX(NOMI,1,JF)         
                     Q0=Q0+S(IND,II)*TRHAR(IMU,JF,1)
                     Q1=Q1+S(IND,II)*TRHAR(IMU,JF,2)
                   ENDDO                       
                   STOT(I0,IMU,1)=WEIGHT*Q0
                   STOT(I0,IMU,2)=WEIGHT*Q1
                 ENDDO
               ENDDO
               DO IMU=1,NMU
                 ZMUI=ZMU(IMU)
                 WEIGHT=WEITF(IL1)*DBLE(WZMU(IMU))
                 DO I=2,NSEG(IL1)-1
                   T2D(I)=HTF(I,IL1)*ZMUI
                 ENDDO
                 CALL SUBFFA(SUBSCH,K,KPN,M,NSEG(IL1),T2D,NOM(1,IL1),
     1                NZON,SIGAL(0,II),STOT(1,IMU,1),STOT(1,IMU,2),
     2                NREG,NMU,NLF,NFUNL,TRHAR,KEYFLX,KEYCUR,IMU,FLUX,B)
               ENDDO
             ENDDO ! ILINE
             PHI(:KPN,II)=PHI(:KPN,II)+FLUX(:KPN)
           ENDIF
         ENDDO ! II
*$OMP END PARALLEL DO
         ENDDO ! IBATCH
*        ---
         DEALLOCATE(T2D)
         ENDIF
         DEALLOCATE(TRHAR,RHARM,B,STOT)
      ELSE IF(NLIN.EQ.3) THEN
*     -----------------------------------------
*     Linear discontinuous source approximation
*     -----------------------------------------
         NDFUNLX=NDIM*NFUNLX
         ALLOCATE(B(6*NMAX))
         ALLOCATE(RHARM(NMU,NFUNLX,NBATCH),TRHAR(NMU,NFUNLX,2))
         ALLOCATE(PHIV(NFUNLX,NREG,NGEFF),DPHIV(NDFUNLX,NREG,NGEFF))
         ALLOCATE(FLUV(NFUNLX,NREG),DFLUV(NDFUNLX,NREG))
         ALLOCATE(STOT(NMAX,NMU,2),DSTOT(NMAX,NMU,2))
         DO II=1,NGEFF            
           IF(NCONV(II)) THEN
             PHIV(:NFUNLX,:NREG,II)=0.0D0
             DPHIV(:NDFUNLX,:NREG,II)=0.0D0
           ENDIF
         ENDDO
         IF(NDIM.EQ.3) THEN
            CALL XABORT('MCGFCF: 3D LDC APPROXIMATION NOT IMPLEMENTED')
         ELSE
*        ---
*        2D calculation -> loop over the polar angle
*        ---
         ALLOCATE(T2D(NMAX))
         IANG0=0
         DO IBATCH=1,(NBTR-1)/NBATCH+1
         DO ILINE=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NBTR)
           IL1=ILINE-(IBATCH-1)*NBATCH
           READ(IFTRAK) NSUB,NSEG(IL1),WEITF(IL1),IANG(IL1),
     1     (NOM(I,IL1),I=1,NSEG(IL1)),(HTF(I,IL1),I=1,NSEG(IL1))
           IF(NSUB.NE.1) CALL XABORT('MCGFCF: NSUB.NE.1.')
           IF(IANG(IL1).NE.IANG0) THEN
             IANG0=IANG(IL1)
             CALL MOCCHR(2,NLFX-1,NFUNLX,NMU,CPO(1),CAZ1(IANG(IL1)),
     1                   CAZ2(IANG(IL1)),RHARM(1,1,IL1))
           ELSE IF(IL1.EQ.1) THEN
             RHARM(:NMU,:NFUNLX,IL1)=RHARM(:NMU,:NFUNLX,NBATCH)
           ELSE
             RHARM(:NMU,:NFUNLX,IL1)=RHARM(:NMU,:NFUNLX,IL1-1)
           ENDIF
         ENDDO
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(IL1,JF,IMU,TRHAR,NOMP,INDP,NOMM,INDM,NOMI,STOT,DSTOT,I0)
*$OMP2 PRIVATE(ZMUI,WEIGHT,T2D,FLUX,FLUV,DFLUV,ILINE,B,Q0,Q1,Q0X,Q1X)
*$OMP3 PRIVATE(Q0Y,Q1Y,IND,INDX,INDY)
         DO II=1,NGEFF            
           IF(NCONV(II)) THEN
             FLUX(:KPN)=0.0D0
             FLUV(:NFUNLX,:NREG)=0.0D0
             DFLUV(:NDFUNLX,:NREG)=0.0D0
             DO ILINE=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NBTR)
               IL1=ILINE-(IBATCH-1)*NBATCH
               DO 35 JF=1,NFUNLX
               DO 30 IMU=1,NMU
                 TRHAR(IMU,JF,1)=ISGNR(1,JF)*RHARM(IMU,JF,IL1)
                 TRHAR(IMU,JF,2)=ISGNR(NMOD,JF)*RHARM(IMU,JF,IL1)
 30            CONTINUE
 35            CONTINUE
               STOT(:NMAX,:NMU,:2)=0.0D0
               DSTOT(:NMAX,:NMU,:2)=0.0D0
*              incoming flux in + direction
               NOMP=NOM(1,IL1)
               INDP=KEYCUR(-NOMP)
*              incoming flux in - direction
               NOMM=NOM(NSEG(IL1),IL1)
               INDM=KEYCUR(-NOMM)
               DO IMU=1,NMU
                  STOT(1,IMU,1)=S(INDP,II)
                  STOT(NSEG(IL1),IMU,2)=S(INDM,II)
               ENDDO
*              regional sources               
               DO I0=2,NSEG(IL1)-1
                 NOMI=NOM(I0,IL1)
                 DO IMU=1,NMU
                   Q0=0.0D0
                   Q1=0.0D0
                   Q0X=0.0D0
                   Q1X=0.0D0
                   Q0Y=0.0D0
                   Q1Y=0.0D0
                   DO JF=1,NFUNL
                     IND=KEYFLX(NOMI,1,JF)         
                     INDX=KEYFLX(NOMI,2,JF)         
                     INDY=KEYFLX(NOMI,3,JF)         
                     Q0=Q0+S(IND,II)*TRHAR(IMU,JF,1)
                     Q1=Q1+S(IND,II)*TRHAR(IMU,JF,2)
                     Q0X=Q0X+S(INDX,II)*TRHAR(IMU,JF,1)
                     Q1X=Q1X+S(INDX,II)*TRHAR(IMU,JF,2)
                     Q0Y=Q0Y+S(INDY,II)*TRHAR(IMU,JF,1)
                     Q1Y=Q1Y+S(INDY,II)*TRHAR(IMU,JF,2)
                   ENDDO                       
                   STOT(I0,IMU,1)=Q0
                   STOT(I0,IMU,2)=Q1
                   DSTOT(I0,IMU,1)=Q0X*CAZ1(IANG(IL1))+Q0Y*
     1             CAZ2(IANG(IL1))
                   DSTOT(I0,IMU,2)=-Q1X*CAZ1(IANG(IL1))-Q1Y*
     1             CAZ2(IANG(IL1))
                 ENDDO
               ENDDO
               DO IMU=1,NMU
                 ZMUI=ZMU(IMU)
                 WEIGHT=WEITF(IL1)*DBLE(WZMU(IMU))
                 DO I0=2,NSEG(IL1)-1
                   T2D(I0)=HTF(I0,IL1)*ZMUI
                 ENDDO
*                MCGFFAL: 'Source Term Isolation' Strategy turned off
                 CALL SUBLDC(SUBSCH,K,KPN,M,NSEG(IL1),T2D,NOM(1,IL1),
     1               NZON,WEIGHT,SIGAL(0,II),STOT(1,IMU,1),
     2               STOT(1,IMU,2),DSTOT(1,IMU,1),DSTOT(1,IMU,2),NREG,
     3               NMU,NLF,NFUNLX,TRHAR,KEYCUR,IMU,B,FLUX,FLUV,DFLUV)
               ENDDO
             ENDDO ! ILINE
             PHI(:KPN,II)=PHI(:KPN,II)+FLUX(:KPN)
             PHIV(:NFUNLX,:NREG,II)=PHIV(:NFUNLX,:NREG,II)+
     1       FLUV(:NFUNLX,:NREG)
             DPHIV(:NDFUNLX,:NREG,II)=DPHIV(:NDFUNLX,:NREG,II)+
     1       DFLUV(:NDFUNLX,:NREG)
           ENDIF
         ENDDO ! II
*$OMP END PARALLEL DO
         ENDDO ! IBATCH
         ALLOCATE(COEFI(2*NFUNLX,2*NFUNLX))
         CALL MCGCOEF(NFUNLX,NMU,ZMU,WZMU,NANGL,CAZ1,CAZ2,COEFI)
         DO II=1,NGEFF            
         IF(NCONV(II)) THEN
           DO IREG=1,NREG
             DPHIV(:,IREG,II)=MATMUL(COEFI,DPHIV(:,IREG,II))
             DO JF=1,NFUNL
               PHI(KEYFLX(IREG,1,JF),II)=PHIV(JF,IREG,II)
               PHI(KEYFLX(IREG,2,JF),II)=DPHIV(JF,IREG,II)
               PHI(KEYFLX(IREG,3,JF),II)=DPHIV(NFUNLX+JF,IREG,II)
             ENDDO
           ENDDO
         ENDIF
         ENDDO
         DEALLOCATE(COEFI)
*        ---
         DEALLOCATE(T2D)
         ENDIF
         DEALLOCATE(DSTOT,STOT,DFLUV,FLUV,DPHIV,PHIV)
         DEALLOCATE(TRHAR,RHARM,B)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FLUX,HTF,NOM,IANG,WEITF,NSEG)
      RETURN
      END
