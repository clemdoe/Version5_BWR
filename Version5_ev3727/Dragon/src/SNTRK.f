*DECK SNTRK
      SUBROUTINE SNTRK(MAXPTS,IPTRK,IPGEOM,IMPX,ISCHM,IELEM,ISPLH,INSB,
     1 NLF,MAXIT,EPSI,ISCAT,IQUAD,LFIXUP,LIVO,ICL1,ICL2,LDSA,NSTART,
     2 NSDSA,IELEMSA,ISOLVSA,LBIHET,LSHOOT,IBFP,MCELL,NMPI,NFOU,
     3 EELEM,ESCHM,IGLK)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover of the geometry and tracking for SN methods.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert and N. Martin
*
*Parameters: input
* MAXPTS  allocated storage for arrays of dimension NEL.
* IPTRK   L_TRACK pointer to the tracking information.
* IPGEOM  L_GEOM pointer to the geometry.
* IMPX    print flag.
* ISCHM   method of spatial discretisation:
*         =1 High-Order Diamond Differencing (HODD) - default;
*         =2 Discontinuous Galerkin finite element method (DG);
*         =3 Adaptive weighted method (AWD).
* IELEM   measure of order of the spatial approximation polynomial:
*         =1 constant - default for HODD;
*         =2 linear - default for DG;
*         >3 higher orders.
* ISPLH   mesh-splitting index for hexagons into lozenges.
* INSB    group vectorization flag (=0/1 group vectorization off/on).
* NLF     SN order for the flux (even number).
* MAXIT   maximum number of inner iterations (default=100).
* EPSI    convergence criterion on inner iterations (default=5.E-5).
* ISCAT   anisotropy of one-speed sources:
*         =1 isotropic sources;
*         =2 linearly anisotropic sources.
* IQUAD   type of SN quadrature (=1 Level symmetric, type IQUAD;
*         =4 Gauss-Legendre and Gauss-Chebyshev; =10 product).
* LFIXUP  flag to enable negative flux fixup.
* LIVO    flag to enable Livolant acceleration.
* ICL1    Number of free iterations with Livolant acceleration.
* ICL2    Number of accelerated iterations with Livolant acceleration.
* LDSA    flag to enable diffusion synthetic acceleration.
* NSTART  restarts the GMRES method every NSTART iterations.
* NSDSA   number of inner flux iterations before enabling SA.
* IELEMSA degree of the Lagrangian finite elements for the SA:
*         <0 order -IELEMSA primal finite elements;
*         >0 order IELEMSA dual finite elements.
* ISOLVSA type of solver to be used for the SA:
*         1 - BIVAC ;
*         2 - TRIVAC.
* LBIHET  flag to enable the double-heterogeneity model.
* LSHOOT  enablig flag for the shooting method.
* IBFP    type of energy proparation relation:
*         =0 no Fokker-Planck term;
*         =1 Galerkin type;
*         =2 heuristic Przybylski and Ligou type.
* MCELL   number of macrocells along each axis (in Cartesian geometry)
*         for the parallelisation using the OpenMP paradigm; OR
*         number of macrocells along the z-axis (in hexagonal geometry)
*         for the parallelisation using the OpenMP paradigm.
* NMPI    number of macrocells along each axis (in Cartesian geometry)
*         or along the z-axis for the hexagonal geometry for the
*         parallelisation using the MPI paradigm, when using WYVERN.
* NFOU    number of Fourier frequencies in the range (2*pi/L) to be
*         investigated.
* ESCHM   method of energy discretisation:
*         =1 High-Order Diamond Differencing (HODD) - default;
*         =2 Discontinuous Galerkin finite element method (DG);
*         =3 Adaptive weighted method (AWD).
* EELEM   measure of order of the energy approximation polynomial:
*         =1 constant - default for HODD;
*         =2 linear - default for DG;
*         >3 higher orders.
* IGLK    angular interpolation type:
*         =0 classical SN method.
*         =1 Galerkin quadrature method (M = inv(D))
*         =2 Galerkin quadrature method (D = inv(M))
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPGEOM
      INTEGER MAXPTS,IMPX,ISCHM,IELEM,ISPLH,INSB,NLF,ISCAT,IQUAD,
     1 MAXIT,ICL1,ICL2,NSTART,NSDSA,IELEMSA,ISOLVSA,MCELL,NMPI,NFOU,
     2 EELEM,ESCHM,IGLK
      REAL EPSI
      LOGICAL LFIXUP,LDSA,LBIHET,LIVO,LSHOOT
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      LOGICAL ILK
      CHARACTER HSMG*131
      INTEGER ISTATE(NSTATE),IGP(NSTATE),NCODE(6),ICODE(6),LOZSWP(3,6),P
      REAL ZCODE(6)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAT,IDL,ISPLX,ISPLY,ISPLZ,
     1 JOP,MRMX,MRMY,MRMZ,IL,IM
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: COORDMAP
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: KEYANI
      REAL, ALLOCATABLE, DIMENSION(:) :: VOL,XXX,YYY,ZZZ,UU,WW,PL,TPQ,
     1 UPQ,VPQ,WPQ,ALPHA,PLZ,SURF,DU,DE,DZ,DC,DB,DA,DAL,WX,WE,CST,MN,DN
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(MAT(MAXPTS),VOL(MAXPTS),IDL(MAXPTS))
*     
      CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
      ITYPE=ISTATE(1)
      ISUB2=ISTATE(9)
*----
*  PARAMETER VALIDATION.
*----
      IF(ISUB2.NE.0) CALL XABORT('SNTRK: DISCRETIZATION NOT AVAILABLE
     1   .(1)')
      IF((ISCHM.NE.1).AND.(ISCHM.NE.2).AND.(ISCHM.NE.3)) THEN
         CALL XABORT('SNTRK: SPATIAL DISCRETIZATION SCHEME NOT AVAILABLE
     1. ONLY VALUES OF 1 (=DIAMOND-DIFFERENCE) OR 2 (=DISCONTINUOUS GALE
     2RKIN) OR 3 (=ADAPTIVE WEIGHTED DIFFERENCE) ARE ALLOWED.')
      ENDIF
      IF((ESCHM.NE.1).AND.(ESCHM.NE.2).AND.(ESCHM.NE.3)) THEN
         CALL XABORT('SNTRK: ENERGY DISCRETIZATION SCHEME NOT AVAILABLE
     1. ONLY VALUES OF 1 (=DIAMOND-DIFFERENCE) OR 2 (=DISCONTINUOUS GALE
     2RKIN) OR 3 (=ADAPTIVE WEIGHTED DIFFERENCE) ARE ALLOWED.')
      ENDIF
      IF(ISCHM.EQ.3.OR.ESCHM.EQ.3) THEN
      IF(.NOT.(ITYPE.EQ.2.OR.ITYPE.EQ.5.OR.ITYPE.EQ.7)) THEN
      CALL XABORT('SNTRK: ADAPTIVE SCHEME ONLY AVAILABLE IN CARTESIAN'
     1 //' GEOMETRIES')
      ELSEIF(IELEM.GT.1.AND.ISCHM.EQ.3.OR.
     1       EELEM.GT.1.AND.ESCHM.EQ.3) THEN
      CALL XABORT('SNTRK: ADAPTIVE SCHEMES ONLY AVAILABLE FOR CONSTANT'
     1 //' ORDER CLOSURE RELATION.') 
      ENDIF        
      ENDIF
      IF((ITYPE.NE.2).AND.(ITYPE.NE.3).AND.(ITYPE.NE.4).AND.
     1   (ITYPE.NE.5).AND.(ITYPE.NE.6).AND.(ITYPE.NE.7).AND.
     2   (ITYPE.NE.8).AND.(ITYPE.NE.9)) THEN
         CALL XABORT('SNTRK: DISCRETIZATION NOT AVAILABLE.(2)')
      ENDIF
      IF((LDSA).AND.(ISOLVSA.EQ.1))THEN
         IF((ITYPE.EQ.7).OR.(ITYPE.EQ.9))
     >      CALL XABORT('SNTRK: SYNTHETIC ACCELERATION WITH BIVAC NOT AV
     1AILABLE IN 3D.')
      ENDIF
*
      ALLOCATE(XXX(MAXPTS+1),YYY(MAXPTS+1),ZZZ(MAXPTS+1))
      ALLOCATE(ISPLX(MAXPTS),ISPLY(MAXPTS),ISPLZ(MAXPTS))
      CALL READ3D(MAXPTS,MAXPTS,MAXPTS,MAXPTS,IPGEOM,IHEX,IR,ILK,SIDE,
     1 XXX,YYY,ZZZ,IMPX,LX,LY,LZ,MAT,NEL,NCODE,ICODE,ZCODE,ISPLX,ISPLY,
     2 ISPLZ,ISPLTH,ISPLTL)
      DEALLOCATE(ISPLZ,ISPLY,ISPLX)
*
      IF(LX*LY*LZ.GT.MAXPTS) THEN
         WRITE (HSMG,'(38HSNTRK: MAXPTS SHOULD BE INCREASED FROM,I7,
     1   3H TO,I7)') MAXPTS,LX*LY*LZ
         CALL XABORT(HSMG)
      ENDIF
      IF((ITYPE.EQ.2).OR.(ITYPE.EQ.3)) THEN
      ! 1-D AND 2-D CYLINDRICAL CASES.
         NCODE(3)=2
         NCODE(4)=5
         ZCODE(3)=1.0
         ZCODE(4)=1.0
         YYY(1)=0.0
         YYY(2)=2.0
      ELSE IF(ITYPE.EQ.6) THEN
         LY=LZ
         DO I=1,LZ+1
            YYY(I)=ZZZ(I)
         ENDDO
         NCODE(3)=NCODE(5)
         NCODE(4)=NCODE(6)
         ZCODE(3)=ZCODE(5)
         ZCODE(4)=ZCODE(6)
         ICODE(3)=ICODE(5)
         ICODE(4)=ICODE(6)
      ENDIF

      IF(IBFP.EQ.0) EELEM=1

*----
*  UNFOLD THE DOMAIN IN DIAGONAL SYMMETRY CASES.
*----
      IF(ITYPE.EQ.7) THEN  
      ! CARTESIAN 3D CASE
      IDIAG=0
      IF((NCODE(2).EQ.3).AND.(NCODE(3).EQ.3)) THEN
         IDIAG=1
         NCODE(3)=NCODE(1)
         NCODE(2)=NCODE(4)
         ZCODE(3)=ZCODE(1)
         ZCODE(2)=ZCODE(4)
         ICODE(3)=ICODE(1)
         ICODE(2)=ICODE(4)
         K=NEL
         DO IZ=LZ,1,-1
            IOFF=(IZ-1)*LX*LY
            DO IY=LY,1,-1
               DO IX=LX,IY+1,-1
                  MAT(IOFF+(IY-1)*LX+IX)=MAT(IOFF+(IX-1)*LY+IY)
               ENDDO
               DO IX=IY,1,-1
                  MAT(IOFF+(IY-1)*LX+IX)=MAT(K)
                  K=K-1
               ENDDO
            ENDDO
         ENDDO
         NEL=LX*LY*LZ
         IF(K.NE.0) THEN
            CALL XABORT('SNTRK: UNABLE TO UNFOLD THE DOMAIN.')
         ENDIF
      ELSE IF((NCODE(1).EQ.3).AND.(NCODE(4).EQ.3)) THEN
         IDIAG=1
         NCODE(1)=NCODE(3)
         NCODE(4)=NCODE(2)
         ZCODE(1)=ZCODE(3)
         ZCODE(4)=ZCODE(2)
         ICODE(1)=ICODE(3)
         ICODE(4)=ICODE(2)
         K=NEL
         DO IZ=LZ,1,-1
            IOFF=(IZ-1)*LX*LY
            DO IY=LY,1,-1
               DO IX=LX,IY,-1
                  MAT(IOFF+(IY-1)*LX+IX)=MAT(K)
                  K=K-1
               ENDDO
            ENDDO
         ENDDO
         DO IZ=1,LZ
            IOFF=(IZ-1)*LX*LY
            DO IY=1,LY
               DO IX=1,IY-1
                  MAT(IOFF+(IY-1)*LX+IX)=MAT(IOFF+(IX-1)*LY+IY)
               ENDDO
            ENDDO
         ENDDO
         NEL=LX*LY*LZ
         IF(K.NE.0) THEN
            CALL XABORT('SNTRK: UNABLE TO UNFOLD THE DOMAIN.')
         ENDIF
      ENDIF

      ELSE 
      ! OTHER CASES
      IF((NCODE(2).EQ.3).AND.(NCODE(3).EQ.3)) THEN
         NCODE(3)=NCODE(1)
         NCODE(2)=NCODE(4)
         ZCODE(3)=ZCODE(1)
         ZCODE(2)=ZCODE(4)
         ICODE(3)=ICODE(1)
         ICODE(2)=ICODE(4)
         K=LX*(LX+1)/2
         DO IY=LY,1,-1
            DO IX=LX,IY+1,-1
               MAT((IY-1)*LX+IX)=MAT((IX-1)*LY+IY)
            ENDDO
            DO IX=IY,1,-1
               MAT((IY-1)*LX+IX)=MAT(K)
               K=K-1
            ENDDO
         ENDDO
         NEL=LX*LY
         IF(K.NE.0) THEN
            CALL XABORT('SNTRK: UNABLE TO UNFOLD THE DOMAIN.')
         ENDIF
      ELSE IF((NCODE(1).EQ.3).AND.(NCODE(4).EQ.3)) THEN
         NCODE(1)=NCODE(3)
         NCODE(4)=NCODE(2)
         ZCODE(1)=ZCODE(3)
         ZCODE(4)=ZCODE(2)
         ICODE(1)=ICODE(3)
         ICODE(4)=ICODE(2)
         K=LX*(LX+1)/2
         DO IY=LY,1,-1
            DO IX=LX,IY,-1
               MAT((IY-1)*LX+IX)=MAT(K)
               K=K-1
            ENDDO
         ENDDO
         DO IY=1,LY
            DO IX=1,IY-1
               MAT((IY-1)*LX+IX)=MAT((IX-1)*LY+IY)
            ENDDO
         ENDDO
         NEL=LX*LY
         IF(K.NE.0) THEN
            CALL XABORT('SNTRK: UNABLE TO UNFOLD THE DOMAIN.')
         ENDIF
      ENDIF
      ENDIF
        
      IF(IMPX.GT.5) THEN
         WRITE(6,600) 'NCODE',(NCODE(I),I=1,6)
         WRITE(6,600) 'ICODE',(ICODE(I),I=1,6)
         WRITE(6,600) 'MAT',(MAT(I),I=1,LX*LY*LZ)
      ENDIF
*---
* CALL TO THE SN TRACKING MODULE RELEVANT TO EACH GEOMETRY
*---
      NDIM=0
      IF(ITYPE.EQ.2) THEN
*        1D SLAB GEOMETRY.
         NDIM=1
         IF(ISCAT.EQ.0) CALL XABORT('SNTRK: SCAT NOT DEFINED.')
         IF(IGLK.NE.0) THEN
            NSCT=NLF
         ELSE
            NSCT=ISCAT
         ENDIF
         ALLOCATE(UU(NLF),WW(NLF),WX(IELEM+1),WE(EELEM+1),
     1   CST(MAX(IELEM,EELEM)),MN(NSCT*NLF),DN(NLF*NSCT),IL(NSCT),
     2   IM(NSCT),PL(NSCT*NLF))
         CALL SNT1DP(IMPX,LX,IELEM,NCODE,ZCODE,XXX,NLF,NSCT,UU,WW,PL,
     1      VOL,IDL,LL4,NUN,LSHOOT,EELEM,WX,WE,CST,IBFP,ISCHM,ESCHM,
     2      IGLK,MN,DN,IL,IM,IQUAD)
         CALL LCMPUT(IPTRK,'U',NLF,2,UU)
         CALL LCMPUT(IPTRK,'W',NLF,2,WW)
         CALL LCMPUT(IPTRK,'WX',IELEM+1,2,WX)
         IF(IBFP.NE.0) CALL LCMPUT(IPTRK,'WE',EELEM+1,2,WE)
         CALL LCMPUT(IPTRK,'CST',MAX(IELEM,EELEM),2,CST)
         CALL LCMPUT(IPTRK,'MN',NSCT*NLF,2,MN)
         CALL LCMPUT(IPTRK,'DN',NLF*NSCT,2,DN)
         CALL LCMPUT(IPTRK,'IL',NSCT,1,IL)
         CALL LCMPUT(IPTRK,'IM',NSCT,1,IM)
         CALL LCMPUT(IPTRK,'PL',NSCT*NLF,2,PL)
         DEALLOCATE(WW,UU,WX,WE,CST,MN,DN,IL,IM,PL)
         ! For Fourier Analysis
         IF(NFOU.GT.0)THEN
            XLEN = XXX(LX+1)
            CALL LCMPUT(IPTRK,'XXX',LX+1,2,XXX)
            CALL LCMPUT(IPTRK,'XLEN',1,2,XLEN)
         ENDIF
      ELSE IF(ITYPE.EQ.3) THEN
*        1D CYLINDRICAL GEOMETRY.
         NDIM=1
         IF(ISCAT.EQ.0) CALL XABORT('SNTRK: SCAT NOT DEFINED.')
         NSCT=(ISCAT/2)*(ISCAT/2+1)+(ISCAT+1)*MOD(ISCAT,2)/2
         IF(IQUAD.GE.10) THEN
*           PRODUCT QUADRATURE.
            NPQ=(NLF**2)/2
         ELSE
            NPQ=NLF*(1+NLF/2)/2
         ENDIF
         ALLOCATE(JOP(NLF/2),UU(NLF/2),WW(NLF/2),TPQ(NPQ),UPQ(NPQ),
     1   VPQ(NPQ),WPQ(NPQ),ALPHA(NPQ),PLZ(NSCT*NLF/2),PL(NSCT*NPQ),
     2   SURF(LX+1),IL(NSCT),IM(NSCT))
         CALL SNT1DC(IMPX,LX,NCODE,ZCODE,XXX,NLF,NPQ,NSCT,IQUAD,JOP,
     1   UU,WW,TPQ,UPQ,VPQ,WPQ,ALPHA,PLZ,PL,VOL,IDL,SURF,IL,IM)
         DEALLOCATE(VPQ,TPQ,WW)
         CALL LCMPUT(IPTRK,'JOP',NLF/2,1,JOP)
         CALL LCMPUT(IPTRK,'U',NLF/2,2,UU)
         CALL LCMPUT(IPTRK,'UPQ',NPQ,2,UPQ)
         CALL LCMPUT(IPTRK,'WPQ',NPQ,2,WPQ)
         CALL LCMPUT(IPTRK,'ALPHA',NPQ,2,ALPHA)
         CALL LCMPUT(IPTRK,'PLZ',NSCT*NLF/2,2,PLZ)
         CALL LCMPUT(IPTRK,'PL',NSCT*NPQ,2,PL)
         CALL LCMPUT(IPTRK,'SURF',LX+1,2,SURF)
         CALL LCMPUT(IPTRK,'IL',NSCT,1,IL)
         CALL LCMPUT(IPTRK,'IM',NSCT,1,IM)
         DEALLOCATE(SURF,PL,PLZ,ALPHA,WPQ,UPQ,UU,JOP,IL,IM)
         LL4=LX*NSCT
         NUN=LL4
      ELSE IF(ITYPE.EQ.4) THEN
*        1D SPHERICAL GEOMETRY.
         NDIM=1
         IF(ISCAT.EQ.0) CALL XABORT('SNTRK: SCAT NOT DEFINED.')
         NSCT=ISCAT
         ALLOCATE(UU(NLF),WW(NLF),ALPHA(NLF),PLZ(NSCT),PL(NSCT*NLF),
     1   SURF(LX+1),IL(NSCT),IM(NSCT))
         CALL SNT1DS(IMPX,LX,NCODE,ZCODE,XXX,NLF,NSCT,UU,WW,ALPHA,
     1   PLZ,PL,VOL,IDL,SURF,IQUAD,IL,IM)
         CALL LCMPUT(IPTRK,'U',NLF,2,UU)
         CALL LCMPUT(IPTRK,'W',NLF,2,WW)
         CALL LCMPUT(IPTRK,'ALPHA',NLF,2,ALPHA)
         CALL LCMPUT(IPTRK,'PLZ',NSCT,2,PLZ)
         CALL LCMPUT(IPTRK,'PL',NSCT*NLF,2,PL)
         CALL LCMPUT(IPTRK,'SURF',LX+1,2,SURF)
         CALL LCMPUT(IPTRK,'XXX',LX+1,2,XXX)
         CALL LCMPUT(IPTRK,'IL',NSCT,1,IL)
         CALL LCMPUT(IPTRK,'IM',NSCT,1,IM)
         DEALLOCATE(SURF,PL,PLZ,ALPHA,WW,UU,IL,IM)
         LL4=LX*NSCT
         NUN=LL4
      ELSE IF((ITYPE.EQ.5).OR.(ITYPE.EQ.6).OR.(ITYPE.EQ.8)) THEN
*        2D GEOMETRIES: CARTESIAN; TUBE; HEXAGONAL
         NDIM=2
         IF(ISCAT.EQ.0) CALL XABORT('SNTRK: SCAT NOT DEFINED.')
         IF(IQUAD.GE.10) THEN
            NPQ=NLF**2  
         ELSE
            NPQ=(NLF+4)*NLF/2
         ENDIF
         IF(IGLK.NE.0) THEN
            NPQ=(NLF+2)*NLF/2
            NSCT=NPQ
         ELSE
            NSCT=ISCAT*(ISCAT+1)/2
         ENDIF
         IGE=0
*
         IF(ITYPE.EQ.5) THEN
            ! 2D Cartesian
            IF(NFOU.GT.0)THEN
               XLEN = XXX(LX+1)
               YLEN = YYY(LY+1)
               CALL LCMPUT(IPTRK,'XXX',LX+1,2,XXX)
               CALL LCMPUT(IPTRK,'XLEN',1,2,XLEN)

               CALL LCMPUT(IPTRK,'YYY',LY+1,2,YYY)
               CALL LCMPUT(IPTRK,'YLEN',1,2,YLEN)
            ENDIF
         ELSEIF(ITYPE.EQ.6) THEN
            ! 2D Tube
            IGE=1
         ELSEIF(ITYPE.EQ.8) THEN
            ! 2D Hexagonal
            IGE=2
            NHEX=LX/(3*ISPLH**2)
            ALLOCATE(COORDMAP(3,NHEX))
            COORDMAP(:,:)=0
            CALL SNTSFH(IMPX,ITYPE,NHEX,LZ,MCELL,ISPLH,MAT,LOZSWP,
     >         COORDMAP)
            CALL LCMPUT(IPTRK,'LOZSWP',3*6,1,LOZSWP)
            CALL LCMPUT(IPTRK,'COORDMAP',3*NHEX,1,COORDMAP)
            CALL LCMPUT(IPTRK,'SIDE',1,2,SIDE)
            DEALLOCATE(COORDMAP)
         ENDIF
*
         ALLOCATE(MRMX(NPQ),MRMY(NPQ))
         ALLOCATE(DU(NPQ),DE(NPQ),WW(NPQ),DB(LX*NPQ),DA(LX*LY*NPQ),
     1   DAL(LX*LY*NPQ),WX(IELEM+1),
     2   WE(EELEM+1),CST(MAX(IELEM,EELEM)),MN(NSCT*NPQ),DN(NPQ*NSCT),
     3   IL(NSCT),IM(NSCT),PL(NSCT*NPQ))
         CALL SNTT2D(IGE,IMPX,LX,LY,SIDE,IELEM,NLF,NPQ,NSCT,IQUAD,
     1   NCODE,ZCODE,MAT,XXX,YYY,VOL,IDL,DU,DE,WW,MRMX,MRMY,DB,DA,DAL,
     2   PL,LL4,NUN,EELEM,WX,WE,CST,IBFP,ISCHM,ESCHM,IGLK,MN,DN,IL,IM,
     3   ISCAT)
         CALL LCMPUT(IPTRK,'DU',NPQ,2,DU)
         CALL LCMPUT(IPTRK,'DE',NPQ,2,DE)
         CALL LCMPUT(IPTRK,'W',NPQ,2,WW)
         CALL LCMPUT(IPTRK,'MRM',NPQ,1,MRMX)
         CALL LCMPUT(IPTRK,'MRMY',NPQ,1,MRMY)
         CALL LCMPUT(IPTRK,'DB',LX*NPQ,2,DB)
         CALL LCMPUT(IPTRK,'DA',LX*LY*NPQ,2,DA)
         IF(IGE.EQ.1) CALL LCMPUT(IPTRK,'DAL',LX*LY*NPQ,2,DAL)
         CALL LCMPUT(IPTRK,'PL',NSCT*NPQ,2,PL)
         CALL LCMPUT(IPTRK,'WX',IELEM+1,2,WX)
         IF(IBFP.NE.0) CALL LCMPUT(IPTRK,'WE',EELEM+1,2,WE)
         CALL LCMPUT(IPTRK,'CST',MAX(IELEM,EELEM),2,CST)
         CALL LCMPUT(IPTRK,'MN',NSCT*NPQ,2,MN)
         CALL LCMPUT(IPTRK,'DN',NPQ*NSCT,2,DN)
         CALL LCMPUT(IPTRK,'IL',NSCT,1,IL)
         CALL LCMPUT(IPTRK,'IM',NSCT,1,IM)
         CALL LCMPUT(IPTRK,'PL',NSCT*NPQ,2,PL)
         DEALLOCATE(DAL,DA,DB,WW,DE,DU,WX,WE,CST,MN,DN,IL,IM,PL)
         DEALLOCATE(MRMY,MRMX)
      ELSE IF((ITYPE.EQ.7).OR.(ITYPE.EQ.9)) THEN
*        3D GEOMETRIES: CARTESIAN; HEXAGONAL
         NDIM=3
         IF(ISCAT.EQ.0) CALL XABORT('SNTRK: SCAT NOT DEFINED.')
         IF(IQUAD.GE.10) THEN
            NPQ=2*NLF**2
         ELSE
            NPQ=(NLF+2)*NLF
         ENDIF
         IF(IGLK.NE.0) THEN
            NSCT=NPQ
         ELSE
            NSCT=ISCAT*(ISCAT+1)
         ENDIF
         IGE=0
*
         IF(ITYPE.EQ.9) THEN
            ! 3D Hexagonal
            IGE=2
            NHEX     =LX/(3*ISPLH**2)
            ALLOCATE(COORDMAP(3,NHEX))
            COORDMAP(:,:)=0
            CALL SNTSFH(IMPX,ITYPE,NHEX,LZ,MCELL,ISPLH,MAT,LOZSWP,
     >         COORDMAP)
            CALL LCMPUT(IPTRK,'LOZSWP',3*6,1,LOZSWP)
            CALL LCMPUT(IPTRK,'COORDMAP',3*NHEX,1,COORDMAP)
            CALL LCMPUT(IPTRK,'SIDE',1,2,SIDE)
            DEALLOCATE(COORDMAP)
         ENDIF
*
         ALLOCATE(MRMX(NPQ),MRMY(NPQ),MRMZ(NPQ))
         ALLOCATE(DU(NPQ),DE(NPQ),DZ(NPQ),WW(NPQ),DC(LX*LY*NPQ),
     1   DB(LX*LZ*NPQ),DA(LY*LZ*NPQ),WX(IELEM+1),
     2   WE(EELEM+1),CST(MAX(IELEM,EELEM)),MN(NSCT*NPQ),DN(NPQ*NSCT),
     3   IL(NSCT),IM(NSCT),PL(NSCT*NPQ))
         CALL SNTT3D(IGE,IMPX,LX,LY,LZ,SIDE,IELEM,NLF,NPQ,NSCT,IQUAD,
     1   NCODE,ZCODE,MAT,XXX,YYY,ZZZ,VOL,IDL,DU,DE,DZ,WW,MRMX,MRMY,
     2   MRMZ,DC,DB,DA,PL,LL4,NUN,EELEM,WX,WE,CST,IBFP,ISCHM,ESCHM,
     3   IGLK,MN,DN,IL,IM,ISCAT)
*
         CALL LCMPUT(IPTRK,'DU',NPQ,2,DU)
         CALL LCMPUT(IPTRK,'DE',NPQ,2,DE)
         CALL LCMPUT(IPTRK,'DZ',NPQ,2,DZ)
         CALL LCMPUT(IPTRK,'W',NPQ,2,WW)
         CALL LCMPUT(IPTRK,'MRMX',NPQ,1,MRMX)
         CALL LCMPUT(IPTRK,'MRMY',NPQ,1,MRMY)
         CALL LCMPUT(IPTRK,'MRMZ',NPQ,1,MRMZ)
         CALL LCMPUT(IPTRK,'DC',LX*LY*NPQ,2,DC)
         CALL LCMPUT(IPTRK,'DB',LX*LZ*NPQ,2,DB)
         CALL LCMPUT(IPTRK,'DA',LY*LZ*NPQ,2,DA)   
         CALL LCMPUT(IPTRK,'PL',NSCT*NPQ,2,PL)
         CALL LCMPUT(IPTRK,'WX',IELEM+1,2,WX)
         IF(IBFP.NE.0) CALL LCMPUT(IPTRK,'WE',EELEM+1,2,WE)
         CALL LCMPUT(IPTRK,'CST',MAX(IELEM,EELEM),2,CST)
         CALL LCMPUT(IPTRK,'MN',NSCT*NPQ,2,MN)
         CALL LCMPUT(IPTRK,'DN',NPQ*NSCT,2,DN)
         CALL LCMPUT(IPTRK,'IL',NSCT,1,IL)
         CALL LCMPUT(IPTRK,'IM',NSCT,1,IM)
         CALL LCMPUT(IPTRK,'PL',NSCT*NPQ,2,PL)
         DEALLOCATE(DA,DB,DC,WW,DZ,DE,DU,WX,WE,CST,MN,DN,IL,IM,PL)
         DEALLOCATE(MRMZ,MRMY,MRMX)
      ELSE
         CALL XABORT('SNTRK: UNKNOWN GEOMETRY.')
      ENDIF
      DEALLOCATE(YYY,ZZZ)
*----
*  THE NUMBER OF UNKNOWNS OF A BOLTZMANN-FOKKER-PLANCK DISCRETIZATION IS
*  INCREASED TO HOLD SLOWING-DOWN ANGULAR FLUXES.
*----
      IF(IBFP.GT.0) THEN
        IF(NDIM.EQ.1) THEN
          NUN=NUN+IELEM*NLF*NEL
        ELSE IF(NDIM.EQ.2) THEN
          NUN=NUN+IELEM**2*NPQ*NEL
        ELSE IF(NDIM.EQ.3) THEN
          NUN=NUN+IELEM**3*NPQ*NEL
        ELSE
         CALL XABORT('SNTRK: FOKKER-PLANCK NOT IMPLEMENTED.')
        ENDIF
      ENDIF
      IF(IMPX.GT.0) WRITE (6,'(/33H SNTRK: ORDER OF LINEAR SYSTEMS =,
     1 I10/8X,37HNUMBER OF UNKNOWNS PER ENERGY GROUP =,I10)') LL4,NUN
*
      IF(IMPX.GT.5) THEN
         I1=1
         DO I=1,(NEL-1)/8+1
            I2=I1+7
            IF(I2.GT.NEL) I2=NEL
            WRITE (6,620) (J,J=I1,I2)
            WRITE (6,630) (MAT(J),J=I1,I2)
            WRITE (6,640) (IDL(J),J=I1,I2)
            WRITE (6,650) (VOL(J),J=I1,I2)
            I1=I1+8
         ENDDO
      ENDIF
*----
*  SYNTHETIC ACCELERATION TRACKING INFORMATION.
*----
      IF(LDSA) THEN
         IF(IMPX.GT.0) WRITE (6,'(/32H SNTRK: SYNTHETIC ACCELERATION I,
     1   19HNFORMATION FOLLOWS:)')
         CALL LCMSIX(IPTRK,'DSA',1)
         ICOL=3  ! Gauss-Legendre quadrature
         NLFSA=2 ! P1 method
         ISPN=1  ! simplified PN method
         ISCSA=1 ! isotropic scattering
         NVD=1   ! SN-type VOID boundary conditions
         NADI=2  ! ADI iteration
         ICHX=2  ! Raviart-Thomas finite elements
         ISEG=0  ! scalar algorithm
         IMPV=0  ! print parameter for vector operations
         IF(MAXPTS.EQ.0) CALL XABORT('SNTRK: MAXPTS NOT DEFINED.')
         IF((ITYPE.EQ.7).OR.(ITYPE.EQ.9)) THEN
            CALL TRITRK(MAXPTS,IPTRK,IPGEOM,IMPX,IELEMSA,ICOL,ICHX,ISEG,
     1      IMPV,NLFSA,NVD,ISPN,ISCSA,NADI)
         ELSE
            IF(ISOLVSA.EQ.1)THEN
               CALL BIVTRK(MAXPTS,IPTRK,IPGEOM,IMPX,IELEMSA,ICOL,
     1            NLFSA,NVD,ISPN,ISCSA)
            ELSEIF(ISOLVSA.EQ.2)THEN
               CALL TRITRK(MAXPTS,IPTRK,IPGEOM,IMPX,IELEMSA,ICOL,ICHX,
     1            ISEG,IMPV,NLFSA,NVD,ISPN,ISCSA,NADI)
            ELSE
               CALL XABORT('SNTRK: UNDEFINED SOLVER OPTION FOR '
     1            //'SYNTHETIC ACCELERATION.')
            ENDIF
         ENDIF
         CALL LCMSIX(IPTRK,' ',2)
      ENDIF
*----
*  SAVE GENERAL AND SN-SPECIFIC TRACKING INFORMATION.
*----
      IGP(:NSTATE)=0
      IGP(1)=NEL
      IGP(2)=NUN
      IF(ILK) THEN
         IGP(3)=0
      ELSE
         IGP(3)=1
      ENDIF
      IGP(4)=ISTATE(7)
      IGP(5)=1
      IGP(6)=ITYPE
      IGP(7)=NSCT
      IGP(8)=IELEM
      IGP(9)=NDIM
      IGP(10)=ISCHM
      IGP(11)=LL4
      IGP(12)=LX
      IGP(13)=LY
      IGP(14)=LZ
      IGP(15)=NLF
      IGP(16)=ISCAT
      IGP(17)=IQUAD
      IGP(18)=0
      IF(LFIXUP) IGP(18)=1
      IGP(19)=0
      IF(LDSA) IGP(19)=1
      IGP(20)=NSTART
      IGP(21)=NSDSA
      IGP(22)=MAXIT
      IF(LIVO) IGP(23)=1
      IGP(24)=ICL1
      IGP(25)=ICL2
      IGP(26)=ISPLH
      IGP(27)=INSB
      IGP(28)=MCELL
      IF((ITYPE.EQ.3).OR.(ITYPE.GE.4)) IGP(29)=1
      IGP(30)=0
      IF(LSHOOT) IGP(30)=1
      IGP(31)=IBFP
      IGP(32)=NMPI
      IGP(33)=ISOLVSA
      IGP(34)=NFOU
      IGP(35)=EELEM
      IGP(36)=ESCHM
      IGP(37)=IGLK
      IF(LBIHET) IGP(40)=1

      CALL LCMPUT(IPTRK,'STATE-VECTOR',NSTATE,1,IGP)
      CALL LCMPUT(IPTRK,'MATCOD',NEL,1,MAT)
      CALL LCMPUT(IPTRK,'VOLUME',NEL,2,VOL)
      CALL LCMPUT(IPTRK,'KEYFLX',NEL,1,IDL)
      CALL LCMPUT(IPTRK,'NCODE',6,1,NCODE)
      CALL LCMPUT(IPTRK,'ICODE',6,1,ICODE)
      CALL LCMPUT(IPTRK,'ZCODE',6,2,ZCODE)
      CALL LCMPUT(IPTRK,'EPSI',1,2,EPSI)
      IF(ITYPE.EQ.4) CALL LCMPUT(IPTRK,'XXX',LX+1,2,XXX)
      DEALLOCATE(XXX)
*----
*  SET KEYFLX$ANIS
*----
      NLIN=IELEM**NDIM*EELEM
      ALLOCATE(KEYANI(NEL,NLIN,NSCT))
      DO IR=1,NEL
         IND=IDL(IR)
         DO IE=1,NLIN
            DO P=1,NSCT
               IF(IND.EQ.0) THEN
                  KEYANI(IR,IE,P)=0
               ELSE
                  KEYANI(IR,IE,P)=IND+(P-1)*NLIN+IE-1
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      CALL LCMPUT(IPTRK,'KEYFLX$ANIS',NEL*NLIN*NSCT,1,KEYANI)
      DEALLOCATE(KEYANI)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IDL,VOL,MAT)
      RETURN
*
  600 FORMAT(/25H SNTRK: VALUES OF VECTOR ,A6,4H ARE/(1X,1P,20I6))
  620 FORMAT (///11H REGION    ,8(I8,6X,1HI))
  630 FORMAT (   11H MIXTURE   ,8(I8,6X,1HI))
  640 FORMAT (   11H POINTER   ,8(I8,6X,1HI))
  650 FORMAT (   11H VOLUME    ,8(1P,E13.6,2H I))
      END
