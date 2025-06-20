SUBROUTINE THMPV(SPEED, POULET, VCOOL, DCOOL, PCOOL, MUT, XFL, HD, NZ, HZ,EPS, RHOL, RHOG, VGJ, IDFM)
!
!-----------------------------------------------------------------------
!
!Purpose:
! Update the pressure and velocity vectors in the THM model to model the
! pressure drop and the velocity of the fluid in the channel
!
!Copyright:
! Copyright (C) 2025 Ecole Polytechnique de Montreal
!
!Author(s): C. Huet
! 02/2025: C. Huet - Creation
!
!Parameters: input
! XFL     quality of the fluid in the channel
! DCOOL   density of the fluid in the channel
! SPEED   inlet velocity of the fluid in the channel
! POULET  pressure drop in the channel
! VCOOL   velocity of the fluid in the channel
! PCOOL   pressure of the fluid in the channel
! MUT     dynamic viscosity of the fluid in the channel
! HD      hydraulic diameter of the channel
! VGJ     drift velocity in the channel
! EPS     coolant void fraction in the channel
! RHOL    density of the liquid fraction
! RHOG    density of the vapour fraction
!
!Parameters: output
! VCOOL   velocity of the fluid in the channel
! PCOOL   pressure of the fluid in the channel
!
!-----------------------------------------------------------------------
!
    USE GANLIB
    IMPLICIT NONE
!----
!   SUBROUTINE ARGUMENTS
!----
    INTEGER NZ, IDFM
    REAL SPEED, POULET, VCOOL(NZ), DCOOL(NZ), PCOOL(NZ), MUT(NZ), XFL(NZ), HD 
    REAL HZ(NZ),VGJ(NZ),RHOL(NZ), RHOG(NZ), EPS(NZ)
!----
!   LOCAL VARIABLES
!----
    REAL g
    REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: A

    INTEGER K, I, J, IER
    REAL PHIL0, TPMULT, TPMULT0
    REAL REY, REY0, FRIC, FRIC0,DELTA

    g = 9.81 !gravity
    ALLOCATE(A(2*NZ,2*NZ+1))
    FORALL (I=1:2*NZ, J=1:2*NZ+1) A(I, J) = 0.0

    PRINT *, 'INSIDE THMPV VGJ', VGJ 

!----
!   MATRIX FILLING FOR THE PRESSURE AND VELOCITY CALCULATION
!----
!   BOTTOM OF THE CHANNEL
!----
    DO K = 1, NZ 
        IF (K .EQ. 1) THEN
            REY0 = ABS(VCOOL(K)*DCOOl(K)) * (1.0 - XFL(K)) * HD / MUT(K)
            REY  = ABS(VCOOL(K+1)*DCOOl(K+1)) * (1.0 - XFL(K+1)) * HD / MUT(K+1)
            CALL THMFRI(REY, FRIC, HD)
            CALL THMFRI(REY0, FRIC0, HD)

            IF (XFL(K) .GT. 0.0) THEN
                CALL THMPLO(PCOOL(K), XFL(K), PHIL0)
                TPMULT0 = PHIL0
                CALL THMPLO(PCOOL(K+1), XFL(K+1), PHIL0)
                TPMULT = PHIL0
            ELSE
                TPMULT = 1.0
                TPMULT0 = 1.0
            ENDIF
            A(1,1) = 1.0
!   MOMENTUM CONSERVATION EQUATION
            IF (IDFM .GT. 0) THEN
                DELTA = ((EPS(K)/1-EPS(K))*RHOL(K)*RHOG(K)/DCOOL(K)*VGJ(K)**2) - &
            ((EPS(K+1)/1-EPS(K+1))*RHOL(K+1)*RHOG(K+1)/DCOOL(K+1)*VGJ(K+1)**2)
            ELSE
                DELTA = 0
            ENDIF
            PRINT *, 'INSIDE THMPV DELTA', DELTA
            A(K+NZ,K) = - (VCOOL(K)*DCOOL(K))*(1.0 - (TPMULT0*FRIC0*HZ(K))/(2.0*HD))
            A(K+NZ,K+1) = (VCOOL(K+1)*DCOOL(K+1))*(1.0 + (TPMULT*FRIC*HZ(K))/(2.0*HD))
            A(K+NZ, 2*NZ+1) =  - ((DCOOL(K+1)* HZ(K+1) + DCOOL(K)* HZ(K)) * g ) /2 + DELTA
            A(K+NZ,K-1+NZ) = 0.0
            A(K+NZ,K+NZ) = -1.0
            A(K+NZ,K+1+NZ) = 1.0

!    MASS CONSERVATION EQUATION
            A(1, 2*NZ+1) = SPEED

!----
!   TOP OF THE CHANNEL
!----
        ELSE IF (K .EQ. NZ) THEN
!   MASS CONSERVATION EQUATION
            A(K,K-1) = - DCOOL(K-1)
            A(K,K) = DCOOL(K)
!   MOMENTUM CONSERVATION EQUATION
            A(K, 2*NZ+1) = 0.0
            A(2*NZ, 2*NZ+1) = POULET
            A(2*NZ, 2*NZ) = 1.0
!----
!   MIDDLE OF THE CHANNEL
!----
        ELSE
            REY = ABS(VCOOL(K+1)*DCOOL(K+1)) * (1.0 - XFL(K+1)) * HD / MUT(K+1)
            REY0 = ABS(VCOOL(K)*DCOOL(K)) * (1.0 - XFL(K)) * HD / MUT(K)

            CALL THMFRI(REY, FRIC, HD)
            CALL THMFRI(REY0, FRIC0,HD)

            IF (XFL(K) .GT. 0.0) THEN
                CALL THMPLO(PCOOL(K+1), XFL(K+1), PHIL0)
                TPMULT = PHIL0
                CALL THMPLO(PCOOL(K), XFL(K), PHIL0)
                TPMULT0 = PHIL0
            ELSE
                TPMULT = 1.0
                TPMULT0 = 1.0
            ENDIF
!   MASS CONSERVATION EQUATION
            A(K,K-1) = - DCOOL(K-1)
            A(K,K) = DCOOL(K)
            A(K,K+1) = 0.0
            A(K, 2*NZ+1) = 0.0 
!----
!   MOMENTUM CONSERVATION EQUATION  
!----
            IF (IDFM .GT. 0) THEN
                DELTA = ((EPS(K)/1-EPS(K))*RHOL(K)*RHOG(K)/DCOOL(K)*VGJ(K)**2) - &
            ((EPS(K+1)/1-EPS(K+1))*RHOL(K+1)*RHOG(K+1)/DCOOL(K+1)*VGJ(K+1)**2)
            ELSE
                DELTA = 0.0
            ENDIF
            PRINT *, 'INSIDE THMPV DELTA', DELTA
            A(K+NZ,K) = - (DCOOL(K)*VCOOL(K))*(1.0 - (TPMULT0*FRIC0*HZ(K))/(2.0*HD))
            A(K+NZ,K+1) = (DCOOL(K+1)*VCOOL(K+1))*(1.0 + (TPMULT*FRIC*HZ(K))/(2.0*HD))
            A(K+NZ, 2*NZ+1) = - ((DCOOL(K+1)* HZ(K+1) + DCOOL(K)* HZ(K)) * g ) /2 + DELTA
            A(K+NZ,K-1+NZ) = 0.0
            A(K+NZ,K+NZ) = -1.0
            A(K+NZ,K+1+NZ) = 1.0
        ENDIF
    END DO
!----
!   SOLVING THE LINEAR SYSTEM
!----
    call ALSBD(2*NZ, 1, A, IER, 2*NZ)

    if (IER /= 0) CALL XABORT('THMPV: SINGULAR MATRIX.')
!----
!   RECOVER THE PRESSURE AND VELOCITY VECTORS
!----
    DO K = 1, NZ
        VCOOL(K) = REAL(A(K, 2*NZ+1))
        PCOOL(K) = REAL(A(K+NZ, 2*NZ+1))
    END DO

    PRINT *, 'INSIDE THMPV VCOOL', VCOOL
    PRINT *, 'INSIDE THMPV PCOOL', PCOOL
    

    DEALLOCATE(A)
    
    RETURN
    END