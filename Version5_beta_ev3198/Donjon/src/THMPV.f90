SUBROUTINE THMPV(MFLXT, SPEED, POULET, VCOOL, DCOOL, &
                DCOOL0, PCOOL, ACOOL, MUT, XFL, HD, DV, NZ, TCOOL, HZ)
!
!-----------------------------------------------------------------------
!
!Purpose:
! Update the pressure and velocity vectors in the THM model to model the pressure drop
! and the velocity of the fluid in the channel
!
!Copyright:
! Copyright (C) 1993 Ecole Polytechnique de Montreal
!
!Author(s): C. Huet
! 03/02/2025: C. Huet - Creation of the model
!
!Parameters: input
! MFLXT   mass flow rate of the fluid in the channel
! XFL     quality of the fluid in the channel
! DCOOL   density of the fluid in the channel
! SPEED   inlet velocity of the fluid in the channel
! POULET  pressure drop in the channel
! VCOOL   velocity of the fluid in the channel
! DCOOL0  density of the fluid at the inlet
! PCOOL   pressure of the fluid in the channel
! ACOOL   acceleration of the fluid in the channel
! MUT     dynamic viscosity of the fluid in the channel
! HD      hydraulic diameter of the channel
! DV      volume of a fluid element in the channel
! NZ      number of nodes in the channel
!
!Parameters: output
! VCOOL   velocity of the fluid in the channel
! PCOOL   pressure of the fluid in the channel
!
!-----------------------------------------------------------------------
    
    USE GANLIB
    IMPLICIT NONE

!----
!   SUBROUTINE ARGUMENTS
!----
    INTEGER NZ
    REAL MFLXT(NZ), SPEED, POULET, VCOOL(NZ), DCOOL(NZ), DCOOL0, PCOOL(NZ), ACOOL, MUT(NZ), XFL(NZ), HD, DV, TCOOL(NZ)
    REAL HZ(NZ)
!----
!   LOCAL VARIABLES
!----

    REAL g
    !REAL A(2*NZ, 2*NZ+1)
    REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: A
    
    INTEGER K, I, J, IER
    REAL PHIL0, TPMULT, TPMULT0
    REAL RET, RET0, FRIC, FRIC0, DZ

    g = - 9.81
    PRINT *, 'HD', HD
    ALLOCATE(A(2*NZ,2*NZ+1))
    FORALL (I=1:2*NZ, J=1:2*NZ+1) A(I, J) = 0.0
    PRINT *, 'DCOOL', DCOOL
    PRINT *, 'VCOOL', VCOOL
    PRINT *, 'PCOOL', PCOOL
    
    PRINT *, "Matrice A :"
    DO I = 1, 2*NZ
        PRINT *, (A(I,J), J=1,2*NZ+1)
    END DO

    DO K = 1, NZ
        IF (K .EQ. 1) THEN

            RET0 = ABS(VCOOL(K)*DCOOl(K)) * (1.0 - XFL(K)) * HD / MUT(K)
            RET  = ABS(VCOOL(K+1)*DCOOl(K+1)) * (1.0 - XFL(K+1)) * HD / MUT(K+1)

            CALL THMFRI(RET, FRIC, HD)
            CALL THMFRI(RET0, FRIC0, HD)



            IF (XFL(K) .GT. 0.0) THEN
                CALL THMPLO(PCOOL(K), XFL(K), PHIL0)
                TPMULT0 = PHIL0
                CALL THMPLO(PCOOL(K+1), XFL(K+1), PHIL0)
                TPMULT = PHIL0
            ELSE
                TPMULT = 1.0
                TPMULT0 = 1.0
            ENDIF

            !FRIC = 0.002
            !FRIC0 = 0.002

            A(1,1) = 1.0
            A(K+NZ,K) = - (VCOOL(K)*DCOOL(K))*(1.0 - (TPMULT0*FRIC0*HZ(K))/(2.0*HD))
            ! Mult par HZ(K) et chg signe (base + et + mtn: + -)
            A(K+NZ,K+1) = (VCOOL(K+1)*DCOOL(K+1))*(1.0 + (TPMULT*FRIC*HZ(K))/(2.0*HD))

            A(1, 2*NZ+1) = SPEED
            A(K+NZ, 2*NZ+1) = 0! - ((DCOOL(K+1) - DCOOL(K)) * g) /2

            A(K+NZ,K-1+NZ) = 0.0
            A(K+NZ,K+NZ) = -1.0
            A(K+NZ,K+1+NZ) = 1.0

        ELSE IF (K .EQ. NZ) THEN
            A(K,K-1) = - DCOOL(K-1)
            A(K,K) = DCOOL(K)
            
            A(K, 2*NZ+1) = 0.0
            A(2*NZ, 2*NZ+1) = POULET
            A(2*NZ, 2*NZ) = 1.0
        
        ELSE 
            RET = ABS(VCOOL(K+1)*DCOOL(K+1)) * (1.0 - XFL(K+1)) * HD / MUT(K+1)
            RET0 = ABS(VCOOL(K)*DCOOL(K)) * (1.0 - XFL(K)) * HD / MUT(K)

            CALL THMFRI(RET, FRIC, HD)
            CALL THMFRI(RET0, FRIC0,HD)

            !FRIC = 0.002
            !FRIC0 = 0.002

            IF (XFL(K) .GT. 0.0) THEN
                CALL THMPLO(PCOOL(K+1), XFL(K+1), PHIL0)
                TPMULT = PHIL0
                CALL THMPLO(PCOOL(K), XFL(K), PHIL0)
                TPMULT0 = PHIL0
            ELSE
                TPMULT = 1.0
                TPMULT0 = 1.0
            ENDIF

            A(K,K-1) = - DCOOL(K-1)
            A(K,K) = DCOOL(K)
            A(K,K+1) = 0.0
            A(K, 2*NZ+1) = 0.0

            PRINT *, "K", K
            PRINT *, "(DCOOL(K+1)- DCOOL(K)) * g", (DCOOL(K+1)- DCOOL(K)) * g
            PRINT *, "A(K+NZ,K)", (DCOOL(K)*VCOOL(K))*(1.0 - (TPMULT0*FRIC0*HZ(K))/(2.0*HD))
            PRINT *, "A(K+NZ,K+1)", (DCOOL(K+1)*VCOOL(K+1))*(1.0 + (TPMULT*FRIC*HZ(K))/(2.0*HD))
            PRINT *, "DCOOL(K+1)", DCOOL(K+1)
            PRINT *, "DCOOL(K)", DCOOL(K)
            PRINT *, "VCOOL(K+1)", VCOOL(K+1)
            PRINT *, "VCOOL(K)", VCOOL(K)
            PRINT *, "TPMULT", TPMULT
            PRINT *, "TPMULT0", TPMULT0
            PRINT *, "FRIC", FRIC
            PRINT *, "FRIC0", FRIC0

            A(K+NZ,K) = - (DCOOL(K)*VCOOL(K))*(1.0 - (TPMULT0*FRIC0*HZ(K))/(2.0*HD))
            ! Mult par HZ(K) et chg signe (base + et + mtn: + -)
            A(K+NZ,K+1) = (DCOOL(K+1)*VCOOL(K+1))*(1.0 + (TPMULT*FRIC*HZ(K))/(2.0*HD))

            A(K+NZ, 2*NZ+1) = 0!- ((DCOOL(K+1) - DCOOL(K)) * g) /2

            A(K+NZ,K-1+NZ) = 0.0
            A(K+NZ,K+NZ) = -1.0
            A(K+NZ,K+1+NZ) = 1.0
        ENDIF
    END DO

    PRINT *, "Matrice A avant ALSBD:"
    DO I = 1, 2*NZ
        PRINT *, (A(I,J), J=1,2*NZ+1)
    END DO

  ! Appel de ALSBD
    call ALSBD(2*NZ, 1, A, IER, 2*NZ)

    ! Affichage de la solution X
    print *, "Solution du système :"
    do I = 1, 2*NZ
        print *, A(I, 2*NZ+1)  ! Les solutions sont dans la dernière colonne
    end do

        ! Vérification d'erreur
    if (IER /= 0) then
        print *, "Erreur : matrice singulière !"
        stop
      end if

    DO K = 1, NZ
        VCOOL(K) = A(K, 2*NZ+1)
        PCOOL(K) = A(K+NZ, 2*NZ+1)
    END DO

    PRINT *, 'INSIDE THMPV'
    PRINT *, 'VCOOL', VCOOL
    PRINT *, 'DCOOL', DCOOL
    PRINT *, 'FRIC', FRIC
    PRINT *, 'TPMULT', TPMULT
    PRINT *, 'HD', HD
    PRINT *, 'HZ', HZ
    PRINT *, 'g', g

    DO K=1, NZ
        PRINT *, 'K', K
        PRINT *, 'a', ((DCOOL(K)*VCOOL(K)**2)*((TPMULT*FRIC*HZ(K))/(HD)) + (DCOOL(K)-DCOOL(K+1))*g/2)
    END DO

    RETURN
    END