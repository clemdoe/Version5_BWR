*DECK THMFRI
      SUBROUTINE THMFRI(REY,FRIC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the value of the friction factor coefficient with the Muller
* Steinhagen correlation formula
*
*Copyright:
* Copyright (C) 2013 Ecole Polytechnique de Montreal.
*
*Author(s): 
* P. Gallet
*
*Parameters: input
* REY     reynolds number 
*
*Parameters: output
* FRIC    friction factor coefficient
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      REAL REY,FRIC
*----
*  COMPUTE VALUE OF THE FRICTION FACTOR COEFFICIENT AS FUNCTION OF THE 
*  REYNOLDS NUMBER
*----

      FRIC = 0.002

      RETURN
      END
