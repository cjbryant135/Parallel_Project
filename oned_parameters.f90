MODULE oned_parameters

CONTAINS


PURE FUNCTION alpha(mu, M)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: M
  DOUBLE PRECISION, DIMENSION(M), INTENT(IN) :: mu
  DOUBLE PRECISION, DIMENSION(M) :: alpha

  alpha(:) = COS(mu(:)) !set alpha boundary condition here!

END FUNCTION alpha

FUNCTION beta(mu, M)
  IMPLICIT NONE
  INTEGER :: M
  DOUBLE PRECISION, DIMENSION(M) :: mu
  DOUBLE PRECISION, DIMENSION(M) :: beta

  beta(:) = SIN(mu(:)) !set beta boundary condition here!

END FUNCTION beta

FUNCTION albedo()
  IMPLICIT NONE
  DOUBLE PRECISION :: albedo 
                !Set albedo here!
  albedo = 0.6d0

END FUNCTION albedo

FUNCTION h(mu1, mu2) 
  IMPLICIT NONE

  DOUBLE PRECISION :: mu1, mu2, h

  h = 0.5d0 !Set h function here

END FUNCTION


FUNCTION tau()
  IMPLICIT NONE
  DOUBLE PRECISION :: tau 
  tau = 0.4d0 !Set tau (the point you want to evaluate at) here!

END FUNCTION tau

FUNCTION taufinal()
  IMPLICIT NONE
  DOUBLE PRECISION :: taufinal 
  taufinal = 1.d0 !Set taufinal!

END FUNCTION taufinal


END MODULE oned_parameters


