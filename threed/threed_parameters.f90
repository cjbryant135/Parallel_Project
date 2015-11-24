MODULE threed_parameters
CONTAINS

FUNCTION phase(mu_a, phi_x, mu_b, phi_y) !DEFINE Phase function here
  IMPLICIT NONE
  DOUBLE PRECISION :: mu_a, phi_x, mu_b, phi_y, phase

  phase = mu_a+phi_x+mu_b+phi_y !This is just a random one I came up with. no physical meaning


END FUNCTION phase

!Parameters
FUNCTION x_max() !Remember we assume x_max = y_max
  IMPLICIT NONE
  DOUBLE PRECISION :: x_max
  
  x_max = 1.d0

END FUNCTION x_max

FUNCTION N() !number of grid points for the x and y grids
  IMPLICIT NONE
  INTEGER :: N
  
  N = 100

END FUNCTION N

FUNCTION pi()
  IMPLICIT NONE
  DOUBLE PRECISION :: pi

  pi = 3.141592653589793d0


END FUNCTION pi

FUNCTION N_mu() !Number of grid points for the mu grid
  IMPLICIT NONE
  INTEGER :: N_mu

  N_mu = 32 !Note: must be even

END FUNCTION N_mu


FUNCTION N_phi() !Number of phi grid points
  IMPLICIT NONE
  INTEGER :: N_phi 
   
  N_phi = 64 !must be even

END FUNCTION N_phi









END MODULE threed_parameters

