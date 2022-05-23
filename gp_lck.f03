!> \\file gp_lck.f03
program gp_lck
  use OMP_LIB
  use FFTW3
  use grid
  use init
  use time

  implicit none

  complex(C_DOUBLE_COMPLEX), allocatable :: psi(:,:,:)

  integer :: T_STEPS = 10000
  integer :: IM_REAL = 0

  integer :: Nx, Ny, Nz
  double precision :: dx, dy, dz 
  double precision :: dt, dt_coef

  double precision :: Nlck = 370.0

  double precision, allocatable :: x(:), y(:), z(:)
  double precision, allocatable :: kx(:), ky(:), kz(:)

  double precision, allocatable :: dk2(:,:,:)

  Nx = 64
  Ny = 64
  Nz = 64

  dx = 0.5
  dy = 0.5
  dz = 0.5  

  dt_coef = 0.1

  ! initialise time-step
  dt = dt_coef*min(dx,dy,dz)**2

  ! set up 3D spatial grid
  write(*,*) "setting up the Cartesian space grid"
  x = space_grid(Nx,dx)
  y = space_grid(Ny,dy)
  z = space_grid(Nz,dz)

  ! set up 3D momentum space grid
  write(*,*) "setting up the momentum space grid"
  kx = mom_grid(Nx,dx)
  ky = mom_grid(Ny,dy)
  kz = mom_grid(Nz,dz)

  ! compute initial profile of wavefunction
  psi = init_wav(x,y,z)
 
  call renorm(psi,dx,dy,dz,Nlck,Nx,Ny,Nz)

  dk2 = exp_lap(kx,ky,kz,dt)

  ! begin time-stepping
  call ssfm(psi,dk2,T_STEPS,dt,dx,dy,dz,Nlck,IM_REAL)

end program
