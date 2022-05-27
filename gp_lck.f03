!> \\file gp_lck.f03
program gp_lck
  use OMP_LIB
  use json_module
  use FFTW3
  use grid
  use init
  use time

  implicit none

  complex(C_DOUBLE_COMPLEX), allocatable :: psi(:,:,:)

  integer :: Nx, Ny, Nz
  double precision :: dx, dy, dz 
  double precision :: dt_coef
  
  integer :: im_t_steps, re_t_steps
  integer :: im_t_save, re_t_save

  double precision :: Nlck

  integer :: im_real
  double complex :: dt
  double precision :: mu

  double precision, allocatable :: x(:), y(:), z(:)
  double precision, allocatable :: kx(:), ky(:), kz(:)

  complex(C_DOUBLE_COMPLEX), allocatable :: dk2(:,:,:)

  type(json_file) :: json
  
  logical :: is_found

  ! initialising the json_file object
  call json%initialize()

  ! loading in the input file
  call json%load_file('config.json'); if (json%failed()) stop

  ! reading in the input data
  ! read in grid size
  call json%get('Nx', Nx, is_found); if (.not. is_found) stop
  call json%get('Ny', Ny, is_found); if (.not. is_found) stop
  call json%get('Nz', Nz, is_found); if (.not. is_found) stop
  ! read in grid spacing
  call json%get('dx', dx, is_found); if (.not. is_found) stop
  call json%get('dy', dy, is_found); if (.not. is_found) stop
  call json%get('dz', dz, is_found); if (.not. is_found) stop
  ! read in time step size
  call json%get('dt_coef', dt_coef, is_found); if (.not. is_found) stop
  ! read in time step numbers
  call json%get('im_t_steps', im_t_steps, is_found); if (.not. is_found) stop
  call json%get('re_t_steps', re_t_steps, is_found); if (.not. is_found) stop
  ! read in time step saving numbers
  call json%get('im_t_save', im_t_save, is_found); if (.not. is_found) stop
  call json%get('re_t_save', re_t_save, is_found); if (.not. is_found) stop
  ! read in theoretical parameters (effective atom number here)
  call json%get('Nlck', Nlck, is_found); if (.not. is_found) stop

  ! initialise time-step
  dt = dt_coef*min(dx,dy,dz)**2

  ! set up 3D spatial grid
  x = space_grid(Nx,dx)
  y = space_grid(Ny,dy)
  z = space_grid(Nz,dz)

  ! set up 3D momentum space grid
  kx = mom_grid(Nx,dx)
  ky = mom_grid(Ny,dy)
  kz = mom_grid(Nz,dz)

  ! compute initial profile of wavefunction
  psi = init_wav(x,y,z)
 
  call renorm(psi,dx,dy,dz,Nlck)

  ! begin time-stepping
  if (im_t_steps > 0) then
    write(*,*) "beginning imaginary time"
    
    ! state that the time-stepping should expect imaginary time
    im_real = 0
    
    ! initialise arrays and parameters
    dk2 = exp_lap(kx,ky,kz,dt)
    mu = -1.0
    
    ! imaginary time function
    call ssfm(psi,dk2,im_t_steps,im_t_save,dt,dx,dy,dz,Nlck,mu,im_real)
  end if
  if (re_t_steps > 0) then
    write(*,*) "beginning real time"
    
    ! state that the time-stepping should expect real time
    im_real = 1

    ! initilise arrays and parameters (in particular a complex time-step)
    dt = complex(0.0,1.0)*dt
    dk2 = exp_lap(kx,ky,kz,dt)
    
    ! real time function
    call ssfm(psi,dk2,re_t_steps,re_t_save,dt,dx,dy,dz,Nlck,mu,im_real)
  end if
  if (im_t_steps == 0 .and. re_t_steps == 0) then
    ! if there are no time-steps for imaginary and real time then stop program
    stop
  end if
end program
