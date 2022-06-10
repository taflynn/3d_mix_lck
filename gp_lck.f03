!> \\file gp_lck.f03
program gp_lck
  use OMP_LIB
  use json_module
  use HDF5
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

  integer :: error
  integer(HID_T) :: file_id
  integer(HID_T) :: dset_id
  integer(HID_T) :: dspace_id
  integer(HSIZE_T), dimension(2) :: dims_r
  character(len=7) :: filename_grid = 'grid.h5'

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
 
  ! initiate the hdf5 environment
  call h5open_f(error)

  ! create the grid.h5 file
  call h5fcreate_f(filename_grid, H5F_ACC_TRUNC_F, file_id, error)
  ! saving x-array
  call h5screate_simple_f(1, dims_r, dspace_id, error); if (Nx .ne. Ny .and. Nx .ne. Nz) stop 
  dims_r = size(x)
  call h5dcreate_f(file_id, 'x', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, dims_r, error)
  call h5dclose_f(dset_id, error)
  ! saving y-array
  dims_r = size(y)
  call h5dcreate_f(file_id, 'y', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, y, dims_r, error)
  call h5dclose_f(dset_id, error)
  ! saving z-array
  dims_r = size(z)
  call h5dcreate_f(file_id, 'z', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, z, dims_r, error)
  call h5dclose_f(dset_id, error)
  ! saving kx-array
  dims_r = size(kx)
  call h5dcreate_f(file_id, 'kx', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, kx, dims_r, error)
  call h5dclose_f(dset_id, error)
  ! saving ky-array
  dims_r = size(ky)
  call h5dcreate_f(file_id, 'ky', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, ky, dims_r, error)
  call h5dclose_f(dset_id, error)
  ! saving kz-array
  dims_r = size(kz)
  call h5dcreate_f(file_id, 'kz', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, kz, dims_r, error)
  call h5dclose_f(dset_id, error)
  call h5sclose_f(dspace_id, error)
  ! close the grid.h5 file
  call h5fclose_f(file_id, error)
  
  call h5close_f(error)
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
