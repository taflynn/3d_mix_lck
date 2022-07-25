!> \\file time.f03
module time
  use OMP_LIB
  use HDF5
  use FFTW3
  use rhs

  implicit none

  contains

  ! main ssfm time-stepping function
  subroutine ssfm(psi,dk2,t_steps,t_save,dt,dx,dy,dz,Nlck,mu,im_real)

    integer, intent(in) :: t_steps, t_save, im_real
    double precision, intent(in) :: dx, dy, dz
    double complex, intent(in) :: dt
    double precision, intent(in) :: Nlck
    complex(C_DOUBLE_COMPLEX), intent(in) :: dk2(:,:,:)
    
    complex(C_DOUBLE_COMPLEX), intent(inout) :: psi(:,:,:)
    double precision, intent(inout) :: mu

    ! local variables 
    integer :: l
    integer :: Nx, Ny, Nz, Nxyz

    type(C_PTR) :: plan_forw, plan_back
    
    complex(C_DOUBLE_COMPLEX), allocatable :: psi_k(:,:,:)
    
    integer(C_INT) :: nthreads
    integer(C_INT) :: error
  
    character(len=4) :: t_current 
    integer :: h5_error 
    integer(HID_T) :: file_id
    integer(HID_T) :: dset_id
    integer(HID_T) :: dspace_id
    integer(HSIZE_T), dimension(3) :: dims
    character(len=19) :: filename_wav

    call h5open_f(h5_error)

    dims = shape(psi)
    Nx = dims(1)
    Ny = dims(2)
    Nz = dims(3)

    Nxyz = Nx*Ny*Nz

    allocate(psi_k(Nx,Ny,Nz))
    
    ! initialising FFTW with threads
    error = fftw_init_threads()
    nthreads = omp_get_max_threads()
    call fftw_plan_with_nthreads(int(nthreads,C_INT))
  
    if (im_real == 0) then 
      write(*,*) "number of threads:" 
      write(*,*) nthreads
    end if 

    ! constructing FFTW plans
    plan_forw = fftw_plan_dft_3d(Nz,Ny,Nx,psi,psi_k,FFTW_FORWARD,FFTW_ESTIMATE)
    plan_back = fftw_plan_dft_3d(Nz,Ny,Nx,psi_k,psi,FFTW_BACKWARD,FFTW_ESTIMATE)

    do l = 1, t_steps

      ! first half-step (non-linear terms)
      call V_rhs(psi,mu,dt,Nx,Ny,Nz)
      
      ! FFT wavefunction to momentum space
      call fftw_execute_dft(plan_forw,psi,psi_k)
      psi_k = psi_k/sqrt(dble(Nxyz))
      
      ! kinetic energy step
      call T_rhs(psi_k,dk2,Nx,Ny,Nz)
      
      ! IFFT wavefunction to real space
      call fftw_execute_dft(plan_back,psi_k,psi)
      psi = psi/sqrt(dble(Nxyz))

      ! second half-step (non-linear terms)
      call V_rhs(psi,mu,dt,Nx,Ny,Nz)

      ! renormalise wavefunction (if in imaginary time)
      if (im_real == 0) then
        call renorm(psi,dx,dy,dz,Nlck)
        !if (mod(l,t_save) == 0) then
          ! FFT wavefunction to real space
          call fftw_execute_dft(plan_forw,psi,psi_k)
          psi_k = psi_k/sqrt(dble(Nxyz))
          mu = chem_pot(psi,psi_k,dk2,plan_back,Nx,Ny,Nz,dt)
        !end if
      end if

      ! data outputting
      if (mod(l,t_save) == 0) then
        !
        write(*,*) "Percentage Completed"
        write(*,*) 100*(dble(l)/dble(T_STEPS))
        write(*,*) "Central Density"
        write(*,*) abs(psi(1+Nx/2,1+Ny/2,1+Nz/2))**2
        write(*,*) "Chemical Potential"
        write(*,*) mu
        ! outputting
        write(t_current,'(I4.4)') 100*l/T_STEPS
        if (im_real == 0) then
          filename_wav = 'psi_im_t_'//trim(t_current)//'.h5'
        elseif (im_real == 1) then
          filename_wav = 'psi_re_t_'//trim(t_current)//'.h5'
        end if
        ! create wavefunction output file
        call h5fcreate_f(filename_wav, H5F_ACC_TRUNC_F, file_id, h5_error)
        ! save real component of wavefunction
        call h5screate_simple_f(3, dims, dspace_id, h5_error)
        call h5dcreate_f(file_id, 'psi_real', H5T_NATIVE_DOUBLE, dspace_id, dset_id, h5_error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, real(psi), dims, h5_error)
        call h5dclose_f(dset_id, h5_error)
        ! save imaginary component of wavefunction
        call h5dcreate_f(file_id, 'psi_imag', H5T_NATIVE_DOUBLE, dspace_id, dset_id, h5_error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, aimag(psi), dims, h5_error)
        call h5dclose_f(dset_id, h5_error)
        call h5sclose_f(dspace_id, h5_error)
        call h5fclose_f(file_id, h5_error)
      end if
    end do
    call h5close_f(h5_error)
  end subroutine ssfm

  subroutine renorm(psi,dx,dy,dz,Nlck)
    
    double precision, intent(in) :: dx, dy, dz
    double precision, intent(in) :: Nlck
   
    complex(C_DOUBLE_COMPLEX), intent(inout) :: psi(:,:,:)

    ! local variables
    double precision :: norm
    
    norm = sum(abs(psi)**2)*dx*dy*dz

    psi = psi*sqrt(Nlck/norm)
  end subroutine renorm

  function chem_pot(psi,psi_k,dk2,plan_back,Nx,Ny,Nz,dt)
    
    type(C_PTR), intent(in) :: plan_back
    
    integer, intent(in) :: Nx, Ny, Nz
    
    double complex, intent(in) :: dk2(:,:,:)

    complex(C_DOUBLE_COMPLEX), intent(in) :: psi(:,:,:)
    complex(C_DOUBLE_COMPLEX) :: psi_k(:,:,:)

    double precision :: chem_pot

    ! local variables
    integer :: i, j, k
   
    double complex :: dt

    double precision, allocatable :: k2(:,:,:)

    complex(C_DOUBLE_COMPLEX), allocatable :: lap_psi(:,:,:), lap_psik(:,:,:)

    allocate(lap_psi(Nx,Ny,Nz))
    allocate(lap_psik(Nx,Ny,Nz))
    allocate(k2(Nx,Ny,Nz))

    ! define the laplacian operator from the exponential laplacian used in the ssfm method
    k2 = (2.0/dt)*log(dk2)

    ! compute second derivative of wavefunction in momentum space
    !$omp parallel do collapse(3)   
    do k = 1, Nz
      do j = 1, Ny
        do i = 1, Nx
          lap_psik(i,j,k) = k2(i,j,k)*psi_k(i,j,k)
        end do
      end do
    end do
    !$omp end parallel do
    
    ! transform kinetic energy term back into real space
    call fftw_execute_dft(plan_back,lap_psik,lap_psi)
    lap_psi = lap_psi/sqrt(dble(Nx*Ny*Nz))
    
    ! compute chemical potential
    chem_pot = sum(-0.5*conjg(psi)*lap_psi - 3.0*abs(psi)**4 + 2.5*abs(psi)**5)/sum(abs(psi)**2)
  
  end function chem_pot

end module time
