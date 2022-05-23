!> \\file time.f03
module time
  use OMP_LIB
  use FFTW3
  use rhs

  implicit none

  contains

  ! main ssfm time-stepping function
  subroutine ssfm(psi,dk2,T_STEPS,dt,dx,dy,dz,Nlck,kx,ky,kz)
    implicit none

    integer, intent(in) :: T_STEPS

    double precision, intent(in) :: dx, dy, dz, dt
    double precision, intent(in) :: Nlck

    double precision, intent(in) :: dk2(:,:,:)
    
    double precision, intent(in) :: kx(:), ky(:), kz(:)

    complex(C_DOUBLE_COMPLEX), intent(inout) :: psi(:,:,:)
   
    ! local variables 
    integer :: i
    integer :: Nx, Ny, Nz

    double precision :: mu
    
    type(C_PTR) :: plan_forw, plan_back
    
    complex(C_DOUBLE_COMPLEX), allocatable :: psi_k(:,:,:)
    
    integer(C_INT) :: nthreads
    integer(C_INT) :: error

    Nx = size(psi, dim = 1)
    Ny = size(psi, dim = 2)
    Nz = size(psi, dim = 3)
 
    allocate(psi_k(Nx,Ny,Nz))

    error = fftw_init_threads()
    nthreads = omp_get_max_threads()
    call fftw_plan_with_nthreads(int(nthreads,C_INT))

    plan_forw = fftw_plan_dft_3d(Nz,Ny,Nx,psi,psi_k,FFTW_FORWARD,FFTW_ESTIMATE)
    plan_back = fftw_plan_dft_3d(Nz,Ny,Nx,psi_k,psi,FFTW_BACKWARD,FFTW_ESTIMATE)
    mu = -1.0

    do i = 1, T_STEPS

      ! first half-step
      call V_rhs(psi,mu,dt)
      
      ! FFT wavefunction to momentum space
      call fftw_execute_dft(plan_forw,psi,psi_k)
      psi_k = psi_k/sqrt(dble(Nx*Ny*Nz))
      
      ! kinetic energy step
      call T_rhs(psi_k,dk2)
      
      ! IFFT wavefunction to real space
      call fftw_execute_dft(plan_back,psi_k,psi)
      psi = psi/sqrt(dble(Nx*Ny*Nz))

      ! second half-step
      call V_rhs(psi,mu,dt)

      ! renormalise wavefunction
      call renorm(psi,dx,dy,dz,Nlck)

      if (mod(i,T_STEPS/10) == 0) then
        ! FFT wavefunction to real space
        call fftw_execute_dft(plan_forw,psi,psi_k)
        psi_k = psi_k/sqrt(dble(Nx*Ny*Nz))
        mu = chem_pot(psi,psi_k,kx,ky,kz,plan_back)
        write(*,*) "Central Density"
        write(*,*) abs(psi(Nx/2,Ny/2,Nz/2))**2
      end if
    end do

  end subroutine ssfm

  subroutine renorm(psi,dx,dy,dz,Nlck)
    implicit none

    double precision, intent(in) :: dx, dy, dz
    double precision, intent(in) :: Nlck
   
    complex(C_DOUBLE_COMPLEX), intent(inout) :: psi(:,:,:)

    ! local variables
    double precision :: norm
    
    norm = sum(abs(psi)**2)*dx*dy*dz
    psi = psi*sqrt(Nlck/norm)
  
  end subroutine renorm

  function chem_pot(psi,psi_k,kx,ky,kz,plan_back)
    implicit none
    
    type(C_PTR), intent(in) :: plan_back
    
    double precision, intent(in) :: kx(:), ky(:), kz(:)

    complex(C_DOUBLE_COMPLEX), intent(in) :: psi(:,:,:), psi_k(:,:,:)

    double precision :: chem_pot

    ! local variables
    integer :: i, j, k
    integer :: Nx, Ny, Nz
    
    complex(C_DOUBLE_COMPLEX), allocatable :: lap_psi(:,:,:), lap_psik(:,:,:)
  
    Nx = size(psi, dim = 1)
    Ny = size(psi, dim = 2)
    Nz = size(psi, dim = 3)

    allocate(lap_psi(Nx,Ny,Nz))
    allocate(lap_psik(Nx,Ny,Nz))
 
    !$omp parallel do collapse(3)   
    do k = 1, Nz
      do j = 1, Ny
        do i = 1, Nx
          lap_psik(i,j,k) = -(kx(i)**2 + ky(j)**2 + kz(k)**2)*psi_k(i,j,k)
        end do
      end do
    end do
    !$omp end parallel do
    
    call fftw_execute_dft(plan_back,lap_psik,lap_psi)
    lap_psi = lap_psi/sqrt(dble(Nx*Ny*Nz))
    
    chem_pot = sum(-0.5*conjg(psi)*lap_psi - 3.0*abs(psi)**4 + 2.5*abs(psi)**5)/sum(abs(psi)**2)
  
  end function chem_pot

end module time
