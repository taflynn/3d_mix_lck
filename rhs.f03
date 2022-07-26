!> \\file rhs.f03
module rhs
 
  use fft

  implicit none

  contains

  ! generate the initial form of the wavefunction 
  subroutine V_rhs(psi,mu,dt,Nx,Ny,Nz)
    implicit none
    
    integer, intent(in) :: Nx, Ny, Nz
    double precision, intent(in) :: mu
    double complex, intent(in) :: dt

    complex(C_DOUBLE_COMPLEX), intent(inout) :: psi(:,:,:)

    ! local variables
    integer :: i, j, k
    
    ! iterate across x, y and z of the two-body and Lee-Huang-Yang terms
    !$omp parallel do collapse(3)
    do k = 1, Nz
      do j = 1, Ny
        do i = 1, Nx
          psi(i,j,k) = psi(i,j,k)*exp(-0.5*dt*(-3.0*abs(psi(i,j,k))**2 + 2.5*abs(psi(i,j,k))**3 - mu))
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine V_rhs

  subroutine T_rhs(psi_k,dk2,Nx,Ny,Nz)
    implicit none

    integer, intent(in) :: Nx, Ny, Nz
    complex(C_DOUBLE_COMPLEX), intent(in) :: dk2(:,:,:)

    complex(C_DOUBLE_COMPLEX), intent(inout) :: psi_k(:,:,:)

    integer :: i, j, k
    
    ! iterate across kx, ky and kz of the kinetic energy term
    !$omp parallel do collapse(3)
    do k = 1, Nz
      do j = 1, Ny
        do i = 1, Nx
          psi_k(i,j,k) = dk2(i,j,k)*psi_k(i,j,k)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine T_rhs
end module rhs
