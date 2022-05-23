!> \\file init.f03
module init

  use FFTW3

  implicit none
  
  contains

  ! generate the initial form of the wavefunction 
  function init_wav(x,y,z) result(psi)
   
    implicit none

    double precision, intent(in) :: x(:), y(:), z(:)
    complex(C_DOUBLE_COMPLEX), allocatable :: psi(:,:,:)
   
    ! Local variables 
    integer :: Nx, Ny, Nz
    integer :: i,j,k
   
    Nx = size(x)
    Ny = size(y)
    Nz = size(z)

    allocate(psi(Nx,Ny,Nz))

    do k = 1, Nz
      do j = 1, Ny
        do i = 1, Nx
          psi(i,j,k) = exp(-(x(i)**2 + y(j)**2 + z(k)**2))
        end do
      end do
    end do

  end function init_wav

end module init
