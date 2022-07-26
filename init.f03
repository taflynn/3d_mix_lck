!> \\file init.f03
module init

  use fft

  implicit none
  
  contains

  ! generate the initial form of the wavefunction 
  function init_wav(x,y,z,init_type,gauss_sig)

    double precision, intent(in) :: x(:), y(:), z(:)
    integer, intent(in) :: init_type    
    double precision, intent(in) :: gauss_sig

    complex(C_DOUBLE_COMPLEX), allocatable :: init_wav(:,:,:)
   
    ! Local variables 
    integer :: Nx, Ny, Nz
    integer :: i, j, k
   
    Nx = size(x)
    Ny = size(y)
    Nz = size(z)

    allocate(init_wav(Nx,Ny,Nz))
    
    ! define initial wavefunction profile as a Gaussian
    if (init_type == 1) then
      write(*,*) "initial input to imaginary time: Gaussian"
      do k = 1, Nz
        do j = 1, Ny
          do i = 1, Nx
            init_wav(i,j,k) = exp(-(x(i)**2.0/(2.0*gauss_sig**2.0) &
                                  + y(j)**2.0/(2.0*gauss_sig**2.0) &
                                  + z(k)**2.0/(2.0*gauss_sig**2.0)))**0.5
          end do
        end do
      end do
    ! define initial wavefunction profile as a Super-Gaussian
    elseif (init_type == 2) then
      write(*,*) "initial input to imaginary time: Super-Gaussian"
      do k = 1, Nz
        do j = 1, Ny
          do i = 1, Nx
            init_wav(i,j,k) = exp(-(x(i)**2.0/(2.0*gauss_sig**2.0) &
                                  + y(j)**2.0/(2.0*gauss_sig**2.0) &
                                  + z(k)**2.0/(2.0*gauss_sig**2.0))**3.0)**0.5
          end do
        end do
      end do
    end if

  end function init_wav

end module init
