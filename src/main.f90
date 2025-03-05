program main

  use iso_fortran_env, only: dp => real64
  use mod_parallel
  implicit none

  ! Parameters
  integer, parameter :: L = 12, Nm = 10**3, Nterm = 500, Nt = 10, Nskip = 10
  real(dp), parameter :: Tmin = 0.1_dp, Tmax = 4.0_dp, DT = Tmax - Tmin
  
  ! Iterators
  integer :: i, isweeps,iskip, it

  ! Arrays
  integer  :: ip(L), im(L)
  real(dp) :: E(Nm), M(Nm)
  real(dp) :: T(Nt), beta(Nt)

  !Parallel variables
  ! Coarrays
  integer :: spin(L,L)[*]
  integer, allocatable :: spin_local(:,:)[:]

  integer :: is, ie, indices(2), tile_size
  integer :: ils, ile, ims, ime, in(2), left, right
  
  
  ! Periodic Boundary conditions
  ip = [(i+1, i = 1, L)] ; ip(L) = 1
  im = [(i-1, i = 1, L)] ; im(1) = L

  ! Temperature array
  T = [(Tmin + DT*i/(Nt - 1), i = 0, Nt-1)]
  beta = 1/T

  ! Open file 
  if(this_image() == 1) open(unit = 10, file = "datosP.dat")

  ! Parallel
  indices = tile_indices(L)
  is = indices(1)
  ie = indices(2)

  in = neighbors()
  left = in(1)
  right = in(2)
  tile_size = L/num_images()
  ils = 1
  ile = tile_size
  ims = ils - 1
  ime = ile + 1
  print*, is, ie, left, right
  allocate(spin_local(ims:ime,L)[*])

  ! Cold Start
  spin_local = 1
  sync all
  
  ! Simulation loop
  temperature : do it = 1, Nt
     ! Thermalization
     do isweeps = 1, Nterm
        call checkerboard_sweeps(spin_local,beta(it))
     end do

     ! Measurements loop
     do isweeps = 1, Nm
        
        do iskip = 1, Nskip
           call checkerboard_sweeps(spin_local,beta(it))
        end do
        
        ! Take measurments
        spin(is:ie,:)[1] = spin_local(1:tile_size,:)
        sync all
        if(this_image() == 1)then
           E(isweeps) = energy_density(spin(:,:))
           M(isweeps) = 1.0_dp*abs(sum(spin(:,:)))/L**2
        end if
     end do
     ! Print measurements
     if(this_image() == 1) write(10,*) T(it), avr(E), avr(M)
     !sync all
  end do temperature
contains

  function avr(x)
    real(dp), intent(in) :: x(:)
    real(dp) :: avr

    avr = sum(x)/size(x)
  end function avr

  subroutine sweeps(spin,beta)
    integer, intent(inout) :: spin(:,:)
    real(dp), intent(in) ::  beta
    integer :: i, j

    do i = 1, L
       do j = 1, L
          !call metropolis(spin,[i,j],beta)
       end do
    end do
  end subroutine sweeps

  subroutine checkerboard_sweeps(spinlocal,beta)
    integer, intent(inout) :: spinlocal(ims:ime,L)[*]
    real(dp), intent(in) ::  beta
    integer :: i, j

    spinlocal(tile_size+1,:)[left]  = spinlocal(1,:)
    spinlocal(0,:)[right] = spinlocal(tile_size,:)
    sync all

    ! Update white squares
    do i = 1, tile_size
       do j = 1, L
          if(mod(i+j,2) == 0) call metropolis(spinlocal,[i,j],beta)
       end do
    end do
    sync all
    spinlocal(tile_size+1,:)[left]  = spinlocal(1,:)
    spinlocal(0,:)[right] = spinlocal(tile_size,:)
    sync all

    ! Update black squares
    do i = 1, tile_size
       do j = 1, L
          if(mod(i+j,2) /= 0) call metropolis(spinlocal,[i,j],beta)
       end do
    end do
    sync all
    
  end subroutine checkerboard_sweeps
  

  subroutine metropolis(spinlocal,x,beta)
    integer, intent(inout) :: spinlocal(0:tile_size+1,L)
    real(dp), intent(in) :: beta
    integer, intent(in) :: x(2)
    integer :: DH
    real(dp) :: r

    DH = DE(spinlocal,x)

    if( DH <= 0 )then
       spinlocal(x(1),x(2)) = -spinlocal(x(1),x(2))
    else
       call random_number(r)
       if( r <= exp(-DH*beta)) spin(x(1),x(2)) = -spin(x(1),x(2))
    end if
    
  end subroutine metropolis

  function DE(spinlocal,x)
    integer, intent(in) :: spinlocal(0:tile_size+1,L)
    integer, intent(in) :: x(2)
    integer :: DE, i,j

    i = x(1)
    j = x(2)

    DE = 2*spin(i,j)*(spin(i+1,j) + spin(i,ip(j))+spin(i-1,j) + spin(i,im(j)))
    
  end function DE

  function energy_density(spin)
    integer, intent(in) :: spin(L,L)
    real(dp) :: energy_density
    integer :: E, i, j

    E = 0
    do i = 1, L
       do j = 1, L
          E = E + spin(i,j) * (spin(ip(i),j) + spin(i,ip(j)))
       end do
    end do
    energy_density = -real(E,dp)/L**2
        
  end function energy_density
  
end program main
