program main

  use iso_fortran_env, only: dp => real64
  use parameters
  use mod_parallel
  implicit none

  ! Iterators
  integer :: i,j, isweeps,iskip, it

  ! Other integers
  integer :: ounit
  
  ! Arrays
  integer, allocatable  :: ip(:), im(:)
  real(dp), allocatable, codimension[:] :: E(:), M(:)
  real(dp), allocatable :: T(:), beta(:)

  !Parallel variables
  ! Coarrays
  integer, allocatable :: spin(:,:)[:]
  
  integer :: is, ie, indices(2), tile_size
  integer :: ils, ile, ims, ime, in(2), left, right
  

  call read_input()

  allocate(ip(L),im(L))
  allocate(E(N_measurements)[*],M(N_measurements)[*])
  allocate(T(Nt),beta(Nt))
  
  ! Periodic Boundary conditions
  ip = [(i+1, i = 1, L)] ; ip(L) = 1
  im = [(i-1, i = 1, L)] ; im(1) = L

  ! Temperature array
  T = [(Tmin + DT*i/(Nt - 1), i = 0, Nt-1)]
  beta = 1/T

  ! Open file 
  if(this_image() == 1) open(newunit = ounit, file = "datosP.dat")

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

  allocate(spin(ims:ime,L)[*])

  ! Cold Start
  spin = 1
  sync all
  
  ! Simulation loop
  temperature : do it = 1, Nt
     ! Thermalization
     do isweeps = 1, N_thermalization
        call checkerboard_sweeps(spin,beta(it))
     end do

     ! Measurements loop
     do isweeps = 1, N_measurements
        do iskip = 1, Nskip
           call checkerboard_sweeps(spin,beta(it))
        end do
        
        ! Take measurements
        E(isweeps) = energy_density(spin(1:ile,:))
        M(isweeps) = 1.0_dp*abs(sum(spin(1:ile,:)))/L**2
        sync all
        call co_sum(E(isweeps),result_image = 1)
        call co_sum(M(isweeps),result_image = 1)
     end do
     
     ! Print measurements
     if(this_image() == 1) then
        write(ounit,*) T(it), avr(E), stderr(E), avr(M), stderr(M)
        flush(ounit)
     end if
     
     sync all
     
  end do temperature
contains

  function avr(x)
    real(dp), intent(in) :: x(:)
    real(dp) :: avr

    avr = sum(x)/size(x)
  end function avr
  
  function var(x)
    real(dp), intent(in) :: x(:)
    real(dp) :: var, avg

    avg = avr(x)
    var = sum((x-avg)**2)/(size(x) - 1)

  end function var

  function stderr(x)
    real(dp), intent(in) :: x(:)
    real(dp) :: stderr

    stderr = sqrt(var(x)/size(x))

  end function stderr
  
  subroutine checkerboard_sweeps(spin,beta)
    integer, intent(inout) :: spin(0:,:)[*]
    real(dp), intent(in) ::  beta
    integer :: i, j

    ! Update white squares
    do i = 1, ile
       do j = 1, L
          if(mod(i+j,2) == 0) call metropolis(spin,[i,j],beta)
       end do
    end do
    spin(ime,:)[left]  = spin(ils,:)
    spin(ims,:)[right] = spin(ile,:)
    sync all

    ! Update black squares
    do i = 1, ile
       do j = 1, L
          if(mod(i+j,2) /= 0) call metropolis(spin,[i,j],beta)
       end do
    end do
    spin(ime,:)[left]  = spin(ils,:)
    spin(ims,:)[right] = spin(ile,:)
    sync all
    
  end subroutine checkerboard_sweeps
  

  subroutine metropolis(spin,x,beta)
    integer, intent(inout) :: spin(0:,:)
    real(dp), intent(in) :: beta
    integer, intent(in) :: x(2)
    integer :: DH
    real(dp) :: r

    DH = DE(spin,x)

    if( DH <= 0 )then
       spin(x(1),x(2)) = -spin(x(1),x(2))
    else
       call random_number(r)
       if( r <= exp(-DH*beta)) spin(x(1),x(2)) = -spin(x(1),x(2))
    end if
    
  end subroutine metropolis

  function DE(spin,x)
    integer, intent(in) :: spin(0:,:)
    integer, intent(in) :: x(2)
    integer :: DE, i,j

    i = x(1)
    j = x(2)

    DE = 2*spin(i,j)*(spin(i+1,j) + spin(i,ip(j)) + spin(i-1,j) + spin(i,im(j)))
    
  end function DE

  function energy_density(spin)
    integer, intent(in) :: spin(:,:)[*]
    real(dp) :: energy_density
    integer :: E, i, j

    E = 0
    do i = 1, ile
       do j = 1, L
          E = E + spin(i,j) * (spin(ip(i),j) + spin(i,ip(j)))
       end do
    end do
    energy_density = -real(E,dp)/L**2
        
  end function energy_density
  
end program main
