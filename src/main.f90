program main

  use iso_fortran_env, only: dp => real64, mp => int32
  use parameters
  use mod_parallel
  implicit none

  ! Iterators
  integer :: i,j,k,k2,isweeps,iskip,it

  ! Other integers
  integer :: ounit, Nw
  
  ! Arrays
  integer, allocatable  :: ip1(:), im1(:), ip2(:), im2(:), p(:,:)
  integer, allocatable  :: black(:,:), white(:,:)  
  real(dp), allocatable :: T(:), beta(:)

  !Parallel variables

  ! Coarrays
  integer(mp), allocatable :: spin(:,:)[:]
  real(dp), allocatable, codimension[:] :: E(:), M(:)
  
  integer :: is, ie, indices(2), tile_sizex, tile_sizey
  integer :: ils, ile, ims, ime, in(2), left, right, up, down, col, row
  

  call read_input()

  call random_init(.false.,.true.)

  ! Allocate variables
  allocate(ip1(d(1)),im1(d(1)))
  allocate(ip2(d(2)),im2(d(2)))
  allocate(p(d(1),d(2)))
  allocate(E(N_measurements)[*],M(N_measurements)[*])
  allocate(T(Nt),beta(Nt))

  !periodic boundary conditions 
  ip1 = [(i+1,i=1,d(1))]; ip1(d(1)) = 1
  ip2 = [(i+1,i=1,d(2))]; ip2(d(2)) = 1

  im1 = [(i-1,i=1,d(1))]; im1(1) = d(1)
  im2 = [(i-1,i=1,d(2))]; im2(1) = d(2)

  if(mod(this_image(),d(1)) == 0) then
     row = d(1)
  else
     row = mod(this_image(),d(1))
  end if

  if(mod(this_image(),d(1)*d(2)) /= 0)then
     col = (mod(this_image(),d(1)*d(2)) - row)/d(1)+1
  else
     col = (d(1)*d(2) - row)/d(1)+1
  end if

  !print*,this_image(), row, col

  do i = 1, d(1)
     do j = 1, d(2)
        p(i,j) = 2*(j-1) + i
     end do
  end do

  left  = p(im1(row),col)
  right = p(ip1(row),col)
  down  = p(row,im2(col))
  up    = p(row,ip2(col))


  ! Temperature array
  T = [(Tmin + DT*i/(Nt - 1), i = 0, Nt-1)]
  beta = 1/T

  ! Open file 
  if(this_image() == 1) open(newunit = ounit, file = "data/datosP.dat")

  ! Parallel
  tile_sizex = L(1)/d(1)
  tile_sizey = L(2)/d(2)

  allocate(spin(0:tile_sizex+1,0:tile_sizey+1)[*])
  ! Number of white squares
  Nw = tile_sizex * tile_sizey / 2
  allocate(white(Nw,2),black(Nw,2))
  
  ! Black and white coordinates
  k = 0; k2 = 0
  do i = 1, tile_sizex
     do j = 1, tile_sizey
        if( mod(i+j,2) == 0) then
           k = k + 1
           white(k,:) = [i,j]
        else
           k2 = k2 + 1
           black(k2,:) = [i,j]
        end if
     end do
  end do
  
  ! Cold Start
  spin = 1
  sync all
  
  !if (this_image() == 4) print*, "Cold start ok", up, down ,left, right
  
  ! Simulation loop
  temperature : do it = 1, Nt
     ! Thermalization
     do isweeps = 1, N_thermalization
        call checkerboard_sweeps(spin,beta(it))
     end do
     !if (this_image() == 1) print*, "Thermalization ok"
  
     ! Measurements loop
     do isweeps = 1, N_measurements
        do iskip = 1, Nskip
           call checkerboard_sweeps(spin,beta(it))
        end do
        
        ! Take measurements
        E(isweeps) = energy_density(spin(:,:))
        M(isweeps) = real(sum(spin(1:tile_sizex,1:tile_sizey)),dp)
        sync all
        call co_sum(E(isweeps),result_image = 1)
        call co_sum(M(isweeps),result_image = 1)
     end do
     !if (this_image() == 1) print*, it," temperature ok"
     
     
     ! Print measurements
     if(this_image() == 1) then
        
        M = abs(M/(L(1)*L(2)))
        !print*, T(it), avr(E)!, stderr(E), avr(M), stderr(M)
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
    integer(mp), intent(inout) :: spin(0:,0:)[*]
    real(dp), intent(in) ::  beta
    integer :: i

    ! Update white squares
    do i = 1, Nw
       call metropolis(spin,[white(i,1),white(i,2)],beta)
    end do
    !if (this_image() == 1) print*, "White squares updated"
  
    spin(tile_sizex+1,:)[left]  = spin(1,:)
    spin(0,:)[right] = spin(tile_sizex,:)
    spin(:,tile_sizey+1)[up]  = spin(:,1)
    spin(:,0)[down] = spin(:,tile_sizey)

    sync all
    
    ! Update black squares
    do i = 1, Nw
       call metropolis(spin,[black(i,1),black(i,2)],beta)
    end do
    spin(tile_sizex+1,:)[left]  = spin(1,:)
    spin(0,:)[right] = spin(tile_sizex,:)
    spin(:,tile_sizey+1)[up]  = spin(:,1)
    spin(:,0)[down] = spin(:,tile_sizey)

    sync all
    
  end subroutine checkerboard_sweeps
  

  subroutine metropolis(spin,x,beta)
    integer(mp), intent(inout) :: spin(0:,0:)
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

  pure function DE(spin,x)
    integer(mp), intent(in) :: spin(0:,0:)
    integer, intent(in) :: x(2)
    integer :: DE, i,j

    i = x(1)
    j = x(2)

    DE = 2*spin(i,j)*(spin(i+1,j) + spin(i,j+1) + spin(i-1,j) + spin(i,j-1))
    
  end function DE

  pure function energy_density(spin)
    integer(mp), intent(in) :: spin(0:,0:)
    real(dp) :: energy_density
    integer :: E, i, j

    E = 0
    do i = 1, tile_sizex
       do j = 1, tile_sizey
          E = E + spin(i,j) * (spin(i+1,j) + spin(i,j+1))
       end do
    end do
    energy_density = -real(E,dp)/(L(1)*L(2))
        
  end function energy_density
  
end program main
