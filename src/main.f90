program main

  use mod_parallel
  use iso_fortran_env, only : dp => real64, i4 => int32
  implicit none

  !parameters
  integer(i4), parameter :: n = 100
  integer(i4), parameter :: N_thermalization = 500, N_measurements = 1000, Nskip = 10
  
  !variables
  integer(i4), allocatable :: spin(:,:)[:], spin_local(:,:)[:]
  integer(i4), dimension(n) :: ip, im
  real(dp) :: E(N_measurements), M(N_measurements)
  
  !iterators
  integer(i4) :: i, j, i_sweeps, j_sweeps, i_beta, iskip
  integer(i4) :: is, ie, indices(2), tile_size, ils, ile,ims, ime, in(2), left, right
  
  !other
  real(dp) :: r
  real(dp), dimension(10) :: beta, T


  do i = 1, n
     ip(i) = i + 1
     im(i) = i - 1
  end do
  ip(n) = 1
  im(1) = n
  
  indices = tile_indices(n)
  is = indices(1)
  ie = indices(2)

  in = neighbors()
  left = in(1)
  right = in(2)
  sync all
  !print*, this_image(), in

  tile_size = n/num_images()
  ils = 1
  ile = tile_size
  ims = ils - 1
  ime = ile + 1
  
  T = [(1.0_dp + 4.0_dp*i_beta/(size(beta) - 1), i_beta = 0, size(beta) - 1)]
  
  ! allocate variables
  allocate(spin_local(ims:ime,n)[*])
  allocate(spin(n,n)[*])

  !Hot start
  do i = ils, ile
     do j = 1, n
        call random_number(r)
        if( r <= 0.5_dp)then
           spin_local(i,j) = 1
        else
           spin_local(i,j) = -1
        end if
     end do
  end do
  spin_local(ime,:)[left] = spin_local(ils,:)
  spin_local(ims,:)[right] = spin_local(ile,:)
  spin(is:ie,:)[1] = spin_local(ils:ile,:)
  sync all
  !if (this_image() == 1) print*, abs(sum(spin))/real(n**2,dp)
  ! Thermalization
  do i_beta = 1, size(beta)
     do i_sweeps = 1, N_thermalization
        do i = 1, ile
           do j = 1, n
              if(mod(i+j,2) == 0) call metropolis(spin_local,[i,j],1/T(i_beta))
           end do
        end do
        spin_local(ime,:)[left] = spin_local(ils,:)
        spin_local(ims,:)[right] = spin_local(ile,:)
        spin(is:ie,:)[1] = spin_local(ils:ile,:)
        sync all
        
        do i = 1, ile
           do j = 1, n
              if (mod(i+j,2) /= 0) call metropolis(spin_local,[i,j],1/T(i_beta))
           end do
        end do
        spin_local(ime,:)[left] = spin_local(ils,:)
        spin_local(ims,:)[right] = spin_local(ile,:)
        spin(is:ie,:)[1] = spin_local(ils:ile,:)
        sync all
        !if (this_image() == 1) print*, abs(sum(spin))/real(n**2,dp)!, this_image()
     end do

     do i_sweeps = 1, N_measurements
        do iskip = 1, Nskip 
           do i = 1, ile
              do j = 1, n
                 if(mod(i+j,2) == 0) call metropolis(spin_local,[i,j],1/T(i_beta))
              end do
           end do
           spin_local(ime,:)[left] = spin_local(ils,:)
           spin_local(ims,:)[right] = spin_local(ile,:)
           spin(is:ie,:)[1] = spin_local(ils:ile,:)
           sync all
           
           do i = 1, ile
              do j = 1, n
                 if (mod(i+j,2) /= 0) call metropolis(spin_local,[i,j],1/T(i_beta))
              end do
           end do
           spin_local(ime,:)[left] = spin_local(ils,:)
           spin_local(ims,:)[right] = spin_local(ile,:)
           spin(is:ie,:)[1] = spin_local(ils:ile,:)
           sync all
        end do
        if(this_image() == 1)then
           M(i_sweeps) = abs(sum(spin))/real(n**2,dp)
           E(i_sweeps) = energy(spin)
        end if
        sync all
     end do
     
     if( this_image() == 1) then
        print*, T(i_beta), avr(E), avr(M)
     end if
     sync all
  end do
        


contains


  function avr(x)
    REAL(dp) :: avr
    real(dp), intent(in) :: x(:)

    avr = sum(x)/size(x)
  end function avr
  
  subroutine metropolis(spin,x,beta)
    integer(i4), intent(inout) :: spin(0:,:)
    integer(i4), intent(in) :: x(2)
    real(dp), intent(in) :: beta
    real(dp) :: r, p

    call random_number(r)
    p = min(1.0_dp,exp(-DS(spin,x,beta)))
    if ( r <= p ) spin(x(1),x(2)) = -spin(x(1),x(2))
  end subroutine metropolis

  function DS(spin,x,beta)
    integer(i4), intent(inout) ::  spin(0:,:)
    integer(i4), intent(in) :: x(2)
    real(dp), intent(in) :: beta
    real(dp) :: DS

    DS = beta*2*spin(x(1),x(2)) * ( spin(x(1)+1,x(2)) + spin(x(1)-1,x(2)) &
                                  + spin(x(1),ip(x(2))) + spin(x(1),im(x(2))) )
  end function DS


  function energy(spin)
    integer(i4), intent(in) :: spin(:,:)
    real(dp) :: energy
    integer(i4) :: i, j

    energy = 0.0_dp
    do i = 1, n
       do j = 1, n
          energy = energy + spin(i,j)*(spin(ip(i),j) + spin(i,ip(j)))
       end do
    end do
    energy = -energy/n**2
  end function energy
end program main
