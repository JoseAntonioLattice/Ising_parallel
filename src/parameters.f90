module parameters
  use iso_fortran_env, only : dp => real64
  implicit none

  ! Parameters

  !Lattice
  integer :: L[*]

  !Simulation
  integer :: N_measurements[*]
  integer :: N_thermalization[*]
  integer :: Nskip[*]

  !Temperature
  integer :: Nt[*]
  real(dp):: Tmin[*], Tmax[*], DT[*]

  namelist/params/ L, N_measurements, N_thermalization, Nskip, Nt, Tmin, Tmax

contains
  
  subroutine read_input()

    character(99) :: filename
    integer :: parameters_unit

    if (this_image() == 1) then
       write(*,'(a)',advance = 'no') 'Type the parameters file: '
       read(*,'(a)') filename
       write(*,'(a)') trim(filename)
    
       
       !filename = 'input/parameters.nml'
       open(newunit = parameters_unit, file = trim(filename))
       read(parameters_unit, nml = params)
       close(parameters_unit)
       write(*,nml = params) 
       DT = Tmax - Tmin
    end if

    call co_broadcast(L,source_image = 1)
    call co_broadcast(N_measurements,source_image = 1)
    call co_broadcast(N_thermalization,source_image = 1)
    call co_broadcast(Nskip,source_image = 1)
    call co_broadcast(Nt,source_image = 1)
    call co_broadcast(Tmax, source_image = 1)
    call co_broadcast(Tmin, source_image = 1)
    call co_broadcast(DT, source_image = 1)
    
  end subroutine read_input

  
end module parameters
