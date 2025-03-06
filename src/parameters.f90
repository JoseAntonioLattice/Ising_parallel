module parameters
  use iso_fortran_env, only : dp => real64
  implicit none

  ! Parameters

  !Lattice
  integer :: L

  !Simulation
  integer :: N_measurements
  integer :: N_thermalization
  integer :: Nskip

  !Temperature
  integer :: Nt
  real(dp):: Tmin, Tmax, DT

  namelist/params/ L, N_measurements, N_thermalization, Nskip, Nt, Tmin, Tmax

contains
  
  subroutine read_input()

    character(99) :: filename
    integer :: parameters_unit

    !write(*,'(a)',advance = 'no') 'Type the parameters file: '
    !read(*,'(a)') filename
    !write(*,'(a)') trim(filename)

    filename = 'input/parameters.nml'
    open(newunit = parameters_unit, file = trim(filename))
    read(parameters_unit, nml = params)
    close(parameters_unit)
    if(this_image() == 1) write(*,nml = params) 

    DT = Tmax - Tmin
    
  end subroutine read_input

  
end module parameters
