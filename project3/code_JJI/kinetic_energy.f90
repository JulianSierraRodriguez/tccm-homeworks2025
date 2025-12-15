! This module contains 1 function that calculates total kinetic energy

module kin_energy

  implicit none

contains

  double precision function T(Natoms, velocity, mass)
    implicit none
    integer, intent(in) :: Natoms
    double precision, intent(in) :: velocity(Natoms,3)
    double precision, intent(in) :: mass(Natoms)
    integer :: i
    double precision :: v2
    
    T = 0.0d0 ! initialize variable

    ! Loop over all atoms
    do i = 1,Natoms
      ! First compute magnitude of velocity
      v2 = velocity(i,1)**2 + velocity(i,2)**2 + velocity(i,3)**2
      ! Then compute kinetic energy T
      T = T + mass(i) * v2
    end do

    ! Finally multiply by the factor 1/2
    T = 0.5d0*T 

  end function T

end module kin_energy
