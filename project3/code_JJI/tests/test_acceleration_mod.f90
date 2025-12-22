! This program tests acceleration calculation in a simple case with in a qualitative way

program test_acceleration_mod

  use acceleration_mod
  implicit none
  integer, parameter :: Natoms = 2
  double precision :: coord(Natoms,3), mass(Natoms), distance(Natoms,Natoms), acceleration(Natoms,3)

  write(*,*) "=== TEST MODULE: ACCELERATION ==="

  ! 1. Initialize variables (2 atoms)

  coord = 0.0d0       ! Positions to 0
  acceleration = 0.0d0     ! Acceleration to 0
  mass = 1.0d0       ! Mass to 1
  distance = 0.0d0       ! Distances to 0


  ! 2. Setup a dummy scenario

  ! Atom 1 at (0,0,0), Atom 2 at (1,0,0)
  coord(2,1) = 1.0d0 
  distance(1,2) = 1.0d0
  distance(2,1) = 1.0d0


  ! 3. Call the subroutine

  ! (Using epsilon=1.0, sigma=1.0)
  call compute_acc(Natoms, coord, mass, distance, 1.0d0, 1.0d0, acceleration)


  ! 4. Check results

  write(*,*) "Result for Atom 1 (X, Y, Z):", acceleration(1,:)

  ! If acceleration in X is NOT zero, the function worked.
  if (abs(acceleration(1,1)) > 1.0d-6) then
     write(*,*) "[OK] Acceleration calculated successfully."
  else
     write(*,*) "[FAIL] Acceleration is zero (something is wrong)."
  end if

  write(*,*) ! Separation line

end program test_acceleration_mod
