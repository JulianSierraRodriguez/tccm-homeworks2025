! This program test independently Kinetic (T), Potential (V), and Total (E) energy calculations

program test_energies

  use energies
  implicit none
  integer, parameter :: Natoms = 2
  double precision :: velocity(Natoms,3), mass(Natoms), distance(Natoms,Natoms)
  double precision :: epsilon=1.0, sigma=1.0
  double precision :: T_tot, V_tot, E_tot

  write(*,*) "=== TEST MODULE: ENERGIES ==="

  ! TEST 1: Function T (Kinetic)

  ! For a single atom:
  ! m=2, v=1 -> T = 0.5*2*1^2 = 1.0

  velocity = 0.0d0; velocity(1,1) = 1.0d0
  mass = 0.0d0; mass(1) = 2.0d0
  T_tot = T(Natoms, velocity, mass)
  
  if (abs(T_tot - 1.0d0) < 1e-6) then
     write(*,*) "[OK] Function T: Correct."
  else
     write(*,*) "[FAIL] Function T: Value", T_tot, "expected 1.0."
  end if


  ! TEST 2: Function V (LJ Potential)

  ! If r = sigma -> V must be 0.0

  distance = 0.0d0; distance(1,2) = 1.0d0 ! equals sigma
  V_tot = V(epsilon, sigma, Natoms, distance)

  if (abs(V_tot) < 1e-6) then
     write(*,*) "[OK] Function V: Correct (Zero at r=sigma)."
  else
     write(*,*) "[FAIL] Function V: Value", V_tot, "expected 0.0."
  end if


  ! TEST 3: Function E (Total)

  ! If V_tot=10.0, T_tot=5.0 -> E_tot must be 15.0

  E_tot = E(10.0d0, 5.0d0)
  if (abs(E_tot - 15.0d0) < 1e-6) then 
     write(*,*) "[OK] Function E: Sum correct."
  else
     write(*,*) "[FAIL] Function E: Value", E_tot, "expected 15.0."
  end if

  write(*,*) ! Separation line

end program test_energies
