! Ths program test independently Kinetic (T), Potential (V), and Total (E) energy calculations

program test_energies
  use energies
  implicit none

  integer, parameter :: Natoms = 2
  double precision :: vel(Natoms,3), mass(Natoms), dist(Natoms,Natoms)
  double precision :: epsilon=1.0, sigma=1.0
  double precision :: res_T, res_V, res_E

  write(*,*) "=== TEST MODULE: ENERGIES ==="

  ! --- TEST 1: Function T (Kinetic) ---
  ! For a single atom:
  ! m=2, v=1 -> T = 0.5*2*1^2 = 1.0
  vel = 0.0d0; vel(1,1) = 1.0d0
  mass = 0.0d0; mass(1) = 2.0d0
  res_T = T(Natoms, vel, mass)
  
  if (abs(res_T - 1.0d0) < 1e-6) then
     write(*,*) "[OK] Function T: Correct."
  else
     write(*,*) "[FAIL] Function T: Value", res_T, "expected 1.0."
  end if

  ! --- TEST 2: Function V (LJ Potential) ---
  ! r = sigma -> V must be 0.0
  dist = 0.0d0; dist(1,2) = 1.0d0 ! equals sigma
  res_V = V(epsilon, sigma, Natoms, dist)

  if (abs(res_V) < 1e-6) then
     write(*,*) "[OK] Function V: Correct (Zero at r=sigma)."
  else
     write(*,*) "[FAIL] Function V: Value", res_V, "expected 0.0."
  end if

  ! --- TEST 3: Function E (Total) ---
  res_E = E(10.0d0, 5.0d0)
  if (abs(res_E - 15.0d0) < 1e-6) then 
     write(*,*) "[OK] Function E: Sum correct."
  end if

end program test_energies
