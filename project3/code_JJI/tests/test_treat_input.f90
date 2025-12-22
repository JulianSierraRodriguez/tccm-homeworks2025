! This program tests all subroutines and functions that are inside treat_input module

program test_treat_input

  use treat_input
  implicit none
  character(len=100) ::  fname = "build/tests/bin/inp_test.txt"
  integer :: Natoms
  double precision, allocatable :: coord(:,:)
  double precision, allocatable :: mass(:)
  double precision, allocatable :: distance(:,:)
  
  write(*,*)
  write(*,*) "=== MODULE TEST: TREAT_INPUT ==="

  ! TEST 1: read_Natoms

  Natoms = read_Natoms(fname)
  
  ! Verification:
  if (Natoms == 2) then
     write(*,*) "[OK] read_Natoms: Read 2 atoms successfully."
  else
     write(*,*) "[FAIL] read_Natoms: Read", Natoms, "and was expecting 2."
  end if


  ! TEST 2: read_molecule (and unit conversion Angstroms - nm)

  allocate(coord(Natoms,3))
  allocate(mass(Natoms))
  
  call read_molecule(fname, Natoms, coord, mass)
   
  ! Verification: Input 10.0 Angstroms -> Code divide by 10 -> It must be 1.0 nm
  if (abs(coord(2,3) - 1.0d0) < 1.0d-6) then
     write(*,*) "[OK] read_molecule: Conversion A -> nm was correct (10.0 -> 1.0)."
  else
     write(*,*) "[FAIL] read_molecule: Conversion failed. Valor:", coord(2,3)
  end if


  ! TEST 3: compute_distances

  allocate(distance(Natoms, Natoms))
  call compute_distances(Natoms, coord, distance)
   
  ! Verification: Distance between (0,0,0) y (0,0,1) must be 1.0
  if (abs(distance(1,2) - 1.0d0) < 1.0d-6) then
     write(*,*) "[OK] compute_distances: Distance calculated successfully."
  else
     write(*,*) "[FAIL] compute_distances: Wrong distance:", distance(1,2)
  end if

  write(*,*) ! Separation line

  ! Memory clean
  deallocate(coord, mass, distance)
  
end program test_treat_input
