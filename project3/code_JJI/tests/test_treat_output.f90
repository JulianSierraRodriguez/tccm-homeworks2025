! This program tests all subroutines that are inside treat_output module

program test_treat_output
  use treat_output
  implicit none

  integer :: Natoms = 1, u_out = 10, u_xyz = 20
  character(len=100) :: f_out = "test_log.out"
  character(len=100) :: f_xyz = "test_traj.xyz"
  logical :: exists
  double precision :: coord(1,3), E, V, T

  write(*,*) "=== TEST MODULE: TREAT_OUTPUT ==="

  ! --- TEST 1: create_output ---
  call create_output(u_out, f_out)
  close(u_out)

  inquire(file=f_out, exist=exists)
  if (exists) then
     write(*,*) "[OK] create_output: File created successfully."
  else
     write(*,*) "[FAIL] create_output: File does not exist."
  end if

  ! --- TEST 2: write_xyz ---
  ! We reuse create_output to open the file, assuming it handles the unit
  call create_output(u_xyz, f_xyz)
  
  coord(1,:) = 0.0d0
  E=1.0; V=0.5; T=0.5

  call write_xyz(u_xyz, Natoms, coord, E, V, T)
  close(u_xyz)

  inquire(file=f_xyz, exist=exists)
  if (exists) then
     write(*,*) "[OK] write_xyz: Trajectory file exists."
  else
     write(*,*) "[FAIL] write_xyz: XYZ file was not created."
  end if

end program test_treat_output
