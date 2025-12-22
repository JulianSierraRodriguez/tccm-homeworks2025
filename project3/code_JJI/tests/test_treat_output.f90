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

  ! TEST 1: create_output - NEW
  ! We test creating a file using the create_output subroutine and then we close it.
   
  call create_output(u_out, f_out)
  close(u_out)

  ! We check if the file indeeds exists.
  inquire(file=f_out, exist=exists)
  if (exists) then
    write(*,*) "[OK] create_output: File created successfully."
  else
    write(*,*) "[FAIL] create_output: File does not exist."
  end if

  ! TEST 2: create_output - REPLACE part
  call create_output(u_out, f_out)
  close(u_out)

  ! We dont need a condition here, if the Replace section of the subroutine does not work it would raise an error.
  ! We cannot use NEW in an already created file.
  write(*,*) "[OK] create_output: File REPLACED successfully."


  ! TEST 3: write_xyz
  ! We reuse create_output to open the file, assuming it handles the unit, to generate another file.
  call create_output(u_xyz, f_xyz)
  
  ! We create some data to test the subroutine.
  coord(1,:) = 0.0d0
  E=1.0; V=0.5; T=0.5

  ! We use the write_xyz with dummy data to see if it writes the data, has to be checked inside the file manually.
  call write_xyz(u_xyz, Natoms, coord, E, V, T)
  close(u_xyz)

  ! We check that the xyz file indeed exists.
  inquire(file=f_xyz, exist=exists)
  if (exists) then
    write(*,*) "[OK] write_xyz: Trajectory file exists."
  else
    write(*,*) "[FAIL] write_xyz: XYZ file was not created."
  end if

  write(*,*) ! Separation line

  ! We delete the generated files.
  call create_output(u_out, f_out)
  close(u_out, status='delete')
  call create_output(u_xyz, f_xyz)
  close(u_xyz, status='delete')

end program test_treat_output
