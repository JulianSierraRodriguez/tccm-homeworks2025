! This module contains 3 subroutines that create and write the output and xyz trajectory

module treat_output

  implicit none

contains

  ! 1 - This subroutine creates or replaces (if it already exists) the output. Not using the correct status creates errors.

  subroutine create_output(u,name_file)
    implicit none
    integer, intent(in) :: u
    character(100), intent(in) :: name_file
    integer :: ios
    logical :: exists

    inquire(file=name_file, exist=exists)
    if (exists) then
    ! if the file already exists we use the status REPLACE, this way we overwrite the alredy existing file.
      open(unit=u, file=name_file, status='REPLACE', action='WRITE', iostat=ios)
    else
    ! if not, we use the status NEW, this way we write on a freshly created file.

      open(unit=u, file=name_file, status='NEW', action='WRITE', iostat=ios)
    end if

    ! Checks that the file was opened correctly using the ios = 0.
    if (ios /= 0) then
      write(10,'(A,I0)') 'Error opening', name_file,' IOSTAT = ', ios
      stop 1
    end if

  end subroutine create_output


  ! Explanation about the formating, after the comma inside write(), used a lot in the the two following subroutines:
  ! A -> string of text
  ! 1X -> one space between what we print, we can change the number for different number of spaces.
  ! I0 -> to write an integer.
  ! FY.Z -> To write a float using a total of Y digits (including the minus sign) with Z digits after the decimal point. 


  ! 2 - Writing a header for the output file with the initial information

  subroutine write_output_head(u_output, input_file, epsilon, sigma, time_step, steps, M, output_file, u_xyz, xyz_file, Natoms)
    implicit none
    integer, intent(in) :: u_output, Natoms, steps, u_xyz, M
    character(100), intent(in) :: input_file, output_file, xyz_file
    double precision :: epsilon, sigma, time_step

    write(u_output,'(A)') "|---------------|"
    write(u_output,'(A)') "|   MD OUTPUT   |"
    write(u_output,'(A)') "|---------------|"
    write(u_output,'()')
    write(u_output,'(A)')"### MD parameters ###"
    write(u_output,'()')
    write(u_output,'(A,1X,A)') "Input Filename =",trim(input_file)
    ! strings too long make the compiler file, we used & to continue writing the same code in the following line
    write(u_output,'(A,1X,I0,1X,A,1X,F5.3,1X,A,1X,I0,1X,A)') "Steps of the simulation =", steps, &  
    "| Time step =", time_step, "ps | Save every",M, "steps"
    write(u_output,'(A,1X,I0,1X,A)') "There are", Natoms, "atoms"
    write(u_output,'()')
    write(u_output,'(A)')"### Output parameters ###"
    write(u_output,'()')
    write(u_output,'(A,1X,A,1X,A,1X,I0)')"Output Filename =", trim(output_file), "| u =", u_output
    write(u_output,'(A,1X,A,1X,A,1X,I0)')"Trajectory Filename =", trim(xyz_file), "| u =", u_xyz
    write(u_output,'()')
    write(u_output,'(A)')"### Lennard-Jones parameters ###"
    write(u_output,'()')
    write(u_output,'(A,1X,F7.4,1X,A,1X,F6.3,1X,A)')"sigma =", sigma, "nm | epsilon =", epsilon, "kJ/mol"
    write(u_output,'()')

  end subroutine write_output_head


  ! 3 - We write the trajectory step in the needed format

  subroutine write_xyz(u_xyz, Natoms, coord, E_total, V_total, T_total)
    implicit none
    integer, intent(in) :: u_xyz, Natoms
    double precision, intent(in) :: coord(Natoms,3)
    double precision, intent(in) :: E_total, V_total, T_total
    integer :: i

    ! We write the number of atoms.
    write(u_xyz,'(I0)') Natoms
    ! We write the Energy results of this step.
    write(u_xyz,'(A,1X,F8.5,1X,A,1X,F8.5,1X,A,1X,F8.5,1X,A)') &
    "E_total = ", E_total,"| V_total = ", V_total,"| T_total = ", T_total, "(kJ/mol)"

    ! Iteration loop that runs over the atoms, to have the coordinates of each atom in a different line
    ! If we write them as write(*,*) coord, we would get everything in a single line.
    do i = 1,Natoms
      write(u_xyz,'(A,1X,F10.7,1X,F10.7,1X,F10.7)')'Ar',coord(i,:)
    end do

  end subroutine write_xyz

end module treat_output

