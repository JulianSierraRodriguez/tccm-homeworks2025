! This module contains 1 function and 2 subroutines to read input file

module read_input

  implicit none

contains

  ! 1 - Function that reads number of atoms in the input file

  integer function read_Natoms(input_file) result(Natoms)
    implicit none
    character(len=*), intent(in) :: input_file
    ! integer :: Natoms !!! esta linea chatGPT dice que hay que quitarla, revisar
    
    open(16,file=input_file, status='old')  ! opens the file
    read(16,*)Natoms ! read the first line of the file 
    close(16)   ! close unit

  end function read_Natoms


  ! 2 - Subroutine that reads coordinates and masses from the input file

  subroutine read_molecule(input_file, Natoms, coord, mass)
    implicit none
    character(len=*), intent(in) :: input_file
    integer, intent(in) :: Natoms
    double precision, intent(out) :: coord(Natoms,3)
    double precision, intent(out) :: mass(Natoms)
    integer :: i

    open(16,file=input_file, status='old')
    read(16,*) ! Pass first line (Natoms)
    ! loop that reads the lines of the file to put the data into the array
    do i = 1,Natoms
      read(16,*)coord(i,1), coord(i,2), coord(i,3), mass(i)
    end do
    close(16)

    ! Convert coordinates from Angstroms to nanometers
    coord = coord/10
    

  end subroutine read_molecule



  ! 3 - Subroutine that calculate internuclear distances

  subroutine compute_distances(Natoms, coord, distance)
    implicit none
    integer, intent(in) :: Natoms
    double precision, intent(in) :: coord(Natoms,3)
    double precision, intent(out) :: distance(Natoms,Natoms)
    integer :: i, j
    

    distance = 0.0d0 ! initialize variable

    ! the code uses a nested loop, that will use an atom, and all the rest
    ! without a repetition of pairs (only upper diagonal part of the matrix is filled)
    do i = 1,Natoms-1
      do j = i+1,Natoms
        ! computes the distance between the pairs of atoms
        distance(i,j) = sqrt(((coord(j,1)-coord(i,1))**2)+((coord(j,2)-coord(i,2))**2)+((coord(j,3)-coord(i,3))**2))
      end do
    end do

  end subroutine compute_distances

end module read_input
