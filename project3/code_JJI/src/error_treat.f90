! This module contains 1 subroutine that controls if the allocation failed or not

module error_treat

  implicit none

contains

  subroutine error_allocate(i_stat)
    implicit none
    integer, intent(in) :: i_stat
    if (i_stat /= 0) then
      write(10,'(A)') "Memory allocation failed!"
      stop 1
    end if
  end subroutine error_allocate

end module error_treat

