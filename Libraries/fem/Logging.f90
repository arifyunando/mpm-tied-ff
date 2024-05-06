module Logging
    INTEGER, PARAMETER, PRIVATE ::iwp=SELECTED_REAL_KIND(15)
    
contains
    
subroutine print_array_real(input_array)
    implicit none
    real(iwp) :: input_array(:, :)
    integer:: counting_index
    integer:: i, columns, page, lpage, starts, ends
        
    lpage = 40
    page = ubound(input_array, 2)/lpage + 1
    do i = 1, page
        starts = (i - 1)*lpage + 1
        ends = starts + lpage - 1
        do counting_index = 1, ubound(input_array, 1)
            if (ends < ubound(input_array, 2)) then
                print "(*(F8.2))", input_array(counting_index, starts:ends)
            else
                print "(*(F8.2))", input_array(counting_index, starts:)
            end if
        end do
        print*
    end do
end subroutine print_array_real
    

subroutine print_array_int(input_array)
    implicit none
    integer:: input_array(:, :)
    integer:: counting_index
    integer:: i, columns, page, lpage, starts, ends
        
    lpage = 40
    page = ubound(input_array, 2)/lpage + 1
    do i = 1, page
        starts = (i - 1)*lpage + 1
        ends = starts + lpage - 1
        do counting_index = 1, ubound(input_array, 1)
            if (ends < ubound(input_array, 2)) then
                print "(*(I8))", input_array(counting_index, starts:ends)
            else
                print "(*(I8))", input_array(counting_index, starts:)
            end if
        end do
        print*
    end do
end subroutine print_array_int
    

subroutine export_2D_array_real(input_array, port)
    implicit none
    real(iwp) :: input_array(:, :)
    integer :: port
    integer :: counting_index
        
    do counting_index = 1, ubound(input_array, 1)
        write (port, '(*(F9.2))') input_array(counting_index, :)
    end do
    write (port, *)     
end subroutine
    

subroutine export_2d_array_int(input_array, port)
    implicit none
    real(iwp) :: input_array(:, :)
    integer :: port
    integer :: counting_index
        
    do counting_index = 1, ubound(input_array, 1)
        write (port, '(*(I10))') input_array(counting_index, :)
    end do
    write (port, *)     
end subroutine


end module Logging
    
    