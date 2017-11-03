    module jaco
    contains
        subroutine jacob(array,topology,coords,iii,ntype)
    implicit none
    ! Variables
    integer,parameter::dd=selected_real_kind(15)
    real, allocatable :: coords(:,:),array(:,:)
    integer, allocatable :: topology(:,:)
    integer :: iii,a,b,jjj,ntype
    
    ! Body of jacob
    !allocate (array(2,ntype))
    do a=1,2
        do b=1,ntype
        array(a,b)=coords(topology(iii,b),a)
    enddo
    enddo
    !array1 = reshape((/ 0, 0, 5, 0, 5, 5, 0, 5 /), shape(array))
    
    end subroutine jacob
    end module jaco