    module assemb
    contains
    
     subroutine assembler(k1,kg,iii,topology,ntype)
    implicit none
!    ! Variables
    integer :: i,j,r,s,iii,ntype
    integer,parameter::dd=selected_real_kind(15)
    real, allocatable :: coords(:,:),k1(:,:)
    real, allocatable :: kg(:,:)
    integer, allocatable :: topology(:,:),topx(:), topy(:)
!    ! Body of assembler
!    
    !kg=0.
    allocate (topx(ntype))
    allocate (topy(ntype))
    do i=1,ntype
        topx(i)=topology(iii,i)
    enddo
    topy=topx
    !topx=(/ 1, 3, 9, 8/)
    !topy=(/ 1, 3, 9, 8/)
!    
    do i=1,ntype
    
        do j=1,ntype
            do s=1,2
        do r=1,2
    kg(2*(topx(i)-1)+r,2*(topy(j)-1)+s)=kg(2*(topx(i)-1)+r,2*(topy(j)-1)+s)+k1(2*(i-1)+r,2*(j-1)+s)
    !kg(2*(topx(i)-1)+r,2*(topy(i)-1)+s)=kg(2*(topx(i)-1)+r,2*(topy(i)-1)+s)+k1(2*(i-1)+r,2*(i-1)+s)
        end do
    end do
    end do
        enddo
    end subroutine assembler
    end module assemb