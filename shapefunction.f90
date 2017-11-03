      subroutine shapefunction(dershapefun,w1,w2,eta,xi)
    implicit none
    ! Variables
    integer,parameter::dd=selected_real_kind(15)
    real :: shapefun(4),dershapefun(4,2)
    integer :: w1, w2
    real :: eta,xi
    ! Body of shapefunction    
    
    shapefun(1)=1./4.*(1-xi)*(1-eta)
    shapefun(2)=1./4.*(1+xi)*(1-eta)
    shapefun(3)=1./4.*(1+xi)*(1+eta)
    shapefun(4)=1./4.*(1-xi)*(1+eta)
    dershapefun(1,1)=-1./4.*(1-eta)
    dershapefun(2,1)=1./4.*(1-eta)
    dershapefun(3,1)=1./4.*(1+eta)
    dershapefun(4,1)=-1./4.*(1+eta)
    dershapefun(1,2)=-1./4.*(1-xi)
    dershapefun(2,2)=-1./4.*(1+xi)
    dershapefun(3,2)=1./4.*(1+xi)
    dershapefun(4,2)=1./4.*(1-xi)
    
    end subroutine shapefunction