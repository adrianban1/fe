    subroutine dmat(dm,young,poisson)
    implicit none
    ! Variables
    integer,parameter::dd=selected_real_kind(15)
    real :: young, poisson,dm(3,3),const
    ! Body of dmat
    const=young/(1.-poisson*poisson)
    
   dm(1,1)=const
   dm(2,2)=const
   dm(1,2)=poisson*const
   dm(2,1)=poisson*const
   dm(3,3)=(const/2.)*(1.-poisson)
    
    end subroutine dmat