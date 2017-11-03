    module interf
    interface
    
    subroutine getname(argv,nlen)
 implicit none
 integer::narg
 integer,intent(out)::nlen
 character(*),intent(out)::argv
 integer::lnblnk,iargc
    end subroutine getname
    end interface
    
    end module interf