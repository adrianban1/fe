subroutine getname(argv,nlen)
!
! this subroutine reads the base name of data file.
!
 implicit none
 integer::narg
 integer,intent(out)::nlen
 integer::lnblnk,iargc
 character(*),intent(out)::argv
 logical found
 narg=iargc()
 if(narg<1)then
   write(*,*)'Please enter the base name of the input file as per the readme: '
   read(*,*)argv
  else
   call getarg(1,argv)
 endif
 nlen=lnblnk(argv)
 inquire(file=argv(1:nlen)//'.inp',exist=found)
 if(.not.found)then
  write(*,*)'Input file not found: ',argv(1:nlen)//'.inp'
  write(*,*)'Please check spelling.'
  stop
 endif
return
end subroutine getname
