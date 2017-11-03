    program q7
    !
    !   Purpose:
    !       To calculate, using the finite element method, the displacements, strains and stresses in a plate
    !       with an elliptical hole. The program provides support for 3-node constant strin triangle elements,
    !       4-node quadratic plane stress elements,with either reduced integration (1 Gauss point) or full 
    !       integration (4 Gauss points). Dynamic dimensioning is used throughout, taking advantage of the 
    !       allocate construct of Fortran-90.The Gauss-Jordan algorithm is used to solve the global system.
    !
    !   Record of revisions:
    !       Date        Programmer          Description of change
    !       ====        ==========          =====================
    !       05.12.15    A. Baniak           Original release.
    !
    use assemb
    use jaco
    !use nrtype
    !use nrutil
    !use gauss
    use interf
    implicit none
    ! Variables
    integer,parameter::dd=selected_real_kind(15)
    real, allocatable :: coords(:,:),array(:,:),k1(:,:),eps1full(:),sig1full(:)
    real, allocatable :: kg(:,:),bk(:,:),del(:,:),bmat2tt(:,:,:),bmat2s(:,:,:),bmat2ss(:,:,:,:)
    character*10 start,endt,difft,time1,time2,date
    character(len=15)::argv
    real :: young, poisson,area,starttime,endtime
    real :: dershapefun(4,2),jay(2,2),jayf(2,2),detm,bmat1(4,2),bmat2(3,8),bmat1t(3,6),bmat2t(3,6)
    real :: x,y,dm(3,3),t,w1,w2,ly(200),eps(3,1),sig(3,1),ket(6,6),lx(200),eps1(3,1),sig1(3,1)
    real, parameter :: augment=1e30
    integer nnodes,nelem,k,l,m,ndof,i,j,iii,iiii,stat,nl(200),ibc,nbc(200),lbc(200,2),s,jj,cc,intpoint,nlen,ntype,intp,jk,ijk
    integer, allocatable :: topology(:,:),ncoord(:),nnelem(:)
    character*80 :: str1,str2,str3,str4,str5,str6,str7,str0,strn1,etype
    logical redint
    
    ! body of q7
    call cpu_time(starttime)
    call date_and_time(time=start)
     call getname(argv,nlen)
    open(1,file=argv(1:nlen)//'.inp') 
    open(2,file=argv(1:nlen)//'.out')
    !Read input file & allocate memory for arrays
    read(1,'(a)')strn1
    read(1,'(a)')etype
    read(1,'(a)')str0
    read(1,*)young,poisson
    read(1,'(a)')str1
    read(1,*)nnodes
    read(1,'(a)')str2
    read(1,*)nelem
    !echo input
    write(2,'(a)')"AGB-Quest, A 2-Dimensional Finite Element Program"
    write(2,'(a)')"Output File"
    write(2,*)
    call DATE_AND_TIME(date)
    write(2,'(a)')"Date:"
    write(2,'(a)')date
    write(2,*)
    write(2,'(a)')"Input File Chosen for Run:"
    write(2,'(a)')argv
    write(2,*)
    write(2,'(a)')"Chosen Element Type:"
    write(2,'(a)')etype
    write(2,*)
    write(2,'(a)')"Youngs Modulus and Poissons Ratio:"
    write(2,*)young,poisson
    write(2,*)
    write(2,'(a)')"Number of Nodes in Model:"
    write(2,*)nnodes
    write(2,*)
    write(2,'(a)')"Number of Elements in Model:"
    write(2,*)nelem
    write(2,*)
    
    read(1,'(a)')str3
    allocate (coords(nnodes,2))
    allocate (ncoord(nnodes))
    coords = 0.
    do i=1,nnodes
        read(1,*)ncoord(i),(coords(i,j),j=1,2)
    enddo
    
    write(2,'(a)')"Nodal Coordinate Data:"
        do i=1,nnodes
        write(2,*)ncoord(i),(coords(i,j),j=1,2)
    enddo
    
    read(1,'(a)')str4
    if(etype=="CST")then
        ntype=3
    else
        ntype=4
    endif
        allocate (k1(ntype*2,ntype*2))
        allocate (topology(nelem,ntype))
        allocate (nnelem(nelem))
        allocate (del(2*ntype,1))
        allocate (bmat2tt(nelem,3,6))
        allocate (bmat2s(nelem,3,8))
        allocate (bmat2ss(nelem,4,3,8))
        allocate (eps1full(nnodes*3))
        allocate (sig1full(nnodes*3))
        eps1full=0.
        sig1full=0.
    do i=1,nelem
    read(1,*)nnelem(i),(topology(i,j),j=1,ntype)
    enddo
    write(2,*)
    write(2,'(a)')"Element Topology Data:"
        do i=1,nelem
        write(2,*)nnelem(i),(topology(i,j),j=1,ntype)
    enddo
    
    read(1,'(a)')str5
    iiii=0
    do
        iiii=iiii+1
   read(1,*, iostat=stat) nl(iiii),lx(iiii),ly(iiii)
   if (stat /= 0) exit
    end do
    backspace(1)
    read(1,'(a)')str6
    ibc=0
    do
        ibc=ibc+1
        read(1,*, iostat=stat) nbc(ibc),(lbc(ibc,j),j=1,2)
           if (stat /= 0) exit
    end do
    
    t=1.
    ndof=nnodes*2
    allocate (kg(ndof,ndof))
    kg=0.
    allocate (array(2,ntype))
    !begin element by element loop. calculate the jacobian, determine derivatives of shape functions
    ! and calculate B matrix.
    do iii=1,nelem
        k1=0.
    call jacob(array,topology,coords,iii,ntype)
    if(etype=="CST")then
        area=(0.5)*(array(1,2)-array(1,1))*(array(2,3)-array(2,1))
        bmat1t(1,1)=array(2,2)-array(2,3)
        bmat1t(1,3)=array(2,3)-array(2,1)
        bmat1t(1,5)=array(2,1)-array(2,2)
        bmat1t(2,2)=array(1,3)-array(1,2)
        bmat1t(2,4)=array(1,1)-array(1,3)
        bmat1t(2,6)=array(1,2)-array(1,1)
        bmat1t(3,1)=bmat1t(2,2)
        bmat1t(3,2)=bmat1t(1,1)
        bmat1t(3,3)=bmat1t(2,4)
        bmat1t(3,4)=bmat1t(1,3)
        bmat1t(3,5)=bmat1t(2,6)
        bmat1t(3,6)=bmat1t(1,5)
        
        bmat2t=((1.)/((2.)*(area)))*bmat1t
        
        bmat2tt(iii,:,:)=bmat2t(:,:)
        
        call dmat(dm,young,poisson)
        
        k1=t*area*matmul(matmul(transpose(bmat2t),dm),bmat2t)
        i=1
        
    ! if reduced integration:
    elseif(etype=="Q4RI")then
    call shapefunction(dershapefun,2,2,0.,0.)
    
    jay=matmul(array,dershapefun)
    
    detm=(jay(1,1)*jay(2,2))-(jay(1,2)*jay(2,1))
    write(3,*)detm
    jayf(1,1)=jay(2,2)
    jayf(1,2)=-jay(1,2)
    jayf(2,1)=-jay(2,1)
    jayf(2,2)=jay(1,1)
    jayf=jayf/detm
    bmat1=matmul(dershapefun,jayf)
    
    
       do m=1,4
     k=2*m
     l=k-1
     x=bmat1(m,1)
     y=bmat1(m,2)
     bmat2(1,l)=x
     bmat2(3,k)=x
     bmat2(2,k)=y
     bmat2(3,l)=y
       end do
       bmat2s(iii,:,:)=bmat2(:,:)
       !calculate D matrix
       call dmat(dm,young,poisson)
       
       !get element stiffness matrix
       w1=2.
       w2=2.
    k1=matmul(matmul(transpose(bmat2),dm),bmat2)*detm*w1*w2*t
    !write(4,*)k1
    i=1
    ! if full integration requested, assemble local stiffness matrix first.
    elseif(etype=="Q4FI")then
        do intpoint=1,4
        if(intpoint==1)call shapefunction(dershapefun,1,1,-1./sqrt(3.),-1./sqrt(3.))
        if(intpoint==2)call shapefunction(dershapefun,1,1,1./sqrt(3.),-1./sqrt(3.))
        if(intpoint==3)call shapefunction(dershapefun,1,1,1./sqrt(3.),1./sqrt(3.))
        if(intpoint==4)call shapefunction(dershapefun,1,1,-1./sqrt(3.),1./sqrt(3.))
        
            jay=matmul(array,dershapefun)
    
    detm=(jay(1,1)*jay(2,2))-(jay(1,2)*jay(2,1))
    write(3,*)detm
    jayf(1,1)=jay(2,2)
    jayf(1,2)=-jay(1,2)
    jayf(2,1)=-jay(2,1)
    jayf(2,2)=jay(1,1)
    jayf=jayf/detm
    bmat1=matmul(dershapefun,jayf)
    
       do m=1,4
     k=2*m
     l=k-1
     x=bmat1(m,1)
     y=bmat1(m,2)
     bmat2(1,l)=x
     bmat2(3,k)=x
     bmat2(2,k)=y
     bmat2(3,l)=y
       end do
       bmat2ss(iii,intpoint,:,:)=bmat2(:,:)
       !calculate D matrix
       call dmat(dm,young,poisson)
       
       !get element stiffness matrix
       w1=1.
       w2=1.
        k1=k1+matmul(matmul(transpose(bmat2),dm),bmat2)*detm*w1*w2*t
        
        enddo
        endif
        
    !Assemble global stiffness matrix
    call assembler(k1,kg,iii,topology,ntype)
    !write(15,*)kg
    enddo

    !Add boundary conditions by augmenting diagonals on global stiffness matrix
    do i=1,ibc-1
        if(lbc(i,1)==1) then
        kg(2*nbc(i)-1,2*nbc(i)-1)=kg(2*nbc(i)-1,2*nbc(i)-1)+augment
        endif
    if(lbc(i,2)==1) then
        kg(2*nbc(i),2*nbc(i))=kg(2*nbc(i),2*nbc(i))+augment
    endif
    enddo
    do i=1,nnodes*2
        write(15,*)kg(i,i)
    enddo
    
    !figure out solving Ax=b by gauss-jordan (NR)
    allocate (bk(nnodes*2,1))
    bk=0.
    do i=1,iiii-1
        bk(2*nl(i),1)=ly(i)
        bk(2*nl(i)-1,1)=lx(i)
        enddo
    !call gaussj(kg,bk)
    call gaussj2(kg,ndof,bk)
    
    write(2,*)
    write(2,'(a)')"Displacements in nodal order "
    write(2,'(a)')"(i.e.: node 1 x, node 1 y, node 2 x, node 2 y etc.):"
    write(2,*)
    do i=1,nnodes*2
    write(2,*)bk(i,1)
    enddo
    write(2,*)
    write(2,'(a)')"Maximum value of displacement:"
    write(2,*)maxval(bk)
    write(2,*)
    write(2,'(a)')"Minimum value of displacement:"
    write(2,*)minval(bk)
    s=0
    !extrapolate the strains and stresses
    write(2,*)
    write(2,*)
    write(2,'(a)')"Strains (left column) and  "
    write(2,'(a)')"Stresses (right column):"
    write(2,*)
    do i=1,nelem
            do jj=1,ntype
            do s=1,0,-1
            del(2*jj-s,1)=bk(2*topology(i,jj)-s,1)
            enddo
            enddo
            
            if(ntype==3)then
            eps=matmul(bmat2tt(i,:,:),del)
            sig=matmul(dm,eps)
                    jk=1
                    do ijk=i*3-2,i*3
                    eps1full(ijk)=eps(jk,1)
                    sig1full(ijk)=sig(jk,1)
                    jk=jk+1
                    enddo
                elseif(etype=="Q4RI")then
                eps=matmul(bmat2s(i,:,:),del)
            sig=matmul(dm,eps)
                    jk=1
                    do ijk=i*3-2,i*3
                    eps1full(ijk)=eps(jk,1)
                    sig1full(ijk)=sig(jk,1)
                    jk=jk+1
                    enddo
                elseif(etype=="Q4FI")then
                    do intpoint=1,4
            eps=matmul(bmat2ss(i,intpoint,:,:),del)
            sig=matmul(dm,eps)
            !jk=1
            !        do ijk=i*3-2,i*3
            !        eps1fullintp(ijk)=eps(jk,1)
            !        sig1fullintp(ijk)=sig(jk,1)
            !        jk=jk+1
            !        enddo
            write(2,*)i
            do cc=1,3
            write(2,*)eps(cc,1), sig(cc,1)
            enddo
            write(2,*)
            eps1=eps1+eps
            sig1=sig1+sig
                    enddo
                    eps1=eps1/4.
                    sig1=sig1/4.
                    jk=1
                    do ijk=i*3-2,i*3
                    eps1full(ijk)=eps1(jk,1)
                    sig1full(ijk)=sig1(jk,1)
                    jk=jk+1
                    enddo
                    
            endif
            if(etype=="Q4FI")then
            write(2,*)i
            do cc=1,3
            write(2,*)eps1(cc,1), sig1(cc,1)
            enddo
            write(2,*)
            
            else
            write(2,*)i
            do cc=1,3
            write(2,*)eps(cc,1), sig(cc,1)
            enddo
            write(2,*)
            endif
            eps1=0.
            sig1=0.
    enddo
    
    write(2,*)
    write(2,'(a)')"Maximum value of stress:"
    write(2,*)maxval(sig1full)
    write(2,*)
    write(2,'(a)')"Minimum value of stress:"
    write(2,*)minval(sig1full)
    
    call date_and_time(TIME=endt)
    ! difft=start-endt
    write(2,*)
    write(2,'(a)')"Run time (seconds):"
    call cpu_time(endtime)
    write(2,*)endtime-starttime
    end program q7