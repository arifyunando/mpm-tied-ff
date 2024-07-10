module sparse_lib
contains
!--------------------------------BASIC INFORMATION------------------------------
!-------------------------------------------------------------------------------
subroutine snfnpn(nf,nn,nodof,ndim,nels,ntot,sn,fn,pn,ebenz)
    ! This subroutine computes the displacement DOFs (sn)
    !     and pore pressure DOFs (fn).
    !     Then it gives an estimation (ebenz).
    implicit none
        integer:: i,nf(:,:),nn,nodof,ntot,ndim,nels,sn,fn,pn,ebenz
        sn=0 ; fn = 0; pn = 0
        fn = sum(nf(3,:)) ;
        do i=1,ndim ; pn = pn + sum(nf(i,:)) ; end do
        do i=4,5; fn = fn + sum(nf(i,:)) ;    end do   
        ebenz = ntot*int(ntot/2)*nels ; ! Estimated number
        return
end subroutine snfnpn
!-------------------------------------------------------------------------------
subroutine form_id(nf,nn,nodof,id)
    ! This subroutine for the identifier array, "id".
    !     nf(:,:) is the original input node freedom array.  
    !     nn - total node number.
    !     nodof - number of freedoms per node.
    implicit none
    integer:: i,nf(:,:),nn,nodof,id(:)
    id(:) = 1; 
    do i = 1,nn ; 
        if(nf(3,i)/=0) id(nf(3,i))=0; 
        if(nf(4,i)/=0) id(nf(4,i))=-1;
        if(nf(5,i)/=0) id(nf(5,i))=-1;
    end do;
    !
    return
end subroutine form_id
!-------------------------------SORTING SUBROUTINES-----------------------------
!-------------------------------------------------------------------------------
subroutine quicksort(uanz,arr,brr,crr)
    ! This subroutine
    ! Quicksort - sorts arr into ascending order, brr and crr change
    !             correspondingly.
    ! quicksort chooses a "pivot" in the set, and explores the array
    ! from both ends, looking for a value > pivot with the increasing
    ! index, for a value <= pivot with the decreasing index, and
    ! swapping them when it has found one of each. !  The array is then
    ! subdivided in 2 ([3]) subsets: { values <= pivot} {pivot}
    ! {values > pivot}. One then call recursively the program to sort
    ! each subset.  When the size of the subarray is small enough, one
    ! uses an insertion sort that is faster for very small sets.
    ! Sorting an array arr(1:n) into ascending order with quicksort,
    ! while making the corresponding rarrangements of arrays brr(1:n)
    ! and crr(1:n).
    ! (Revised from ORDERPACK codes)
    !     uanz: the nonzero number of arr (or brr, crr).
    !--------------------------------
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    !real(8):: crr(:)
    REAL(iwp),INTENT(IN OUT)::crr(:)
    INTEGER,INTENT(IN OUT)::uanz,arr(:),brr(:)
    !integer::uanz,arr(:),brr(:)
    !
    call subsort(arr,brr,crr,1, uanz) ;
    call inssor(arr,brr,crr,uanz) ;
    !
    return
end subroutine quicksort
!-------------------------------------------------------------------------------
Recursive subroutine subsort(arr,brr,crr, ideb1, ifin1)
    !  This subroutine sorts arr from ideb1 to ifin1
        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        REAL(iwp),INTENT(IN OUT)::crr(:)
        REAL(iwp)::zwrk
        !real(8)::  crr(:), zwrk
        integer, intent (in) :: ideb1, ifin1
        INTEGER,INTENT(IN OUT)::arr(:), brr(:)
        !integer :: arr(:), brr(:), icrs, ideb, idcr, ifin, imil, xpiv, &
        !           xwrk, ywrk, nins = 16  ! Max for insertion sort
        INTEGER::icrs, ideb, idcr, ifin, imil,xpiv,xwrk, ywrk, nins = 16  ! Max for insertion sort
        ideb = ideb1
        ifin = ifin1
    !  if we don't have enough values to make it worth while, we leave
    !  them unsorted, and the final insertion sort will take care of them.
        if ((ifin - ideb) > nins) Then
        imil = (ideb+ifin) / 2
    !  One chooses a pivot, median of 1st, last, and middle values
    !
        if (arr(imil) < arr(ideb)) Then
    xwrk = arr(ideb) ;      ywrk=brr(ideb);         zwrk=crr(ideb)
    arr(ideb) = arr(imil); brr(ideb) = brr(imil); crr(ideb) = crr(imil)
    arr(imil) = xwrk;      brr(imil) = ywrk;      crr(imil) = zwrk
        end if
        !
        if (arr(imil) > arr(ifin)) Then
    xwrk = arr(ifin);       ywrk = brr(ifin);       zwrk = crr(ifin)
    arr(ifin) = arr(imil); brr(ifin) = brr(imil); crr(ifin) = crr(imil)
    arr(imil) = xwrk;      brr(imil) = ywrk;      crr(imil) = zwrk
            if (arr(imil) < arr(ideb)) Then
    xwrk = arr(ideb);       ywrk = brr(ideb);       zwrk = crr(ideb)
    arr(ideb) = arr(imil); brr(ideb) = brr(imil); crr(ideb) = crr(imil)
    arr(imil) = xwrk;      brr(imil) = ywrk;      crr(imil) = zwrk
            end if
        end if
        xpiv = arr(imil)
    !
    !  One exchanges values to put those > pivot in the end and
    !  those <= pivot at the beginning
    !
        icrs = ideb
        idcr = ifin
            ech2: do         !---------------------------
            do
                icrs = icrs + 1
                if (icrs >= idcr) Then
    !
    !  the first > pivot is idcr
    !  the last <= pivot is icrs-1
    !  Note: if one arrives here on the first iteration, then
    !        the pivot is the maximum of the set, the last value is equal
    !        to it, and one can reduce by one the size of the set to
    !        process, as if arr (ifin) > xpiv
    !
                    exit ech2
    !
                end if
                if (arr(icrs) > xpiv) exit
            end do
            !
            do
                if (arr(idcr) <= xpiv) exit
                idcr = idcr - 1
                if (icrs >= idcr) Then
    !
    !  The last value < pivot is always icrs-1
    !
                    exit ech2
                end if
            end do
    !
            xwrk = arr(idcr);      ywrk = brr(idcr);  zwrk = crr(idcr)
            arr(idcr)=arr(icrs);brr(idcr)=brr(icrs);crr(idcr)=crr(icrs)
            arr(icrs)=xwrk;      brr(icrs)=ywrk;    crr(icrs)=zwrk
        end do ech2         !---------------------------
    !
    !  One now sorts each of the two sub-intervals
    !
        call subsort (arr,brr,crr, ideb1, icrs-1)
        call subsort (arr,brr,crr, idcr, ifin1)
        end if
        !
        return
end subroutine subsort
!-------------------------------------------------------------------------------
subroutine inssor(arr,brr,crr,uanz)
    ! This subroutine sorts arr into increasing order (Insertion sort)
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    INTEGER,INTENT(IN OUT)::arr(:),brr(:),uanz
    INTEGER:: icrs, idcr, xwrk,ywrk
    REAL(iwp),INTENT(IN OUT):: crr(:)
    REAL(iwp):: zwrk
    !integer :: arr(:),brr(:),uanz, icrs, idcr, xwrk,ywrk
    !real(8) :: crr(:),zwrk
    !
    do icrs = 2, uanz
        xwrk = arr(icrs);     ywrk = brr(icrs);   zwrk = crr(icrs);
        if (xwrk >= arr(icrs-1)) cycle
        arr(icrs)=arr(icrs-1);brr(icrs)=brr(icrs-1);crr(icrs)=crr(icrs-1)
        do idcr = icrs - 2, 1, - 1
        if (xwrk >= arr(idcr)) exit
        arr (idcr+1) = arr (idcr)
        brr (idcr+1) = brr (idcr)
        crr (idcr+1) = crr (idcr)
        end do
        arr(idcr+1) = xwrk;   brr(idcr+1) = ywrk;    crr(idcr+1) = zwrk
        end do
        !
        return
end subroutine inssor
!-------------------------------------------------------------------------------
subroutine sortadd(uanz,arr,brr,crr,ni,nnz,penpos)
    ! For the same arr index, subsort brr, and at the same time, crr
    ! changes correspondingly with brr. After this work, adding up all crr
    ! components with the same (arr, brr) or (brr, arr) index, and the
    ! zero-value crr entry will be removed. Finally forming the Compressed
    ! Sparse Row (CSR) format or Compressed Sparse Column (CSC) format
    ! to overwrite arr,brr,crr.
    !         uanz: the nonzero number of arr (or brr, crr).
    !  arr,brr,crr: three vectors required to be sorted.
    !           ni: = n + 1 (n is dimension of A)
    !          nnz: the nonzero number of crr.
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    !integer:: i,j,k,k1,k2,m,arr(:),brr(:),uanz,nnz,ni 
    integer, allocatable:: itep(:)
    !real(8):: crr(:), aa 
    INTEGER::i,j,k,k1,k2,m,acum
    INTEGER,INTENT(IN OUT)::arr(:),brr(:),uanz,nnz,penpos(:)!,ni
    INTEGER,INTENT(IN)::ni
    REAL(iwp),INTENT(IN OUT)::crr(:)
    REAL(iwp)::aa
    allocate (itep(ni))
    call quicksort(uanz,arr,brr,crr) ; ! sorting three vectors
    k=1;  itep(1)=1
    do i=2, uanz
        if(arr(i)/=arr(i-1)) then
        k=k+1 ;  itep(k)=i
        end if
    end do
    itep(k+1)=uanz+1
    !----------------------------
    do i=1, k
        k1=itep(i);  k2=itep(i+1)-1
        j=k2-k1+1
        if(j<=16) then    ! sub-brr sorting by Insertion sort if j <= 16.
        call subbrr2(brr(k1:k2),crr(k1:k2),j)
        else               ! quick sorting when j is larger (>16).
            call quicksort2(j,brr(k1:k2),crr(k1:k2))
        end if
    end do
    !----------------------------
    m = 0 ;  
    acum=1
    do i=1, k
    k1=itep(i);  k2=itep(i+1)-1  ;  m=m+1;   
    arr(i) = m ;  brr(m) = brr(k1) ; aa = .0
    do j=k1, k2-1
        aa = aa + crr(j) ; 
        if(brr(j+1)/=brr(j) ) then
            if(aa /=.0) then
            crr(m) = aa
            IF(crr(m)>10.0E10)THEN
                penpos(acum)=m
                acum=acum+1
            END IF
            m=m+1 ; 
            brr(m)= brr(j+1)  
            aa = .0
        else            !  aa is removed when it is zero.
            brr(m)= brr(j+1)
        end if
        end if
    end do
        crr(m) = aa + crr(k2)
        if(crr(m)==.0) m=m-1
        
    end do
    arr(k+1)=m+1;  nnz=m
    !
    return
end subroutine sortadd
!-------------------------------------------------------------------------------
subroutine quicksort2(uanz,arr,crr)
    ! This subroutine Quicksort2
    !       - sorts arr into ascending order, crr changes correspondingly.
    !  Sorts arr into ascending order - Quicksort
    !  Quicksort chooses a "pivot" in the set, and explores the
    !  array from both ends, looking for a value > pivot with the
    !  increasing index, for a value <= pivot with the decreasing
    !  index, and swapping them when it has found one of each.
    !  The array is then subdivided in 2 ([3]) subsets:
    !  { values <= pivot} {pivot} {values > pivot}
    !  One then call recursively the program to sort each subset.
    !  When the size of the subarray is small enough, one uses an
    !  insertion sort that is faster for very small sets.
    !  Sorting an array arr(1:n) into ascending order with quicksort,
    !  while making the corresponding rarrangements of arrays crr(1:n).
    ! (Revised from ORDERPACK codes)
    !-------------------------------
    real(8):: crr(:)
    integer::uanz,arr(:)
    !
    call subsort2(arr,crr,1, uanz)
    call inssor2(arr,crr,uanz)
    !
    return
end subroutine quicksort2
!-------------------------------------------------------------------------------
Recursive subroutine subsort2(arr,crr, ideb1, ifin1)
    ! This subroutine sorts arr from ideb1 to ifin1
    real(8)::  crr(:), zwrk
    integer, intent (in) :: ideb1, ifin1
    integer :: arr(:),  icrs, ideb, idcr, ifin, imil, xpiv, &
                xwrk,  nins = 16  ! Max for insertion sort
        ideb = ideb1
        ifin = ifin1
    ! if we don't have enough values to make it worth while, we leave
    ! them unsorted, and the final insertion sort will take care of them
    if ((ifin - ideb) > nins) Then
        imil = (ideb+ifin) / 2
    ! One chooses a pivot, median of 1st, last, and middle values
    !
        if (arr(imil) < arr(ideb)) Then
    xwrk = arr(ideb) ;              zwrk=crr(ideb)
    arr(ideb) = arr(imil);  crr(ideb) = crr(imil)
    arr(imil) = xwrk;        crr(imil) = zwrk
        end if
        if (arr(imil) > arr(ifin)) Then
    xwrk = arr(ifin);              zwrk = crr(ifin)
    arr(ifin) = arr(imil);  crr(ifin) = crr(imil)
    arr(imil) = xwrk;       crr(imil) = zwrk
            if (arr(imil) < arr(ideb)) Then
    xwrk = arr(ideb);              zwrk = crr(ideb)
    arr(ideb) = arr(imil);  crr(ideb) = crr(imil)
    arr(imil) = xwrk;       crr(imil) = zwrk
            end if
        end if
        xpiv = arr(imil)
    !
    !  One exchanges values to put those > pivot in the end and
    !  those <= pivot at the beginning
    !
        icrs = ideb
        idcr = ifin
        ech2: do        !-------------------------
            do
                icrs = icrs + 1
                if (icrs >= idcr) Then
    !
    !  the first  >  pivot is idcr
    !  the last   <= pivot is icrs-1
    !  Note: if one arrives here on the first iteration, then
    !        the pivot is the maximum of the set, the last value is equal
    !        to it, and one can reduce by one the size of the set to
    !        process, as if arr(ifin) > xpiv
    !
                    exit ech2
    !
                end if
                if (arr(icrs) > xpiv) exit
            end do
            !
            do
                if (arr(idcr) <= xpiv) exit
                idcr = idcr - 1
                if (icrs >= idcr) Then
    !
    !  The last value < pivot is always icrs-1
    !
                    exit ech2
                end if
            end do
            !
            xwrk = arr(idcr);               zwrk = crr(idcr)
            arr(idcr) = arr(icrs);   crr(idcr) = crr(icrs)
            arr(icrs) = xwrk;        crr(icrs) = zwrk
        end do ech2       !---------------------------
    !
    !  One now sorts each of the two sub-intervals
    !
        call subsort2 (arr,crr, ideb1, icrs-1)
        call subsort2 (arr,crr, idcr, ifin1)
        end if
        !
        return
end subroutine subsort2
!-------------------------------------------------------------------------------
subroutine inssor2(arr,crr,uanz)
    ! This subroutine sorts arr into increasing order (Insertion sort)
    integer :: arr(:),uanz, icrs, idcr, xwrk
    real(8) :: crr(:),zwrk
    !
    do icrs = 2, uanz
        xwrk = arr (icrs);         zwrk = crr (icrs);
        if (xwrk >= arr(icrs-1)) cycle
            arr (icrs) = arr (icrs-1);  crr (icrs) = crr (icrs-1)
            do idcr = icrs - 2, 1, - 1
            if (xwrk >= arr(idcr)) exit
            arr (idcr+1) = arr (idcr)
            crr (idcr+1) = crr (idcr)
            end do
            arr (idcr+1) = xwrk;       crr (idcr+1) = zwrk
    end do
    !
    return
end subroutine inssor2
!-------------------------------------------------------------------------------
subroutine subbrr2(br,cr,n)
    ! For the same arr index, subsort brr, and at the same time, crr
    ! changes correspondingly with brr.Because of small number of sub-brr,
    ! Insertion sort should be faster.
    integer:: br(n), icrs,idcr,n,ywrk
    real(8):: cr(n), zwrk
    !
    do icrs = 2, n
        ywrk = br (icrs);       zwrk = cr (icrs)
        if (ywrk >= br(icrs-1)) Cycle
            br (icrs) = br (icrs-1); cr (icrs) = cr (icrs-1)
            do idcr = icrs - 2, 1, - 1
            if (ywrk >= br(idcr)) exit
            br (idcr+1) = br (idcr)
            cr (idcr+1) = cr (idcr)
            end do
            br (idcr+1) = ywrk;    cr (idcr+1) = zwrk
    end do
end subroutine subbrr2
!-------------------------------------------------------------------------------
subroutine formspars(ntot,g,ke,iebea,jebea,ebea,ebeanz)
    ! This subroutine collect non-zero entries for each new generated    &
    !                 element stiffness matrix, forming the element-level&
    !                 three vectors which store nonzero entries of upper &
    !                 triangular part of A.
    !           ntot: total freedoms per element;
    !              g: element steering vector;
    !             ke: element "stiffness" matrix;
    !          iebea: global row index;
    !          jebea: global column index;
    !           ebea: correspondent value of the nonzero element stiffness
    !                 entry;
    !         ebeanz: true total number of nonzero element-level entries,
    !                 not estimated number any more when returned.
    implicit none
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    !real(8):: ke(:,:),ebea(:)
    REAL(iwp),INTENT(IN) ::ke(:,:)
    REAL(iwp),INTENT(IN OUT) ::ebea(:)
    INTEGER,INTENT(IN OUT)  ::ntot,ebeanz,g(:),iebea(:),jebea(:)
    INTEGER ::i,j
    REAL(iwp) ::zero=0.0_iwp
    !integer::i,j,ntot,ebeanz,g(:),iebea(:),jebea(:)
    !--- Storing upper triangle of element stiffness column by column ---
        do j=1, ntot
            CONTINUE
            do i=1, j
            if(g(i)/=0.and.g(j)/=0) then
                if(ke(i,j)/=zero)then
                if(g(i)<=g(j))then
                    ebeanz=ebeanz+1 ;  iebea(ebeanz)=g(i)
                    jebea(ebeanz)=g(j) ; ebea(ebeanz)=ke(i,j)
                else
                    ebeanz=ebeanz+1 ;  iebea(ebeanz)=g(j)
                    jebea(ebeanz)=g(i) ; ebea(ebeanz)=ke(i,j)
                end if
                end if
            end if
            end do
        end do
        !
        return
end subroutine formspars
!-------------------------------------------------------------------------------
subroutine formspars_unsym(ntot,g,ke,iebea,jebea,ebea,ebeanz)
    ! This subroutine collect non-zero entries for each new generated    &
    !                 element stiffness matrix, forming the element-level&
    !                 three vectors which store nonzero entries of upper &
    !                 triangular part of A.
    !           ntot: total freedoms per element;
    !              g: element steering vector;
    !             ke: element "stiffness" matrix;
    !          iebea: global row index;
    !          jebea: global column index;
    !           ebea: correspondent value of the nonzero element stiffness
    !                 entry;
    !         ebeanz: true total number of nonzero element-level entries,
    !                 not estimated number any more when returned.
    implicit none
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp)::ke(:,:),ebea(:),zero=0.0_iwp
    !real(8):: ke(:,:),ebea(:)
    integer::i,j,ntot,ebeanz,g(:),iebea(:),jebea(:)
    !--- Storing upper triangle of element stiffness column by column ---
        do j=1, ntot
            do i=1, ntot
            if(g(i)/=0.and.g(j)/=0) then
                !if(ke(i,j)/=.0)then !Old rule
                if(ABS(ke(i,j))/=zero)then
                !if(g(i)<=g(j) )then
                    ebeanz=ebeanz+1 ;  iebea(ebeanz)=g(i)
                    jebea(ebeanz)=g(j) ; ebea(ebeanz)=ke(i,j)
                !else
                !  ebeanz=ebeanz+1 ;  iebea(ebeanz)=g(j)
                !  jebea(ebeanz)=g(i) ; ebea(ebeanz)=ke(i,j)
                !end if
                end if
            end if
            end do
        end do
        !
        return
end subroutine formspars_unsym
!-------------------------------------------------------------------------------
!-------- SUBROUTINES FOR MATRIX-VECTOR PRODUCTS AND TRIANGULAR SOLVERS --------
!-------------------------------------------------------------------------------
subroutine cscbx(icscb,jcscb,cscb,x,y)
    ! Compute y=B*x, and B is stored in CSC (icscb,jcscb,cscb) format.
    ! x: input vector
    ! y: output vector
    implicit none
    real(8):: cscb(:),x(:),y(:)
    integer::i,j,n,k1,k2,icscb(:),jcscb(:)
        n=ubound(jcscb,1)-1 ; y=.0 ;
        do i=1, n
        if(x(i)/=.0) then
            k1=jcscb(i); k2=jcscb(i+1)-1
            do j=k1, k2
            y(icscb(j))=y(icscb(j))+cscb(j)*x(i)
            end do
        end if
        end do
end subroutine cscbx
!-------------------------------------------------------------------------------
subroutine cscbtx(icscb,jcscb,cscb,x,y)
    ! Compute y=B'*x, and B is stored in CSC (icscb,jcscb,cscb) format.
    ! x: input vector
    ! y: output vector
    implicit none
    real(8):: cscb(:),x(:),y(:)
    integer::i,j,n,k1,k2,icscb(:),jcscb(:)
        n=ubound(jcscb,1)-1 ; y=.0 ;
        do j=1, n
        k1=jcscb(j); k2=jcscb(j+1)-1
        do i=k1, k2
            y(j)=y(j)+cscb(i)*x(icscb(i))
        end do
        end do
end subroutine cscbtx
!-------------------------------------------------------------------------------
subroutine csrbx(icsrb,jcsrb,csrb,x,y)
    ! Compute y=B*x, and B is stored in CSR (icsrb,jcsrb,csrb) format.
    ! x: input vector
    ! y: output vector
    implicit none
    real(8):: csrb(:),x(:),y(:)
    integer::i,j,k1,k2,n,icsrb(:),jcsrb(:)
        n=ubound(icsrb,1)-1 ; y=.0 ;
        do i=1, n
        k1=icsrb(i); k2=icsrb(i+1)-1
        do j=k1, k2
            y(i)=y(i)+csrb(j)*x(jcsrb(j))
        end do
        end do
end subroutine csrbx
!-------------------------------------------------------------------------------
subroutine csrbtx(icsrb,jcsrb,csrb,x,y)
    ! Compute y=B'*x, and B is stored in CSR (icsrb,jcsrb,csrb) format.
    ! x: input vector
    ! y: output vector
    implicit none
    real(8):: csrb(:),x(:),y(:)
    integer::i,j,n,k1,k2,icsrb(:),jcsrb(:)
        n=ubound(icsrb,1)-1 ;  y=.0 ;
        do i=1, n
        if(x(i)/=.0) then
            k1=icsrb(i); k2=icsrb(i+1)-1
            do j=k1, k2
            y(jcsrb(j))=y(jcsrb(j))+csrb(j)*x(i)
            end do
        end if
        end do
end subroutine csrbtx
!-------------------------------------------------------------------------------
subroutine csrax(icsr,jcsr,csra,x,y)
    ! Compute y = Ax with A is symmetric and square, only upper triangular
    !                 part is stored.
    !              n: dimension of coefficient matrix A;
    ! icsr,jcsr,csra: CSR storage of upper triangular part of matrix A;
    !              x: is input vector;
    !              y: is output vector.
    implicit none
    real(8):: csra(:),x(:),y(:)
    integer::i,j,k1,k2,n,icsr(:),jcsr(:)
    n = ubound(icsr,1)-1 ; y=.0 ;
    do i=1, n
        k1=icsr(i); k2=icsr(i+1)-1
        do j=k1, k2
        y(i)=y(i)+csra(j)*x(jcsr(j))
        end do
        !
        if(x(i)/=.0) then
        do j=k1+1, k2
        y(jcsr(j))=y(jcsr(j))+csra(j)*x(i)
        end do
        end if
    end do
end subroutine csrax
!-------------------------------------------------------------------------------
subroutine cscax(icsc,jcsc,csca,x,y)
    ! Compute y = Ax with A is symmetric and square, only upper triangular
    !                 part is stored.
    !              n: dimension of coefficient matrix A;
    ! icsc,jcsc,csca: CSC storage of upper triangular part of matrix A;
    !              x: is input vector;
    !              y: is output vector.
    implicit none
    real(8):: csca(:),x(:),y(:),tmp
    integer::j,k,r,k1,k2,n,jcsc(:),icsc(:)
    n = ubound(jcsc,1)-1 ; y=.0
    do j=1, n
        if(x(j)/=.0)then
        k1=jcsc(j); k2=jcsc(j+1)-1 ;
        do k=k1, k2
            r=icsc(k) ;
            y(r)=y(r)+csca(k)*x(j) ;
        end do
        end if
        !
        tmp=.0 ; k1=jcsc(j); k2=jcsc(j+1)-2 ;
        do k=k1, k2
        r=icsc(k) ;
        tmp = tmp+x(r)*csca(k) ;
        end do
        y(j) = y(j)+tmp ;
    end do
    !
    return
end subroutine cscax
!-------------------------------------------------------------------------------
subroutine lsolve(n, da1,icsc,jcsc,csca,b, x)
    !  This subroutine performs forward solve of MSSOR, that is,
    !                 (L+DA)x = b (provided for sqmrmssor subroutine).
    !              n: dimension of coefficient matrix A;
    ! icsc,jcsc,csca: CSC storage of upper triangular part of matrix A;
    !            da1: inverse of da (da: modified diagonal for MSSOR);
    !              b: it is right hand vector b;
    implicit none
    real(8):: da1(:), csca(:), b(:), x(:), tmp
    integer::i, j, k1, k2, n, icsc(:), jcsc(:)
    ! -------- forward substitution --------
    x(1)=b(1)*da1(1);
    do j=2, n
        k1=jcsc(j); k2=jcsc(j+1)-1 ; tmp=.0
        do i=k1, k2-1
            tmp = tmp + csca(i) * x(icsc(i)) ;
        end do
            tmp = b(j) - tmp
            x(j) = tmp*da1(j) ;
    end do
    !
    return
end subroutine lsolve
!-------------------------------------------------------------------------------
subroutine usolve(n, da1,icsc,jcsc,csca,b, x)
    !  This subroutine performs backward solve of MSSOR, that is,
    !                 (DA+U)x = b (provided for sqmrmssor subroutine).
    !              n: dimension of coefficient matrix A;
    ! icsc,jcsc,csca: CSC storage of upper triangular part of matrix A;
    !            da1: inverse of da (da: modified diagonal for MSSOR);
    !              b: it is right hand vector b;
    implicit none
    real(8):: da1(:),csca(:),b(:),x(:)
    real(8),allocatable:: tmp(:)
    integer::j,r,k,k1,k2,n,icsc(:), jcsc(:)
    allocate(tmp(n) )
    ! ----- backward substitution -----
    tmp = b ;
    do k = n, 2, -1
        x(k) = tmp(k)*da1(k)
        do j = jcsc(k), jcsc(k+1)-2
            r = icsc(j)
            tmp(r) = tmp(r) - x(k)*csca(j)
        end do
    end do
    x(1) = tmp(1)*da1(1)
    !
    return
end subroutine usolve
!-------------------------------------------------------------------------------
subroutine gtor(n,icsc,jcsc,csca,da,g,r)
    !  This subroutine performs backward solve of MSSOR, that is,
    !                 (L+DA)g = r (provided for sqmrmssor subroutine).
    !              n: dimension of coefficient matrix A;
    ! icsc,jcsc,csca: CSC storage of upper triangular part of matrix A;
    !             da: the modified diagonal for MSSOR;
    !              g: it is input vector(preconditioned residual);
    !              r: output vector (returned true residual).
    implicit none
    real(8):: da(:),csca(:),g(:),r(:)
    integer::i,j,n,k1,k2,icsc(:), jcsc(:)
    r=.0 ; r(1)=r(1)+da(1)*g(1);
    do j=2, n
        k1=jcsc(j); k2=jcsc(j+1)-2
        do i= k1, k2
        r(j)=r(j)+csca(i)*g(icsc(i))
        end do
        r(j)=r(j)+da(j)*g(j);
    end do
    !
    return
end subroutine gtor
!-------------------------------------------------------------------------------
subroutine kpu(icsc,jcsc,csca,theta,id,x,y)
    ! this subroutine performs the KP*u product for Biot's incremental
    !                 formula.
    ! icsc,jcsc,csca: CSC storage of upper triangular part of A.
    !          theta: implicit time-stepping parameter [0.5, 1].
    !          id: identifier array for displacement or porepressure DOF;
    !              id(i) = 1, displacement DOF;
    !              id(i) = 0, pore pressure DOF.
    !              x: x = u is the current excess pore pressure.
    implicit none
    real(8):: theta, csca(:), x(:), y(:)
    integer::j,k,r,k1,k2,n,jcsc(:),icsc(:),id(:)
        !
        n=ubound(jcsc,1)-1 ;  y =.0
        do j=1, n
        if( id(j)==0 ) then             ! pore pressure DOF
            if( x(j)/=.0)then
            k1=jcsc(j);  k2=jcsc(j+1)-1
            do k = k1, k2
                if(id(icsc(k))==0) y(icsc(k))=y(icsc(k)) + csca(k)*x(j)
            end do
            end if
            !
            k1=jcsc(j);  k2=jcsc(j+1)-2
            do k = k1, k2
            if(id(icsc(k))==0) y(j)=y(j) + csca(k)*x(icsc(k))
            end do
        end if
        end do
        !
        do j=1, n;
        if(id(j)==0) y(j)=-y(j)/theta ;
        end do
        !
        return
end subroutine kpu
!-------------------------------------------------------------------------------
!------------------ PRECONDITIONED SPARSE SOLVER SUBROUTINES--------------------
!-------------------------------------------------------------------------------
subroutine formda(n,icsc,jcsc,csca,icho,ipre,coef,omega,id,d,da,da1)
    ! This subroutine forms GJ diagonal vector - da (d and da1);
    ! In this routine:
    !              n: dimension of coefficient matrix A;
    ! icsc,jcsc,csca: CSC storage of upper triangular part of matrix A;
    !           icho: choose standard or modified preconditioenr.
    !               =1: standard preconditioner.
    !               =2: generalized or modified preconditioner. 
    !          ipre: choose preconditioner,
    !               =1: Jacobi preconditioner.
    !               =2: SSOR preconditioner (omega should be applied).
    !           coef: the scaling factor for GJ diagonal vector.
    !          omega: relaxation parameter, which is applied to MSSOR.
    !             id: a vector to indicate the type of current DOF,
    !                 id(j)= 0 for pore water pressure DOF;
    !                 id(j)= 1 for displacement DOF.
    !              d: diagonal of A;
    !             da: modified diagonal for MSSOR preconditioner;
    !            da1: inverse of da;
    implicit none
    real(8):: coef,omega,d(:),da(:),da1(:),csca(:),absv, maxabs,minabs
    integer::n,j,r,k,k1,k2,icho,ipre,id(:),icsc(:), jcsc(:)
        !
        do j=1, n;
        r=jcsc(j+1)-1 ; da(j) = csca(r);
        end do
        !
        if(ipre==2) d = da ; ! Transfer diagonal of A from da to d;
        if(icho == 2)then     ! For generalized or modified preconditioner
        !
        do j=2, n
            k1=jcsc(j) ;  k2=jcsc(j+1)-2
            if(id(j)==1) then
            do k=k1 , k2
                if(id(icsc(k))==0 ) then
                da(icsc(k))=da(icsc(k))-csca(k)**2/da(j) ;
                end if
            end do
            else                     ! id(j)==0
            do k=k1 , k2
                if(id(icsc(k))== 1 ) then
                da(j)=da(j)-csca(k)**2/da(icsc(k)) ;
                end if
            end do
            end if
        end do
        !
        coef = coef/omega ;
        do j=1, n          ! coef-scaling factor (negative is preferred)
            if(id(j)==0)then ! modified diagonal with relaxation parameter.
            da(j)=coef*abs(da(j))
            else             ! id(j)==0
            da(j)=da(j)/omega
            end if
        end do
        end if
        da1 = 1./da ;
        !
        return
end subroutine formda
!-------------------------------------------------------------------------------
subroutine sqmrmssor(n,icsc,jcsc,csca,d,da,da1,rhs,maxit,tol,icc, &
                        iinc,qmriters,relres)
    ! Modified SSOR preconditioned SQMR for symmetric Ax=b linear system.
    ! Combining with Eisenstat trick.
    ! In this routine:
    !              n: dimension of coefficient matrix A;
    ! icsc,jcsc,csca: CSC storage of upper triangular part of matrix A;
    !              d: diagonal of A;
    !             da: modified diagonal for MSSOR preconditioner;
    !            da1: inverse of da;
    !            rhs: at input, it is right hand vector b;
    !                 at output,it is returned approximate solution x;
    !          maxit: user-defined maximum iteration count;
    !            tol: it is the user-defined stopping tolerance;
    !            icc: choice for convergence criterion;
    !                 = 1, relative improvement norm criterion.
    !                 = 2, relative residual norm criterion (x0=.0)
    !           iinc: Check convergence every 'iinc' iteration.
    !       qmriters: the iterative count when SQMR converges;
    !         relres: the relative residual when SQMR converges.
    implicit none
    integer::i,n,maxit,icc,iinc,qmriters,ic,icsc(:),jcsc(:)
    real(8)::tol,relres,tao,theta0,theta,rho0,rho,phi,nrmb,nrmr,sigma,   &
            kesa,beta,rhs(:),csca(:),d(:),da(:),da1(:)
    real(8),allocatable::x(:), xold(:), z(:), r(:), g(:),v(:),w(:),c(:), &
                        t1(:),t2(:),t(:),p(:),q(:)
    allocate(x(n),xold(n),z(n),r(n),g(n),v(n),w(n),c(n),t1(n),t2(n),t(n),&
            p(n),q(n) )
    !---------- initialize vectors ------------
        x = .0 ; z = .0                           ! initial guess
        r=rhs;
        call lsolve(n,da1,icsc,jcsc,csca,r, g); ! preconditioned residual
        v =g;   w=da*v;
        tao = dot_product(g, g);
        rho0 = dot_product(g, w)
        c=.0 ; theta0=.0 ;
        nrmb=sqrt(dot_product(r, r))*tol ;
        p=d-2.0*da  ;                       ! p is used in each iteration
        ic = 0;
    !---------- SQMR iteration ----------
    iteration: do i=1, maxit
        !----- matrix-vector product with Eisenstat trick -----
                call usolve(n, da1,icsc,jcsc,csca,w, t1);
                t2 = p*t1+w ;
                call lsolve(n, da1,icsc,jcsc,csca,t2, t);
                t =  t1 + t;
        !-----------------------------------------------------
        sigma=dot_product(w, t);
        if(sigma==.0) then
            write(11,*) 'SQMR stops due to Sigma=0 ';  stop
        end if
        kesa=rho0/sigma ; g = g - kesa*t ;
        theta=dot_product(g, g)/tao ;
        phi=1./(1.+theta) ;
        tao=tao*theta*phi ;
        c=phi*(theta0*c+kesa*v) ;
        z = z + c ;
        select case (icc)
            case (1)
            call ccri(n,i,icsc,jcsc,csca,da,da1,z,iinc,tol, xold, &
                        rhs,ic,qmriters,relres)             
            case (2)
            call ccrb(n,i,icsc,jcsc,csca,da,da1,g,z,iinc,tol,nrmb,&
                        rhs,ic,qmriters,nrmr,relres)
        end select
        if(ic==1) return ;
        if(rho0==.0)then
            write(11,*) 'SQMR stops due to rho0=0 ' ; stop
        end if
        q = da * g ;
        rho=dot_product(g, q) ;
        beta = rho/rho0 ;  v= g + beta*v;
        w= da*v;
        theta0=theta;  rho0= rho
    end do iteration
    !--------- End iteration --------
    write(11,*)'********************************************** '
    write(11,*) 'SQMR does not converge to user-defined tolerance. '
        z = da * z;
        call usolve(n, da1,icsc,jcsc,csca,z, x)
        if(icc==1) relres = nrmr*tol/nrmb ;
        qmriters = maxit ;  rhs = x ;
    return
end subroutine sqmrmssor
!-------------------------------------------------------------------------------
subroutine ccrb(n,i,icsc,jcsc,csca,da,da1,g,z,iinc,tol,nrmb,rhs,ic, &
                iters,nrmr,relres)
    ! This subroutine performs convergence check in terms of relative
    !               residual with x0=.0 is chosen,
    !            n: i.e. neq - number of total DOFs or equations;
    !            i: current iteration # ;
    ! icsc,jcsc,csca: CSC storage of upper triangular part of matrix A;
    !           da: modified diagonal for MSSOR preconditioner;
    !          da1: inverse of da;
    !            g: preconditioned residual;
    !            z: "preconditioned" solution;
    !         iinc: Check convergence every 'iinc' iteration.
    !          tol: it is the user-defined stopping tolerance;
    !         nrmb: computed initial residual (b) norm multiplied by tol;
    !          rhs: at input, it is right hand vector b;
    !               at convcergence,it is returned approximate solution x;
    !           ic: = 1 converged (ic is a identifier);
    !               = 0 doesn't satisfy the convergence criterion;
    !        iters: returned iteration count when converged;
    !         nrmr: norm of true residual.
    !       relres: returned relative residual when converged.
    implicit none
    integer:: i,j,n,iinc,iters,ic,icsc(:),jcsc(:)
    real(8):: tol,nrmb,nrmr,relres,csca(:),da(:),da1(:),g(:),z(:),rhs(:)
    real(8),allocatable:: r(:)
    allocate( r(n) )
    !
        if(mod(i, iinc)==0)then       !  per iinc steps, check convergence
        call gtor(n,icsc,jcsc,csca,da,g,r) ; ! true residual is computed
        nrmr=sqrt(dot_product(r, r)) ;		   
        if(nrmr < nrmb)then          !  solver converged
            write(11,*)'********************************************** '
            write(11,*) 'SQMR converges to user-defined tolerance. '
            z = da * z ;
            call usolve(n,da1,icsc,jcsc,csca,z,rhs) ! rhs is the solution.
            ic =1 ; iters = i ; relres = nrmr*tol/nrmb ;
            return 
        end if
        end if
        !
        return
end subroutine ccrb
!-------------------------------------------------------------------------------
subroutine ccri(n,i,icsc,jcsc,csca,da,da1,z,iinc,tol,xold,rhs,ic,iters,&
                relres)
    ! This subroutine performs convergence check in terms of rleative
    !                 'improvement' norm criterion ( when set icc = 2 )
    !              n: i.e. neq - number of total DOFs or equations;
    !              i: current iteration # ;
    ! icsc,jcsc,csca: CSC storage of upper triangular part of matrix A;
    !             da: modified diagonal for MSSOR preconditioner;
    !            da1: inverse of da;
    !              z: "preconditioned" solution;
    !           iinc: Check convergence every 'iinc' iteration.(iinc > 1)
    !            tol: it is the user-defined stopping tolerance;
    !           xold: solution of the previous iterative step;
    !            rhs: at input, it is right hand vector b;
    !                 at convcergence, it is returned solution;
    !             ic: = 1 converged (ic is a identifier);
    !                 = 0 doesn't satisfy the convergence criterion;
    !          iters: returned iteration count when converged;
    !         relres: returned relative residual when converged.
    implicit none
    integer:: i,j,n,iinc,iters,ic,icsc(:),jcsc(:)
    real(8):: big,tol,ratio,relres,csca(:),da(:),da1(:),rhs(:),z(:),&
                xold(:)
    real(8),allocatable:: t(:)
    allocate(t(n) )
    !
        if( mod(i, iinc)==0 )then      !  per iinc steps, compute a xold;
        t = da * z ;
        call usolve(n, da1,icsc,jcsc,csca,t, rhs) ;
        xold = rhs ;
        !
        else if( mod(i, iinc)==1 )then ! check convergence, closely next
                                        !             per iinc iterations;
        t = da * z ;
        call usolve(n, da1,icsc,jcsc,csca,t, rhs) ;
        big=.0 ;
        do j = 1, n; if( abs(rhs(j)) > big ) big=abs(rhs(j));  end do
        relres =.0
        do j = 1, n;
            ratio = abs(rhs(j)-xold(j))/big;
            if(ratio > relres ) relres = ratio;
        end do
        if(relres < tol) ic=1;       !  Solver converged
        if(ic==1)then
            write(11,*)'********************************************** '
            write(11,*) 'SQMR converges to user-defined tolerance. '
            iters = i
            !
            return
        end if
        !
        end if
        !
        return
end subroutine ccri
!-------------------------------------------------------------------------------
subroutine psqmr(n,icsc,jcsc,csca,pr,rhs,maxit,tol,icc,qmriters,relres)
    ! This subroutine uses SQMR to solve Ax=b linear system with a right
    !                 diagonal preconditioner.
    ! In this routine:
    !              n: dimension of coefficient matrix A;
    ! icsc,jcsc,csca: CSC storage of coefficient matrix A;
    !             pr: right preconditioner (which is inverted at input);
    !            rhs: at input, it is right hand vector b;
    !                 at output,it is returned approximate solution x;
    !          maxit: user-defined maximum iteration count;
    !            tol: it is the user-defined stopping tolerance;
    !            icc: choice for convergence criterion;
    !                 = 1, relative improvement norm criterion
    !                 = 2, relative residual norm criterion (x0=.0)
    !             ic: indentifier of convergence;
    !                 = 1, solver converged;
    !                 = 0, not converge.
    !       qmriters: the iterative count when SQMR converges;
    !         relres: the relative residual when SQMR converges.
    implicit none
    integer::i,n,maxit,qmriters,icc,ic,icsc(:),jcsc(:)
    real(8)::tol,relres,tao,theta0,rho0,rho,nrmb,nrmr,sigma,alpha,beta,&
            theta,cj,rhs(:),csca(:),pr(:)
    real(8),allocatable::x(:),xold(:),r(:),t(:),q(:),d(:),u(:)
    allocate(x(n),xold(n),r(n),t(n),q(n),d(n),u(n)  )
    !------ Initial vectors of SQMR iterations ------
        x=.0                           ! assumed initial guess
        r=rhs;
        t=r;                           ! left preconditioning
        q=pr*t;                        ! right preconditioning
        tao=sqrt(dot_product(t, t) );
        theta0=.0;
        rho0=dot_product(r,q);
        nrmb=sqrt(dot_product(r, r))*tol;
        d=.0; ic=0 ; xold = x ;
    !------------ Sart SQMR iterations -------------
    iteration: do i=1, maxit
    call cscax(icsc,jcsc,csca,q,t)   !  t=A*q, Matrix-vector product.
    sigma=dot_product(q,t);
    alpha=rho0/sigma ;
    r = r-alpha*t ;
    t= r;                            ! left preconditioning
    theta=sqrt(dot_product(t, t) )/tao ;
    cj=1./sqrt(1+theta*theta);
    tao=tao*theta*cj;
    d=(cj*theta0)**2*d+(cj*cj)*alpha*q ;
    x = x + d;
    nrmr=sqrt(dot_product(r, r)) ;
    select case (icc)
        case (1) 
        call pccri(n,i,xold,x,rhs,tol,ic,qmriters,relres)
        xold = x ;
        case (2)
        call pccrb(n,i,r,x,rhs,tol,nrmb,ic,qmriters,nrmr,relres)                          
    end select
    if(ic==1) return ;          ! SQMR converged
    u=pr*t ;                    ! right preconditioning
    rho=dot_product(r,u);
    beta=rho/rho0;  q = u + beta*q ;
    !
    rho0=rho; theta0=theta; 
    end do iteration
    write(11,*)'********************************************** '
    write(11,*) 'SQMR does not converge to user-defined tolerance. '
    relres=nrmr*tol/nrmb ;  qmriters = maxit ;  rhs=x
    !
    return
end subroutine psqmr
!-------------------------------------------------------------------------------
subroutine pccrb(n,i,r,x,rhs,tol,nrmb,ic,iters,nrmr,relres)
    ! This subroutine performs convergence check in terms of relative
    !             residual with x0=.0 is chosen,
    !          n: i.e. neq - number of total DOFs or equations;
    !          i: current iteration # ;
    !          r: current residual;
    !        rhs: at input, it is right hand vector b;
    !             at convcergence,it is returned approximate solution x;
    !        tol: it is the user-defined stopping tolerance;
    !       nrmb: computed initial residual (b) norm multiplied by tol;
    !         ic: = 1 converged (ic is a identifier);
    !             = 0 doesn't satisfy the convergence criterion;
    !      iters: returned iteration count when converged;
    !       nrmr: norm of true residual.
    !     relres: returned relative residual when converged.
    implicit none
    integer:: i,j,n,iters,ic
    real(8):: tol,nrmb,nrmr,relres,r(:),x(:),rhs(:)
    !
        nrmr=sqrt(dot_product(r, r)) ;
        if(nrmr < nrmb)then          !  solver converged
        write(11,*)'********************************************** '
        write(11,*) 'Psolver converges to user-defined tolerance. '
        ic =1 ; iters = i ; relres = nrmr*tol/nrmb ; rhs = x;
        !
        return
        end if
        !
        return
end subroutine pccrb
!-------------------------------------------------------------------------------
subroutine pccri(n,i,xold,x,rhs,tol,ic,iters,relres)
    ! This subroutine performs convergence check in terms of rleative
    !                 'improvement' norm criterion ( when set icc = 2 )
    !              n: i.e. neq - number of total DOFs or equations;
    !              i: current iteration # ;
    !           xold: approximate solution of the previous iterative step;
    !              x: approximate solution of current iterative step;
    !            rhs: at input, it is right hand vector b;
    !                 at convcergence, it is returned solution;
    !            tol: it is the user-defined stopping tolerance;
    !             ic: = 1 converged (ic is a identifier);
    !                 = 0 doesn't satisfy the convergence criterion;
    !          iters: returned iteration count when converged;
    !         relres: returned relative residual when converged.
    implicit none
    integer:: i,j,n,iters,ic
    real(8):: big,tol,ratio,relres,rhs(:),xold(:),x(:)
    !
        big=.0 ;
        do j = 1, n; if( abs(x(j)) > big ) big=abs(x(j));  end do
        relres =.0
        do j = 1, n;
            ratio = abs(x(j)-xold(j))/big;
            if(ratio > relres ) relres = ratio;
        end do
        if(relres < tol) ic=1;       !  Solver converged
        if(ic==1)then
            write(11,*)'********************************************** '
            write(11,*) 'Psolver converges to user-defined tolerance. '
            iters = i ; rhs = x ;
            !
            return
        end if
        !
        return
end subroutine pccri
! ----------------------------------------------------------------------
! ------------- Preconditioned CG Method for Linear System Ax = b  -------------
! ----------------------------------------------------------------------
subroutine pcg(n,icsc,jcsc,csca,pr,rhs,maxit,tol,icc,iters,relres)
    ! This subroutine uses PCG to solve Ax=b linear system with a right
    !                 diagonal preconditioner.
    ! In this routine:
    !              n: dimension of coefficient matrix A;
    ! icsc,jcsc,csca: CSC storage of coefficient matrix A;
    !             pr: right preconditioner (which is inverted at input);
    !            rhs: at input, it is right hand vector b;
    !                 at output,it is returned approximate solution x;
    !          maxit: user-defined maximum iteration count;
    !            tol: it is the user-defined stopping tolerance;
    !            icc: choice for convergence criterion;
    !                 = 1, relative improvement norm criterion
    !                 = 2, relative residual norm criterion (x0=.0)
    !             ic: indentifier of convergence;
    !                 = 1, solver converged;
    !                 = 0, not converge.
    !          iters: the iterative count when PCG converges;
    !         relres: the relative residual when PCG converges.
    ! This subroutine is for preconditioned PCG iterative method.
    ! Refer to the PCG Algorithm by van der Vorst's book or Template Book.
    integer:: i,n,ic,icc,icsc(:),jcsc(:),maxit,iters
    real(8):: pr(:),csca(:),tol,rhs(:),nrmb,nrmr,rho,rho0,alpha,beta,relres
    real(8),allocatable:: r(:),z(:),p(:),q(:),x(:),xold(:)
    allocate (r(n),z(n),p(n),q(n),x(n),xold(n) )  
    !  b -- is RHS vector when inputting, while it is the Solution Vector 
    !  when returning.
        x = .0 ;  r = rhs ;             ! x0=0 is the initial solution guess
        nrmb=sqrt(dot_product(rhs, rhs))*tol; 
        ic=0 ; xold = x ;
    pcg_iter: do i=1, maxit
        z = pr*r ; rho = dot_product(r, z) ;
        if ( i > 1 )then               ! direction vector
            beta = rho/rho0;
            p = z + beta*p;
        else
            p = z; 
        end if
        call cscax(icsc,jcsc,csca,p,q)   !  q=Ap, Matrix-vector product.
        alpha = rho/dot_product(p, q) ;
        x = x + alpha * p ;
        r = r - alpha * q ;
        nrmr=sqrt(dot_product(r, r)) ;
        select case (icc)
            case (1) 
            call pccri(n,i,xold,x,rhs,tol,ic,iters,relres)
            xold = x ;
            case (2)
            call pccrb(n,i,r,x,rhs,tol,nrmb,ic,iters, nrmr,relres)
                        
        end select
        if(ic==1)return ;          ! PCG converged
        rho0 = rho
    end do pcg_iter
    write(11,*)'********************************************** '
    write(11,*) 'PCG does not converge to user-defined tolerance. '
    relres=nrmr*tol/nrmb ;  iters = maxit ;  rhs=x
    return
end subroutine pcg
! ----------------------------------------------------------------------
subroutine pcgmssor(n,icsc,jcsc,csca,d,da,da1,rhs,maxit,tol,icc,iinc, &
                    iters,relres)
    ! Modified SSOR preconditioned PCG for symmetric Ax=b linear system.
    ! Combining with Eisenstat trick.
    ! In this routine:
    !              n: dimension of coefficient matrix A;
    ! icsc,jcsc,csca: CSC storage of upper triangular part of matrix A;
    !              d: diagonal of A;
    !             da: modified diagonal for MSSOR preconditioner;
    !            da1: inverse of da;
    !            rhs: at input, it is right hand vector b;
    !                 at output,it is returned approximate solution x;
    !          maxit: user-defined maximum iteration count;
    !            tol: it is the user-defined stopping tolerance;
    !            icc: choice for convergence criterion;
    !                 = 1, relative improvement norm criterion.
    !                 = 2, relative residual norm criterion (x0=.0)
    !           iinc: Check convergence every 'iinc' iteration.
    !          iters: the iterative count when PCG converges;
    !         relres: the relative residual when PCG converges.
    implicit none
    integer:: i,n,icc,ic,iinc,icsc(:),jcsc(:),maxit,iters
    real(8):: csca(:),d(:),da(:),da1(:),tol,rhs(:), nrmb, nrmr, rho,rho0,&
            alpha,beta,relres
    real(8),allocatable:: tmp(:),x(:),xold(:),z(:),t(:),t1(:),t2(:),r(:),&
            g(:),p(:),q(:)          
    allocate(tmp(n),x(n),xold(n),z(n),t(n),t1(n),t2(n),r(n),g(n),p(n),q(n))
    !---------- initialize vectors ------------
    !  b is RHS vector when inputting, while it is the Solution Vector  &
    !    when returning.
    x = 0 ; z =.0;  r = rhs ;       ! x0=0 is the initial solution guess
    nrmb=sqrt(dot_product(r, r))*tol;  
    call lsolve(n,da1,icsc,jcsc,csca,r, g);    ! preconditioned residual
    tmp = d - 2.0*da ; ic=0 ;  
    pcg_iter: do i=1, maxit
        t = da*g ; rho = dot_product(g,t) ;        !(t = z ; g =r)
        if(i==1)then
        p = t;
        else
        beta = rho/rho0 ; p = t + beta*p ;
        end if
        !  q=Ap, Matrix-vector product.
        !----- matrix-vector product with Eisenstat trick -----
        call usolve(n, da1,icsc,jcsc,csca,p, t1);
        t2 = tmp*t1 + p ;
        call lsolve(n, da1,icsc,jcsc,csca,t2, q);
        q =  t1 + q;
        !-----------------------------------------------------
        alpha = rho/dot_product(p, q) ;
        z = z + alpha * p ;
        g = g - alpha * q ;
        select case (icc)
            case (1)
            call cgccri(n,i,icsc,jcsc,csca,da,da1,z,iinc,tol, xold, &
                        rhs,ic,iters,relres)             
            case (2)
            call cgccrb(n,i,icsc,jcsc,csca,da,da1,g,z,iinc,tol,nrmb,&
                        rhs,ic,iters,nrmr,relres)
        end select
        if(ic==1) return ;
        rho0 = rho
    end do pcg_iter
    !
    write(11,*)'********************************************** '
    write(11,*) 'PCG does not converge to user-defined tolerance. '
        call usolve(n, da1,icsc,jcsc,csca,z, x)
        if(icc==1) relres = nrmr*tol/nrmb ;
        iters = maxit ;  rhs = x ;
    return
end subroutine pcgmssor
!-------------------------------------------------------------------------------
subroutine cgccrb(n,i,icsc,jcsc,csca,da,da1,g,z,iinc,tol,nrmb,rhs,ic, &
                iters,nrmr,relres)
    ! This subroutine performs convergence check in terms of relative
    !             residual with x0=.0 is chosen,
    !          n: i.e. neq - number of total DOFs or equations;
    !          i: current iteration # ;
    ! icsc,jcsc,csca: CSC storage of upper triangular part of matrix A;
    !         da: modified diagonal for MSSOR preconditioner;
    !        da1: inverse of da;
    !          g: preconditioned residual;
    !          z: "preconditioned" solution;
    !       iinc: Check convergence every 'iinc' iteration.
    !        tol: it is the user-defined stopping tolerance;
    !       nrmb: computed initial residual (b) norm multiplied by tol;
    !        rhs: at input, it is right hand vector b;
    !             at convcergence,it is returned approximate solution x;
    !         ic: = 1 converged (ic is a identifier);
    !             = 0 doesn't satisfy the convergence criterion;
    !      iters: returned iteration count when converged;
    !       nrmr: norm of true residual.
    !     relres: returned relative residual when converged.
    implicit none
    integer:: i,j,n,iinc,iters,ic,icsc(:),jcsc(:)
    real(8):: tol,nrmb,nrmr,relres,csca(:),da(:),da1(:),g(:),z(:),rhs(:)
    real(8),allocatable:: r(:)
    allocate( r(n) )
    !
    if(mod(i, iinc)==0)then        !  per iinc steps, check convergence
        call gtor(n,icsc,jcsc,csca,da,g,r) ; ! true residual is computed
        nrmr=sqrt(dot_product(r, r)) ;
        if(nrmr < nrmb)then          !  solver converged
            write(11,*)'********************************************** '
            write(11,*) 'PCG converges to user-defined tolerance. '
            ! z = da * z ;
            call usolve(n,da1,icsc,jcsc,csca,z,rhs) ! rhs is the solution.
            ic =1 ; iters = i ; relres = nrmr*tol/nrmb ;
            return 
        end if
    end if
    !
    return
end subroutine cgccrb
!-------------------------------------------------------------------------------
subroutine cgccri(n,i,icsc,jcsc,csca,da,da1,z,iinc,tol,xold,rhs,ic, &
                    iters,relres)
    ! This subroutine performs convergence check in terms of rleative
    !                 'improvement' norm criterion ( when set icc = 2 )
    !              n: i.e. neq - number of total DOFs or equations;
    !              i: current iteration # ;
    ! icsc,jcsc,csca: CSC storage of upper triangular part of matrix A;
    !             da: modified diagonal for MSSOR preconditioner;
    !            da1: inverse of da;
    !              z: "preconditioned" solution;
    !           iinc: Check convergence every 'iinc' iteration.(iinc>1)
    !            tol: it is the user-defined stopping tolerance;
    !           xold: solution of the previous iterative step;
    !            rhs: at input, it is right hand vector b;
    !                 at convcergence, it is returned solution;
    !             ic: = 1 converged (ic is a identifier);
    !                 = 0 doesn't satisfy the convergence criterion;
    !          iters: returned iteration count when converged;
    !         relres: returned relative residual when converged.
    implicit none
    integer:: i,j,n,iinc,iters,ic,icsc(:),jcsc(:)
    real(8):: big,tol,ratio,relres,csca(:),da(:),da1(:),rhs(:),z(:), &
                xold(:)
    real(8),allocatable:: t(:)
    allocate(t(n) )
    !
        if( mod(i, iinc)==0 )then      !  per iinc steps, compute a xold;
        ! t = da * z ;
        call usolve(n, da1,icsc,jcsc,csca,z, rhs) ;
        xold = rhs ;
        !
        else if( mod(i, iinc)==1 )then ! check convergence, closely next
                                        !             per iinc iterations;
        !t = da * z ;
        call usolve(n, da1,icsc,jcsc,csca,z, rhs) ;
        big=.0 ;
        do j = 1, n; if( abs(rhs(j)) > big ) big=abs(rhs(j));  end do
        relres =.0
        do j = 1, n;
            ratio = abs(rhs(j)-xold(j))/big;
            if(ratio > relres ) relres = ratio;
        end do
        if(relres < tol) ic=1;       !  Solver converged
        if(ic==1)then
            write(11,*)'********************************************** '
            write(11,*) 'PCG converges to user-defined tolerance. '
            iters = i
            !
            return
        end if
        !
        end if
        !
        return
end subroutine cgccri
! ------------------------------------------------------------------------------
! ---------- MINRES method for symmetric and possible indefinite Ax=b ----------
! ------------------------------------------------------------------------------
subroutine minres(n,icsc,jcsc,csca,pr,rhs,maxit,tol,icc,iters,relres)
    ! This subroutine uses MINRES to solve Ax=b linear system with a right
    !                 diagonal preconditioner.
    ! In this routine:
    !              n: dimension of coefficient matrix A;
    ! icsc,jcsc,csca: CSC storage of coefficient matrix A;
    !             pr: right preconditioner (which is inverted at input);
    !            rhs: at input, it is right hand vector b;
    !                 at output,it is returned approximate solution x;
    !          maxit: user-defined maximum iteration count;
    !            tol: it is the user-defined stopping tolerance;
    !            icc: choice for convergence criterion;
    !                 = 1, relative improvement norm criterion
    !                 = 2, relative residual norm criterion (x0=.0)
    !             ic: indentifier of convergence;
    !                 = 1, solver converged;
    !                 = 0, not converge.
    !          iters: the iterative count when MINRES converges;
    !         relres: the relative residual when MINRES converges.
    ! This subroutine is Jacobi preconditioned MINRES iterative method.
    ! Refer to the MINRES Algorithm. (Wang,2004)
    integer:: i,n,icc,ic,icsc(:),jcsc(:),maxit,iters
    real(8):: pr(:),csca(:),tol,rhs(:),nrmb,nrmr,rho,rho0,alpha,beta,relres
    real(8),allocatable::r(:),u(:),p(:),q(:),x(:),xold(:),z0(:),z(:),t0(:),&
                t(:)
    allocate (r(n),u(n),p(n),q(n),x(n),xold(n),z0(n),z(n),t0(n),t(n) )  
    x = 0;  r = rhs ! x0=0 is the initial guess.               
    p = pr*r ; z0 = p ; 
    call cscax(icsc,jcsc,csca,p,q)   !  q=Ap, Matrix-vector product.
    t0 = q; nrmb = sqrt(dot_product(r, r))*tol ;
    ic=0 ; xold = x ;
    mr_iter: do i=1, maxit
        u = pr * q ; 
        alpha = dot_product(z0,q)/dot_product(q,u);
        x = x + alpha * p ;
        r = r - alpha * q ;
        z = z0 - alpha * u ;
        call cscax(icsc,jcsc,csca,z,t)   !  t = Az, Matrix-vector product.
        beta = dot_product(z, t)/dot_product(z0, t0) ;
        p = z + beta*p ; 
        q = t + beta*q ;
        nrmr = sqrt(dot_product(r, r));
        select case (icc)
            case (1) 
            call pccri(n,i,xold,x,rhs,tol,ic,iters,relres)
            xold = x ;
            case (2)
            call pccrb(n,i,r,x,rhs,tol,nrmb,ic,iters,nrmr,relres)       
        end select
        if(ic==1) return ;          ! MINRES converged
        z0 = z ; t0 = t ;
    end do mr_iter
    write(11,*)'********************************************** '
    write(11,*) 'MINRES does not converge to user-defined tolerance. '
    relres=nrmr*tol/nrmb ;  iters = maxit ;  rhs=x
    return
end subroutine minres
!-------------------------------------------------------------------------------
subroutine mrmssor(n,icsc,jcsc,csca,d,da,da1,rhs, maxit, tol, icc,iinc,&
                    iters,relres)
    ! Modified SSOR preconditioned MINRES for symmetric Ax=b linear system.
    ! Combining with Eisenstat trick.
    ! In this routine:
    !              n: dimension of coefficient matrix A;
    ! icsc,jcsc,csca: CSC storage of coefficient matrix A;
    !             pr: right preconditioner (which is inverted at input);
    !            rhs: at input, it is right hand vector b;
    !                 at output,it is returned approximate solution x;
    !          maxit: user-defined maximum iteration count;
    !            tol: it is the user-defined stopping tolerance;
    !            icc: choice for convergence criterion;
    !                 = 1, relative improvement norm criterion
    !                 = 2, relative residual norm criterion (x0=.0)
    !             ic: indentifier of convergence;
    !                 = 1, solver converged;
    !                 = 0, not converge.
    !          iters: the iterative count when MINRES converges;
    !         relres: the relative residual when MINRES converges.
    ! This subroutine is Jacobi preconditioned MINRES iterative method.
    ! Refer to the MINRES Algorithm. (Wang,2004)
    integer::i,n,icc,ic,iinc,icsc(:),jcsc(:),maxit,iters
    real(8)::d(:),da(:),da1(:),csca(:),tol,rhs(:),nrmb,nrmr,rho,rho0,alpha,&
                beta,relres
    real(8),allocatable::tmp(:),r(:),g(:),u(:),p(:),q(:),x(:),xold(:),y(:),&
                y0(:),z(:),t0(:),t1(:),t2(:),t(:)
    allocate (tmp(n),r(n),g(n),u(n),p(n),q(n),x(n),xold(n),y(n),y0(n),z(n),&
                t0(n),t1(n),t2(n),t(n) )  
    x = .0 ;  z = .0 ;  r = rhs ;                 
    ! x0 = 0 is the initial guess, y is the preconditioned solution.
    call lsolve(n,da1,icsc,jcsc,csca,r,g); ! g -preconditioned residual
    p = da*g ; y0 = p ; tmp = d - 2.0*da ;
    !----- matrix-vector product q=A*p with Eisenstat trick ------
    call usolve(n, da1,icsc,jcsc,csca,p,t1);
    t2 = tmp*t1+p ;
    call lsolve(n, da1,icsc,jcsc,csca,t2,q);
    q =  t1 + q; 
    !-------------------------------------------------------------
    t0 = q ;  ic=0 ; nrmb = sqrt(dot_product(r, r))*tol ;
    mr_iter: do i=1, maxit
        u = da * q ; 
        alpha = dot_product(y0,q)/dot_product(q,u);
        z = z + alpha * p ;
        g = g - alpha * q ;
        y = y0 - alpha * u ;
        !----- matrix-vector product t=A*y with Eisenstat trick -----
        call usolve(n, da1,icsc,jcsc,csca,y, t1);
        t2 = tmp*t1+y ;
        call lsolve(n, da1,icsc,jcsc,csca,t2, t);
        t = t + t1; 
        !-------------------------------------------------------------
        beta = dot_product(y, t)/dot_product(y0, t0) ;
        p = y + beta*p ; 
        q = t + beta*q ;
        select case (icc)
            case (1)
                call cgccri(n,i,icsc,jcsc,csca,da,da1,z,iinc,tol, xold, &
                        rhs,ic,iters,relres)             
            case (2)
                call cgccrb(n,i,icsc,jcsc,csca,da,da1,g,z,iinc,tol,nrmb,&
                        rhs,ic,iters,nrmr,relres)
        end select
        if(ic==1) return ;          ! MINRES converged
        y0 = y ; t0 = t ;
    end do mr_iter
    write(11,*)'********************************************** '
    write(11,*) 'MINRES does not converge to user-defined tolerance. '
        call usolve(n, da1,icsc,jcsc,csca,z, x)
        if(icc==1) relres = nrmr*tol/nrmb ;
        iters = maxit ;  rhs = x ;
    return
end subroutine mrmssor
!-------------------------------------------------------------------------------
subroutine psolver(n,icsc,jcsc,csca,d,da,da1,b,maxit,tol,isolver,icho, &
                    ipre,icc,iinc,iters,relres)
    ! Choose preconditioned iterative methods:
    !              n: i.e. neq - number of total DOFs or equations;
    ! icsc,jcsc,csca: CSC storage of upper triangular part of matrix A;
    !              d: true diagonal of A;
    !             da: modified diagonal basing on GJ;
    !            da1: inverse of da;
    !              b: right hand side vector;
    !          maxit: user-defined maximal iteration number;
    !            tol: user-defined stopping tolerance;
    !        isolver: Iterative solver selection;
    !                 = 1, SQMR iterative solver;
    !                 = 2, PCG iterative solver;
    !                 = 3, MINRES iterative solver;
    !           icho: choose standard or modified preconditioenr.
    !               =1: standard preconditioner.
    !               =2: generalized or modified preconditioner. 
    !          ipre: choose preconditioned iterative solver;
    !                 = 1, GJ preconditioned iterative method;
    !                 = 2, MSSOR preconditioned iterative method;
    !            icc: choose convergence criterion;
    !                 = 1, relative improvement norm criterion;
    !                 = 2, relative residual norm criterion;
    !           iinc: check convergence every "iinc" step when ipre = 2;
    !          iters: returned iteration number when converged;
    !         relres: returned relative residual when converged;
    implicit none
    integer::n,icsc(:),jcsc(:),maxit,isolver,icho,ipre,icc,iinc,iters
    real(8)::csca(:),d(:),da(:),da1(:),b(1:),tol, relres
    !real(8),allocatable::d1(:)
    !allocate(d1(n))
    !d1(1:n) = 1./d(1:n) ;
    if(isolver==1)then         ! (SQMR Iterative Solver)
        select case (ipre)
        case(1)
        if(icho==1) write(11,*) '  ---> SJ preconditioned SQMR solver'
        if(icho==2) write(11,*) '  ---> GJ preconditioned SQMR solver'
            select case (icc)
            case (1)
                write(11,*) '   with relative improvement norm criterion!'
            case (2)
                write(11,*) '   with relative residual norm criterion!'
            end select
            !
            call psqmr(n,icsc,jcsc,csca,da1,b,maxit,tol,icc,iters,relres)
        case(2)
        if(icho==1)write(11,*) ' ---> SSOR preconditioned SQMR solver'
        if(icho==2)write(11,*) ' ---> MSSOR preconditioned SQMR solver'
            select case (icc)
            case (1)
                write(11,*) '   with relative improvement norm criterion!'
            case (2)
                write(11,*) '   with relative residual norm criterion!'
            end select
            call sqmrmssor(n,icsc,jcsc,csca,d,da,da1,b,maxit,tol,icc, &
                            iinc,iters,relres)
        case default
            write(*,*) ' No preconditioned solver is chosen, stop here! '
            Stop
        end select
    else if(isolver==2)then  ! (PCG Iterative Solver)
        select case (ipre)
        case(1)
            if(icho==1)write(11,*) ' ---> SJ preconditioned PCG solver'
            if(icho==2)write(11,*) ' ---> GJ preconditioned PCG solver'
            select case (icc)
            case (1)
                write(11,*) '   with relative improvement norm criterion!'
            case (2)
                write(11,*) '   with relative residual norm criterion!'
            end select
            !		      
            call pcg(n,icsc,jcsc,csca,da1,b,maxit,tol,icc,iters,relres)
        case(2)
            if(icho==1)write(11,*) ' ---> SSOR preconditioned PCG solver'
            if(icho==2)write(11,*) ' ---> MSSOR preconditioned PCG solver'
            select case (icc)
            case (1)
                write(11,*) '   with relative improvement norm criterion!'
            case (2)
                write(11,*) '   with relative residual norm criterion!'
            end select
            call pcgmssor(n,icsc,jcsc,csca,d,da,da1,b,maxit,tol,icc, &
                            iinc,iters,relres)
        case default
            write(*,*) ' No preconditioned solver is chosen, stop here! '
            Stop
        end select
    else if(isolver==3)then   ! isolver =3  ! (MINRES Iterative Solver)
        select case (ipre)
        case(1)
            if(icho==1)write(11,*) ' ---> SJ preconditioned MINRES solver'
            if(icho==2)write(11,*) ' ---> GJ preconditioned MINRES solver'
            select case (icc)
            case (1)
                write(11,*) '   with relative improvement norm criterion!'
            case (2)
                write(11,*) '   with relative residual norm criterion!'
            end select
            !
            call minres(n,icsc,jcsc,csca,da1,b,maxit,tol,icc,iters,relres)
        case(2)
        if(icho==1)write(11,*) ' ---> SSOR preconditioned MINRES solver'
        if(icho==2)write(11,*) ' --->MSSOR preconditioned MINRES solver'
            select case (icc)
            case (1)
                write(11,*) '   with relative improvement norm criterion!'
            case (2)
                write(11,*) '   with relative residual norm criterion!'
            end select
            call mrmssor(n,icsc,jcsc,csca,d,da,da1,b,maxit,tol,icc, &
                            iinc,iters,relres)
        case default
            write(*,*) '  No preconditioned solver is chosen, stop here! '
            Stop
        end select
        !
    else
        write(11,*) '  No Iterative Solver is selected, and STOP! ' ;
        stop
    end if
    !
    return
end subroutine psolver
!
end module sparse_lib