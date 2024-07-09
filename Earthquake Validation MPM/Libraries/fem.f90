module fem
  save
contains
  SUBROUTINE ecmat2(ecm,fun,ndof,nodof)
    !
    ! This subroutine forms the element consistent mass matrix.
    !
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::fun(:)
    REAL(iwp),INTENT(OUT)::ecm(:,:)
    INTEGER,INTENT(IN)::nodof,ndof
    INTEGER::nod,i,j
    REAL::nt(ndof,nodof),tn(nodof,ndof),zero=0.0_iwp
    ecm=zero
    nod=ndof/nodof
    nt=zero
    tn=zero
    DO i=1,nod
      DO j=1,nodof
        nt((i-1)*nodof+j,j)=fun(i)
        tn(j,(i-1)*nodof+j)=fun(i)
      END DO
    END DO
    
    ecm=MATMUL(nt,tn)
    RETURN
  END SUBROUTINE ecmat2

  SUBROUTINE bantmul(kb,loads,ans)
    !
    ! This subroutine multiplies an unsymmetrical band kb by the vector loads.
    ! Could be much improved for vector processors.
    !
     IMPLICIT NONE
     INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
     REAL(iwp),INTENT(IN)::kb(:,:),loads(0:)
     REAL(iwp),INTENT(OUT)::ans(0:)
     INTEGER::i,j,k,l,m,n,iw
     REAL(iwp)::x,zero=0.0_iwp
     n=SIZE(kb,1)
     l=SIZE(kb,2)
     iw=(l-1)/2
     DO i=1,n
       x=zero
       k=iw+2
       DO j=1,l
         k=k-1
         m=i-k+1
         IF(m<=n.AND.m>=1)x=x+kb(i,j)*loads(m)
       END DO
       ans(i)=x
     END DO
    RETURN
  END SUBROUTINE bantmul
  
  SUBROUTINE formlump(diag,emm,g)
    !
    ! This subroutine forms the lumped global mass matrix as a vector diag.
    !
     IMPLICIT NONE
     INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
     REAL(iwp),INTENT(IN)::emm(:,:)
     REAL(iwp),INTENT(OUT)::diag(0:)
     INTEGER,INTENT(IN)::g(:)
     INTEGER::i,ndof
     ndof=UBOUND(emm,1)
     DO i=1,ndof
       diag(g(i))=diag(g(i))+emm(i,i)
     END DO
    RETURN
  END SUBROUTINE formlump
  
  SUBROUTINE formtb(pb,km,g)
    !
    ! This subroutine assembles an unsymmetrical band matrix pb from
    ! element constituent matrices km.
    !
     IMPLICIT NONE
     INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
     REAL(iwp),INTENT(IN)::km(:,:)
     INTEGER,INTENT(IN)::g(:)
     REAL(iwp),INTENT(OUT)::pb(:,:)
     INTEGER::i,j,idof,icd,iw
     idof=SIZE(km,1)
     iw=(SIZE(pb,2)-1)/2
     DO i=1,idof
       IF(g(i)/=0)THEN
         DO j=1,idof
           IF(g(j)/=0)THEN
             icd=g(j)-g(i)+iw+1
             pb(g(i),icd)=pb(g(i),icd)+km(i,j)
           END IF
         END DO
       END IF
     END DO
    RETURN
  END SUBROUTINE formtb
end module fem