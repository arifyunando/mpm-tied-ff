MODULE fem
  SAVE
CONTAINS
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
  
  
  SUBROUTINE get_ff_displacement(mpm_disp, mpm_counter, iel_boundary, mpm_g_num, mpm_nf, ff_disp, ff_g_num, ff_nf)
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(OUT)::mpm_disp(0:)
    INTEGER,INTENT(OUT)::mpm_counter(0:)
    REAL(iwp),INTENT(IN)::ff_disp(0:)
    INTEGER,INTENT(IN)::iel_boundary(:,:),mpm_g_num(:,:),ff_g_num(:,:),mpm_nf(:,:),ff_nf(:,:)
    INTEGER::i, iel_bc_mpm, iel_bc_ff
    INTEGER,ALLOCATABLE::g_mp(:), g_ff(:)
    mpm_disp = 0.0_iwp
    mpm_counter = 0
    DO i=1,size(iel_boundary, 2)
      ! get boundary element index in ff & mpm
      iel_bc_ff  = iel_boundary(1,i)*2 ! 2 column for each row
      iel_bc_mpm = iel_boundary(2,i)
      
      ! if there is no element number, skip!
      IF(iel_bc_ff==0) CYCLE
      
      ! get g in mpm & ff from element number
      g_mp = mpm_nf(1,mpm_g_num(:,iel_bc_mpm))
      g_ff = ff_nf(1,ff_g_num(:,iel_bc_ff))
    
      ! take displacement with from ff to mpm and add counter
      mpm_disp(g_mp) = mpm_disp(g_mp) + ff_disp(g_ff)
      mpm_counter(g_mp) = mpm_counter(g_mp) + 1

      ! get g in mpm & ff from element number
      g_mp = mpm_nf(2,mpm_g_num(:,iel_bc_mpm))
      g_ff = ff_nf(2,ff_g_num(:,iel_bc_ff))
    
      ! take displacement with from ff to mpm and add counter
      mpm_disp(g_mp) = mpm_disp(g_mp) + ff_disp(g_ff)
      mpm_counter(g_mp) = mpm_counter(g_mp) + 1
    END DO
    ! cleanup memory
    mpm_disp(0)=0.0_iwp; mpm_counter(0)=0
    
    ! normalize displacement
    DO i=0, size(mpm_disp)-1
      IF(mpm_counter(i) > 0) mpm_disp(i) = mpm_disp(i) / mpm_counter(i)
    END DO
    mpm_disp(0)=0.0_iwp; mpm_counter(0)=0
  END SUBROUTINE get_ff_displacement
  
   
END MODULE fem