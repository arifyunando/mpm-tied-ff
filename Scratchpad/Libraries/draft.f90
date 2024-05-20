mbod(bod)%kv=zero
    mbod(bod)%mv=zero
    mbod(bod)%kp=zero
    KM_MV:DO k=1,mbod(bod)%nmps
      ! This loop goes until 4 since 4 is the maximum number of 
      ! elements affected by a material point
      El_stiff:DO i=1,4 
        iel = mbod(bod)%elemmpoints(k,i)
        ! If is true, an element is affected by a material point 
        ! (material point can be outside the element)
        Full_el:IF(mbod(bod)%elemmpoints(k,i)>0)THEN 
          num = g_num(:,iel)
          nod_num(:,1) = num
          coord = TRANSPOSE(g_coord(:,num))
          g = mbod(bod)%g_g(:,iel)

          km_gauss = zero   
          waverage = 4.0/nip

          !-Double mapping technique
          Double_map:IF(smethod>1)THEN 

            CALL iGIMP_funder3(                                 &
              k,mbod(bod)%mpoints,mbod(bod)%lp_mp,nip,coord,  &
              cellsize,mbod(bod)%gm_coord,mbod(bod)%gimptol,  &
              nod_num,g_coord,mbod(bod)%a_ele,mbod(bod)%c_ele,&
              iel,der,fun                                     &
            )

            IF(fun(1)<=0.or.fun(2)<=0.or.fun(3)<=0.or.fun(4)<=0) THEN
              PRINT*,'Stiffness problem'
              PRINT*,'shape functions',fun(:)
              PRINT*,''
              PRINT*,'material point',k
              PRINT*,''
              PRINT*,'element',iel
              PRINT*,''
              PRINT*,'lp',iel
              PAUSE
            END IF    

            CALL sample(element,points,weights)
            DO s=1,nip
              CALL shape_der(der_gauss,points,s)
              CALL shape_fun(fun_gauss,points,s)
              scalar_dee = fun*fun_gauss*waverage
              sum_scalar = SUM(scalar_dee)
              dee_scal = mbod(bod)%dee*sum_scalar
              jac = MATMUL(der_gauss,coord)
              det = determinant(jac)
              CALL invert(jac)
              deriv_gauss = MATMUL(jac,der_gauss)
              CALL beemat(bee_gauss,deriv_gauss)
              IF(mbod(bod)%c_ele(iel)<1) waverage=1.5_iwp
              km_gauss = km_gauss + MATMUL(MATMUL(TRANSPOSE(bee_gauss),dee_scal),bee_gauss)*det*weights(i)
            END DO
            CALL fsparv(mbod(bod)%kv,km_gauss,g,mbod(bod)%kdiag)
          END IF Double_map
        END IF Full_el
      END DO El_stiff
      !-----------END STIFNESS MATRICES---------------

      !---Diagonal mass matrix---  
      values=mbod(bod)%valuesg(k)
      ALLOCATE(                                                           &
          derextend(nodof,values),funextend2(values*2,2),                 &
          eldddylds(values*2),mm_gimp(values*2,values*2)                  &
      )
      CALL GIMP_funder2(                                                  &
          k,nip,g_coord,cellsize,mbod(bod)%gm_coord,mbod(bod)%lp_mp,      &
          mbod(bod)%GIMP_nodes,mbod(bod)%gimptol,derextend,funextend2     &
      )

      mm_gimp = zero
      mvval   = 2
      CALL gimpfunform(                                                   &
          k,iel,eldddylds,mbod(bod)%nf,mbod(bod)%GIMP_nodes,              &
          values,mbod(bod)%g_g,mvval                                      &
      ) 

      mvval=1
      a=1
      DO j=1,values
          DO q=1,2
              mm_gimp(a,a)=funextend2(a,q)*mbod(bod)%m_mass(k)
              a=a+1
          END DO
      END DO
      CALL fsparv(mbod(bod)%mv,mm_gimp,eldddylds,mbod(bod)%kdiag)
      DEALLOCATE(derextend,funextend2,eldddylds,mm_gimp)
    END DO KM_MV