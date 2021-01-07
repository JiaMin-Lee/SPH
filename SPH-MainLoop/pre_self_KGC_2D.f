        subroutine pre_self(j1,kind_p1,ini_kind_p2)
        include 'common.2D'
        
        do kind_p2=ini_kind_p2,2

           jj_start = 1

           if (kind_p1.gt.kind_p2) then
              jj_start = nplink_max + 1      
           endif

           do ii = 1,nc(j1,kind_p1)

              i = ibox(j1,kind_p1,ii)

              if (kind_p1.eq.kind_p2) then
                 jj_start = jj_start + 1
              end if
            
              do jj=jj_start,nc(j1,kind_p2)
            
                 j = ibox(j1,kind_p2,jj)
            
                 drx = xp(i) - xp(j)
                 drz = zp(i) - zp(j)
                 rr2 = drx*drx + drz*drz
              
                 if (rr2.lt.rrk2h2.and.rr2.gt.1.e-18) then
                 
                    call kernel(drx,drz,i,j,j1,j1,rr2)             
                    
                    V_j = pVol(j) 
                    frxh = V_j*frx
                    frzh = V_j*frz
            
                    !frxhi * (xp(j) - xp(i)) = frxh * -(drx)
                    aM_a11(i) = aM_a11(i) - frxh*drx
                    aM_a12(i) = aM_a12(i) - frxh*drz
                    aM_a21(i) = aM_a21(i) - frzh*drx
                    aM_a22(i) = aM_a22(i) - frzh*drz
            
                    V_i = pVol(i)
                    frxh = V_i*frx
                    frzh = V_i*frz
            
                    !frxhj * (xp(i) - xp(j)) = -frxh * drx
                    aM_a11(j) = aM_a11(j) - frxh*drx
                    aM_a12(j) = aM_a12(j) - frxh*drz
                    aM_a21(j) = aM_a21(j) - frzh*drx
                    aM_a22(j) = aM_a22(j) - frzh*drz 

                 endif ! if q<2
              enddo  ! Box jj
           enddo ! Box ii
        enddo ! Kind of particle
 
        return

        end subroutine pre_self

