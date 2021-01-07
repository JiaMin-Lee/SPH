        subroutine pre_celij(j1,j2,kind_p1,ini_kind_p2,lx2)

        include 'common.2D' 

        do kind_p2=ini_kind_p2,2

           if (nc(j2,kind_p2).ne.0) then

              do ii=1,nc(j1,kind_p1)

                 i = ibox(j1,kind_p1,ii)
    
                 do jj=1,nc(j2,kind_p2)

                    j = ibox(j2,kind_p2,jj)
         
                    drx = xp(i) - xp(j)
                    drz = zp(i) - zp(j)
   
                    call periodicityCorrection(i,j,drx,drz,lx2)
 
                    rr2 = drx*drx + drz*drz
              
                    if (rr2.lt.rrk2h2.and.rr2.gt.1.e-18) then
         
                       call kernel(drx,drz,i,j,j1,j2,rr2)  
                  
                       V_j = pVol(j)
                       frxh = V_j*frx
                       frzh = V_j*frz
                    
                       aM_a11(i) = aM_a11(i) - frxh*drx
                       aM_a12(i) = aM_a12(i) - frxh*drz
                       aM_a21(i) = aM_a21(i) - frzh*drx
                       aM_a22(i) = aM_a22(i) - frzh*drz
                    
                       V_i = pVol(i)  
                       frxh = V_i*frx
                       frzh = V_i*frz
                    
                       aM_a11(j) = aM_a11(j) - frxh*drx
                       aM_a12(j) = aM_a12(j) - frxh*drz 
                       aM_a21(j) = aM_a21(j) - frzh*drx 
                       aM_a22(j) = aM_a22(j) - frzh*drz
 
                    endif ! if q<2 
                 enddo ! Box jj
              enddo ! Box ii
           endif  ! Box jj is not empty
        enddo   ! Kind of particle
 
        return

        end subroutine pre_celij
