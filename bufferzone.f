        subroutine bufferzone(np_initial)
        include 'common.2D'

        nb_start = nb+1
        do i_num_Buff = 1,num_Buff
           if (ioptBuff(i_num_Buff).eq.2) then
              nbavoids = nBuff_ini(i_num_Buff)
              nbavoide = nBuff_end(i_num_Buff)
           endif 
        enddo

        do i_num_Buff = 1, num_Buff

           do ii = nb_start,np_initial          !np_initial is np from
                                                !previous timestep

              if (itype(ii).eq.4) then

                  up(ii) = uo(ii)+grx*iflag(ii)*dt
                  wp(ii) = wo(ii)+grz*iflag(ii)*dt

                 if (ioptBuff(i_num_Buff).eq.1.and.
     &               (ii.lt.nbavoids.or.ii.gt.nbavoide)) then

                     zp(ii) = zp(ii) + wo(ii)*dt
                     xp(ii) = xp(ii) + uo(ii)*dt
     

                     if (xp(ii).ge.Buffint(i_num_Buff)) then
             
                        itype(ii) = 2

                        !!! Create new particle !!!
                        np = np + 1
              
                        itype(np) = 4
                        rhop(np) = rhop(ii)
                        pp(np) = cs0*cs0*(rhop(np)-rho0)
                        pm(np) = pm(ii)
                        iflag(np) = 1
           
                        xp(np) = xp(ii) - (2*ikappa+1)*dx0
                        zp(np) = zp(ii) 
!                        up(np) =(-7.7295*zp(ii)*zp(ii)+0.0116*zp(ii))
!     &                          *1000
                        up(np) = VXbuff(i_num_Buff)
                        wp(np) = VZbuff(i_num_Buff)

                     end if

                 elseif (ioptBuff(i_num_Buff).eq.2.and.
     &                   (ii.ge.nbavoids.or.ii.le.nbavoide)) then

                     xp(ii) = xp(ii)
                     zp(ii) = zp(ii)

!                     rrmin = rrk2h2

!                     do kk = nb_start,np_initial

!                        if (itype(kk).eq.2) then

!                           xxdiff = (xp(ii) - xp(kk))*(xp(ii) - xp(kk))
!                           zzdiff = (zp(ii) - zp(kk))*(zp(ii) - zp(kk))
!                           rrdiff = xxdiff + zzdiff

!                           if (rrdiff.lt.rrmin) then

!                              iipart(ii) = kk
!                              rrmin = rrdiff
        
!                           endif

!                        endif

!                     enddo

!                     iipart(ii) = kk
                     rhop(ii) = rho0!rhop(kk)
                     pp(ii) = p0!pp(kk)
!                     pp(ii) = cs0*cs0*rhop(ii)
                     nbavoids = nBuff_ini(i_num_Buff)
                     nbavoide = nBuff_end(i_num_Buff)

                 endif

               endif

           enddo

        enddo

        end subroutine bufferzone
