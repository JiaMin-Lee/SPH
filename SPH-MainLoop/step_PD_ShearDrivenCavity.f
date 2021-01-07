        subroutine step
        include 'common.2D'

        character supp*4,name*40,detsupp*4
        character tname1*40,tname2*40

!!!     ----- PREDICTOR-CORRECTOR METHOD -----  !!!

        dt2 = 0.5*dt

        call check_limits_2D
        call ini_divide(2)
        call divide(nvmb+1,np,2)

        if (nbf+1.lt.nbfm) then
           call recover_list
           call divide(nbf+1,nvmb,1)
        end if
        
        call ac

        call correct

!!!     -----   PREDICTOR   -----    !!!

        do ip = nstart,np

           xo(ip) = xp(ip)
           zo(ip) = zp(ip)
           uo(ip) = up(ip)
           wo(ip) = wp(ip)      
           po(ip) = pp(ip)
           rhoo(ip) = rhop(ip)
           TEo(ip) = TEp(ip)

        enddo 

        if (iBC.eq.3.or.iBC.eq.1) nrho_start = nvmb+1
        if (iBC.eq.4) nrho_start = nstart

        do ip = nrho_start,np

           rhop(ip) = 1000!rhoo(ip) + dt2*rdot(ip)
           pVol(ip) = pm(ip)/rhop(ip)
           TEp(ip) = TEo(ip) + dt2*TEdot(ip)
           call equation_of_state(rhop(ip), TEp(ip), pp(ip), cs(ip))

        enddo

        do ip = 1,np
!        do ip = nb+1,np
           
           if (itype(ip).eq.2) then

              xp(ip) = xo(ip) + dt2*xdot(ip)
              zp(ip) = zo(ip) + dt2*zdot(ip)

              up(ip) = uo(ip) + dt2*udot(ip)
              wp(ip) = wo(ip) + dt2*wdot(ip)

           elseif (iBC.eq.4.and.itype(ip).lt.0) then
!           elseif (iBC.eq.4.and.abs(itype(ip)).eq.1) then

              up(ip) = uo(ip) + dt2*udot(ip)
              wp(ip) = wo(ip) + dt2*wdot(ip)

           end if

        enddo

!        if (nvmb+1.lt.nb) then

!           call virtualmovingboundary

!        endif 

!!!     -----   CORRECTOR   -----    !!!     


        call check_limits_2D
        call ini_divide(2)
        call divide(nvmb+1,np,2)

        if (nbf+1.lt.nbfm) then
           call recover_list
           call divide(nbf+1,nvmb,1)
        end if

        call ac

        call correct

        do ip = nrho_start,np

           rhop(ip) = 1000!rhoo(ip) + dt2*rdot(ip)
           TEp(ip) = TEo(ip) + dt2*TEdot(ip)

        enddo

        do ip = 1,np
!        do ip = nb+1,np

           if (itype(ip).eq.2) then

              xp(ip) = xo(ip) + dt2*xdot(ip)
              zp(ip) = zo(ip) + dt2*zdot(ip)

              up(ip) = uo(ip) + dt2*udot(ip)
              wp(ip) = wo(ip) + dt2*wdot(ip)

           elseif (iBC.eq.4.and.itype(ip).lt.0) then
!           elseif (iBC.eq.4.and.abs(itype(ip)).eq.1) then

              up(ip) = uo(ip) + dt2*udot(ip)
              wp(ip) = wo(ip) + dt2*wdot(ip)

           endif

        enddo

!!!     -----   FINAL INTEGRATION   -----    !!!

        do ip = nrho_start,np

!           rhop(ip) = 2.*rhop(ip) - rhoo(ip)
           rhop(ip) = 1000
           pVol(ip) = pm(ip) / rhop(ip)
           TEp(ip) = 2.*TEp(ip) - TEo(ip)
           call equation_of_state(rhop(ip), TEp(ip), pp(ip), cs(ip))

        enddo

        do ip = 1,np
!        do ip = nb+1,np
        
           if (itype(ip).eq.2) then

              xp(ip) = 2.*xp(ip) - xo(ip)
              zp(ip) = 2.*zp(ip) - zo(ip)

              up(ip) = 2.*up(ip) - uo(ip)
              wp(ip) = 2.*wp(ip) - wo(ip)

              if (iVelocity.eq.1.and.xp(ip).le.(ikappa*hsml).and.
     &            xp(ip).ge.0.) then

!                 if (zp(ip).gt.0..and.zp(ip).le.0.2e-3) then
!                    up(ip) = 29./2.*zp(ip)
!                 elseif(zp(ip).lt.1.5e-3.and.zp(ip).ge.1.3e-3) then
!                    up(ip) = 29./2.*abs(zp(ip)-1.5e-3)
!                 else
!                    up(ip) = VXInlet
!                 endif

                up(ip) = (-7.7295*zp(ip)*zp(ip)+0.0116*zp(ip))*1000
!                 wp(ip) = VZInlet


!                 up(ip) = VXInlet
!                 wp(ip) = VZInlet

              end if
        
           elseif (iBC.eq.4.and.itype(ip).lt.0) then
!           elseif (iBC.eq.4.and.abs(itype(ip)).eq.1) then

              up(ip) = uo(ip) + dt2*udot(ip)
              wp(ip) = wo(ip) + dt2*wdot(ip)

           end if 

        enddo 

        if (nvmb+1.lt.nb) then

           call virtualmovingboundary

        endif

        if (num_Buff.gt.0) then

           call bufferzone(np)

        end if 

        end subroutine step
