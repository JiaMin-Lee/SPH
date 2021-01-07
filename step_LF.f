
        subroutine step
        include 'common.2D'

        character supp*4,name*40,detsupp*4
        character tname1*40,tname2*40

!!!     ----- LEAP-FROG METHOD -----  !!!

        dt2 = 0.5*dt

        call check_limits_2D
        call ini_divide(2)
        call divide(nvmb+1,np,2)

        if (nbf+1.lt.nbfm) then
           call recover_list
           call divide(nbf+1,nvmb,1)
        end if
        
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

        if (itime.ne.0) then            ! nitime .ne. 1

           do ip = nrho_start,np

              rhop(ip) = rhop(ip) + dt2 * rdot(ip)
              pVol(ip) = pm(ip)/rhop(ip)
              TEp(ip) = TEp(ip) + dt2 * TEdot(ip)
              call equation_of_state(rhop(ip), TEp(ip), pp(ip), cs(ip))

              if (itype(ip).eq.2) then

                 up(ip) = up(ip) + dt2 * udot(ip)
                 wp(ip) = wp(ip) + dt2 * wdot(ip)
                
              endif

           enddo

        end if

        call ac

        call correct


        if (itime.eq.0) then    ! nitime .eq. 1

           do ip = nrho_start,np

              rhop(ip) = rhop(ip) + dt2 * rdot(ip)
              pVol(ip) = pm(ip)/rhop(ip)
              TEp(ip) = TEp(ip) + dt2 * TEdot(ip)
              call equation_of_state(rhop(ip), TEp(ip), pp(ip), cs(ip))

              if (itype(ip).eq.2) then

                 up(ip) = up(ip) + dt2 * udot(ip)
                 wp(ip) = wp(ip) + dt2 * wdot(ip)
                
                 xp(ip) = xp(ip) + dt * xdot(ip)
                 zp(ip) = zp(ip) + dt * zdot(ip)

              end if

           enddo

        else    ! nitime .gt. 1

           do ip = nrho_start,np

              rhop(ip) = rhoo(ip) + dt2 * rdot(ip)
              pVol(ip) = pm(ip)/rhop(ip)
              TEp(ip) = TEo(ip) + dt2 * TEdot(ip)
              call equation_of_state(rhop(ip), TEp(ip), pp(ip), cs(ip)) 

              if (itype(ip).eq.2) then

                 up(ip) = uo(ip) + dt * udot(ip)
                 wp(ip) = wo(ip) + dt * wdot(ip)
                
                 xp(ip) = xo(ip) + dt * xdot(ip)
                 zp(ip) = zo(ip) + dt * zdot(ip)

              end if

           enddo

        end if


        if (nvmb+1.lt.nb) then

           call virtualmovingboundary

        endif

        if (num_Buff.gt.0) then

           call bufferzone(np)

        end if

        end subroutine step

