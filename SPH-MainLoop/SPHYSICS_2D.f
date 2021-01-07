        program SPHYSICS

        include 'common.2D'

        character supp*4,name*40,detsupp*4

        REAL time_begin, time_end
        CALL CPU_TIME (time_begin)

!!!     -- OPEN FILES

        open(80,file='sph.out')
        write(80,*) 'Opening sph.out'

!        open(19,file='DT')
!        open(50,file='ENERGY')
!        open(54,file='MovingBodyHistory.out')

        pi=4.*atan(1.)
        g = grz

        call getdata
        call energy(g)

!!!     -- BINNED FIXED PARTICLES FIRST --> SO NO NEED TO KEEP BINNING

        call ini_divide(1)
        call divide(1,nbf,1)
        call keep_list

        if (i_restartRun.eq.1.or.i_restartRun.eq.3) then
           open(19,file='DT',status='old',access='append')
        else
           open(19,file='DT')
           write(19,*) 'time dt1 dt2 dt_new'
        endif
        close(19)

        tdetail=0
        ngrabdet=0
        grab_P=0
        grab_E=0

        write(*,*) ' '
        write(*,*) 'g= ',g
        write(*,*) 'grx= ',grx
        write(*,*) 'grz= ',grz
        write(80,*) ' '
        write(80,*) 'g= ',g
        write(80,*) 'grx= ',grx
        write(80,*) 'grz= ',grz
        write(80,*) ''

        itimestart = itime
        nmaxtimestep = int(tmax/dt)
        nprintstep = int(tout/dt)

        write(*,*) ' nmaxtimestep = ', nmaxtimestep 
        write(*,*) ' nprintstep = ', nprintstep

        do nitime = itimestart+1,nmaxtimestep

           if (mod(nitime,nprintstep).eq.0) then
              ipoute = 1
           else
              ipoute = 0
           endif

           call step

           itime = nitime
           time = itime * dt

!           write(*,*) itime, time
!           write(*,*) ''

           if (ipoute.eq.1) then
              ngrab = ngrab+1
              write(supp,'(i4.4)') ngrab
              name = 'PART_'//supp
              write(80,*) name
              write(*,*) name
              open(23, file=name)

              call poute(23)

              close(23)

              open(44,file='RESTART')
              write(44,100) nitime, nitime*dt, ngrab, dt, np
              close(44)

              grap_P = 0

           else 
              grap_P = grap_P+dt
           endif

100     format(i8,e16.8,i8,e16.8,i8)

        if (i_densityFilter.gt.0.and.
     &        mod(itime,ndt_FilterPerform).eq.0)then
!     &        itime/ndt_FilterPerform*ndt_FilterPerform.eq.itime)then
            call densityFilter
        endif

        enddo

        open(33, file='EPART')
        call poute(33)
        close(33)

        open(333, file='IFLAG')
        do ip = 1, np
           write(333,101) iflag(ip)
        enddo
        close(333)
101     format(i8)

!!!     -- CLOSE FILES
        close(19)
        close(50)
!        close(54)

        close(80)

        end program
