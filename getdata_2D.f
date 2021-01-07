      subroutine getdata
      character*10  ch1,ch2,ch3,ch4,ch5,ch6,ch7,ch8,ch9,ch10,ch11,ch12

      include 'common.2D'

      character supp_ini*4,name_ini*40,normals_ini*40
      character (LEN=40) :: paddle_fileName, MovingObjectsFile 


      open(11, file='INDAT')

      read(11,*) i_restartRun
      read(11,*) i_FlowDir
      read(11,*) iVelocity
      read(11,*) i_kernel
      read(11,*) i_kernelcorrection
      read(11,*) i_densityFilter
      read(11,*) ndt_FilterPerform
      read(11,*) i_EoS
      read(11,*) iBC
      read(11,*) coef
      read(11,*) p0
      read(11,*) cs0
      read(11,*) rho0
      read(11,*) igamma
      read(11,*) ikappa
      read(11,*) hsml
      read(11,*) grx
      read(11,*) grz
      read(11,*) i_periodicOBs(1)
      read(11,*) i_periodicOBs(2)
      read(11,*) VXInlet
      read(11,*) VZInlet
      read(11,*) dt
      read(11,*) tmax
      read(11,*) tout
      read(11,*) viscos_val
      read(11,*) vlx
      read(11,*) vly
      read(11,*) dx
      read(11,*) dz
      read(11,*) nbf
      read(11,*) nbfm
      read(11,*) nvmb
      read(11,*) nb
      read(11,*) np

      dx0 = dx
      dz0 = dz
      eps = 0.5

      write(*,*) ''
      if (i_FlowDir.eq.1) write(*,*)'Flowing Direction is towards RIGHT'
      if (i_periodicOBs(1).eq.1) write(*,*) 'X-Direction Periodic BC'
      if (i_FlowDir.eq.2) write(*,*)'Flowing Direction is UPWARDS'
      if (i_periodicOBs(2).eq.1) write(*,*) 'Z-Direction Periodic BC'

      if (np.ge.npar) then
         write(80,*) 'Number of particles exceeds npar in common!'
         write(80,*) ' '
         write(*,*) 'Number of particles exceeds npar in common!'
         write(*,*) ' '
         stop
      endif

      if (nb.ge.nb_max) then
         write(80,*) 'Number of boundary particles exceeds nb_max in'
         write(80,*) 'common'
         write(80,*) ' '
         write(*,*) 'Number of boundary particles exceeds nb_max in'
         write(*,*) 'common'
         write(*,*) ' '
         stop
      endif

      write(80,*) 'nbf = ',nbf,'nbfm = ', nbfm,'nvmb = ', nvmb
      write(80,*) 'nb = ', nb, 'np = ', np
      write(80,*) 'hsml = ', hsml
      write(80,*) 'dx0 = ',dx0, 'dz0 = ', dz0
      write(80,*) 'p0 = ', p0, 'rho0 = ', rho0, 'mu = ',viscos_val
      write(80,*) 'dt = ', dt, 'tmax = ',tmax

      write(*,*) 'nbf = ',nbf,'nbfm = ', nbfm,'nvmb = ', nvmb
      write(*,*) 'nb = ', nb, 'np = ', np
      write(*,*) 'hsml = ', hsml
      write(*,*) 'dx0 = ',dx0, 'dz0 = ', dz0
      write(*,*) 'p0 = ', p0, 'rho0 = ', rho0, 'mu = ',viscos_val
      write(*,*) 'dt = ', dt, 'tmax = ',tmax

      close(11)


      !!! Accounting for PARTICLES NUMBERING
      nbfp = nbfp + 1   !Non-fixed boundary particles onwards

      nbp1 = nb + 1

!!!   BOUNDARY CONDITION

      if (iBC.eq.1) then
         beta_coef2 = (coef/10.)**2
         beta_coef = (coef/10.)
         one_over_beta_coef2 = 1.0/beta_coef2
         one_over_beta_coef = 1.0/beta_coef
      else
         one_over_beta_coef2 = 0.0
         one_over_beta_coef = 0.0
      endif


!!!   RESTART RUN

      if (i_restartRUN.ge.1) then

         i_RESTART_exists = 1
         open(44,file='RESTART', err=45)
         read(44,100,err=45) itime,time,ngrab,dt,np
         close(44)

         i_RESTART_exists = i_RESTART_exists + 1
45       i_RESTART_exists = i_RESTART_exists - 1

         if (i_RESTART_exists.eq.1.and.itime.gt.0) then

            write(supp_ini,'(i4.4)') ngrab
            name_ini='PART_'//supp_ini
            write(*,*) name_ini
            normals_ini = 'NORMALS.init'
            MovingObjectsFile = 'Floating_Bodies.RESTART'
         else
            itime = 0
            ngrab = 0
            grab_P = 0.0
            grab_E = 0.0
            name_ini = 'IPART'
            write(*,*) name_ini
            normals_ini = 'NORMALS.init'
            MovingObjectsFile = 'FloatingBodies.txt'
         endif
     
      else

         itime = 0
         ngrab = 0
         grab_P = 0.0
         grab_E = 0.0
         name_ini = 'IPART'
         write(*,*) name_ini
         normals_ini = 'NORMALS.init'
         MovingObjectsFile = 'FloatingBodies.txt'

      endif
100   format(i8,e16.8,i8,e16.8,i8)

      write(80,*) ' '
      write(80,*) 'i_restartRUN ',i_restartRUN
      write(80,*) 'Reading in ', name_ini
      write(80,*) 'Using dt ', dt
      write(*,*) ' '
      write(*,*) 'i_restartRUN ',i_restartRUN
      write(*,*) 'Reading in ', name_ini
      write(*,*) 'Using dt ', dt



!!!   READ INITIAL FILE TO LOAD PARTICLE DATA

      ! Initialise thermal energy
      do ii = 1,np
         TEp(ii) = 0.
      enddo

      ! Read initial file
      open(13,file=name_ini)

196   format(7e16.8,i5,i5)

      do ii = 1,np

         read(13,196) xp(ii),zp(ii),up(ii),wp(ii),rhop(ii),pp(ii),
     &                pm(ii),itype(ii),iflag(ii)

      enddo

!!!   READ NORMALS FILE
      if (iBC.eq.1) then

!         open(15,file = normals_ini)

!         do ii = 1,nb

!            read(15,*)
!     &      xnb(ii), znb(ii)
!     &      iBP_Pointer_Info(ii,1),
!     &      iBP_Pointer_Info(ii,2),
!     &      iBP_Pointer_Info(ii,3),
!     &      iBP_Pointer_Info(ii,4),
!     &      BP_xz_Data(ii,1),
!     &      BP_xz_Data(ii,2)

!         enddo

!         close(15)
!         call updatesNormals_2D(1,nb)

         !- Angular Momentum Correction
         two_alpha_h = 2.0*1.3
         rNum_h = 3.0*hsml
         one_over_rNum_h = 1.0/rNum_h
         one_over_r_ij_dot_n = 1.3/hsml

      endif

!!! READ MOVING PARTICLES INFORMATION
      open(77, file = MovingObjectsFile)
      
      read(77,*) iopt_FloatingBodies
      i_FB_Pointer_Info = 0
      
      if (iopt_FloatingBodies.ge.1) then
         read(77,*) nbfm
         print*
         print*, ' --- READING MOVING PARTICLES INFORMATION --- '
         i_FB_finish = 0
         i_ini = nbfm + 1
         num_FB_counter = 0
        
         do while (i_FB_finish.eq.0)

            read(77,*,END = 76) i_num_FB
            
            if (i_num_FB.gt.num_FB_max) then
               print*, 'Number of moving particles exceeds allowed max '
     &                  ,num_FB_max
               print*
               stop
            endif

            read(77,*)iopt_FB_type(i_num_FB)
            read(77,*)bigMass(i_num_FB)
            read(77,*)bigInertiaYY(i_num_FB)
            read(77,*)XcylinderDimension(i_num_FB),
     &                ZcylinderDimension(i_num_FB)
            read(77,*)cylinderDensity(i_num_FB)
            read(77,*)FB_specificDensity(i_num_FB)
            read(77,*)friction_coeff(i_num_FB)
            read(77,*)Box_XC(i_num_FB),Box_ZC(i_num_FB)
            read(77,*)bigU(i_num_FB),bigW(i_num_FB),bigOmegaY(i_num_FB)
            read(77,*)nb_FB(i_num_FB)

            print*,' Particle ', i_num_FB, ' Information '
            if (iopt_FB_type(i_num_FB).eq.1) print*, ' SQUARE Particle' 
            if (iopt_FB_type(i_num_FB).eq.2) print*, ' CIRCLE Particle'
            print*,' bigMass ', bigMass(i_num_FB)
            print*,' bigInertiaYY ', bigInertiaYY(i_num_FB)
            print*,' XcylinderDimension ', XcylinderDimension(i_num_FB)
            print*,' ZcylinderDimension ', ZcylinderDimension(i_num_FB)
            print*,' Particle Density ', cylinderDensity(i_num_FB)
            print*,' Specific Densit  ', FB_specificDensity(i_num_FB)
            print*,' Friction Coefficient ', friction_coeff(i_num_FB)
            print*
            print*,' INITIAL CONDITIONS: '
            print*,' X- and Z- CG ', Box_XC(i_num_FB),Box_ZC(i_num_FB)
            print*,' X- and Z- Linear Velocity ',
     &               bigU(i_num_FB),bigW(i_num_FB)
            print*,' Angular Velocity, OmegaY ', bigOmegaY(i_num_FB)
            print*,' nn_start, nn_end ',i_ini, nb_FB(i_num_FB)

            do ip = i_ini, nb_FB(i_num_FB)

               i_FB_Pointer_Info(ip) = i_num_FB
               xp_minus_R0(ip) = xp(ip) - Box_XC(i_num_FB)
               zp_minus_R0(ip) = zp(ip) - Box_ZC(i_num_FB)

            enddo

            ! --- INITIALISE TO ZERO FORCE --- !
            BodyForce_x = 0
            BodyForce_z = 0
            
            X_nonFriction = 0
            Z_nonFriction = 0

            i_ini = nb_FB(i_num_FB) + 1
            num_FB_counter = num_FB_counter + 1
            u_parallel_max = cs0*1.0e-4
            one_over_u_parallel_max = 1.0/u_parallel_max

            if (i_FB_finish.eq.0) then
               i_FB_finish = i_FB_finish - 1
            end if
76             i_FB_finish = i_FB_finish + 1
            print*

         enddo

         num_FB = num_FB_counter
         print*, ' Number of Moving Particle in Fluid Domain ', num_FB
         print*, ' Particles Indices ', nbfm+1, nb_FB(num_FB)
         print*
         print*,'--- Moving Particle Data Loaded ---'

      endif
      close(77)

      if (iopt_FloatingBodies.ge.1) then
         iopt_movingObject = 1
      elseif (iopt_FloatingBodies.eq.0) then
         iopt_movingObject = 0
      endif
      
      write(80,*) ' '
      write(80,*) 'Moving Object ',iopt_movingObject
      write(*,*) ' '
      write(*,*) 'Moving Object ',iopt_movingObject

!!!   Particles incorporated in interactions
      if (iBC.eq.1.and.nb.eq.nvmb) then
         nstart = nb+1
      elseif (iBC.eq.1.and.nvmb.lt.nb) then
         nstart = nvmb + 1
      elseif (iBC.eq.2) then
         nstart = 1
      elseif (iBC.eq.3.or.iBC.eq.4) then
         nstart = 1
      endif

      write(80,*) ' '
      write(80,*) 'nstart = ',nstart
      write(*,*) ' '
      write(*,*) 'nstart = ',nstart


!!!   Create initial corresponding internal particle list for Virtual
!!!   Moving Boundary Condition

      open(12,file = 'VirtualMovingBoundary')
      read(12,*) iopt_VMB

      if (iopt_VMB.gt.0) then
        
         print*
         print*, '--- Reading Virtual Moving Boundary Input File ---'

         i_VMB_finish = 0
         i_ini = nbfm + 1
         num_MB = 0

        do while(i_VMB_finish.eq.0)

           num_MB = num_MB + 1

           read(12,*,END=296)iVMBwall(num_MB)
           read(12,*) nVMB_ini(num_MB)
           read(12,*) nVMB_l1(num_MB)
           read(12,*) nVMB_end(num_MB)
           read(12,*) nVMB_lt(num_MB)
           read(12,*) VXvel(num_MB), VZvel(num_MB)
!           if (iVMBwall(num_MVB).eq.1) then
              
!           endif 

           print*, 'num_MB ', num_MB
           print*, 'VMB Boundary Position ', iVMBwall(num_MB)
           print*, 'Initial Particle Index of num_MB ', nVMB_ini(num_MB)
           print*, 'Last Particle Index of FIRST Layer', nVMB_l1(num_MB)
           print*, 'Last Particle Index of num_MB ', nVMB_end(num_MB)
           print*, 'Total number of layer ', nVMB_lt(num_MB)

           !- Exit loop
           if (i_VMB_finish.eq.0) then
              i_VMB_finish = i_VMB_finish - 1
           endif
296        i_VMB_finish = i_VMB_finish + 1
           print*

        enddo 

        num_MB = num_MB - 1

        print*,'TOTAL Number of Virtual Moving Boundary ', num_MB
        print*
        print*,'--- Virtual Moving Boundary Data Loaded ---'

      endif 
      close(12)

!!!   Create initial corresponding internal particle list for Virtual
!!!   Moving Boundary Condition

      open(14,file = 'BufferZone')
      read(14,*) iopt_Buff

      if (iopt_Buff.eq.1) then
        
         print*
         print*, '--- Reading Buffer Zone Input File ---'

         i_Buff_finish = 0

        do while(i_Buff_finish.eq.0)

           read(14,*,END=297) num_Buff
           read(14,*) ioptBuff(num_Buff)
           read(14,*) nBuff_ini(num_Buff)
           read(14,*) nBuff_end(num_Buff)
           read(14,*) VXbuff(num_Buff), VZbuff(num_Buff)
           read(14,*) Buffint(num_Buff)

           print*, 'num_Buff ', num_Buff
           print*, 'Buffer Zone Inlet/Outlet ', ioptBuff(num_Buff)
           print*, 'Initial Particle of num_Buff ', nBuff_ini(num_Buff)
           print*, 'Last Particle of num_Buff ', nBuff_end(num_Buff)
           print*, 'Velocity Inlet/Out, VX, VZ ', VXbuff(num_Buff),
     &                                            VZbuff(num_Buff)
           print*, 'Interface X-coordinate ', Buffint(num_Buff)

           !- Exit loop
           if (i_Buff_finish.eq.0) then
              i_Buff_finish = i_Buff_finish - 1
           endif
297        i_Buff_finish = i_Buff_finish + 1
           print*

        enddo 

        print*,'TOTAL Number of Buffer Zone ', num_Buff
        print*
        print*,'--- Buffer Zone Data Loaded ---'

      else

        num_Buff = 0

      endif 
      close(14)

      write(80,*)''
      write(80,*)'-----------------------------------------------------'
      write(80,*)'All particles information loaded! '
      write(*,*)''
      write(*,*)'-----------------------------------------------------'
      write(*,*)'All particles information loaded! '

!!!   VISCOSITY TERM
      eta = 0.1 * hsml
      eta2 = eta * eta      

!!!   DEFINE DOMAIN

      xmin = xp(1)
      xmax = xp(1)
      zmin = zp(1)
      zmax = zp(1)
      
      do ii=2,nbf
        
         if (xp(ii).lt.xmin) xmin = xp(ii)
         if (xp(ii).gt.xmax) xmax = xp(ii)
         if (zp(ii).lt.zmin) zmin = zp(ii)
         if (zp(ii).gt.zmax) zmax = zp(ii)

      enddo

      ! -- Domain
      xmax_container = xmax !+ ikappa*hsml
      xmin_container = xmin !- ikappa*hsml
      zmin_container = zmin !- ikappa*hsml
      zmax_container = zmax !+ ikappa*hsml

      if (i_FlowDir.eq.1.and.i_periodicOBs(1).eq.0) then
          xmin_container = xmin_container - ikappa*hsml
      elseif (i_FlowDir.eq.2.and.i_periodicOBs(2).eq.0) then
          zmin_container = zmin_container - ikappa*hsml
      endif

      if(i_periodicOBs(1).eq.1)then
         xmax_container = xmax + dx0   !BDR
      endif

      write(80,*) ' '
      write(80,*) 'xmin_container ', xmin_container
      write(80,*) 'xmax_container ', xmax_container
      write(80,*) 'zmin_container ', zmin_container
      write(80,*) 'zmax_container ', zmax_container
      write(80,*) ' '
      write(*,*) ' '
      write(*,*) 'xmin_container ', xmin_container
      write(*,*) 'xmax_container ', xmax_container
      write(*,*) 'zmin_container ', zmin_container
      write(*,*) 'zmax_container ', zmax_container
      write(*,*) ' '

      do ii=nbf+1,np

         if (xp(ii).lt.xmin) xmin = xp(ii)
         if (xp(ii).gt.xmax) xmax = xp(ii)
         if (zp(ii).lt.zmin) zmin = zp(ii)
         if (zp(ii).gt.zmax) zmax = zp(ii)

      enddo


      ! -- Determine number of mesh
      one_over_h = 1.0/hsml
      one_over_kh = 1.0/(ikappa*hsml)

      xmax = xmax + 0.1*hsml
      xmin = xmin - 0.1*hsml
      zmin = zmin - 0.1*hsml
      zmax = zmax + 0.1*hsml
      zmax_ini = zmax           ! used to keep moving obstacles inside
      xmax_ini = xmax           ! the initial boundaries

      ncx = int( (xmax-xmin) * one_over_kh ) + 1
      ncz = int( (zmax-zmin) * one_over_kh ) + 1

      xmax_container_double = dble(xmax_container)
      xmin_container_double = dble(xmin_container)

      if(i_periodicOBs(1).eq.1)then
         ncx_min = int((xmax_container-xmin)*one_over_kh)
         if(ncx.gt.ncx_min) ncx = ncx_min
         write(80,*)'Corrected ncx for X-Periodicity'
         write(80,*) 'ncx, ncx_min', ncx, ncx_min
         write(80,*) 'xmin, xmax_container', xmin, xmax_container
         write(80,*) ''
         write(*,*) 'Corrected ncx for X-Periodicity'
         write(*,*) 'ncx, ncx_min', ncx, ncx_min
         write(*,*) 'xmin, xmax_container', xmin, xmax_container
         write(*,*) ''
       endif

      nct = ncx * ncz
      ncx_ini = ncx
      ncz_ini = ncz
      nct_ini = nct

      if (i_FlowDir.eq.2) then
         zcontrol = zmin + ncz*ikappa*hsml
         write(80,*) 'zcontrol', zcontrol, 'zmax_ini', zmax_ini
      elseif (i_FlowDir.eq.1) then
         xcontrol = xmin + ncx*ikappa*hsml
         write(80,*) 'xcontrol', xcontrol, 'xmax_ini', xmax_ini
      endif

      if (nct.gt.nct_max) then

         write(80,*) ' '
         write(80,*) 'No. of mesh exceeds nct_max in common.2D'
         write(80,*) 'nct', nct, 'nct_max', nct_max
         write(80,*) ' '
         write(*,*) ' '
         write(*,*) 'No. of mesh exceeds nct_max in common.2D'
         write(*,*) 'nct', nct, 'nct_max', nct_max
         write(*,*) ' '
         stop

      end if
      
      write(80,*) 'xmin = ', xmin, 'xmax = ', xmax
      write(80,*) 'zmin = ', zmin, 'zmax = ', zmax
      write(80,*) 'ncx = ',ncx, 'ncz =  ', ncz
      write(80,*) 'nct = ',nct      
      write(*,*) 'xmin = ', xmin, 'xmax = ', xmax
      write(*,*) 'zmin = ', zmin, 'zmax = ', zmax
      write(*,*) 'ncx = ',ncx, 'ncz =  ', ncz
      write(*,*) 'nct = ',nct      


      ! -- Marching 
      
      if (i_FlowDir.eq.1) then

         ncall1 = ncz           ! EAST
         ncall2 = 1             ! NORTH
         ncall3 = ncz + 1       ! NORTH-EAST
         ncall4 = 1 - ncz       ! NORTH-WEST

         if (i_periodicOBs(1).eq.1) then

            ncall11 = ncz - nct           ! EAST
            ncall12 = ncz + 1 - nct       ! NORTH-EAST
            ncall13 = 1 - ncz + nct       ! NORTH-WEST

         endif

      else if (i_FlowDir.eq.2) then

         ncall1 = 1             ! EAST
         ncall2 = ncx           ! NORTH
         ncall3 = ncx + 1       ! NORTH-EAST
         ncall4 = ncx - 1       ! NORTH-WEST

      endif

      write(80,*) ncall1, ncall2, ncall3, ncall4
      write(*,*) ncall1, ncall2, ncall3, ncall4
      if (i_periodicOBs(1).eq.1) then
         write(80,*) ncall11, ncall12, ncall13
         write(*,*) ncall11, ncall12, ncall13
      endif

      ! -- Miscellaneous Constants
      rrk2h2 = ikappa * ikappa * hsml * hsml
      write(80,*) ''
      write(*,*) ''    
      if (i_kernel.eq.1) then
         write(80,*) ' Cubic Spline Kernel'      
         write(*,*) ' Cubic Spline Kernel'      
      elseif (i_kernel.eq.2) then
         write(80,*) ' Quintic Spline Kernel'      
         write(*,*) ' Quintic Spline Kernel'      
      endif
      write(80,*) 'ikappa = ',ikappa,'rrk2h2 = ',rrk2h2
      write(*,*) 'ikappa = ',ikappa,'rrk2h2 = ',rrk2h2

      end subroutine getdata
