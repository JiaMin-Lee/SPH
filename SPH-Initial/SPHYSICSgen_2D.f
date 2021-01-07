        program SPHYSICSgen_2D

        include 'common.gen2D'

        print*,'Making input file'
        print*


        write(*,*) 'Choose Starting options:   Start new RUN = 0 '
        write(*,*) '                       : Restart old RUN = 1 '
!        write(*,*) '     with CheckPointing:   Start new RUN = 2 '
!        write(*,*) '                       : Restart old RUN = 3 '
        read(*,*) i_restartRun
        write(*,*) i_restartRun

        write(*,*) 'Flow Direction 1.Right 2.Upwards'
        read(*,*) i_FlowDir
        write(*,*) i_FlowDir
        print*

        write(*,*) 'Periodic Boundary Condition? 1.Yes, 2.No.'
        read(*,*) iperiodic
        if (iperiodic.eq.1.and.i_FlowDir.eq.1) then
            write(*,*) 'X-Periodic Boundary at Left & Right Boundaries'
            i_periodicOBs(1) = 1
            i_periodicOBs(2) = 0
        else if (iperiodic.eq.1.and.i_FlowDir.eq.2) then
            write(*,*) 'Z-Periodic Boundary at Left & Right Boundaries'
            i_periodicOBs(1) = 0
            i_periodicOBs(2) = 1
        else if (iperiodic.eq.0) then
            i_periodicOBs(1) = 0
            i_periodicOBs(2) = 0
        endif
        write(*,*) i_periodicOBs(1), i_periodicOBs(2)
        print*

        write(*,*) 'Velocity Inlet/Out Boundary Condition'
        read(*,*) iVelocity
        write(*,*) iVelocity
        print*
        
        if (iVelocity.eq.1) then

           write(*,*) 'VXInlet, VZInlet'
           read(*,*) VXInlet, VZInlet
           write(*,*) VXInlet, VZInlet
           print*

        else if (iVelocity.eq.0) then

           VXInlet = 0.
           VZInlet = 0.
           write(*,*) VXInlet, VZInlet
           print*

        endif 

        write(*,*) 'Body Force grx, grz'
        read(*,*) grx, grz
        write(*,*) grx, grz
        print*

        write(*,*) 'Kinetic Viscosity, Laminar'
        read(*,*) viscos_val
        write(*,*) viscos_val
        print*

        write(*,*) 'Equation of State, 1.Morris, 2.Tait'
        read(*,*) i_EoS
        write(*,*) i_EoS
        print*

        write(*,*) 'Density'
        read(*,*) rho0
        write(*,*) rho0
        print*

        write(*,*) 'Speed of sound'
        read(*,*) cs0
        write(*,*) cs0
        print*

        write(*,*) 'Initial pressure'
        read(*,*) p0
        write(*,*) p0
        print*

        write(*,*) 'Gamma'
        read(*,*) igamma
        write(*,*) igamma
        print*

        write(*,*) 'i_kernel: (1)cubic spline, (2)quintic spline '
        read(*,*) i_kernel
        write(*,*) i_kernel
        print*

        if (i_kernel.eq.1)then
           ikappa = 2
        else if (i_kernel.eq.2) then       
           ikappa = 3
        endif

        write(*,*) 'Smothing length, hsml = a0*dx, a0 = ?'
        read(*,*) a0
        write(*,*) a0
        print*

        write(*,*) 'Kernel Gradient Correction (KGC)? '
        read(*,*) i_kernelcorrection
        write(*,*) i_kernelcorrection
        print*

        write(*,*) 'Density Filtering? (0)None (1)MLS'
        read(*,*) i_densityFilter
        write(*,*) i_densityFilter
        print*

        if (i_densityFilter.eq.1) then
           write(*,*) 'Time Step Interval to apply Density Filtering '
           read(*,*) ndt_FilterPerform
           write(*,*) ndt_FilterPerform
           print*
        else if (i_densityFilter.eq.0) then
           ndt_FilterPerform = 0
        endif
            
        call specgeom(dx, dz)

        if (dx.eq.dz) then
            hsml = a0*dx
        endif

        write(*,*) 'hsml = ',hsml
        print*

        write(*,*) 'dt'
        read(*,*) dt
        write(*,*) dt 
        print*

        write(*,*) 'tmax, tout'
        read(*,*) tmax, tout
        write(*,*) tmax, tout
        print*

!!!!!   WRITING IPART
        open(13,file='IPART')
        print*, 'Writing IPART'
        do ii = 1,np

           iflag(ii) = 1
           write(13,122)xp(ii),zp(ii),up(ii),wp(ii),rhop(ii),pp(ii),
     +                  pm(ii),itype(ii),iflag(ii)
        enddo

122     format(7e16.8, I5, I5)
        close(13)

        print*, 'Closed IPART'        

        open(11, file='INDAT')

        write(11,*) i_restartRun
        write(11,*) i_FlowDir
        write(11,*) iVelocity
        write(11,*) i_kernel
        write(11,*) i_kernelcorrection
        write(11,*) i_densityFilter
        write(11,*) ndt_FilterPerform
        write(11,*) i_EoS
        write(11,*) iBC
        write(11,*) coef
        write(11,*) p0
        write(11,*) cs0
        write(11,*) rho0
        write(11,*) igamma
        write(11,*) ikappa
        write(11,*) hsml
        write(11,*) grx
        write(11,*) grz
        write(11,*) i_periodicOBs(1)
        write(11,*) i_periodicOBs(2)
        write(11,*) VXInlet
        write(11,*) VZInlet
        write(11,188) dt
        write(11,188) tmax
        write(11,188) tout
        write(11,*) viscos_val
        write(11,*) vlx
        write(11,*) vly
        write(11,188) dx
        write(11,188) dz
        write(11,*) nbf
        write(11,*) nbfm
        write(11,*) nvmb
        write(11,*) nb
        write(11,*) np        

188     format(es13.7)

        close(11)

        open(17, file = 'matlabin')

        write(17,*) nbf
        write(17,*) nbfm
        write(17,*) nvmb
        write(17,*) nb
        write(17,*) np

        close(17)

        xb_min = 0.0            !x-boundary min
        xb_max = vlx
        zb_min = 0.0
        zb_max = vlz

        call position_check(dx,dz)

        call tocompile_gfortran 

        end program SPHYSICSgen_2D

!!!!!   BUILD GEOMETRY
        
        subroutine specgeom(dx, dz)

        include 'common.gen2D'

        pi=4.*atan(1.)
        nn = 0

        write(*,*) 'Computational Domain vlx, vlz'
        read(*,*) vlx, vlz
        write(*,*) vlx, vlz
        print*

        write(*,*) 'dx, dz'
        read(*,*)  dx, dz
        write(*,*) dx, dz
        print*

        write(*,*) 'Boundary Condition'
        write(*,*) '1. Monaghan-Repulsive from Cretan Beach (1999)'
        write(*,*) '2. Dalrymple'
        write(*,*) '3. Monaghan-Lennard Jones (1994)'
        write(*,*) '4. Liu CD-SBT (2011)'
        read(*,*) iBC
        write(*,*) iBC
        print*

        if (iBC.eq.1) then
           read(*,*) coef
           print*, 'coef ', coef
           print*
        else
           coef = 0.0
           print*,'coef ', coef
           print*
        endif

        ! Fixed boundary (1 --> nbf)
        write(*,*) 'Generate a Fixed Wall Boundary (L=1, R=2, B=3, T=4)'
        read(*,*)  ifixedwall
        write(*,*) ifixedwall
        print*

        do while (ifixedwall.ge.1)

          vlxstart = 0.0
          vlzstart = 0.0
          vlxend = 0.0
          vlzend = 0.0

           if (ifixedwall.eq.1.or.ifixedwall.eq.2) then

               write(*,*) 'Starting X-Coord'
               read(*,*) vlxstart
               write(*,*) vlxstart
               print*

               write(*,*) 'Starting Z-Coord, Ending Z-Coord'
               read(*,*) vlzstart, vlzend
               write(*,*) vlzstart, vlzend
               print*

               write(*,*) ' Slope (deg) of the inclined plane (beta)?? '
               read(*,*)  beta_deg
               write(*,*) 'beta_deg ',beta_deg
               beta = beta_deg*pi/180.0
               print*
  
               write(*,*) ' How many layers of boundary particles ?? '
               read(*,*)  nlayers
               write(*,*) 'No. of boundary layers  ',nlayers
               print*
           
               if (ifixedwall.eq.1) then

                  call boundaryleft(nn, vlxstart, vlzstart, vlzend,
     +                                0., 0., dx, dz, beta, nlayers)

               else

                  call boundaryright(nn, vlxstart, vlzstart, vlzend,
     +                                0., 0., dx, dz, beta, nlayers)
               endif


           elseif (ifixedwall.eq.3.or.ifixedwall.eq.4) then

               write(*,*) 'Starting Z-Coord'
               read(*,*) vlzstart
               write(*,*) vlzstart
               print*

               write(*,*) 'Starting X-Coord, Ending X-Coord'
               read(*,*) vlxstart, vlxend
               write(*,*) vlxstart, vlxend
               print*

               write(*,*) ' Slope (deg) of the inclined plane (beta)?? '
               read(*,*)  beta_deg
               write(*,*) 'beta_deg ',beta_deg
               beta = beta_deg*pi/180.0
               print*

               write(*,*) ' How many layers of boundary particles ?? '
               read(*,*)  nlayers
               write(*,*) 'No. of boundary layers  ',nlayers
               print*

               if (ifixedwall.eq.3) then

                  call boundarybottom(nn, vlxstart, vlxstart, vlxend,
     +                                0., 0., dx, dz, beta, nlayers)

               else

                  call boundarytop(nn, vlzstart, vlxstart, vlxend,
     +                                0., 0., dx, dz, beta, nlayers)
               endif
           endif

          write(*,*) 'Generate a Fixed Wall Boudary (L=1, R=2, B=3,T=4)'
          read(*,*)  ifixedwall
          write(*,*) ifixedwall
          print*

        enddo

        nbf = nn        

        ! Moving boundary (example gate here) (nbf+1 --> nbfm)

        nbfm = nn       

        ! Floating bodies here (nbfm+1 --> nvmb)

        num_FB = 0
        open(18, file = 'FloatingBodies.txt')
        write(*,*) ' Add Moving Particle in Fluid Domain? '
        write(*,*) ' 1. Square, 2.Circle'
        read(*,*) iopt_FloatingBodies
        write(*,*) iopt_FloatingBodies
        write(18,*) iopt_FloatingBodies
        write(18,*) nbfm       
 
        do while (iopt_FloatingBodies.ge.1)
 
           if (iopt_FloatingBodies.eq.1) then
              call FloatingSquareParticle(nn, dx, dz)
           elseif (iopt_FloatingBodies.eq.2) then
              call FloatingCircleParticle(nn, dx, dz)
           endif             

           write(*,*) 'Add another Moving Particle in Fluid Domain'
           write(*,*) ' 1. Square, 2.Circle'
           read(*,*) iopt_FloatingBodies
           write(*,*) iopt_FloatingBodies

        enddo 

        close(18)

        nvmb = nn

        ! Virtual Moving boundary here (nvmb+1 --> nb)
        write(*,*) ' Generate a Virtual Moving Boundary'
        write(*,*) ' (L=1, R=2, B=3, T=4)'
        read(*,*)  iVMBwall
        write(*,*) iVMBwall
        print*

        open(12,file = 'VirtualMovingBoundary')
        write(12,*) iVMBwall

        do while (iVMBwall.ge.1)

           call VMBoundary(nn, dx, dz, iVMBwall)

           write(*,*) 'Generate another VMB (L=1, R=2, B=3,T=4)'
           read(*,*)  iVMBwall
           write(*,*) iVMBwall
           print*

        enddo

        close(12)
        
        nb = nn

        ! Fluid here (nb+1 --> np)
        write(*,*) 'Generate fluid domain '
        read(*,*) i_fluid
        write(*,*) i_fluid
        print*

        open(14,file = 'BufferZone')
        num_Buff = 0

        do while (i_fluid.eq.1)

          call fill_part(nn, XXmin, XXmax, ZZmin, ZZmax, dx, dz, beta)

          write(*,*) 'Generate another fluid domain '
          read(*,*) i_fluid
          write(*,*) i_fluid
          print*
      
        enddo

        close(14)

        np = nn

        end subroutine specgeom


!------------------------------------------------------------------------------------------!
!!!!!   FIXED WALL BOUNDARY
        subroutine boundaryleft(nn, vlxstart, vlzstart, vlzend,
     +                          uvel, wvel, dx, dz, beta, nlayers)
        
        include 'common.gen2D'

        print*, 'Left wall start at', nn + 1

        N_start = nint(vlxstart/dx)
        L_start = nint(vlzstart/dz)
        L_end = nint(vlzend/dz)

        if (iBC.eq.1) then

            do nly=0,(nlayers-1)

               do kk = L_start,L_end

                  i_type=1
                  nn = nn+1
                  call pos_veloc(nn, (N_start-nly)*dx, kk*dz, uvel,wvel,
     +                           i_type)                  
                  call property(nn, dx, dz)

                  if (kk.eq.L_start) then
                      i_minus1 = nn + 1
                      i_plus1 = nn
                  elseif(kk.eq.L_end) then
                      i_minus1 = nn
                      i_plus1 = nn -1
                  else
                      i_minus1 = nn + 1
                      i_plus1 = nn -1
                  end if

                  iBP_Pointer_Info(nn,3) = i_minus1
                  iBP_Pointer_Info(nn,4) = i_plus1
               enddo
            enddo

        elseif (iBC.eq.3.or.iBC.eq.4) then

            if (iBC.eq.4) nlayers = ikappa + 1

            do nly=0,(nlayers-1)

               do kk = L_start,L_end

                  i_type=1
                  if (iBC.eq.4.and.nly.gt.0) i_type = -1

                  nn = nn+1
                  call pos_veloc(nn, (N_start-nly)*dx, kk*dz, uvel,wvel,
     +                           i_type)
                  call property(nn, dx, dz)

                enddo

            enddo

        endif

        end subroutine boundaryleft

        subroutine boundaryright(nn, vlxstart, vlzstart, vlzend,
     +                          uvel, wvel, dx, dz, beta, nlayers)
        
        include 'common.gen2D'

        print*, 'Right wall start at', nn + 1

        N_start = nint(vlxstart/dx)
        L_start = nint(vlzstart/dz)
        L_end = nint(vlzend/dz)

        if (iBC.eq.1) then

            do nly=0,(nlayers-1)

               do kk = L_start,L_end

                  i_type=1
                  nn = nn+1
                  call pos_veloc(nn, (N_start+nly)*dx, kk*dz, uvel,wvel,
     +                           i_type)                  
                  call property(nn, dx, dz)

                  if (kk.eq.L_start) then
                      i_minus1 = nn
                      i_plus1 = nn + 1
                  elseif(kk.eq.L_end) then
                      i_minus1 = nn - 1
                      i_plus1 = nn
                  else
                      i_minus1 = nn - 1
                      i_plus1 = nn + 1
                  end if

                  iBP_Pointer_Info(nn,3) = i_minus1
                  iBP_Pointer_Info(nn,4) = i_plus1
               enddo
            enddo

        elseif (iBC.eq.3.or.iBC.eq.4) then

            if (iBC.eq.4) nlayers = ikappa + 1

            do nly=0,(nlayers-1)

               do kk = L_start,L_end

                  i_type=1
                  if (iBC.eq.4.and.nly.gt.0) i_type = -1

                  nn = nn+1
                  call pos_veloc(nn, (N_start+nly)*dx, kk*dz, uvel,wvel,
     +                           i_type)
                  call property(nn, dx, dz)

                enddo

            enddo

        endif

        end subroutine boundaryright

        subroutine boundarybottom(nn, vlzstart, vlxstart, vlxend,
     +                          uvel, wvel, dx, dz, beta, nlayers)
        
        include 'common.gen2D'

        print*, 'Bottom wall start at', nn + 1

        L_start = nint(vlzstart/dx)
        N_start = nint(vlxstart/dz)
        N_end = nint(vlxend/dz)

        if (iBC.eq.1) then

            do nly=0,(nlayers-1)

               do ii = N_start,N_end

                  i_type=1
                  nn = nn+1
                  call pos_veloc(nn, ii*dx, (L_start-nly)*dz, uvel,wvel,
     +                           i_type)                  
                  call property(nn, dx, dz)

                  if (ii.eq.N_start) then
                      i_minus1 = nn
                      i_plus1 = nn + 1
                  elseif(ii.eq.N_end) then
                      i_minus1 = nn - 1
                      i_plus1 = nn
                  else
                      i_minus1 = nn - 1
                      i_plus1 = nn + 1
                  end if

                  iBP_Pointer_Info(nn,3) = i_minus1
                  iBP_Pointer_Info(nn,4) = i_plus1
               enddo
            enddo

        elseif (iBC.eq.3.or.iBC.eq.4) then

            if (iBC.eq.4) nlayers = ikappa + 1

            do nly=0,(nlayers-1)
               
               do ii = N_start,N_end

                  i_type = 1
                  if (iBC.eq.4.and.nly.gt.0) i_type = -1

                  nn = nn+1
                  call pos_veloc(nn, ii*dx, (L_start-nly)*dz, uvel,wvel,
     +                           i_type)
                  call property(nn, dx, dz)

                enddo

            enddo

        endif

        end subroutine boundarybottom
 
        subroutine boundarytop(nn, vlzstart, vlxstart, vlxend,
     +                          uvel, wvel, dx, dz, beta, nlayers)
        
        include 'common.gen2D'

        print*, 'Top wall start at', nn + 1

        L_start = nint(vlzstart/dx)
        N_start = nint(vlxstart/dz)
        N_end = nint(vlxend/dz)

        if (iBC.eq.1) then

            do nly=0,(nlayers-1)

               do ii = N_start,N_end

                  i_type=1
                  nn = nn+1
                  call pos_veloc(nn, ii*dx, (L_start+nly)*dz,uvel, wvel,
     +                           i_type)                  
                  call property(nn, dx, dz)

                  if (ii.eq.N_start) then
                      i_minus1 = nn + 1
                      i_plus1 = nn 
                  elseif(ii.eq.N_end) then
                      i_minus1 = nn
                      i_plus1 = nn - 1
                  else
                      i_minus1 = nn + 1
                      i_plus1 = nn - 1
                  end if

                  iBP_Pointer_Info(nn,3) = i_minus1
                  iBP_Pointer_Info(nn,4) = i_plus1
               enddo
            enddo

        elseif (iBC.eq.3.or.iBC.eq.4) then

            if (iBC.eq.4) nlayers = ikappa + 1

            do nly=0,(nlayers-1)

               do ii = N_start,N_end

                  i_type=1
                  if (iBC.eq.4.and.nly.gt.0) i_type = -1

                  nn = nn+1
                  call pos_veloc(nn, ii*dx, (L_start+nly)*dz, uvel,wvel,
     +                           i_type)
                  call property(nn, dx, dz)

                enddo

            enddo

        endif

        end subroutine boundarytop

!------------------------------------------------------------------------------------------!
!!!!!   MOVING PARTICLES IN FLUID DOMAIN

        subroutine FloatingSquareParticle(nn, dx, dz)

        include 'common.gen2D'

        double precision x1, x2, z1, z2

        pi=4.*atan(1.)

        num_FB = num_FB + 1

        if (num_FB.gt.num_FB_max) then
           print*, 'Number of Moving Particles ', num_FB
           print*, 'More than maximum ',num_FB_max
           stop
        endif

        iopt_FloatingObject_type(num_FB) = iopt_FloatingBodies

        !------------------ READING INPUT INFORMATION OF PARTICLE ------------------!
        write(*,*) 'SQUARE PARTICLE ',iopt_FloatingObject_type(num_FB)
        print*

        write(*,*) 'X and Z Lengths of Square/Rectangle '
        read(*,*) XcylinderDimension(num_FB),ZcylinderDimension(num_FB)
        write(*,*) XcylinderDimension(num_FB),ZcylinderDimension(num_FB)
        print*

        write(*,*) 'Specific Density of Particle, rhoP = SD x rho0'
        read(*,*) FB_specificDensity(num_FB)
        write(*,*) FB_specificDensity(num_FB)
        print*

        write(*,*) 'CG X,Z coordinates of Particle '
        read(*,*) Box_XC(num_FB), Box_ZC(num_FB)
        write(*,*) Box_XC(num_FB), Box_ZC(num_FB)
        print*

        write(*,*) 'Shift of CG for Parallel Axis Theorem'
        write(*,*) 'Shift will make Particle slightly unbalanced'
        read(*,*) xdiff_CofG, zdiff_CofG
        print*

        write(*,*) 'Initial Velocity of Particle'
        read(*,*) bigU(num_FB), bigW(num_FB)
        write(*,*) bigU(num_FB), bigW(num_FB)
        print*

        write(*,*) 'Initial Body Angle (counterclockwise)' 
        write(*,*) 'and Angular Velocity'
        read(*,*) body_Angle(num_FB), bigOmega(num_FB)
        write(*,*) body_Angle(num_FB), bigOmega(num_FB)
        print*

        if (body_Angle(num_FB).gt.180) then
            body_Angle(num_FB) = body_Angle(num_FB) - 360.0
        endif
        body_Angle(num_FB) = body_Angle(num_FB)*pi/180.0

        write(*,*) 'Coefficient of Friction'
        read(*,*) friction_coeff(num_FB) 
        write(*,*) friction_coeff(num_FB)
        print*

        
        cylinderDensity(num_FB) = rho0*FB_specificDensity(num_FB)
        bigMass(num_FB) = cylinderDensity(num_FB) 
     &          *(XcylinderDimension(num_FB)*ZcylinderDimension(num_FB))
        bigInertiaYY(num_FB) = bigMass(num_FB)
     &     *(XcylinderDimension(num_FB)*XcylinderDimension(num_FB) 
     &     + ZcylinderDimension(num_FB)*ZcylinderDimension(num_FB))/12.0
        x_BottomLeft = Box_XC(num_FB) - 0.5*XcylinderDimension(num_FB)
        z_BottomLeft = Box_ZC(num_FB) - 0.5*ZcylinderDimension(num_FB)
        Box_XC(num_FB) = Box_XC(num_FB) + xdiff_CofG
        Box_ZC(num_FB) = Box_ZC(num_FB) + zdiff_CofG
        bigInertiaYY(num_FB) = bigInertiaYY(num_FB) 
     &  + bigMass(num_FB)*(xdiff_CofG*xdiff_CofG+zdiff_CofG*zdiff_CofG)
        
        print*
        print*, '-- Generate Moving Particle --'
        print*, 'Moving Particle Index ',num_FB
        print*, ' INITIAL CONDITION '
        print*, 'Density of Particle ', cylinderDensity(num_FB)
        print*, 'Mass of Particle ', bigMass(num_FB)
        print*, 'Moment of Inertia of Particle ', bigInertiaYY(num_FB)
        print*, 'X-Length of Particle ', XcylinderDimension(num_FB)
        print*, 'Z-Length of Particle ', ZcylinderDimension(num_FB)
        print*, 'X- and Z-CG of Particle ',Box_XC(num_FB),Box_ZC(num_FB)
        print*, 'Velocity ', bigU(num_FB), bigW(num_FB)
        print*, 'Angular Velocity and Body Angle ', bigOmega(num_FB),
     &                                              body_Angle(num_FB)

        nn_initial = nn

        x_cyl_min(num_FB)=Box_XC(num_FB)-0.5*XcylinderDimension(num_FB)
        x_cyl_max(num_FB)=Box_XC(num_FB)+0.5*XcylinderDimension(num_FB)
        z_cyl_min(num_FB)=Box_ZC(num_FB)-0.5*ZcylinderDimension(num_FB)
        z_cyl_max(num_FB)=Box_ZC(num_FB)+0.5*ZcylinderDimension(num_FB)

        print*,'x_cyl_min ',x_cyl_min(num_FB),' x_cyl_max ',
     &          x_cyl_max(num_FB)
        print*,'z_cyl_min ',z_cyl_min(num_FB),' z_cyl_max ',
     &          z_cyl_max(num_FB)

        N_start = nint(x_cyl_min(num_FB)/dx)
        N_end = nint((x_cyl_min(num_FB)+XcylinderDimension(num_FB))/dx)
        
        L_start = nint(z_cyl_min(num_FB)/dz)
        L_end = nint((z_cyl_min(num_FB)+ZcylinderDimension(num_FB))/dz)

        xx = N_start*dx
        zz = L_start*dz 

        if (iBC.eq.3) then

           N_end = N_end - 1
           L_end = L_end - 1

        endif

        print*, 'Bottom Left Corner X- and Z-Coordindates', xx,zz
        print*
        ! --- Left Side (Bottom to Up) --- !
        zz = zz - dz
        ntemp = nn

        do kk = L_start, L_end

           i_type = 1
           zz = zz + dz
           nn = nn + 1

           dx_Box = xx - Box_XC(num_FB)
           dz_Box = zz - Box_ZC(num_FB)
           dx_temp = dx_Box*cos(body_Angle(num_FB))
     &             - dz_Box*sin(body_Angle(num_FB))
           dz_temp = dx_Box*sin(body_Angle(num_FB))
     &             + dz_Box*cos(body_Angle(num_FB))

           x1 = Box_XC(num_FB) + dx_temp
           z1 = Box_ZC(num_FB) + dz_temp
           u1 = bigU(num_FB) - bigOmega(num_FB)*(z1-Box_ZC(num_FB))
           w1 = bigW(num_FB) + bigOmega(num_FB)*(x1-Box_XC(num_FB))

           call pos_veloc(nn, x1, z1, u1, w1, i_type)
           call property(nn,dx,dz)

        enddo

        ! --- Top Side (Left to Right) --- !
        xx = xx - dx 
        if (iBC.eq.3) zz = zz + dz ! Since L_end - 1

        do ii = N_start, N_end

           i_type = 1
           xx = xx + dx
           nn = nn + 1

           dx_Box = xx - Box_XC(num_FB)
           dz_Box = zz - Box_ZC(num_FB)
           dx_temp = dx_Box*cos(body_Angle(num_FB))
     &             - dz_Box*sin(body_Angle(num_FB))
           dz_temp = dx_Box*sin(body_Angle(num_FB))
     &             + dz_Box*cos(body_Angle(num_FB))

           x1 = Box_XC(num_FB) + dx_temp
           z1 = Box_ZC(num_FB) + dz_temp
           u1 = bigU(num_FB) - bigOmega(num_FB)*(z1-Box_ZC(num_FB))
           w1 = bigW(num_FB) + bigOmega(num_FB)*(x1-Box_XC(num_FB))

           call pos_veloc(nn, x1, z1, u1, w1, i_type)           
           call property(nn,dx,dz)

        enddo

        ! --- Right Side (Top to Bottom) --- !
        zz = zz + dz
        if (iBC.eq.3) xx = xx + dx

        do kk = L_start, L_end

           i_type = 1
           zz = zz - dz
           nn = nn + 1

           dx_Box = xx - Box_XC(num_FB)
           dz_Box = zz - Box_ZC(num_FB)
           dx_temp = dx_Box*cos(body_Angle(num_FB))
     &             - dz_Box*sin(body_Angle(num_FB))
           dz_temp = dx_Box*sin(body_Angle(num_FB))
     &             + dz_Box*cos(body_Angle(num_FB))

           x1 = Box_XC(num_FB) + dx_temp
           z1 = Box_ZC(num_FB) + dz_temp
           u1 = bigU(num_FB) - bigOmega(num_FB)*(z1-Box_ZC(num_FB))
           w1 = bigW(num_FB) + bigOmega(num_FB)*(x1-Box_XC(num_FB))

           call pos_veloc(nn, x1, z1, u1, w1, i_type)
           call property(nn,dx,dz)

        enddo        

        ! --- Bottom Side (Right to Left) --- !
        xx = xx + dx
        if (iBC.eq.3) zz = zz - dz

        do ii = N_start, N_end

           i_type = 1
           xx = xx - dx
           nn = nn + 1

           dx_Box = xx - Box_XC(num_FB)
           dz_Box = zz - Box_ZC(num_FB)
           dx_temp = dx_Box*cos(body_Angle(num_FB))
     &             - dz_Box*sin(body_Angle(num_FB))
           dz_temp = dx_Box*sin(body_Angle(num_FB))
     &             + dz_Box*cos(body_Angle(num_FB))

           x1 = Box_XC(num_FB) + dx_temp
           z1 = Box_ZC(num_FB) + dz_temp
           u1 = bigU(num_FB) - bigOmega(num_FB)*(z1-Box_ZC(num_FB))
           w1 = bigW(num_FB) + bigOmega(num_FB)*(x1-Box_XC(num_FB))

           call pos_veloc(nn, x1, z1, u1, w1, i_type)
           call property(nn,dx,dz)

        enddo

        nb_FB(num_FB) = nn

        write(18,*)num_FB
        write(18,*)iopt_FloatingObject_type(num_FB)
        write(18,*)bigMass(num_FB)
        write(18,*)bigInertiaYY(num_FB)
        write(18,*)XcylinderDimension(num_FB),ZcylinderDimension(num_FB)
        write(18,*)cylinderDensity(num_FB)
        write(18,*)FB_specificDensity(num_FB)
        write(18,*)friction_coeff(num_FB)
        write(18,*)Box_XC(num_FB),Box_ZC(num_FB)
        write(18,*)bigU(num_FB),bigW(num_FB),bigOmega(num_FB)
        write(18,*)nb_FB(num_FB)

        return

        end subroutine FloatingSquareParticle




        subroutine FloatingCircleParticle(nn, dx, dz)

        include 'common.gen2D'

        double precision x1, x2, z1, z2
        
        pi=4.*atan(1.)

        num_FB = num_FB + 1

        if (num_FB.gt.num_FB_max) then
           print*, 'Number of Moving Particles ', num_FB
           print*, 'More than maximum ',num_FB_max
           stop
        endif

        iopt_FloatingObject_type(num_FB) = iopt_FloatingBodies

        !------------------ READING INPUT INFORMATION OF PARTICLE ------------------!        
        write(*,*) 'CIRCLE PARTICLE ', iopt_FloatingObject_type(num_FB)
        print*
        
        write(*,*) 'X and Z Length  of Circle (Radius) '
        read(*,*) XcylinderDimension(num_FB), ZcylinderDimension(num_FB)
        write(*,*) XcylinderDimension(num_FB),ZcylinderDimension(num_FB)
        print*

        write(*,*) 'Specific Density of Particle, rhoP = SD x rho0'
        read(*,*) FB_specificDensity(num_FB)
        write(*,*) FB_specificDensity(num_FB)
        print*

        write(*,*) 'CG X,Z coordinates of Particle '
        read(*,*) Box_XC(num_FB), Box_ZC(num_FB)
        write(*,*) Box_XC(num_FB), Box_ZC(num_FB)
        print*

        write(*,*) 'Shift of CG for Parallel Axis Theorem'
        write(*,*) 'Shift will make Particle slightly unbalanced'
        read(*,*) xdiff_CofG, zdiff_CofG
        print*

        write(*,*) 'Initial Velocity of Particle'
        read(*,*) bigU(num_FB), bigW(num_FB)
        write(*,*) bigU(num_FB), bigW(num_FB)
        print*

        write(*,*) 'Initial Body Angle (counterclockwise) about CG' 
        write(*,*) 'and Angular Velocity'
        read(*,*) body_Angle(num_FB), bigOmega(num_FB)
        write(*,*) body_Angle(num_FB), bigOmega(num_FB)
        print*

        if (body_Angle(num_FB).gt.180) then
            body_Angle(num_FB) = body_Angle(num_FB) -360.0
        endif
        body_Angle(num_FB) = body_Angle(num_FB)*pi/180.0

        write(*,*) 'Coefficient of Friction'
        read(*,*) friction_coeff(num_FB) 
        write(*,*) friction_coeff(num_FB)
        print*

     
        cylinderDensity(num_FB) = rho0*FB_specificDensity(num_FB)
        bigMass(num_FB) = cylinderDensity(num_FB)*pi 
     &          *(XcylinderDimension(num_FB)*ZcylinderDimension(num_FB))
        bigInertiaYY(num_FB) = bigMass(num_FB)
     &     *(XcylinderDimension(num_FB)*XcylinderDimension(num_FB) 
     &     + ZcylinderDimension(num_FB)*ZcylinderDimension(num_FB))/4.0
        x_BottomLeft = Box_XC(num_FB) - 0.5*XcylinderDimension(num_FB)
        z_BottomLeft = Box_ZC(num_FB) - 0.5*ZcylinderDimension(num_FB)
        Box_XC(num_FB) = Box_XC(num_FB) + xdiff_CofG
        Box_ZC(num_FB) = Box_ZC(num_FB) + zdiff_CofG
        bigInertiaYY(num_FB) = bigInertiaYY(num_FB) 
     &  + bigMass(num_FB)*(xdiff_CofG*xdiff_CofG+zdiff_CofG*zdiff_CofG)
        
        print*
        print*, '-- Generate Moving Particle --'
        print*, 'Moving Particle Index ',num_FB
        print*, ' INITIAL CONDITION '
        print*, 'Density of Particle ', cylinderDensity(num_FB)
        print*, 'Mass of Particle ', bigMass(num_FB)
        print*, 'Moment of Inertia of Particle ', bigInertiaYY(num_FB)
        print*, 'X-Radius of Particle ', XcylinderDimension(num_FB)
        print*, 'Z-Radius of Particle ', ZcylinderDimension(num_FB)
        print*, 'X- and Z-CG of Particle ',Box_XC(num_FB),Box_ZC(num_FB)
        print*, 'Velocity ', bigU(num_FB), bigW(num_FB)
        print*, 'Angular Velocity and Body Angle ', bigOmega(num_FB),
     &                                              body_Angle(num_FB)
        
        nn_initial = nn

        ! Rough estimation of number of particles
        radiusRC = max(ZcylinderDimension(num_FB),
     &                XcylinderDimension(num_FB))

        N_start = 1
        N_end = nint(2*pi*radiusRC/(0.5*(dx+dz)))
        darc(num_FB) = 360.0/N_end

        xCenter(num_FB) = Box_XC(num_FB)-xdiff_CofG
        zCenter(num_FB) = Box_ZC(num_FB)-zdiff_CofG

        print*, 'Arc length ', 0.5*(dx+dz)
        print*, 'Arc angle ', darc(num_FB)
        print*
        
        do ii = N_start,N_end
        
           tt = (ii-1)*darc(num_FB)*pi/180
           nn = nn + 1
           i_type = 1
           
           xx = xCenter(num_FB)
     &     +  XcylinderDimension(num_FB)*cos(tt)!*cos(body_Angle(num_FB))
!     &     -  ZcylinderDimension(num_FB)*sin(tt)*sin(body_Angle(num_FB))
           zz = zCenter(num_FB)
!     &     + XcylinderDimension(num_FB)*cos(tt)*sin(body_Angle(num_FB))
     &     + ZcylinderDimension(num_FB)*sin(tt)!*cos(body_Angle(num_FB))

           dx_Box = xx - Box_XC(num_FB)
           dz_Box = zz - Box_ZC(num_FB)

           ! Body Angle about CG
           xp_temp = dx_Box*cos(body_Angle(num_FB))
     &             - dz_Box*sin(body_Angle(num_FB))
           zp_temp = dx_Box*sin(body_Angle(num_FB))
     &             + dz_Box*cos(body_Angle(num_FB))
           x1=Box_XC(num_FB) + xp_temp
           z1=Box_ZC(num_FB) + zp_temp


           u1 = bigU(num_FB) - bigOmega(num_FB)*(z1-Box_ZC(num_FB))
           w1 = bigW(num_FB) + bigOmega(num_FB)*(x1-Box_XC(num_FB))

           call pos_veloc(nn, x1, z1, u1, w1, i_type)
           call property(nn,dx,dz)

        enddo

        nb_FB(num_FB) = nn

        write(18,*)num_FB
        write(18,*)iopt_FloatingObject_type(num_FB)
        write(18,*)bigMass(num_FB)
        write(18,*)bigInertiaYY(num_FB)
        write(18,*)XcylinderDimension(num_FB),ZcylinderDimension(num_FB)
        write(18,*)cylinderDensity(num_FB)
        write(18,*)FB_specificDensity(num_FB)
        write(18,*)friction_coeff(num_FB)
        write(18,*)Box_XC(num_FB),Box_ZC(num_FB)
        write(18,*)bigU(num_FB),bigW(num_FB),bigOmega(num_FB)
        write(18,*)nb_FB(num_FB)

        return


        end subroutine FloatingCircleParticle

!------------------------------------------------------------------------------------------!
!!!!!   VIRTUAL MOVING BOUNDARY

        subroutine VMBoundary(nn, dx, dz, iVMBwall)

        include 'common.gen2D'

        write(12,*) iVMBwall
        write(12,*) nn+1        

        vlxstart = 0.
        vlzstart = 0.
        vlxend = 0.
        vlzend = 0.
        
        N_start = 0
        N_end = 0
        L_start = 0       
        L_end = 0

        if (iVMBwall.eq.1.or.iVMBwall.eq.2) then

            write(*,*) 'Starting X-Coord'
            read(*,*) vlxstart
            write(*,*) vlxstart
            print*

            write(*,*) 'Starting Z-Coord, Ending Z-Coord'
            read(*,*) vlzstart, vlzend
            write(*,*) vlzstart, vlzend
            print*

            write(*,*) ' Velocity X, Velocity Z '
            read(*,*)  VXwall, VZwall
            write(*,*) VXwall, VZwall       
            print*
  
            write(*,*) ' How many layers of boundary particles ?? '
            read(*,*)  nlayers
            write(*,*) 'No. of boundary layers  ',nlayers
            print*
             
            N_start = nint(vlxstart/dx)
            L_start = nint(vlzstart/dz)
            L_end = nint(vlzend/dz)
        
            do nly=0,(nlayers-1)

               if (iVMBwall.eq.1) then
                   nly1 = -nly
               else
                   nly1 = nly
               endif

               do kk = L_start,L_end

                  i_type=3
                  nn = nn+1
                  call pos_veloc(nn, (N_start+nly1)*dx, kk*dz,VXwall,
     +                               VZwall,i_type)             
                  call property(nn, dx, dz)

                  if (nly.eq.0) then

                     nlayer1 = nn

                  endif
                
                enddo

             enddo                

        elseif (iVMBwall.eq.3.or.iVMBwall.eq.4) then

            write(*,*) 'Starting Z-Coord'
            read(*,*) vlzstart
            write(*,*) vlzstart
            print*

            write(*,*) 'Starting X-Coord, Ending X-Coord'
            read(*,*) vlzstart, vlzend
            write(*,*) vlzstart, vlzend
            print*

            write(*,*) ' Velocity '
            read(*,*)  VXwall, VZwall
            write(*,*) VXwall, VZwall       
            print*

            write(*,*) ' How many layers of boundary particles ?? '
            read(*,*)  nlayers
            write(*,*) 'No. of boundary layers  ',nlayers
            print*

            L_start = nint(vlzstart/dz)
            N_start = nint(vlxstart/dx)
            N_end = nint(vlxend/dz)
        
            do nly=0,(nlayers-1)

               if (iVMBwall.eq.3) then
                   nly1 = -nly
               else
                   nly1 = nly
               endif

               do ii = N_start,N_end

                  i_type=3
                  nn = nn+1
                  call pos_veloc(nn, ii*dx, (L_start+nly1)*dz,VXwall,
     +                               VZwall,i_type)             
                  call property(nn, dx, dz)

               enddo

            enddo
        endif

        write(12,*) nlayer1
        write(12,*) nn
        write(12,*) nlayers
        write(12,*) VXwall, VZwall

        end subroutine VMBoundary

!------------------------------------------------------------------------------------------!
!!!!!   FILL DOMAIN WITH FLUID

        subroutine fill_part(nn, XXmin, XXmax, ZZmin, ZZmax, dx, dz,
     +                       beta)

        include 'common.gen2D'

        double precision x1, z1

        XXmin = 0.
        XXmax = 0.
        ZZmin = 0.
        ZZmax = 0.

        N_start = 0
        N_end = 0 
        L_start = 0
        L_end = 0

        write(*,*) 'Cube containing fluid'
        write(*,*) 'XXmin, XXmax'
        read(*,*) XXmin, XXmax
        write(*,*) XXmin, XXmax
        print*

        write(*,*) 'ZZmin, ZZmax'
        read(*,*) ZZmin, ZZmax
        write(*,*) ZZmin, ZZmax
        print*

        write(*,*) 'VXwater, VZwater'
        read(*,*) VXwater, VZwater
        write(*,*) VXwater, VZwater
        print*

        N_start = nint(XXmin/dx)
        N_end = nint(XXmax/dx)
        L_start = nint(ZZmin/dz)
        L_end = nint(ZZmax/dz)

!!!     ----- BUFFER ZONE -----     !!!                         

        write(*,*) 'Include Buffer Zone Interface?'
        write(*,*) '1. Inlet 2. Outlet'
        read(*,*) ibuffer
        write(*,*) ibuffer
        print*

        if (ibuffer.eq.2.and.i_periodicOBs(1).eq.1) then
            ibuffer = 0
            write(*,*) 'Buffer Zone at Outlet and X-Periodic cannot'
            write(*,*) 'coexist!!!'
            stop
        endif
 
        if (ibuffer.eq.0) then
           if (num_Buff.eq.0) write(14,*) ibuffer
        else
           if (num_Buff.eq.0) write(14,*) 1
            num_Buff = num_Buff + 1
        endif

        do while (ibuffer.ge.1)

           write(14,*) num_Buff
           write(14,*) ibuffer 

           write(*,*) 'Buffer Zone Particles Velocity'
           read(*,*) VXbuffer, VZbuffer
           write(*,*) VXbuffer, VZbuffer

           if (ibuffer.eq.1) then

              write(14,*) nn + 1

              Nxbuffint = N_start + 2*ikappa

              do ii = N_start, Nxbuffint
                 do kk = L_start, L_end

                    nn = nn + 1
                    i_type = 4

!                    VXbuffer =(-7.7295*kk*dz*kk*dz+0.0116*kk*dz)
!     &                         *1000

                    call pos_veloc(nn, ii*dx, kk*dz,VXbuffer,VZbuffer,
     &                             i_type)
                    call property(nn, dx, dz)                  

                 enddo
              enddo 

              N_start =  Nxbuffint + 1
              Vxbuffint = dx*(Nxbuffint+0.5)

              write(14,*) nn              
              write(14,*) VXbuffer, VZbuffer
              write(14,*) Vxbuffint

           elseif (ibuffer.eq.2) then

              write(14,*) nn + 1

              Nxbuffint = N_end - 2*ikappa
             
              do ii = Nxbuffint,N_end
                 do kk = L_start, L_end

                    nn = nn + 1
                    i_type = 4

!                    VXbuffer =(-7.7295*kk*dz*kk*dz+0.0116*kk*dz)
!     &                         *1000


                    call pos_veloc(nn, ii*dx, kk*dz,VXbuffer,VZbuffer,
     &                             i_type)
                    call property(nn, dx, dz)                  

                 enddo
              enddo
        
              N_end = Nxbuffint - 1
              Vxbuffint = dx*(Nxbuffint-0.5)
               
              write(14,*) nn
              write(14,*) VXbuffer, VZbuffer
              write(14,*) Vxbuffint

           end if

           write(*,*) 'Another Buffer Zone? 1.Inlet, 2.Outlet'
           read(*,*) ibuffer
           write(*,*) ibuffer    

           if (ibuffer.gt.0) num_Buff = num_Buff + 1
           
        enddo
        
!!!!!   FILL DOMAIN WITH FLUID PARTICLES
        XXmin = N_start*dx
        XXmax = N_end*dx
        ZZmin = L_start*dz
        ZZmax = L_end*dz

!        write(*,*) XXmin, XXmax
!        write(*,*) ZZmin, ZZmax

        do ii = N_start, N_end
           do kk = L_start, L_end

              x1 = XXmin + (ii - N_start)*dx
              z1 = ZZmin + (kk - L_start)*dz

              i_createFluid = 1

              if (num_FB.gt.0.and.i_createFluid.eq.1) then
                 i_num_FB = 0
                 i_outsideFB = 1

                 do while (i_num_FB.lt.num_FB.and.i_outsideFB.eq.1)

                    i_num_FB = i_num_FB + 1
                    
                    if (iopt_FloatingObject_type(i_num_FB).eq.1) then

                       clearanceF = 0.99
                        
                       ! Fluid particle distance with respect to FB CG 
                       ! in global coordinate
                       dx_newParticle = x1 - Box_XC(i_num_FB)
                       dz_newParticle = z1 - Box_ZC(i_num_FB)

                       ! Transform distance with respect to FB body
                       ! angle (local coordinate)
                       x1_Tr = dx_newParticle*cos(body_Angle(i_num_FB))
     &                       - dz_newParticle*sin(body_Angle(i_num_FB))
                       z1_Tr = dx_newParticle*sin(body_Angle(i_num_FB))
     &                       + dz_newParticle*cos(body_Angle(i_num_FB))

                       ! Local bounded perimenter of rectangle from CG
                       xtmin = - 0.5*XcylinderDimension(i_num_FB)
     &                         - clearanceF*dx
                       xtmax =   0.5*XcylinderDimension(i_num_FB)
     &                         + clearanceF*dx  
                       ztmin = - 0.5*ZcylinderDimension(i_num_FB)
     &                         - clearanceF*dz
                       ztmax =   0.5*ZcylinderDimension(i_num_FB)
     &                         + clearanceF*dz

                       if (x1_Tr.gt.xtmin.and.x1_Tr.lt.xtmax) then
                          if (z1_Tr.gt.ztmin.and.z1_Tr.lt.ztmax) then
                             i_outsideFB = 0
                          else
                             i_outsideFB = 1 
                          endif
                       else
                          i_outsideFB = 1
                       endif

                    elseif(iopt_FloatingObject_type(i_num_FB).eq.2) then
                        
                       clearanceF = 0.99
                        
                       ! Fluid particle distance with respect to FB CG 
                       ! in global coordinate
                       dx_newParticle = x1 - Box_XC(i_num_FB)
                       dz_newParticle = z1 - Box_ZC(i_num_FB)

                       ! Transform distance with respect to FB body
                       ! angle (local coordinate)
                       x1_Tr = dx_newParticle*cos(body_Angle(i_num_FB))
     &                       - dz_newParticle*sin(body_Angle(i_num_FB))
                       z1_Tr = dx_newParticle*sin(body_Angle(i_num_FB))
     &                       + dz_newParticle*cos(body_Angle(i_num_FB))

                       ! Local bounded perimenter of rectangle from CG
                       if (x1_Tr.eq.0.) then
                          anglett = pi/4.
                       else 
                          anglett = atan(z1_Tr/x1_Tr)
                       endif

                       xtmin = (XcylinderDimension(i_num_FB)
     &                       +  clearanceF*dx)   
     &                       * cos(anglett)
                             
                       ztmin = (ZcylinderDimension(i_num_FB)
     &                       +  clearanceF*dz)
     &                       *  sin(anglett)

                       dist_Tr = z1_Tr*z1_Tr+x1_Tr*x1_Tr
                       distmin = xtmin*xtmin + ztmin*ztmin

                       if (dist_Tr.lt.distmin) then
                          i_outsideFB = 0
                       else
                          i_outsideFB = 1
                       endif

                    endif

                 enddo !do while(i_num_FB.lt.num_FB.and.i_outsideFB.eq.1)

                 if (i_outsideFB.eq.1) then
                    i_createFluid = 1
                 else
                    i_createFluid = 0
                 endif

              endif ! end if(num_FB.gt.0.and.i_createFluid.eq.1)

              if (i_createFluid.eq.1) then

                 nn=nn+1
                 i_type = 2

!                 if (z1.gt.0..and.z1.le.0.2e-3) then
!                    VXwater = 29./2.*z1
!                 elseif(z1.lt.1.5e-3.and.z1.ge.1.3e-3) then
!                    VXwater = 29./2.*abs(z1-1.5e-3)
!                 else
!                    VXwater = 2.9e-3
!                 endif

                 VXwater=(-7.7295*z1*z1+0.0116*z1)*1000

                 call pos_veloc(nn,x1,z1,VXwater,VZwater,i_type)
                 call property(nn, dx, dz)

              endif

           enddo
        enddo

        end subroutine fill_part


!------------------------------------------------------------------------------------------!
!!!!!   POSITION AND PROPERTIES OF PARTICLES

        subroutine pos_veloc(nn, xpos, zpos, uvel, wvel, i_type)
        include 'common.gen2D'
         
        double precision xpos, zpos
 
        xp(nn) = xpos
        zp(nn) = zpos
        up(nn) = uvel
        wp(nn) = wvel
        itype(nn) = i_type
       
        end subroutine pos_veloc


        subroutine property(nn, dx, dz)
        include 'common.gen2D'

        if (itype(nn).gt.0.or.itype(nn).eq.-1) then

           pp(nn)=p0
           if(itype(nn).eq.1) pp(nn) = 0.
           rhop(nn) = rho0
           pm(nn) = rhop(nn)*dx*dz

        endif

        end subroutine property


!------------------------------------------------------------------------------------------!
!!!!!   CHECK PARTICLES ARE NOT OVERLAPPING

        subroutine position_check(dx, dz)
        include 'common.gen2D'

        dx_thresh = 0.10*dx     !! tolerance
        dz_thresh = 0.10*dz

        !- Generate checking mesh -

        one_over_kh = 1./(ikappa*hsml)
       
        xGC_min = 0.0 - 2*ikappa*hsml
        xGC_max = vlx
        zGC_min = 0.0 - 2*ikappa*hsml
        zGC_max = vlz

        ncx = int((xGC_max - xGC_min) * one_over_kh) + 1
        ncz = int((zGC_max - zGC_min) * one_over_kh) + 1

        nsheet = ncx * ncz
        nct = nsheet            !! total number of cells
        nc = 0
       
        if (nct.ge.nct_max) then

           print*
           print*, 'ERROR in position_check'
           print*, 'No. of cells exceed nct_max in common.gen2D'
           print*, 'nct = ',nct
           print*, 'nct_max = ',nct_max
           stop

        endif

        print*
        print*, '--- Checking particles initial positions ---'
        print*, 'dx ', dx, 'dz', dz
        print*, 'dx_thresh ', dx_thresh, 'dz_thresh', dz_thresh
        print*, 'vlx', vlx, 'vlz', vlz
        print*, 'ncx ', ncx, 'ncz', ncz

        call divide_dr_check(1, np)
        call ac_dr_check(dx_thresh,dz_thresh) 

        end subroutine position_check


!!!!!   BIN PARTICLES TO CHECKING MESH

        subroutine divide_dr_check(n_start,n_end)
        include 'common.gen2D'

        do kk = n_start,n_end
        
           dxp = xp(kk) - xGC_min
           dzp = zp(kk) - zGC_min

           icell = int( dxp * one_over_kh ) + 1
           kcell = int( dzp * one_over_kh ) + 1

           ! -- Mesh number (Left to right, bottom to top)
           ii = icell + (kcell - 1) * ncx

           if (ii.lt.1) then

              print*
              print*, 'ii < 1 '
              print*, 'particle no. ', kk
              stop

           endif

           !!! No. of particle in cell ii
           nc(ii) = nc(ii) + 1
         
           !!! Bin particle kk (index) to cell ii and nc (index in cell
           !!! ii array)

           ibox(ii, nc(ii)) = kk

        enddo

        end subroutine divide_dr_check


        subroutine ac_dr_check(dx_thresh,dz_thresh)
        include 'common.gen2D'

        ncx_start = 1
        ncx_finish = ncx
        ncz_start = 1
        ncz_finish = ncz
        
        do lz = ncz_start,ncz_finish

           do lx = ncx_start,ncx_finish

              j1 = lx + (lz - 1)*ncx

                if (nc(j1).gt.0) then      ! mesh contains particles

                   lxR = lx + 1         ! Right of current mesh
                   lzT = lz + 1         ! Top of cirrent mesh


                   ! -- EAST MESH
                   if (lxR.le.ncx) then

                      call celij(j1, j1+1, dx_thresh, dz_thresh)

                   endif
                
                   ! -- NORTH MESH
                   if (lzT.le.ncz) then

                      call celij(j1, j1+ncx, dx_thresh, dz_thresh)

                      ! -- NORTH-EAST MESH

                      if (lxR.le.ncx) then

                         call celij(j1, j1+ncx+1, dx_thresh, dz_thresh)

                      endif

                   endif

                  if (lxR.le.ncx.and.lzT.gt.2) then

                     call celij(j1, j1-ncx+1, dx_thresh, dz_thresh)

                  endif

                end if

           enddo

        enddo 

        do jllx = ncx_start, ncx_finish
           do jllz = ncz_start, ncz_finish

              j11 = jllx + (jllz-1)*ncx

              if (nc(j11).gt.0) then
                 call self(j11, dx_thresh, dz_thresh)
              endif
           enddo
        enddo 

        end subroutine ac_dr_check

        subroutine celij(j1, j2, dx_thresh, dz_thresh)
        include 'common.gen2D'
        
        do ii = 1, nc(j1)

           ip = ibox(j1, ii)

           do jj = 1, nc(j2)                   

              jp = ibox(j2, jj)

              drx = xp(ip) - xp(jp)
              drz = zp(ip) - zp(jp)

              if (abs(drx).lt.dx_thresh.and.abs(drz).lt.dz_thresh) then

                 print*
                 write(*,*) ' WARNING: The following  particles are'
                 write(*,*) ' closer than dx_thresh and dz_thresh ',i,j
                 
                 if(abs(drx).lt.0.1*dx_thresh.and.
     +              abs(drz).lt.0.1*dz_thresh) then
                 if ((iBC.eq.1.and.(itype(ip).eq.2.and.itype(jp).eq.2))
     +              .or.iBC.eq.2.or.iBC.eq.3) then

                    print*, 'Particle are too close ', ip, jp
                    stop

                 endif
                 endif    

              endif
           
           enddo

        enddo

        end subroutine celij

        subroutine self(j1, dx_thresh, dz_thresh)
        include 'common.gen2D'

        jj_start = 1
        do ii = 1,nc(j1)

           ip = ibox(j1,ii)

           jj_start = jj_start + 1

           do jj = jj_start,nc(j1)
        
              jp = ibox(j1,jj)

              drx = xp(ip) - xp(jp)
              drz = zp(ip) - zp(jp)

              if (abs(drx).lt.dx_thresh.and.abs(drz).lt.dz_thresh) then

                 print*
                 write(*,*) ' WARNING: The following  particles are'
                 write(*,*) ' closer than dx_thresh and dz_thresh ',i,j

                 if(abs(drx).lt.0.1*dx_thresh.and.
     +              abs(drz).lt.0.1*dz_thresh) then
                 if ((iBC.eq.1.and.(itype(ip).eq.2.and.itype(jp).eq.2))
     +              .or.iBC.eq.2.or.iBC.eq.3) then

                    print*, 'Particle are too close ', ip, jp
                    stop

                 endif
                 endif

              endif

           enddo
           
        enddo        

        end subroutine self

        subroutine tocompile_gfortran
        include 'common.gen2D'
        CHARACTER(LEN=10) :: FMT,FMT1
        CHARACTER(LEN=1)  :: TAB,LS

        TAB=CHAR(9)
        LS =CHAR(92)

        FMT="(A)"
        FMT1="(2A)"
        open(22,file='SPHYSICS.mak')
        write(22,FMT) 'FC=gfortran'
        write(22,FMT) 'OPTIONS= -O3'
     &   //' -Wall'
     &   //' -Wextra'
     &   //' -Wimplicit-interface'
     &   //' -fPIC'
     &   //' -fmax-errors=1'
     &   //' -g'
     &   //' -pedantic'
     &   //' -fcheck=all'
     &   //' -fbacktrace'
c     !- Uncomment as required -
c    &   //' -ftree-vectorize'
c    &   //' -ffast-math -funroll-loops'
c    &   //' -fbounds-check'
c    &   //' -frange-check'
c    &   //' -Wunderflow'
c    &   //' -Wuninitialized'
c    &   //' -ffpe-trap=invalid,zero,overflow'
c    &   //' -fstack-check'
c    &   //' -fstack-protector'
c    &   //' -ftrapv'
c    &   //' -ftracer'
c    &   //' -g'
c    &   //' -fimplicit-none'
c     &   //' -fno-automatic'
      write(22,FMT) 'srcdir=.'
      write(22,FMT) 'idir=../../execs'
      write(22,FMT) 'bakdir=../../execs.bak'
      write(22,FMT) 'objects=energy_2D.o recover_list_2D.o ini_divide_2D
     &.o \'
      write(22,FMT1)TAB,'keep_list_2D.o SPHYSICS_2D.o  \'
      write(22,FMT1)TAB,'getdata_2D.o check_limits.o \'
      write(22,FMT1)TAB,'step.o divide.o \'
      write(22,FMT1)TAB,'ac_main.o poute.o \'
      write(22,FMT1)TAB,'ac.o \' !self.o celij.o \'
      write(22,FMT1)TAB,'kernel.o kernel_correction.o \'
      write(22,FMT1)TAB,'gradients_calc_basic_2D.o viscosity.o \'
      write(22,FMT1)TAB,'LennardJones.o correct_2D.o \'
      write(22,FMT1)TAB,'equation_of_state.o \'
      write(22,FMT1)TAB,'virtualmovingboundary.o \'
      write(22,FMT1)TAB,'bufferzone.o \'
      write(22,FMT1)TAB,'pre_celij_KGC_2D.o pre_self_KGC_2D.o \'
      write(22,FMT1)TAB,'periodicityCorrection.o \'

      if(iBC.eq.1.or.iBC.eq.3)then
        write(22,FMT1)TAB,'self.o celij.o \'
      elseif(iBC.eq.4) then
        write(22,FMT1)TAB,'self_CDSBT.o celij_CDSBT.o \'
        write(22,FMT1)TAB,'viscosityCDSBT.o \'
        write(22,FMT1)TAB,'CDSBT.o \'
      endif

      if(i_densityFilter.eq.1) then
        write(22,FMT1)TAB,'densityFilter_MLS_2D.o \'
        write(22,FMT1)TAB,'ac_MLS_2D.o \'
        write(22,FMT1)TAB,'LU_decomposition_2D.o \'
        write(22,FMT1)TAB,'pre_celij_MLS_2D.o \'
        write(22,FMT1)TAB,'pre_self_MLS_2D.o \'
      elseif(i_densityFilter.eq.0) then
        write(22,FMT1)TAB,'densityFilter_NONE_2D.o \'
      endif 

      write(22,FMT)'#'
      write(22,FMT)'SPHYSICS_2D: $(objects)'
      write(22,FMT1)TAB,'$(FC) $(OPTIONS) -o SPHYSICS_2D $(objects)'
      write(22,FMT)'#'
      write(22,FMT1)TAB,'if [ -d $(bakdir) ]; then \'
      write(22,FMT1)TAB,'echo "execs.bak Directory Exists"; else \'
      write(22,FMT1)TAB,'mkdir $(bakdir); \'
      write(22,FMT1)TAB,'echo "execs.bak Directory Created"; \'
      write(22,FMT1)TAB,'fi'
      write(22,FMT)'#'
      write(22,FMT1)TAB,'if [ -d $(idir) ]; then \'
      write(22,FMT1)TAB,'echo "execs Directory Exists"; else \'
      write(22,FMT1)TAB,'mkdir $(idir); \'
      write(22,FMT1)TAB,'echo "execs Directory Created"; \'
      write(22,FMT1)TAB,'fi'
      write(22,FMT)'#'
      write(22,FMT1)TAB,'-if [ -f $(idir)/SPHYSICS_2D ]; then \'
      write(22,FMT1)TAB,'mv -f $(idir)/SPHYSICS_2D $(bakdir)/; \'
      write(22,FMT1)TAB,'echo Old SPHYSICS_2D moved to execs.bak'
     &//' from execs; \'
      write(22,FMT1)TAB,'fi'
      write(22,FMT)'#'
      write(22,FMT1)TAB,'mv SPHYSICS_2D $(idir)'
      write(22,FMT1)TAB,'echo New SPHYSICS_2D moved to execs'
      write(22,FMT)'#'
      write(22,FMT)'clean:'
      write(22,FMT1)TAB,'rm *.o'
      write(22,FMT1)TAB,'rm *~'
      write(22,FMT)'#'
      write(22,FMT)'%.o: %.f'
      write(22,FMT1)TAB,'$(FC) $(OPTIONS) -c -o $@ $<'


      close(22)

      return
      end subroutine tocompile_gfortran

!------------------------------------------------------------------------------------------!

      subroutine expression
      include 'common.gen2D'


      
      end subroutine expression
