c    "Copyright 2009 Prof. Robert Dalrymple, Prof. M. Gomez Gesteira, Dr Benedict Rogers, 
c     Dr Alejandro Crespo, Dr Muthukumar Narayanaswamy, Dr Shan Zou, Dr Andrea Panizzo "
c
c    This file is part of SPHYSICS.
c
c    SPHYSICS is free software; you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation; either version 3 of the License, or
c    (at your option) any later version.
c
c    SPHYSICS is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.
c
c    You should have received a copy of the GNU General Public License
c    along with this program.  If not, see <http://www.gnu.org/licenses/>.


      subroutine check_limits_2D
      include 'common.2D'
      
      double precision dxp_dble

      if (i_FlowDir.eq.1) xmax=0
      if (i_FlowDir.eq.2) zmax=0
      ncases=0

      xlimit = 100000

      do i_num_Buff = 1, num_Buff
         if (ioptBuff(i_num_Buff).eq.2) then
            xlimit = Buffint(i_num_Buff)
            nlimits = nBuff_ini(i_num_Buff)
            nlimite = nBuff_end(i_num_Buff)
         endif
      enddo

      do i=nb+1,np ! Check if any fluid particle is over the box

         if (i_FlowDir.eq.2.and.zp(i).gt.zmax) zmax=zp(i)

         if (i_FlowDir.eq.1.and.xp(i).gt.xmax) xmax=xp(i)
 
!         if (iflag(i).eq.1.and.(xp(i).gt.xmax_container.or.
!     +         xp(i).lt.xmin_container.or.zp(i).lt.zmin_container)) then
         
         if (iflag(i).eq.1.and.itype(i).eq.2.and.
     +       (xp(i).lt.xmin_container.or.
     +        zp(i).lt.zmin_container.or.
     +        (i_FlowDir.eq.2.and.xp(i).gt.xmax_container).or.
     +        (i_FlowDir.eq.1.and.zp(i).gt.zmax_container).or.
     +        (i_FlowDir.eq.1.and.i_periodicOBs(1).eq.0.and.
     +         xp(i).ge.xlimit))) then

            write(80,*) 'Particle out of limits: ',i, 
     +      '  last X position: ',xp(i), '  last Z-position: ',zp(i)
            write(*,*) 'Particle out of limits: ',i, 
     +      '  last X position: ',xp(i), '  last Z-position: ',zp(i)
       
            iflag(i)=0
            ncases=ncases+1

            up(i)=0.
            wp(i)=0.
            um1(i)=0.
            uo(i)=0.
            wm1(i)=0.
            wo(i)=0.

         endif       
      enddo

      if (i_periodicOBs(1).eq.1) then

         ! if particle is outside boundary face, re-enters from
         ! OPPOSITE boundary face
         do ip = nb+1,np

            if (iflag(ip).eq.1.and.itype(ip).eq.2.and.
     &          xp(ip).gt.xmax_container) then

               dxp_dble = dble(xp(ip) - xmax_container_double)
               xp(ip) = real(dxp_dble + xmin_container_double)
               
            end if

            if (xp(ip).le.(ikappa*hsml).and.xp(ip).gt.0..and.
     &          itype(ip).eq.2.and.iVelocity.eq.1) then

!                  write(*,*) 'ip = ',ip
!                if (zp(ip).gt.0..and.zp(ip).le.0.2e-3) then
!                   up(ip) = 29./2.*zp(ip)
!                elseif(zp(ip).lt.1.5e-3.and.zp(ip).ge.1.3e-3) then
!                   up(ip) = 29./2.*abs(zp(ip)-1.5e-3)
!                else
!                   up(ip) = VXInlet
!                endif

                up(ip) = (-7.7295*zp(ip)*zp(ip)+0.0116*zp(ip))*1000
!                wp(ip) = VZInlet               

            endif

         enddo

      endif

      if (i_FlowDir.eq.2) then      
         differ=(zmax-zcontrol)
   
         if (differ.gt.0) then     
             nn=int(differ*one_over_kh)+1
             zmax = zcontrol + ikappa*hsml*nn
             ncz = int( (zmax-zmin)*one_over_kh ) + 1
             nct = ncx*ncz !2D

             do i=nct_ini+1,nct
                nc(i,1)  = 0 !No Fixed Boundary Particles in these cells               
             enddo
          else
             zmax=zcontrol
             nct=nct_ini
             ncz=ncz_ini
          endif 

       elseif (i_FlowDir.eq.1.and.i_periodicOBs(1).eq.0) then
          differ = xmax-xcontrol

          if (differ.gt.0) then
             nn = int(differ*one_over_kh) + 1
             xmax = xcontrol + ikappa*hsml*nn
             ncx = int( (xmax-xmin)*one_over_kh) + 1
             nct = ncx*ncz

             write(*,*) nct

             do ii = nct_ini+1,nct
                nc(ii,1) = 0
             enddo
          else
             xmax = xcontrol
             nct = nct_ini
             ncx = ncx_ini
          endif   
 
       endif


      if (ncases.ne.0) write(80,*) 'Number ',ncases,' time',time
      if (ncases.ne.0) write(*,*) 'Number ',ncases,' time',time

      if (ncases.ge.100000) stop

       if(nct.gt.nct_max.or.ncz.lt.1)then
         write(80,*) ' '
         write(80,*)'ERROR in check_limits_2D.f'
         write(80,*)'nct.gt.nct_max.or.ncz.lt.1'
         write(80,*)'ncx/z    = ',ncx,ncz
         write(80,*)'nct= ',nct
         write(80,*)'nct_max ',nct_max
         write(80,*)'itime  ', itime
!!!      -- screen printout ---              
         write(*,*) ' '
         write(*,*)'ERROR in check_limits_2D.f'
         write(*,*)'nct.gt.nct_max.or.ncz.lt.1'
         write(*,*)'ncx/z    = ',ncx,ncz
         write(*,*)'nct= ',nct
         write(*,*)'nct_max ',nct_max
         write(*,*)'itime  ', itime
         stop
       endif

      return
	end
