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

      subroutine self(j1,kind_p1,ini_kind_p2)
c
      include 'common.2D'
      
!     For kind_p1 = 1, ini_kind_p2 = 2
!     For kind_p1 = 2, ini_kind_p2 = 1

!     1->1
!     1->2
!     2->2

      if (kind_p1.eq.1) then
         kind_p2_start = 1
      else
         kind_p2_start = ini_kind_p2
      endif


      do kind_p2=kind_p2_start,2
    
         jj_start = 1

         if (kind_p1.gt.kind_p2) then
                                            ! if kind_p1 is FP, kind_p2 is BP
            jj_start = nplink_max + 1       ! jj_start will be more than
                                            ! maximum, loop CANNOT START
         endif
        
         do ii=1,nc(j1,kind_p1)

            ip = ibox(j1,kind_p1,ii)       ! ip is particle itself
         
            if (kind_p1.eq.kind_p2)then

               jj_start = jj_start + 1       ! If fluid-fluid interaction jj_start = 2

            end if
         
            do jj = jj_start,nc(j1,kind_p2)  
        
               jp = ibox(j1,kind_p2,jj)      ! jp is neighbour particle

               drx = xp(ip) - xp(jp)
               drz = zp(ip) - zp(jp)
               rr2 = drx*drx + drz*drz

               if (rr2.lt.rrk2h2.and.rr2.gt.1.e-18) then

                  dux = up(ip) - up(jp)
                  duz = wp(ip) - wp(jp)

!!!               Calculating kernel & Normalized Kernel Gradient

                  call kernel(drx,drz,ip,jp,j1,j1,rr2) 
                  call kernel_correction(ip,jp)

!!!               Average Density and Speed of Sound

                  robar  = 0.5*( rhop(ip) + rhop(jp) )
                  one_over_rhobar = 2.0/(rhop(ip) + rhop(jp))

                  cbar = 0.5*( cs(ip) + cs(jp) ) 

!!!               Inner product of Relative Distance . Relative Velocity

                  dot = drx*dux + drz*duz

!!!               Used to calculate the time step due to viscosity

                  visc_dt=max(dot/(rr2 + eta2),visc_dt) 


!!!               PRESSURE  FORCES (1992; 3.3)

                  p_v = pr(ip) + pr(jp)  !+  pi_visc

 
!                 For TENSILE CORRECTION (Monaghan , JCP.  159 (2000) 290- 311)
!                 Only to be activated with CUBIC SPLINE kernel 

                  if (index_tensile.eq.1) then   

                     fab=Wab*od_Wdeltap                 !NOTE: We'll use a non-normalized
                     fab=fab*fab                        !kernel to calculate tensile correction
                     fab=fab*fab                        !  It's the (Wab/Wdeltap)**4  of Monaghan's paper

                     if (pp(ip).gt.0) then
                        Ra= 0.01 *pr(ip)
                     else
                        Ra= 0.2 *abs(pr(ip))
                     endif

                     if (pp(jp).gt.0) then
                        Rb= 0.01 *pr(jp)
                     else
                        Rb= 0.2 *abs(pr(jp))
                     endif
 
                     p_v = p_v+ (Ra+Rb)*fab

                  endif

!                  if (ip.gt.nvmb.and.jp.gt.nvmb) then   !both fluid particles and virtual moving boundary particles
                  if (itype(ip).ge.2.and.itype(ip).ge.2) then  
 
                     ax(ip) = ax(ip) - pm(jp) * p_v * frxi
                     az(ip) = az(ip) - pm(jp) * p_v * frzi

                     ax(jp) = ax(jp) + pm(ip) * p_v * frxj
                     az(jp) = az(jp) + pm(ip) * p_v * frzj

                     call gradients_calc(ip,jp,dux,duz)               

                     call viscosity(dot,drx,drz,dux,duz,rr2,cbar,
     +                robar,one_over_rhobar,ip,jp,j1,j1,term2i,term2j)

!!!                  DENSITY ACCELERATION (1992; 3.9)

                     dot2i = dux*frxi + duz*frzi
                     dot2j = dux*frxj + duz*frzj
                     ar(ip) = ar(ip) + pm(jp)*dot2i
                     ar(jp) = ar(jp) + pm(ip)*dot2j       

!!!                  THERMAL ENERGY
!                    (Monaghan, JCP 110 (1994) 399- 406)

                     term1i=0.5 * p_v *( frxi*dux+frzi*duz)
                     term1j=0.5 * p_v *( frxj*dux+frzj*duz)
  
                     aTE(ip)=aTE(ip)+pm(jp) * (term1i+term2i)
                     aTE(jp)=aTE(jp)+pm(ip) * (term1j+term2j)

!!!                  XSPH correction

                     pmj_Wab_over_rhobar = pm(jp)*Wab*one_over_rhobar
                     ux(ip) = ux(ip) - dux*pmj_Wab_over_rhobar  !pm(j) * dux * Wab / robar
                     wx(ip) = wx(ip) - duz*pmj_Wab_over_rhobar  !pm(j) * duz * Wab / robar

        
                     pmi_Wab_over_rhobar = pm(ip)*Wab*one_over_rhobar
                     ux(jp) = ux(jp) + dux*pmi_Wab_over_rhobar   !pm(i) * dux * Wab / robar
                     wx(jp) = wx(jp) + duz*pmi_Wab_over_rhobar   !pm(i) * duz * Wab / robar

                  ! ip is FLUID and jp is BOUNDARY

!                  else if (ip.gt.nvmb.and.jp.le.nvmb) then             
                  else if (itype(ip).ge.2.and.abs(itype(jp)).le.1) then 

                     write(*,*) ' NOT ALLOWED HERE'

                  ! ip is BOUNDARY and jp is FLUID

!                  else if (jp.gt.nvmb.and.ip.le.nvmb) then 
                  else if (itype(jp).ge.2.and.abs(itype(ip)).le.1) then 

                     ax(ip) = ax(ip) - pm(jp) * p_v * frxi
                     az(ip) = az(ip) - pm(jp) * p_v * frzi

                     ax(jp) = ax(jp) + pm(ip) * p_v * frxj
                     az(jp) = az(jp) + pm(ip) * p_v * frzj

                     call gradients_calc(ip,jp,dux,duz)

                     call viscosity(dot,drx,drz,dux,duz,rr2,cbar,
     +              robar,one_over_rhobar,ip,jp,j1,j1,term2i,term2j)

!!!                  REPULSIVE FORCE only on FLUID PARTICLES
                     if (itype(ip).eq.1.and.itype(jp).ge.2) then

                        call CDSBT(jp,ip,rr2,-drx,-drz,
     &                                      fxbp,fzbp)
                        ax(jp) = ax(jp) + fxbp
                        az(jp) = az(jp) + fzbp

                     endif

                     dot2i = dux*frxi + duz*frzi
                     dot2j = dux*frxj + duz*frzj
                     ar(ip) = ar(ip) + pm(jp)*dot2i
                     ar(jp) = ar(jp) + pm(ip)*dot2j

                  ! Both BOUNDARY
                  else if (abs(itype(ip)).le.1.and.
     &                     abs(itype(jp)).le.1) then    
                     
                     ! ip is GHOST and jp is REPULSIVE
                     if (itype(ip).lt.0.and.itype(jp).eq.1) then

                        write(*,*) 'CANNOT BE HERE'

                        ax(ip) = ax(ip) - pm(jp) * p_v * frxi
                        az(ip) = az(ip) - pm(jp) * p_v * frzi

                        call gradients_calc(ip,jp,dux,duz)

                        call viscosityCDSBT(dot,drx,drz,dux,duz,rr2,
     +                                  ip,jp,term2i,term2j)
                        
!                        call viscosity(dot,drx,drz,dux,duz,rr2,cbar,
!     +                 robar,one_over_rhobar,ip,jp,j1,j1,term2i,term2j)

                        call CDSBT(ip,jp,rr2,drx,drz,
     &                                      fxbp,fzbp)
                        ax(ip) = ax(ip) + fxbp
                        az(ip) = az(ip) + fzbp

                        dot2i = dux*frxi + duz*frzi
                        ar(ip) = ar(ip) + pm(jp)*dot2i

!                        dot2j = dux*frxj + duz*frzj
!                        ar(jp) = ar(jp) + pm(ip)*dot2j


                     ! jp is GHOST and ip is REPULSIVE
                     elseif (itype(ip).eq.1.and.itype(jp).lt.0) then

                        ax(jp) = ax(jp) + pm(ip) * p_v * frxj
                        az(jp) = az(jp) + pm(ip) * p_v * frzj

                        call gradients_calc(ip,jp,dux,duz)

                        call viscosityCDSBT(dot,drx,drz,-dux,-duz,rr2,
     +                                  jp,ip,term2i,term2j)

!                        call viscosity(dot,drx,drz,dux,duz,rr2,cbar,
!     +                 robar,one_over_rhobar,ip,jp,j1,j1,term2i,term2j)


                        call CDSBT(jp,ip,rr2,-drx,-drz,
     &                                      fxbp,fzbp)
                        ax(jp) = ax(jp) + fxbp
                        az(jp) = az(jp) + fzbp

                        dot2j = dux*frxj + duz*frzj
                        ar(jp) = ar(jp) + pm(ip)*dot2j

!                        dot2i = dux*frxi + duz*frzi
!                        ar(ip) = ar(ip) + pm(jp)*dot2i

                     endif
 
                  endif ! Interaction fluid-fluid or fluid-boundary
	        endif ! if q<2
             enddo ! Box jj
          enddo   ! Box ii
      enddo   ! Kind of particle

      end
c
