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

      subroutine ac_main
c
      include 'common.2D'

c  ...  store useful arrays
c
      !- Need to zero each object for multiobjects -
      bigUdot   = 0.0
      bigWdot   = 0.0
        
      X_Friction   = 0.0
      Z_Friction   = 0.0
        
      nb_inFriction = 0

      do i=nbfm+1,nb
         pr(i) = pp(i)/(rhop(i)*rhop(i))

         one_over_rhop(i) = 1.0 / rhop(i)

c        -- Zeroing Variables for Free-Moving Objects --         
         ax(i) = 0.
         az(i) = 0.

         ar(i) = 0.

         ux(i) = 0.
         wx(i) = 0.


         aTE(i)=0.
      enddo 

      do ii = 1,np!nstart,np

         pr(ii) = pp(ii)/(rhop(ii)*rhop(ii))

         one_over_rhop(ii) = 1.0 / rhop(ii)

c        -- Zeroing Variables --         
         ax(ii) = 0.
         az(ii) = 0.
              
         ar(ii) = 0.
              
         ux(ii) = 0.
         wx(ii) = 0.

         
         aTE(ii)=0.

         sum_wab(ii) = 0.
         rhop_sum(ii) = 0.

      enddo

           

       do kind_p1=1,2

          ini_kind_p2=mod(kind_p1,2)+1   ! when kind_p1=1, ini_kind_p2=2
                                         ! when kind_p1=2, ini_kind_p2=1 
       
          do lz = 1,ncz

            do lx = 1,ncx

               if (i_FlowDir.eq.1) then
                  j1 = lz + (lx-1)*ncz
               elseif (i_FlowDir.eq.2) then    
                  j1 = lx + (lz-1)*ncx
               endif
         
               if (nc(j1,kind_p1).gt.0) then       ! if the cell is not empty, then
                                                   ! loop over it and over neighboring
                                                   ! cells

                  lx2 = lx + 1

                  if (lx2.le.ncx) then             ! EAST CELL EXISTS
                     call celij(j1,j1+ncall1,kind_p1,ini_kind_p2,lx2) !East
                  endif      
      
                  lz2 = lz + 1
                  if (lz2.le.ncz) then             ! NORTH CELL EXISTS
       
                     call celij(j1,j1+ncall2,kind_p1,ini_kind_p2,lx2)   !North
       
                     lx2=lx+1

                     if(lx2.gt.2) call celij(j1,j1+ncall4,
     &                                       kind_p1,ini_kind_p2,lx2)   !North-West

                     lx2=lx+1
                     if(lx2.le.ncx) call celij(j1,j1+ncall3,
     &                                         kind_p1,ini_kind_p2,lx2) !North-East                 
       
                  endif    !End of:   if(lz2.le.ncz)then
               endif    !End of:  if(nc(j1,kind_p1).gt.0) then
            enddo
          enddo
       
          do j1=1,nct
             if (nc(j1,kind_p1).gt.0) 
     +          call self(j1,kind_p1,ini_kind_p2)
          enddo

          if (i_periodicOBs(1).eq.1) then
         !- Special Treatment for Right Column Cells (lx=ncx) -
             lx = ncx
             lx2 = lx+1

             do lz = 1,ncz-1

                if (i_FlowDir.eq.1) then      ! Flow Right
                   j1 = lz + (lx-1)*ncz
                elseif (i_FlowDir.eq.2) then  ! Flow Up
                   j1 = lx + (lz-1)*ncx
                endif

                call celij(j1,j1+ncall11,kind_p1,ini_kind_p2,lx2) !East X-Periodic
                call celij(j1,j1+ncall12,kind_p1,ini_kind_p2,lx2) !North-East X-Periodic
             end do    !End of:  do lz = 1,ncz-1
 
             !- Special Treatment for Corner Cell (lx=ncx, lz=ncz) -
             !- Note this is for X-Periodicity Only!               -
             !Easts of lx = 1
             lz = ncz

             if (i_FlowDir.eq.1) then      ! Flow Right
                j1 = lz + (lx-1)*ncz
             elseif (i_FlowDir.eq.2) then  ! Flow Up
                j1 = lx + (lz-1)*ncx
             endif

             call celij(j1,j1+ncall11,kind_p1,ini_kind_p2,lx2) !East X-Periodic
 
             !- Special Treatment for Left Column Cells (lx=1) -
             lx = 1
             lx2 = lx+1
             do lz = 1,ncz-1

                if (i_FlowDir.eq.1) then      ! Flow Right
                   j1 = lz + (lx-1)*ncz
                elseif (i_FlowDir.eq.2) then  ! Flow Up
                   j1 = lx + (lz-1)*ncx
                endif

                call celij(j1,j1+ncall13,kind_p1,ini_kind_p2,lx2) !North-West X-Periodic
             end do    !End of:  do lz = 1,ncz-1
          endif   !End of:  if(i_periodicOBs(1).eq.1)then
                    
       enddo   !End of:  do kind_p1=1,2
       
       do ii = nstart,np
    
!          if (itype(ii).ge.2) then
     
             udot(ii) = ax(ii)
             wdot(ii) = az(ii)
             rdot(ii) = ar(ii)     
      
             xcor(ii) = eps*ux(ii)
             zcor(ii) = eps*wx(ii)
       
             TEdot(ii)=aTE(ii)

!           endif 
         
       enddo
             
       return
       end

