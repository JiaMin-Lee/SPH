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

      subroutine ac  ! ac_NoCorrection
c
      include 'common.2D'

      ! avoids compiling kernel corrections subroutines when not required 
      icall = 0

      do i=nstart,np

        aM_a11(i)=0.    !Mi(1,1)
        aM_a12(i)=0.    !Mi(1,2)
        aM_a21(i)=0.    !Mi(2,1)
        aM_a22(i)=0.    !Mi(2,2)
        aL_a11(i)=0.    !Li(1,1)
        aL_a12(i)=0.    !Li(1,2)
        aL_a21(i)=0.    !Li(2,1)
        aL_a22(i)=0.    !Li(2,2)


      enddo

      if (i_kernelcorrection.eq.1) then

         do kind_p1=1,2
      
            ini_kind_p2=mod(kind_p1,2)+1
      
            do lz=1,ncz
      
               do lx=1,ncx
      
                  if (i_FlowDir.eq.1) then
                     j1 = lz + (lx-1)*ncz
                  elseif (i_FlowDir.eq.2) then
                     j1 = lx + (lz-1)*ncx
                  endif
      
                  if (nc(j1,kind_p1).gt.0) then     ! if the cell is not empty, then
                                                    ! loop over it and over neighboring
                                                    ! cells
      
                     lx2 = lx+1
      
                     if(lx2.le.ncx)then             ! EAST CELL EXIST
                       call pre_celij(j1,j1+ncall1,kind_p1,ini_kind_p2,
     &                                lx2)
                     endif
      
                     lz2=lz+1
      
                     if(lz2.le.ncz)then
      
                       ! NORTH CELL 
                       call pre_celij(j1,j1+ncall2,kind_p1,ini_kind_p2,
     &                                lx2)
      
                       lx2=lx+1
                       if(lx2.gt.2) call pre_celij(j1,j1+ncall4,
     &                                       kind_p1,ini_kind_p2,lx2) !North-West        
      
                       lx2=lx+1
                       if(lx2.le.ncx) call pre_celij(j1,j1+ncall3,
     &                                       kind_p1,ini_kind_p2,lx2) !North-East        
                      
                     endif    !End of:   if(lz2.le.ncz)then
      
                   endif    !End of:  if(nc(j1,kind_p1).gt.0) then
               enddo
            enddo
      
            do j1=1,nct
               if (nc(j1,kind_p1).gt.0) then
                  call pre_self(j1,kind_p1,ini_kind_p2)
               endif
            enddo
      
            if (i_periodicOBs(1).eq.1) then

               if (i_FlowDir.eq.1) then
                  j1 = lz + (lx-1)*ncz
               elseif (i_FlowDir.eq.2) then
                  j1 = lx + (lz-1)*ncx
               endif

               !- Special Treatment for Right Column Cells (lx=ncx) -
               lx = ncx
               lx2 = lx+1
               do lz = 1, ncz-1

                  if (i_FlowDir.eq.1) then      ! Flow Right
                     j1 = lz + (lx-1)*ncz
                  elseif (i_FlowDir.eq.2) then  ! Flow Up
                     j1 = lx + (lz-1)*ncx
                  endif

                  call pre_celij(j1,j1+ncall11,kind_p1,ini_kind_p2,lx2) !East X-Periodic
                  call pre_celij(j1,j1+ncall12,kind_p1,ini_kind_p2,lx2) !North-East X-Periodic
               end do    !End of:  do lz = 1,ncz-1

               !- Special Treatment for Corner Cell (lx=ncx, lz=ncz) -
               !- Note this is for X-Periodicity Only!               -
               lz = ncz
               if (i_FlowDir.eq.1) then      ! Flow Right
                  j1 = lz + (lx-1)*ncz
               elseif (i_FlowDir.eq.2) then  ! Flow Up
                  j1 = lx + (lz-1)*ncx
               endif
               call pre_celij(j1,j1+ncall11,kind_p1,ini_kind_p2,lx2) !East X-Periodic

               !- Special Treatment for Left Column Cells (lx=1) -
               lx = 1
               lx2 = lx+1
               do lz = 1,ncz-1

                  if (i_FlowDir.eq.1) then      ! Flow Right
                     j1 = lz + (lx-1)*ncz
                  elseif (i_FlowDir.eq.2) then  ! Flow Up
                     j1 = lx + (lz-1)*ncx
                  endif

                  call pre_celij(j1,j1+ncall13,kind_p1,ini_kind_p2,lx2) !North-West X-Periodic
               end do    !End of:  do lz = 1,ncz-1

            endif   !End of:  if(i_periodicOBs(1).eq.1)then
         enddo

         !- KGC inverse matrix L = M^(-1)

         do i=nstart,np
      
!            if (itype(i).ge.2.and.iflag(i).gt.0) then !Only for Fluid particles
            if (iflag(i).gt.0) then                                    !inside the limits

               aM_a12(i)=0.5*(aM_a12(i)+aM_a21(i)) ! M should be symmetric
               aM_a21(i)=aM_a12(i)
      
               detM=aM_a11(i)*aM_a22(i)-aM_a12(i)*aM_a21(i)  !Determinant of matrix M det(M)
               
               if (abs(detM).gt.0.01.and.abs(aM_a11(i)).gt.0.25
     +             .and.abs(aM_a22(i)).gt.0.25) then  !Only for "well conditioned" matrices

                  one_over_detM = 1.0/detM
                  aL_a11(i) =  aM_a22(i)*one_over_detM !Li(1,1)=Mi(1,1)/det(M)
                  aL_a22(i) =  aM_a11(i)*one_over_detM !Li(2,2)=Mi(2,2)/det(M)
                  aL_a12(i) = -aM_a12(i)*one_over_detM !Li(1,2)=-Mi(1,2)/det(M)
                  aL_a21(i) = -aM_a21(i)*one_over_detM !Li(2,1)=-Mi(2,1)/det(M)
      
               else    !No correction
                  aL_a11(i) = 1.
                  aL_a12(i) = 0.
                  aL_a21(i) = 0.
                  aL_a22(i) = 1.
      
               endif
      
            else !No correction
               aL_a11(i) = 1.
               aL_a12(i) = 0.
               aL_a21(i) = 0.
               aL_a22(i) = 1.
            endif

         enddo

      endif  

      call ac_main
        
      return

      end subroutine ac
