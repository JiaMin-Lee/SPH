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

      subroutine kernel_correction(ip,jp)

      include 'common.2D'   

      if (i_kernelcorrection.eq.0) then

         frxi=frx
         frzi=frz

         frxj=frx
         frzj=frz

      elseif (i_kernelcorrection.eq.1) then       

         frxi=aL_a11(ip)*frx+aL_a12(ip)*frz
         frzi=aL_a21(ip)*frx+aL_a22(ip)*frz
       
         frxj=aL_a11(jp)*frx+aL_a12(jp)*frz
         frzj=aL_a21(jp)*frx+aL_a22(jp)*frz

      endif

      return
      end subroutine kernel_correction
