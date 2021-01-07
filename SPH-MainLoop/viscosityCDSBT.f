      subroutine viscosityCDSBT(dot,drx,drz,dux,duz,rr2,
     +           ip,jp,term2i,term2j)

      include 'common.2D'

!         --- Viscous Diffusion terms (Lo & Shao 2002) ---
!     For CDSBT particles interaction, ghost particles only
!     ip = GP, jp = RP

      tempi =2.0*viscos_val*one_over_rhobar*
     &     ((drx*frxi+drz*frzi)/(rr2 + eta2))

!      tempj =2.0*viscos_val*one_over_rhobar*
!     &     ((drx*frxj+drz*frzj)/(rr2 + eta2))

      ax(ip) = ax(ip) + pm(jp)*tempi*dux
      az(ip) = az(ip) + pm(jp)*tempi*duz

!     ax(j) = ax(j) - pm(i)*tempj*dux
!     az(j) = az(j) - pm(i)*tempj*duz

      term2i=  -0.5 * tempi *( dux*dux+duz*duz)
      term2j=  -0.5 * tempj *( dux*dux+duz*duz)

      return

      end subroutine viscosityCDSBT

