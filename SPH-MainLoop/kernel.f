        subroutine kernel(drx,drz,ii,jj,j1,j2,rr2)

        include 'common.2D'

!!!     -- CUBLIC SPLINE KERNEL

        rad = sqrt(rr2)
        qq = rad/hsml         

        if (i_kernel.eq.1) then
          index_tensile = 0!1

          alphad = 10./(7.*pi)
          alphadh2 = alphad/(hsml*hsml)

          if (qq.ge.0.and.qq.lt.1.) then
             
             Wab = alphadh2 * (1. - 1.5*qq*qq + 0.75*qq*qq*qq)

             fab = alphadh2 * (-3.*qq - 2.25*qq*qq)
             frx = fab * drx / (hsml * rad)
             frz = fab * drz / (hsml * rad)

          else if (qq.ge.1..and.qq.lt.2.) then

             Wab = alphadh2 * 0.25 * ((2.-qq)*(2.-qq)*(2.-qq))

             fab = - alphadh2 * 0.75 *((2.-qq)*(2.-qq))
             frx = fab * drx / (hsml * rad)
             frz = fab * drz / (hsml * rad)
        
          else

             Wab = 0.
             frx = 0.
             frz = 0.

          endif

        else if (i_kernel.eq.2) then

          index_tensile = 0

          alphad = 7./(478.*pi)
          alphadh2 = alphad/(hsml*hsml)

          termm1 = (3.-qq)*(3.-qq)*(3.-qq)*(3.-qq)*(3.-qq)
          termm2 = (2.-qq)*(2.-qq)*(2.-qq)*(2.-qq)*(2.-qq)
          termm3 = (1.-qq)*(1.-qq)*(1.-qq)*(1.-qq)*(1.-qq)

          ftermm1 = (3.-qq)*(3.-qq)*(3.-qq)*(3.-qq)
          ftermm2 = (2.-qq)*(2.-qq)*(2.-qq)*(2.-qq)
          ftermm3 = (1.-qq)*(1.-qq)*(1.-qq)*(1.-qq)

          if (qq.ge.0.and.qq.lt.1) then
     
             Wab = alphadh2 * (termm1 - 6.*termm2 + 15.*termm3)

             fab = alphadh2 * (-5.*ftermm1 + 30.*ftermm2 - 75.*ftermm3)
             frx = fab * drx / (hsml * rad)
             frz = fab * drz / (hsml * rad)

          else if (qq.ge.1..and.qq.lt.2.) then

             Wab = alphadh2 * (termm1 - 6.*termm2)

             fab = alphadh2 * (-5.*ftermm1 + 30.*ftermm2)
             frx = fab * drx / (hsml * rad)
             frz = fab * drz / (hsml * rad)


          else if (qq.ge.2..and.qq.lt.3.) then

             Wab = alphadh2 * termm1
        
             fab = alphadh2 * (-5.*ftermm1)
             frx = fab * drx / (hsml * rad)
             frz = fab * drz / (hsml * rad)

          else

             Wab = 0.
             frx = 0.
             frz = 0.

          endif

        endif        

        end subroutine kernel
