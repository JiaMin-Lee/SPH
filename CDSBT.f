        subroutine CDSBT(iip,jjp,rr2,drx,drz,fxbp,fzbp)
        include 'common.2D'

        rr = sqrt(rr2)
        etaa = rr / (0.75*hsml)
        coD = 1.0

        ! drx must be xp(FP/GP) - xp(BP) 
        ! drz must be zp(FP/GP) - zp(BP)

        if (rr.lt.(coD*dx0).and.rr.gt.0.) then

           chii = 1. - rr/(coD*dx0)
            
           if (etaa.le.(2./3.).and.etaa.gt.0.) then 

              ffeta = 2./3. 
                
           elseif (etaa.gt.(2./3.).and.etaa.le.1.) then

              ffeta = 2.*etaa - 1.5*etaa*etaa
           
           elseif (etaa.gt.(1.).and.etaa.lt.2.) then

              ffeta = 0.5*(2.-etaa)*(2.-etaa)

           else

              ffeta = 0.

           endif

           fbp = 0.01*cs0*cs0*chii*ffeta
           fxbp = fbp*drx/rr2
           fzbp = fbp*drz/rr2

        else

           chii = 0.

           fxbp = 0.
           fzbp = 0.

        endif

        end subroutine CDSBT
