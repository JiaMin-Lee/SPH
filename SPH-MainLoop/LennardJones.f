        subroutine LennardJonesBC(iip,jjp,rr2,drx,drz,fxbp,fzbp)
        include 'common.2D'

        rr0 = 0.5*min(dx0,dz0)
        dd = 1e-2
        p1 = 12
        p2 =  4
        rr = sqrt(rr2)

        ! drx must be xp(FP) - xp(BP) 
        ! drz must be zp(FP) - zp(BP)

        if (rr.lt.rr0) then 

           fbp = ((rr0/rr)**p1 - (rr0/rr)**p2)/rr2
           fxbp = fbp * drx * dd
           fzbp = fbp * drz * dd

        else

           fxbp = 0
           fzbp = 0

        endif

        end subroutine LennardJonesBC
