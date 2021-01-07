        subroutine equation_of_state(rho_EoS,TE_EoS,press_EoS,cs_EoS)
        include 'common.2D'

        if (i_EoS.eq.1) then

           press_EoS =  cs0 * cs0 *(rho_EoS-rho0)
           cs_EoS = cs0

        elseif (i_EoS.eq.2) then
           
           psicoeff = rho0 * cs0 * cs0 / igamma
           press_EoS = psicoeff * ((rho_EoS/rho0)**igamma-1.0)
           cs_EoS = cs0

        endif

        end subroutine equation_of_state
