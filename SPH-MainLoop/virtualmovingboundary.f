        subroutine virtualmovingboundary
        include 'common.2D'
        
        do iVMB = 1, num_MB

           ilayer = 1

           do iip = nVMB_ini(iVMB), nVMB_end(iVMB)

              nperlayer = nVMB_l1(iVMB)-nVMB_ini(iVMB)+1

              ndiff = iip - nVMB_ini(iVMB) + 1

              if (iVMBwall(iVMB).eq.1) then !build from right to left

                 npart = nb + ndiff - nperlayer*(ilayer-1)

                 up(iip) = up(npart)
                 wp(iip) = wp(npart)

                 xp(iip) = xp(npart) - dx0*ilayer
                 zp(iip) = zp(npart)

              elseif (iVMBwall(iVMB).eq.2) then !build from left to right

                 npart = np + ndiff - nperlayer*ilayer
        
                 up(iip) = up(npart)
                 wp(iip) = wp(npart)

                 xp(iip) = xp(npart) + dx0*ilayer
                 zp(iip) = zp(npart)

              endif

              if (mod(ndiff,nperlayer).eq.0) then

                 ilayer = ilayer + 1

              endif

           enddo

        enddo

        end subroutine virtualmovingboundary
