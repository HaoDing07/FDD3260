       zz = 0.5*dz*fdz0(1)
       w0(1) = 0.
       do k = 1, nz
       
        w0(k)  = -Ddiv*zz
        if (k <= kpert+1) then
         k0(k) = 1.
        else
         k0(k) = 0.
        endif
!      
       if (k/=nz) zz = zz + 0.5*dz*(fdz0(k)+fdz0(k+1))      
       enddo
