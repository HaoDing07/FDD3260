!
       zz = 0.5*dz*fdz0(1)
       u00 = 6.
       v00 = 0.
       w0 = 0.
       k0 = 0.01
       do k = 1, nz
!
         if (zz <= 12200.) then
           t0(k) = 291.65 - 0.007*zz
         else
           t0(k) = 291.65 - 0.007*12200.
         endif
!
         if (zz <= 2000.) then
           rh0(k) = 0.85
         else if ( zz > 2000. .and. zz <= 8000.) then
           rh0(k) = 0.85 - 0.15*(1. + tanh((zz - 5000.)/750.))
         else
           rh0(k) = 0.55
         endif
!
         if (k < nz) zz = zz + 0.5*dz*(fdz0(k)+fdz0(k+1))	   
       enddo
!
       with_t = .true.
       with_rh = .true.

