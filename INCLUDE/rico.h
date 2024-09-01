!
!  Core values
!
       zz = 0.5*dz*fdz0(1)
       zw = 0.
       w0(1) = 0.
       do k = 1, nz
         u00(k) = -9.9 + 8.*zz/4000.
         v00(k) = -3.8
         w0(k)  = max(-0.5,-0.5*zw/2260.)*1.e-2
!
         if (zz <= 740.) then
          pt0(k) = 297.9
          qv0(k) = (16. - 2.2*zz/740.)*1.e-3
         else if (zz <= 3260.) then
          pt0(k) = 297.9 + 19.1*(zz - 740.)/3260.
          qv0(k) = (13.8 - 11.4*(zz - 740.)/2520.)*1.e-3	 
	 else if (zz > 3260.) then
          pt0(k) = 297.9 + 19.1*(zz - 740.)/3260.
          qv0(k) = (2.4 - 0.6*(zz - 3260.)/740.)*1.e-3
         endif
!
         k0(k) = 0.01
!
         if (k < nz) zz = zz + 0.5*dz*(fdz0(k)+fdz0(k+1))	   
         if (k < nz) zw = zw + dz*fdz0(k)
       enddo
