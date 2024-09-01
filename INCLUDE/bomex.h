!
!  Core values
!
       zz = 0.5*dz*fdz0(1)
       zw = 0.
       w0(1) = 0.
       do k = 1, nz
	 if (zz <= 700.) then
           u00(k) = -8.75
	 else
	   u00(k) = -8.75 + 4.14 * (zz - 700.) / 2300.
	 endif
         v00(k) = 0.0
	 if (zz <= 1500.) then
           w0(k) = -0.65 * zw / 1500.
	 else if (zz <= 2100.) then
	   w0(k) = -0.65 + 0.65 * (zw - 1500.) / 600.
	 else
	   w0(k) = 0.
	 endif
!
         if (zz <= 520.) then
          pt0(k) = 298.7
          qv0(k) = (17. - 0.7*zz/520.)*1.e-3
         else if (zz <= 1480.) then
          pt0(k) = 298.7 + 3.7*(zz - 520.)/960.
          qv0(k) = (16.3 - 5.6*(zz - 520.)/960.)*1.e-3	 
	 else if (zz <= 2000.) then
          pt0(k) = 302.4 + 5.8*(zz - 1480.)/520.
          qv0(k) = (10.7 - 6.5*(zz - 1480.)/520.)*1.e-3
	 else
          pt0(k) = 308.2 + 3.65*(zz - 2000.)/1000.
          qv0(k) = (4.2 -  1.2*(zz - 2000.)/1000.)*1.e-3
         endif
!
         k0(k) = 0.01
!
         if (k < nz) zz = zz + 0.5*dz*(fdz0(k)+fdz0(k+1))	   
         if (k < nz) zw = zw + dz*fdz0(k)
       enddo
