       xrbubble = 250.
       zrbubble = 250.
       xbubble = 1800.
       zbubble = 800.
       dtbubble = 0.
       
       zz = 0.
       do k = 1, nz
	 pt0(k) = 296.5*exp(1.3e-5*zz)

         u00(k) = 0.
         v00(k) = 0.
         w0(k)  = 0.
	 k0(k)  = 0.
	 qt0(k) = 0.
	 zz = zz + dz
       enddo
