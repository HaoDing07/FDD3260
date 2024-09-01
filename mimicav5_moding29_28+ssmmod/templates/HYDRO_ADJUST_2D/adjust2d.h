       zz = 0.5*dz
       do k = 1, nz
	 p0(k)  = 100000.
	 pt0(k) = 300.*exp(1.3e-5*zz)
	 
         u00(k) = 0.
         v00(k) = 0.
         w0(k)  = 0.
	 k0(k)  = 0.
	 qt0(k) = 0.
	 zz=zz+dz
       enddo
