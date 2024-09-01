       with_t=.true.
       with_rh=.true.

       u00 = 0.0
       v00 = 0.0
       w0  = 0.0
       k0  = 0.01
   
       rh00 = 0.85
       rh01 = 0.65
       rh02 = rh00

       zz = 0.5*dz*fdz0(1)
       do k = 1, nz
         if (zz < 10000.) then
           rh0(k) = rh00 - (rh00 - rh01)*(zz/10000.)**(1.2)
         else if (z0(k) < 15000.) then
           rh0(k) = rh01 - (rh01 - rh02)*((zz-10000.)/5000.)**(3.)
         else
           rh0(k) = rh02*exp(-(zz-15000.)/1500.)
         endif

         if (zz < 15000.) then
	   t0(k) = (300. - 0.0067*zz) 
         else
	   t0(k) = (300. - 0.0067*15000.) 
         endif
!
         if (k < nz) zz = zz + 0.5*dz*(fdz0(k)+fdz0(k+1))	   
       enddo       
