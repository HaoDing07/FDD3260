       with_t=.true.

       u00 = 0.0
       v00 = 0.0
       w0  = 0.0

       do k = 1, nz
         if (z0(k) < 15000.) then
           qv0(k) = 0.01865*exp(-z0(k)/4000.)*exp(-(z0(k)/7500.)**2.)
	   t0(k) = (300.*(1. + 0.608*0.01865) - 0.0067*z0(k)) / (1. + 0.608*qv0(k))
         else
           qv0(k) = 1.e-14
	   t0(k) = (300.*(1. + 0.608*0.01865) - 0.0067*15000.) / (1. + 0.608*qv0(k))
         endif
         scal0(k,:) = 0.0
       enddo       