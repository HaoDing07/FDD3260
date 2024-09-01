      u00(:) = 20.
      v00(:) = 0.
      w0(:)  = 0.
      p0(:)  = 100000.
      pt0(:) = 300.
      qt0(:) = 0.
      ql0(:) = 0.
      qi0(:) = 0.
      k0(:)  = 1.e-16
          
      dsurf = psurf / ((cp_a-cv_a)*300.)

      zz = 0.5*dz
      pt0(1) = 300.*exp(1.e-4*zz/g)
      p0(1) = 100000.*((psurf/100000.)**((cp_a-cv_a)/cp_a) - g/(2.*cp_a)*(1./pt0(1) + 1./300.)*zz)**(cp_a/(cp_a-cv_a))
      den0(1) = - dsurf - 2./g*(p0(1) - psurf)/zz 
      
      zz = zz + dz
      do k = 2, nz
        pt0(k) = 300.*exp(1.e-4*zz/g)
        p0(k) = 100000.*((p0(k-1)/100000.)**((cp_a-cv_a)/cp_a) - g/(2.*cp_a)*(1./pt0(k) + 1./pt0(k-1))*dz)**(cp_a/(cp_a-cv_a))
        den0(k) = - den0(k-1) - 2./g*(p0(k) - p0(k-1))/dz 
        zz = zz + dz
      enddo

      avden0(:) = den0(:)
      scal0(:,:) = 0.
      
