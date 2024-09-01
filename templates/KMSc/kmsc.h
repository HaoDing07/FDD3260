      u00(:) = 0.
      v00(:) = 0.
      w0(:)  = 0.
      k0(:)  = 0.
      ql0(:) = 0.
      qt0(:) = 9.e-3
      pt0(:) = 288.  
     
      ptsurf = 288. - flv00/cp_a*(100000./psurf)**((cp_a-cv_a)/cp_a)*qt0(1)
      dsurf = psurf / ((cp_a-cv_a)*ptsurf*(psurf/100000.)**((cp_a-cv_a)/cp_a))

      do l = 1, 5 
        p0(1) = 100000.*((psurf/100000.)**((cp_a-cv_a)/cp_a) - g/(2.*cp_a)*(1./pt0(1) + 1./ptsurf)*0.5*dz)**(cp_a/(cp_a-cv_a))
        den0(1) = - dsurf - 2./g*(p0(1) - psurf)/(0.5*dz)
        pt0(1) = 288. - flv00/cp_a*(100000./p0(1))**((cp_a-cv_a)/cp_a)*qt0(1) 
      enddo

      do k = 2, nz
        pt0(k) = pt0(k-1) 
        do l = 1, 5
          p0(k) = 100000.*((p0(k-1)/100000.)**((cp_a-cv_a)/cp_a) - g/(2.*cp_a)*(1./pt0(k) + 1./pt0(k-1))*dz)**(cp_a/(cp_a-cv_a))
          den0(k) = - den0(k-1) - 2./g*(p0(k) - p0(k-1))/dz 
          pt0(k) = 288. - flv00/cp_a*(100000./p0(k))**((cp_a-cv_a)/cp_a)*qt0(k) 
        enddo
      enddo


