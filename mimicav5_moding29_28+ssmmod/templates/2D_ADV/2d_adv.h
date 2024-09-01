      u00(:) = 1.
      v00(:) = 0.
      w0(:)  = 1.
      p0(:)  = 100000.
      pt0(:) = 300.
      qt0(:) = 0.
      k0(:)  = 1.e-16
      
      den0(1:nz) = p0(1:nz) / (Ra*Pt0(:)*(p0(1:nz)/1.E+05)**(Ra/cp))
      den0(nwz) = p0(nwz) / (Ra*Pt0(nz)*(p0(nwz)/1.E+05)**(Ra/cp))
      avden0(:) = den0(:)
      scal0(:,:) = 0.
      
