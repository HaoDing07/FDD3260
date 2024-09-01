      zz = dz*fdz0(1)*0.5
      Pt0(1) = 263.42
      qv0(1) = 1.79625e-3
      u00(1) = -7.
      v00(1) = -2.
      w0(1) = 0.
      do k = 2, nz
!
!  qt0 is given in kg/kg
!  T0 is temperature
!
      if (zz < 400.) then
        Pt0(k) = 265. + 0.004*(zz - 400.) 
        qv0(k) = (1.5 - 0.00075*(zz - 400.)) * 1.e-3
      else if (zz >= 400. .and. zz < 825.) then
        Pt0(k) = 265.
        qv0(k) = 1.5e-3
      else if (zz >= 825. .and. zz < 2045.) then
        Pt0(k) = 266. + (zz - 825.)**0.3
        qv0(k) = 1.2e-3
      else if (zz >= 2045.) then
        Pt0(k) = 271. + (zz - 2000.)**0.33
        qv0(k) = (0.5 - 0.000075*(zz - 2045.)) * 1.e-3
      endif
!
      u00(k) = -7.
      v00(k) = -2. + 0.003*zz
      if (zz < 825.) then
        w0(k)  = -Ddiv * (zz + dz*fdz0(k)*0.5)
      else
        w0(k)  = w0(k-1)
      endif
!
      if (k <= kpert+1) then
        k0(k) = 0.1
      else
        k0(k) = 0.
      endif
!      
      if (k /= nz) zz = zz + 0.5*dz*(fdz0(k)+fdz0(k+1))      
      enddo
!
!  Surface values
!
      k0(1) = k0(2) 
