      zz = dz*fdz0(2)*0.5
      p0(1) = (psurf/100000.)**(Ra/cp)
      Pt0(1) = 264.5*(psurf/100000.)**(Ra/cp)
      do k = 2, nz
      p0(k) = (psurf/100000.)**(Ra/cp) - g/(cp*Pt0(1)) * zz
!
!  qt0 is given in kg/kg
!  T0 is temperature
!
      if (zz < 235.) then
        Pt0(k) = (261.9 - 0.011064*(zz - 235.)) * p0(k)
        qt0(k) = (1.642 - 0.000647*(zz - 235.)) * 1.e-3
      else if (zz >= 235. .and. zz < 430.) then
        Pt0(k) = (261.2 - 0.00359*(zz - 430.)) * p0(k)
        qt0(k) = 1.642e-3
      else if (zz >= 430. .and. zz < 465.) then
        Pt0(k) = (267.5 + 0.0315*(zz - 630.)) * p0(k)
        qt0(k) = 1.642e-3
      else if (zz >= 465. .and. zz < 630.) then
        Pt0(k) = (267.5 + 0.0315*(zz - 630.)) * p0(k)
        qt0(k) = (1.433 - 0.000173*(zz - 1670.)) * 1.e-3
      else if (zz >= 630.) then
        Pt0(k) = (263.15 - 0.004182*(zz - 1670.)) * p0(k)
        qt0(k) = (1.433 - 0.000173*(zz - 1670.)) * 1.e-3
      endif
!
      u00(k) =  0.
      v00(k) = -0.01144*zz
      if (zz < 630.) then
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
      if (k/=nz) zz = zz + 0.5*dz*(fdz0(k)+fdz0(k+1))      
      enddo
