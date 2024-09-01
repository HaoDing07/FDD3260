      zz = dz*fdz0(2)*0.5
      p0(1) = (psurf/100000.)**(Ra/cp)
      Pt0(1) = 264.75*(psurf/100000.)**(Ra/cp)
      do k = 2, nz
      p0(k) = (psurf/100000.)**(Ra/cp) - g/(cp*Pt0(1)) * zz
!
!  qt0 is given in kg/kg
!  T0 is temperature
!
      if (zz < 292.) then
        Pt0(k) = (260.15 - 0.00729*(zz - 631.)) * p0(k)
        qt0(k) = (1.78 - 0.000246*(zz - 292.)) * 1.e-3
      else if (zz >= 292. .and. zz < 631.) then
        Pt0(k) = (260.15 - 0.00729*(zz - 631.)) * p0(k)
        qt0(k) = qt0(k-1)
      else if (zz >= 631. .and. zz < 755.) then
        Pt0(k) = (261.85 + 0.01371*(zz - 755.)) * p0(k)
        qt0(k) = qt0(k-1)
      else if (zz >= 755.) then
        Pt0(k) = (260.45 - 0.005833*(zz - 995.)) * p0(k)
        qt0(k) = (1.42 - 0.00075*(zz - 995.)) * 1.e-3
      endif
!
      u00(k) = 0.0013*zz
      v00(k) = 0.
      if (zz < 755.) then
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
