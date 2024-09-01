      zz = dz*fdz0(2)*0.5
      p0(1) = (psurf/100000.)**(Ra/cp)
      Pt0(1) = 264.7*(psurf/100000.)**(Ra/cp)
      do k = 2, nz
      p0(k) = (psurf/100000.)**(Ra/cp) - g/(cp*Pt0(1)) * zz
!
!  qt0 is given in kg/kg
!  T0 is temperature
!
      if (zz < 540.) then
        Pt0(k) = (261.5 - 0.005926*(zz - 540.)) * p0(k)
        qt0(k) = (1.35 - 0.000741*(zz - 540.)) * 1.e-3
      else if (zz >= 540. .and. zz < 1100.) then
        Pt0(k) = (256.8 - 0.008393*(zz - 1100.)) * p0(k)
        qt0(k) = 1.35e-3
      else if (zz >= 1100. .and. zz < 1150.) then
        Pt0(k) = (261. + 0.086*(zz - 1150.)) * p0(k)
        qt0(k) = 1.35e-3
      else if (zz >= 1150.) then
        Pt0(k) = (258.7 - 0.004364*(zz - 1700.)) * p0(k)
        qt0(k) = (1.05 - 0.000909*(zz - 1700.)) * 1.e-3
      endif
!
      u00(k) = -11. - 0.003529*(zz - 1700.)
      v00(k) = 3.
      if (zz < 1150.) then
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
!
!  Surface values
!
      u00(1) = u00(2)
      v00(1) = v00(2)
      pt0(1) = pt0(2)
      w0(1) = 0.
      qt0(1) = qt0(2)
      k0(1) = k0(2) 
