       zz = 0.5*dz*fdz0(2)
       do k = 2, nz
        u00(k) = 0.
        v00(k) = 0.
        w0(k)  = 0.
!
!  qt0 is given in kg/kg
!  pt0 is the liquid potential temperature
!
      if (zz < 687.5) then
        pt0(k) = 288.
        qt0(k) = 0.
      else if (zz >= 687.5 .and. zz <= 712.5) then
        pt0(k) = 288. + 0.28*(zz - 687.5)
      qt0(k) = 0.
      else
        pt0(k) = 295. + 0.0001*(zz - 712.5)
        qt0(k) = 0.
      endif
!      
      if (k/=nz) zz = zz + 0.5*dz*(fdz0(k)+fdz0(k+1))      
      enddo
      w0(nwz) = 0.
!
!  Surface values
!
      u00(1) = u00(2)
      v00(1) = v00(2)
      pt0(1) = pt0(2)
      w0(1)  = 0.
      qt0(1) = qt0(2)
      k0     = 0.
