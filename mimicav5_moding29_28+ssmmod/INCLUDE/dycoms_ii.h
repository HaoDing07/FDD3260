       zz = 0.5*dz*fdz0(1)
       w0(1) = 0.
       do k = 1, nz
        u00(k) = (3. + 4.3*zz/1000.) 
        v00(k) = (-9. + 5.6*zz/1000.) 
!      u00(k) = (1.5 + 2.*zz/1000.) 
!      v00(k) = (-3.5 + 5.*zz/1000.) 
!      u00(k) = (6.5 + 2.*zz/1000.) 
!      v00(k) = (-8.5 + 5.*zz/1000.) 
!
!  qv0 is given in kg/kg
!  pt0 is the liquid potential temperature
!
        if (zz < 795.) then
         pt0(k) = 288.3
         qv0(k) = 9.45e-3
        else
         pt0(k) = 295. + (zz - 795.)**(1./3.)
         qv0(k) = (5. - 3.*(1.-exp((795. - zz)/500.))) * 1.e-3
        endif

        if (zz>1200. .AND. zz<2000.) then
        qv0(k) = 6.0e-3  
        endif          
!
        w0(k)  = -Ddiv*zz
        if (k <= kpert+1) then
         k0(k) = 1.
        else
         k0(k) = 0.
        endif
!      
        if (k/=nz) zz = zz + 0.5*dz*(fdz0(k)+fdz0(k+1))      
       enddo
