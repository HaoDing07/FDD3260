      
      ! --- 090295
      !  rktable1 is the main table storing photochemical
      !       reaction rate or associated formulae as
      !       FUNCTION OF TEMPERATURE with interval of 0.5
      !       degree STARTED FROM 200.5 K to 300 K (300 elements)
      !       for the nract-th reaction
      !
      !  rktable2 is a additional table storing rk13, 15, & 20's
      !       rt00 formulae indexed as:
      !       1 for rk(13)
      !       2 for rk(15)
      !       3 for rk(20)
      !
      ! --- 042596
      !  rktable3 includes photorates parameter for rk(1),rk(8),
      !      rk(17), rk(24), and rk(26), as well as cos znith angle:
      !      zangle = 0, 10, 20, 30, 40, 50, 60, 70, 78, 86
      !      
      real, dimension(34,300)  :: rktable1
      real, dimension( 3,300)  :: rktable2
      ! --- rktable 3
      real, dimension(10) :: cosza4rk = (/                            &
           0.100000E+01,  0.984808E+00,  0.939693E+00,  0.866025E+00,   &
           0.766044E+00,  0.642788E+00,  0.500000E+00,  0.342020E+00,   &
           0.207912E+00,  0.697565E-01/)

      real, dimension(10) :: rk08gama = (/                            &
           0.958260E-05,  0.956540E-05,  0.947786E-05,  0.933214E-05,   &
           0.904837E-05,  0.853986E-05,  0.790570E-05,  0.644646E-05,   &
           0.540693E-05,  0.640000E-05/)

      real, dimension(10) :: rk08aaa = (/                             &
          -0.213357E-08, -0.212833E-08, -0.210322E-08, -0.206239E-08,   &
          -0.198520E-08, -0.185147E-08, -0.169861E-08, -0.133796E-08,   &
          -0.109912E-08, -0.131294E-08/)

      real, dimension(10) :: rk01table = (/                           &
           0.553973E-07,  0.541612E-07,  0.501300E-07,  0.436050E-07,   &
           0.346763E-07,  0.241176E-07,  0.141663E-07,  0.593145E-08,   &
           0.238321E-08,  0.148548E-08/)

      real, dimension(10) :: rk17table = (/                           &
           0.761430E-08,  0.755268E-08,  0.730677E-08,  0.691454E-08,   &
           0.626346E-08,  0.536943E-08,  0.430112E-08,  0.288139E-08,   &
           0.199044E-08,  0.209279E-08/)

      real, dimension(10) :: rk24table = (/                           &
           0.535338E-08,  0.531281E-08,  0.515418E-08,  0.489785E-08,   &
           0.447421E-08,  0.387878E-08,  0.318768E-08,  0.217960E-08,   &
           0.155503E-08,  0.170283E-08/)

      real, dimension(10) :: rk26table = (/                           &
           0.109106E-06,  0.108473E-06,  0.105771E-06,  0.101587E-06,   &
           0.941011E-07,  0.831775E-07,  0.701839E-07,  0.494384E-07,   &
           0.364853E-07,  0.410051E-07/)!
