subroutine pwl_basis_1d ( nd, xd, ni, xi, bk )

!*****************************************************************************80
!
!! PWL_BASIS_1D evaluates a 1D piecewise linear basis function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!
!    Input, real ( kind = 8 ) XD(ND), the data points.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(NI), the interpolation points.
!
!    Output, real ( kind = 8 ) BK(NI,ND), the basis functions at the 
!    interpolation points.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  real ( kind = 8 ) bk(ni,nd)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) t
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi(ni)

  bk(1:ni,1:nd) = 0.0D+00

  if ( nd == 1 ) then
    bk(1:ni,1:nd) = 1.0D+00
    return
  end if

  do i = 1, ni

    do j = 1, nd

      if ( j == 1 .and. xi(i) <= xd(j) ) then

        t = ( xi(i) - xd(j) ) / ( xd(j+1) - xd(j) )
        bk(i,j) = 1.0D+00 - t

      else if ( j == nd .and. xd(j) <= xi(i) ) then

        t = ( xi(i) - xd(j-1) ) / ( xd(j) - xd(j-1) )
        bk(i,j) = t

      else if ( xd(j-1) < xi(i) .and. xi(i) <= xd(j) ) then

        t = ( xi(i) - xd(j-1) ) / ( xd(j) - xd(j-1) )
        bk(i,j) = t

      else if ( xd(j) <= xi(i) .and. xi(i) < xd(j+1) ) then

        t = ( xi(i) - xd(j) ) / ( xd(j+1) - xd(j) )
        bk(i,j) = 1.0D+00 - t

      end if

    end do

  end do

  return
end

subroutine pwl_interp_2d ( nxd, nyd, xd, yd, zd, ni, xi, yi, zi )

!*****************************************************************************80
!
!! PWL_INTERP_2D: piecewise linear interpolant to data defined on a 2D grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NXD, NYD, the number of X and Y data values.
!
!    Input, real ( kind = 8 ) XD(NXD), YD(NYD), the sorted X and Y data.
!
!    Input, real ( kind = 8 ) ZD(NXD,NYD), the Z data.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(NI), YI(NI), the coordinates of the
!    interpolation points.
!
!    Output, real ( kind = 8 ) ZI(NI), the value of the interpolant.
!
  implicit none

  integer ni
  integer nxd
  integer nyd

  integer i
  integer j
  integer k
  integer r8vec_bracket5
  
  real alpha
  real beta
  real det
  real dxa
  real dxb
  real dxi
  real dya
  real dyb
  real dyi
  real gamma
  real xd(nxd)
  real xi(ni)
  real yd(nyd)
  real yi(ni)
  real zd(nxd,nyd)
  real zi(ni)

  do k = 1, ni

    i = r8vec_bracket5 ( nxd, xd, xi(k) )
    if ( i == -1 ) then
      zi(k) = huge (alpha)
      cycle
    end if

    j = r8vec_bracket5 ( nyd, yd, yi(k) )
    if ( j == -1 ) then
      zi(k) = huge (alpha)
      cycle
    end if

    if ( yi(k) < yd(j+1) &
      + ( yd(j) - yd(j+1) ) * ( xi(k) - xd(i) ) / ( xd(i+1) - xd(i) ) ) then

      dxa = xd(i+1) - xd(i)
      dya = yd(j)   - yd(j)

      dxb = xd(i)   - xd(i)
      dyb = yd(j+1) - yd(j)

      dxi = xi(k)   - xd(i)
      dyi = yi(k)   - yd(j)

      det = dxa * dyb - dya * dxb

      alpha = ( dxi * dyb - dyi * dxb ) / det
      beta =  ( dxa * dyi - dya * dxi ) / det
      gamma = 1.0D+00 - alpha - beta

      zi(k) = alpha * zd(i+1,j) + beta * zd(i,j+1) + gamma * zd(i,j)

    else

      dxa = xd(i)   - xd(i+1)
      dya = yd(j+1) - yd(j+1)

      dxb = xd(i+1) - xd(i+1)
      dyb = yd(j)   - yd(j+1)

      dxi = xi(k)   - xd(i+1)
      dyi = yi(k)   - yd(j+1)

      det = dxa * dyb - dya * dxb

      alpha = ( dxi * dyb - dyi * dxb ) / det
      beta =  ( dxa * dyi - dya * dxi ) / det
      gamma = 1.0D+00 - alpha - beta

      zi(k) = alpha * zd(i,j+1) + beta * zd(i+1,j) + gamma * zd(i+1,j+1)

    end if

  end do

  return
end

subroutine pwl_interp_1d ( nd, xd, yd, ni, xi, yi )

!*****************************************************************************80
!
!! PWL_VALUE_1D evaluates the piecewise linear interpolant.
!
!  Discussion:
!
!    The piecewise linear interpolant L(ND,XD,YD)(X) is the piecewise
!    linear function which interpolates the data (XD(I),YD(I)) for I = 1
!    to ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!    ND must be at least 1.
!
!    Input, real ( kind = 8 ) XD(ND), the data points.
!
!    Input, real ( kind = 8 ) YD(ND), the data values.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(NI), the interpolation points.
!
!    Output, real ( kind = 8 ) YI(NI), the interpolated values.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ) t
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) yi(ni)

  yi(1:ni) = 0.0D+00

  if ( nd == 1 ) then
    yi(1:ni) = yd(1)
    return
  end if

  do i = 1, ni

    if ( xi(i) <= xd(1) ) then

      t = ( xi(i) - xd(1) ) / ( xd(2) - xd(1) )
      yi(i) = ( 1.0D+00 - t ) * yd(1) + t * yd(2)

    else if ( xd(nd) <= xi(i) ) then

      t = ( xi(i) - xd(nd-1) ) / ( xd(nd) - xd(nd-1) )
      yi(i) = ( 1.0D+00 - t ) * yd(nd-1) + t * yd(nd)

    else

      do k = 2, nd

        if ( xd(k-1) <= xi(i) .and. xi(i) <= xd(k) ) then

          t = ( xi(i) - xd(k-1) ) / ( xd(k) - xd(k-1) )
          yi(i) = ( 1.0D+00 - t ) * yd(k-1) + t * yd(k)
          exit

        end if

      end do

    end if

  end do
  
  return
end

function r8vec_bracket5 ( nd, xd, xi )

!*****************************************************************************80
!
!! R8VEC_BRACKET5 brackets data between successive entries of a sorted R8VEC.
!
!  Discussion:
!
!    We assume XD is sorted.
!
!    If XI is contained in the interval [XD(1),XD(N)], then the returned 
!    value B indicates that XI is contained in [ XD(B), XD(B+1) ].
!
!    If XI is not contained in the interval [XD(1),XD(N)], then B = -1.
!
!    This code implements a version of binary search which is perhaps more
!    understandable than the usual ones.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data values.
!
!    Input, real ( kind = 8 ) XD(N), the sorted data.
!
!    Input, real ( kind = 8 ) XD, the query value.
!
!    Output, integer ( kind = 4 ) R8VEC_BRACKET5, the bracket information.
!
  implicit none

  integer nd

  integer r8vec_bracket5
  integer b
  integer l
  integer m
  integer r
  
  real xd(nd)
  real xi

  if ( xi < xd(1) .or. xd(nd) < xi ) then

    b = -1

  else

    l = 1
    r = nd

    do while ( l + 1 < r )
      m = ( l + r ) / 2
      if ( xi < xd(m) ) then
        r = m
      else
        l = m
      end if
    end do

    b = l

  end if

  r8vec_bracket5 = b

  return
end
