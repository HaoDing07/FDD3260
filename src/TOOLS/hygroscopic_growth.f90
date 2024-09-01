#include "ctrparam.h"

module hygro_growth

 !General modules

    USE gridno 
    USE shared_data    
    USE shared_all
    USE shared_state
    USE averages
    USE gradients
  
    IMPLICIT NONE

    public :: aerosol_swelling
    
    contains

    !  ===================================================
     subroutine aerosol_swelling
    !  ===================================================
    !
      implicit none
      integer :: i, j, k
      if (verbose > 0) write(7,*) 'Starting subroutine aerosol_swelling'      

      if (with_aero_swelling) then
        write(7,*) 'Before aerosol swelling', 'qv(k=50) = ', thermo%qv(ip_start,jp_start,50)
        do k = 1,100
          do j = jp_start, jp_end
            do i = ip_start, ip_end
              ! thermo%qv(i,j,k) = thermo%qv(i,j,k)*0.5
              thermo%qv(i,j,k) = thermo%qv(i,j,k) - (aeroi(1)%init%n0-hydromtr(drop)%n(i,j,k))*rhow*pi/6.0*(aeroi(1)%size%rmean*2.0*10.0**(-6.0))**3.0*(hgf**3.0-1.0)
             
            enddo
          enddo
        enddo
        write(7,*) 'With aerosol swelling', 'qv(k=50) = ' , thermo%qv(ip_start,jp_start,50)
        write(7,*) 'With aerosol swelling', 'Na = ' , aeroi(1)%init%n0*10.0**(-6.0), hydromtr(drop)%n(ip_start,jp_start,50)*10.0**(-6.0)
        write(7,*) 'With aerosol swelling', 'Rg = ' , (aeroi(1)%size%rmean*2.0*10.0**(-6.0)), hgf
        write(7,*) 'aerosol_sewlling is working at this timestep'
        
      endif 
    

    return
    end

end module hygro_growth