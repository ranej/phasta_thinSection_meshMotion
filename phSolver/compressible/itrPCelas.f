c-----------------------------------------------------------------------
c
c    Predict solution at time n+1
c
c-----------------------------------------------------------------------
      subroutine itrPredictElas (disp,umesh,dt)
      
      include "common.h"
      
      real*8   disp(nshg,nelas),umesh(numnp,nsd),dt
c
      disp = zero
c      disp = umesh*dt
c
      return
      end

c-----------------------------------------------------------------------
c
c    Correct solution at time n+1
c
c-----------------------------------------------------------------------
      subroutine itrCorrectElas (disp, Dy)
      
      include "common.h"
      
      real*8   disp(nshg,nelas), Dy(nshg,nelas)
c      
      disp = disp + Dy
c
      return
      end

c-----------------------------------------------------------------------
c
c    Update solution at end of time step
c
c-----------------------------------------------------------------------
      subroutine itrUpdateElas (xold, x)
      
      include "common.h"
      
      real*8   xold(numnp,nsd), x(numnp,nsd)
c      
      xold    = x
c
      return
      end

