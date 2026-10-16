!-----------------------------------------------------------------------
! Massless tracers: dx/dt = u(x), with u interpolated from nekRS.
!-----------------------------------------------------------------------
      subroutine ppiclf_user_SetYdot
!
      implicit none
!
      include "PPICLF"
!
      integer*4 i
!
      do i=1,ppiclf_npart
         ppiclf_ydot(PPICLF_JX,i) = ppiclf_rprop(PPICLF_R_JUX,i)
         ppiclf_ydot(PPICLF_JY,i) = ppiclf_rprop(PPICLF_R_JUY,i)
         ppiclf_ydot(PPICLF_JZ,i) = ppiclf_rprop(PPICLF_R_JUZ,i)
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_user_MapProjPart(map,y,ydot,ydotc,rprop)
!
      implicit none
!
      real*8 map(*), y(*), ydot(*), ydotc(*), rprop(*)
!
      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_user_EvalNearestNeighbor
     >                                        (i,j,yi,rpropi,yj,rpropj)
!
      implicit none
!
      integer*4 i, j
      real*8 yi(*), rpropi(*), yj(*), rpropj(*)
!
      return
      end
!-----------------------------------------------------------------------
