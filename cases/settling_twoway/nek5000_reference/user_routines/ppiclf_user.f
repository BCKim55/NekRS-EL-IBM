!-----------------------------------------------------------------------
! settling_twoway on the original ppiclF + Nek5000: Schiller-Naumann drag,
! gravity g = -1 in z, fluid rho = 1, mu = 0.01; the force on the fluid is
! the drag reaction, projected with the particle volume.
!-----------------------------------------------------------------------
      subroutine ppiclf_user_SetYdot
      implicit none
      include "PPICLF"
      real*8 rpi, rmu, rhof, rmass, dp, urx, ury, urz, rep, fd
      integer*4 i
      rpi  = 4.0d0*atan(1.0d0)
      rmu  = 0.01d0
      rhof = 1.0d0
      do i=1,ppiclf_npart
         rmass = ppiclf_rprop(PPICLF_R_JVOLP,i)
     >          *ppiclf_rprop(PPICLF_R_JRHOP,i)
         dp    = ppiclf_rprop(PPICLF_R_JDP,i)
         urx = ppiclf_rprop(PPICLF_R_JUX,i) - ppiclf_y(PPICLF_JVX,i)
         ury = ppiclf_rprop(PPICLF_R_JUY,i) - ppiclf_y(PPICLF_JVY,i)
         urz = ppiclf_rprop(PPICLF_R_JUZ,i) - ppiclf_y(PPICLF_JVZ,i)
         rep = sqrt(urx**2 + ury**2 + urz**2)*dp*rhof/rmu
         fd  = 3.0d0*rpi*rmu*dp*(1.0d0 + 0.15d0*rep**0.687d0)
         ppiclf_ydot(PPICLF_JX ,i) = ppiclf_y(PPICLF_JVX,i)
         ppiclf_ydot(PPICLF_JY ,i) = ppiclf_y(PPICLF_JVY,i)
         ppiclf_ydot(PPICLF_JZ ,i) = ppiclf_y(PPICLF_JVZ,i)
         ppiclf_ydot(PPICLF_JVX,i) = fd*urx/rmass
         ppiclf_ydot(PPICLF_JVY,i) = fd*ury/rmass
         ppiclf_ydot(PPICLF_JVZ,i) = fd*urz/rmass - 1.0d0
         ppiclf_ydotc(PPICLF_JVX,i) = -fd*urx
         ppiclf_ydotc(PPICLF_JVY,i) = -fd*ury
         ppiclf_ydotc(PPICLF_JVZ,i) = -fd*urz
      enddo
      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_user_MapProjPart(map,y,ydot,ydotc,rprop)
      implicit none
      include "PPICLF"
      real*8 map(PPICLF_LRP_PRO), y(PPICLF_LRS), ydot(PPICLF_LRS),
     >       ydotc(PPICLF_LRS), rprop(PPICLF_LRP)
      map(PPICLF_P_JFX)   = ydotc(PPICLF_JVX)
      map(PPICLF_P_JFY)   = ydotc(PPICLF_JVY)
      map(PPICLF_P_JFZ)   = ydotc(PPICLF_JVZ)
      map(PPICLF_P_JPHIP) = rprop(PPICLF_R_JVOLP)
      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_user_EvalNearestNeighbor
     >                                        (i,j,yi,rpropi,yj,rpropj)
      implicit none
      integer*4 i, j
      real*8 yi(*), rpropi(*), yj(*), rpropj(*)
      return
      end
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine ppiclf_user_budget(pp,sf,rm)
! rank-local particle z momentum, drag reaction sum and mass, so that the
! .usr need not include PPICLF (duplicate DATA-initialised commons on macOS)
      implicit none
      include "PPICLF"
      real*8 pp, sf, rm
      integer*4 i
      pp = 0.0d0
      sf = 0.0d0
      rm = 0.0d0
      do i=1,ppiclf_npart
         rm = rm + ppiclf_rprop(PPICLF_R_JRHOP,i)
     >            *ppiclf_rprop(PPICLF_R_JVOLP,i)
         pp = pp + ppiclf_rprop(PPICLF_R_JRHOP,i)
     >            *ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_y(PPICLF_JVZ,i)
         sf = sf + ppiclf_ydotc(PPICLF_JVZ,i)
      enddo
      return
      end
!-----------------------------------------------------------------------
