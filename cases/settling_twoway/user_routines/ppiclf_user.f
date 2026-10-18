!-----------------------------------------------------------------------
! Point particles under gravity with Schiller-Naumann drag, two-way
! coupled. The force each particle exerts on the fluid is kept in
! rprop(PPICLF_R_JF*) and projected with the particle volume.
!-----------------------------------------------------------------------
      subroutine ppiclc_user_SetParams(rmu_in,rhof_in,grav_in)
     > bind(C, name="ppiclc_user_SetParams")
!
      implicit none
!
      real*8 rmu_in, rhof_in, grav_in(3)
!
      real*8 rmu, rhof, grav(3)
      common /ppiclf_user_fluid/ rmu, rhof, grav
!
      rmu     = rmu_in
      rhof    = rhof_in
      grav(1) = grav_in(1)
      grav(2) = grav_in(2)
      grav(3) = grav_in(3)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_user_SetYdot
!
      implicit none
!
      include "PPICLF"
!
      real*8 rmu, rhof, grav(3)
      common /ppiclf_user_fluid/ rmu, rhof, grav
!
      real*8 rpi, rmass, dp, urx, ury, urz, rep, fd
      integer*4 i
!
      rpi = 4.0d0*atan(1.0d0)

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
         ppiclf_ydot(PPICLF_JVX,i) = fd*urx/rmass + grav(1)
         ppiclf_ydot(PPICLF_JVY,i) = fd*ury/rmass + grav(2)
         ppiclf_ydot(PPICLF_JVZ,i) = fd*urz/rmass + grav(3)

         ! reaction of the drag on the fluid
         ppiclf_rprop(PPICLF_R_JFX,i) = -fd*urx
         ppiclf_rprop(PPICLF_R_JFY,i) = -fd*ury
         ppiclf_rprop(PPICLF_R_JFZ,i) = -fd*urz
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_user_MapProjPart(map,y,ydot,ydotc,rprop)
!
      implicit none
!
      include "PPICLF"
!
      real*8 map  (PPICLF_LRP_PRO)
      real*8 y    (PPICLF_LRS)
      real*8 ydot (PPICLF_LRS)
      real*8 ydotc(PPICLF_LRS)
      real*8 rprop(PPICLF_LRP)
!
      map(PPICLF_P_JFX)   = rprop(PPICLF_R_JFX)
      map(PPICLF_P_JFY)   = rprop(PPICLF_R_JFY)
      map(PPICLF_P_JFZ)   = rprop(PPICLF_R_JFZ)
      map(PPICLF_P_JPHIP) = rprop(PPICLF_R_JVOLP)

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
