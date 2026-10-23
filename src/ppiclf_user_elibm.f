!-----------------------------------------------------------------------
! ppiclF user routines for coupled EL point particles and IBM rigid
! bodies (layout in PPICLF_ELIBM.h, body state in ELIBM.inc).
!
! EL particle : Schiller-Naumann drag, gravity, soft-sphere contact
!               with other EL particles (ppiclF neighbour search) and
!               with IBM bodies (direct loop over the few bodies).
! IBM marker  : moves with its body; explicit direct forcing as in
!               NekIBM (ibm_f_marker), f = rho_f dV (u_body - u)/dt,
!               exerted on the fluid. The body integrates the reaction
!               in C++ (src/elibm.hpp).
! Contact     : linear spring-dashpot, normal only; the dashpot gives
!               restitution coefficient e for an isolated binary hit.
!-----------------------------------------------------------------------
      subroutine ppiclc_elibm_SetFluid(rmu,rhof,grav,dtf)
     > bind(C, name="ppiclc_elibm_SetFluid")
!
      implicit none
!
      include "PPICLF"
      include "ELIBM.inc"
!
      real*8 rmu, rhof, grav(3), dtf
!
      elibm_mu      = rmu
      elibm_rhof    = rhof
      elibm_grav(1) = grav(1)
      elibm_grav(2) = grav(2)
      elibm_grav(3) = grav(3)
      elibm_dtf     = dtf

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclc_elibm_SetContact(rkn,erest,ifnb)
     > bind(C, name="ppiclc_elibm_SetContact")
!
! ifnb = 1 enables EL-EL contact (needs ppiclc_solve_InitNeighborBin)
!
      implicit none
!
      include "PPICLF"
      include "ELIBM.inc"
!
      real*8 rkn, erest
      integer*4 ifnb
!
      elibm_kn    = rkn
      elibm_erest = erest
      elibm_ifnb  = ifnb

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclc_elibm_SetBodies(nb,bx,bv,bw,brot,brad,bmass)
     > bind(C, name="ppiclc_elibm_SetBodies")
!
! Rigid-body state, one entry per body: centre, velocity, angular
! velocity, rotation matrix (body -> lab, column-major), radius, mass.
!
      implicit none
!
      include "PPICLF"
      include "ELIBM.inc"
!
      integer*4 nb
      real*8 bx(3,*), bv(3,*), bw(3,*), brot(3,3,*), brad(*), bmass(*)
!
      integer*4 b, i, j
!
      if (nb .gt. elibm_nbmax)
     >   call ppiclf_exittr('Increase elibm_nbmax$',0.0d0,nb)

      elibm_nbody = nb
      do b=1,nb
         do i=1,3
            elibm_bx(i,b) = bx(i,b)
            elibm_bv(i,b) = bv(i,b)
            elibm_bw(i,b) = bw(i,b)
            do j=1,3
               elibm_brot(i,j,b) = brot(i,j,b)
            enddo
         enddo
         elibm_brad(b)  = brad(b)
         elibm_bmass(b) = bmass(b)
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine elibm_MinImage(d)
!
! Shortest periodic image of the separation vector d.
!
      implicit none
!
      include "PPICLF"
!
      real*8 d(3)
!
      real*8 rlen
      integer*4 k
!
      do k=1,ppiclf_ndim
         if (ppiclf_iperiodic(k) .eq. 0) then
            rlen = ppiclf_xdrange(2,k) - ppiclf_xdrange(1,k)
            d(k) = d(k) - rlen*anint(d(k)/rlen)
         endif
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclc_elibm_BeginStep()
     > bind(C, name="ppiclc_elibm_BeginStep")
!
! Zero the EL-IBM contact impulses and the RK stage counter; call once
! per fluid step before the particle sub-steps.
!
      implicit none
!
      include "PPICLF"
      include "ELIBM.inc"
!
      integer*4 i, k
!
      elibm_istage = 0
      do i=1,ppiclf_npart
         do k=1,3
            ppiclf_rprop(PPICLF_R_JIX-1+k,i) = 0.0d0
         enddo
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclc_elibm_PlaceMarkers()
     > bind(C, name="ppiclc_elibm_PlaceMarkers")
!
! Put every local marker at its rigid position for the current body
! state and give it the rigid-body velocity there.
!
      implicit none
!
      include "PPICLF"
      include "ELIBM.inc"
!
      integer*4 i, b, k
      real*8 r(3), rlen
!
      do i=1,ppiclf_npart
         if (nint(ppiclf_rprop(PPICLF_R_JTYPE,i)) .ne.
     >       PPICLF_TYPE_MARKER) cycle
         b = nint(ppiclf_rprop(PPICLF_R_JID,i)) + 1
         do k=1,3
            r(k) = elibm_brot(k,1,b)*ppiclf_rprop(PPICLF_R_JRX,i)
     >           + elibm_brot(k,2,b)*ppiclf_rprop(PPICLF_R_JRY,i)
     >           + elibm_brot(k,3,b)*ppiclf_rprop(PPICLF_R_JRZ,i)
            ppiclf_y(PPICLF_JX-1+k,i) = elibm_bx(k,b) + r(k)
            if (ppiclf_iperiodic(k) .eq. 0) then
               rlen = ppiclf_xdrange(2,k) - ppiclf_xdrange(1,k)
               ppiclf_y(PPICLF_JX-1+k,i) = ppiclf_xdrange(1,k) +
     >            modulo(ppiclf_y(PPICLF_JX-1+k,i)
     >                   - ppiclf_xdrange(1,k), rlen)
            endif
         enddo
         ppiclf_y(PPICLF_JVX,i) = elibm_bv(1,b)
     >      + elibm_bw(2,b)*r(3) - elibm_bw(3,b)*r(2)
         ppiclf_y(PPICLF_JVY,i) = elibm_bv(2,b)
     >      + elibm_bw(3,b)*r(1) - elibm_bw(1,b)*r(3)
         ppiclf_y(PPICLF_JVZ,i) = elibm_bv(3,b)
     >      + elibm_bw(1,b)*r(2) - elibm_bw(2,b)*r(1)
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine elibm_ContactForce(xi,vi,ri,rmi,xj,vj,rj,rmj,f)
!
! Normal spring-dashpot force on i from j (zero without overlap).
!
      implicit none
!
      include "PPICLF"
      include "ELIBM.inc"
!
      real*8 xi(3), vi(3), ri, rmi, xj(3), vj(3), rj, rmj, f(3)
!
      real*8 d(3), dist, delta, rn(3), vn, rmeff, eta, rpi2, rle
      integer*4 k
!
      f(1) = 0.0d0
      f(2) = 0.0d0
      f(3) = 0.0d0

      do k=1,3
         d(k) = xi(k) - xj(k)
      enddo
      call elibm_MinImage(d)
      dist  = sqrt(d(1)**2 + d(2)**2 + d(3)**2)
      delta = ri + rj - dist
      if (delta .le. 0.0d0 .or. dist .le. 0.0d0) return

      vn = 0.0d0
      do k=1,3
         rn(k) = d(k)/dist
         vn    = vn + (vi(k) - vj(k))*rn(k)
      enddo

      rpi2  = 9.869604401089358d0
      rle   = log(elibm_erest)
      rmeff = 1.0d0/(1.0d0/rmi + 1.0d0/rmj)
      eta   = -2.0d0*rle*sqrt(rmeff*elibm_kn)/sqrt(rpi2 + rle**2)

      do k=1,3
         f(k) = (elibm_kn*delta - eta*vn)*rn(k)
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_user_SetYdot
!
      implicit none
!
      include "PPICLF"
      include "ELIBM.inc"
!
      real*8 rpi, rmass, dp, ur(3), rep, fd, fb(3), fc(3), r(3),
     >       ud(3), vol, wimp, rnp
      integer*4 i, b, k
!
      rpi = 4.0d0*atan(1.0d0)

      ! weight of this RK3 stage's RHS in the step (ppiclF's SSP RK3:
      ! stages contribute dt/6, dt/6, 2dt/3), used to integrate the
      ! contact impulse the IBM bodies receive
      elibm_istage = elibm_istage + 1
      wimp = ppiclf_dt/6.0d0
      if (mod(elibm_istage,3) .eq. 0) wimp = 2.0d0*ppiclf_dt/3.0d0

      do i=1,ppiclf_npart

         if (nint(ppiclf_rprop(PPICLF_R_JTYPE,i)) .eq.
     >       PPICLF_TYPE_MARKER) then
            b = nint(ppiclf_rprop(PPICLF_R_JID,i)) + 1
            do k=1,3
               r(k) = ppiclf_y(PPICLF_JX-1+k,i) - elibm_bx(k,b)
            enddo
            call elibm_MinImage(r)
            ud(1) = elibm_bv(1,b) + elibm_bw(2,b)*r(3)
     >                            - elibm_bw(3,b)*r(2)
            ud(2) = elibm_bv(2,b) + elibm_bw(3,b)*r(1)
     >                            - elibm_bw(1,b)*r(3)
            ud(3) = elibm_bv(3,b) + elibm_bw(1,b)*r(2)
     >                            - elibm_bw(2,b)*r(1)
            vol = ppiclf_rprop(PPICLF_R_JVOLP,i)
            do k=1,3
               ppiclf_ydot(PPICLF_JX -1+k,i) = ud(k)
               ppiclf_ydot(PPICLF_JVX-1+k,i) = 0.0d0
               ppiclf_rprop(PPICLF_R_JFX-1+k,i) = elibm_rhof*vol
     >            *(ud(k) - ppiclf_rprop(PPICLF_R_JUX-1+k,i))
     >            /elibm_dtf
            enddo
            cycle
         endif

         ! an EL particle is a parcel of rnp droplets of diameter dp
         ! with total liquid volume JVOLP; the equation of motion is
         ! per droplet, the feedback on the fluid per parcel
         dp    = ppiclf_rprop(PPICLF_R_JDP,i)
         rmass = rpi/6.0d0*dp**3*ppiclf_rprop(PPICLF_R_JRHOP,i)
         rnp   = ppiclf_rprop(PPICLF_R_JVOLP,i)/(rpi/6.0d0*dp**3)

         do k=1,3
            ur(k) = ppiclf_rprop(PPICLF_R_JUX-1+k,i)
     >            - ppiclf_y(PPICLF_JVX-1+k,i)
         enddo
         rep = sqrt(ur(1)**2 + ur(2)**2 + ur(3)**2)*dp*elibm_rhof
     >        /elibm_mu
         fd  = 3.0d0*rpi*elibm_mu*dp*(1.0d0 + 0.15d0*rep**0.687d0)

         ! EL-EL contact, accumulated by ppiclF into ydotc
         do k=1,3
            fc(k) = 0.0d0
         enddo
         if (elibm_ifnb .eq. 1) then
            call ppiclf_solve_NearestNeighbor(i)
            do k=1,3
               fc(k) = ppiclf_ydotc(PPICLF_JVX-1+k,i)
            enddo
         endif

         ! EL-IBM contact
         do k=1,3
            ppiclf_rprop(PPICLF_R_JCX-1+k,i) = 0.0d0
         enddo
         do b=1,elibm_nbody
            call elibm_ContactForce(ppiclf_y(PPICLF_JX,i),
     >         ppiclf_y(PPICLF_JVX,i),0.5d0*dp,rmass,
     >         elibm_bx(1,b),elibm_bv(1,b),elibm_brad(b),
     >         elibm_bmass(b),fb)
            do k=1,3
               ppiclf_rprop(PPICLF_R_JCX-1+k,i) =
     >            ppiclf_rprop(PPICLF_R_JCX-1+k,i) + fb(k)
               ppiclf_rprop(PPICLF_R_JIX-1+k,i) =
     >            ppiclf_rprop(PPICLF_R_JIX-1+k,i) + wimp*fb(k)
            enddo
         enddo

         do k=1,3
            ppiclf_ydot(PPICLF_JX -1+k,i) = ppiclf_y(PPICLF_JVX-1+k,i)
            ppiclf_ydot(PPICLF_JVX-1+k,i) = (fd*ur(k) + fc(k)
     >         + ppiclf_rprop(PPICLF_R_JCX-1+k,i))/rmass
     >         + elibm_grav(k)
            ! reaction of the drag on the fluid, for the whole parcel
            ppiclf_rprop(PPICLF_R_JFX-1+k,i) = -rnp*fd*ur(k)
         enddo
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
      integer*4 k
!
      do k=1,PPICLF_LRP_PRO
         map(k) = 0.0d0
      enddo

      if (nint(rprop(PPICLF_R_JTYPE)) .eq. PPICLF_TYPE_MARKER) then
         do k=1,3
            map(PPICLF_P_JFIBX-1+k) = rprop(PPICLF_R_JFX-1+k)
         enddo
      else
         do k=1,3
            map(PPICLF_P_JFELX-1+k) = rprop(PPICLF_R_JFX-1+k)
         enddo
         map(PPICLF_P_JPHIP) = rprop(PPICLF_R_JVOLP)
      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_user_EvalNearestNeighbor
     >                                        (i,j,yi,rpropi,yj,rpropj)
!
! EL-EL contact; markers and walls (j = 0) take no part.
!
      implicit none
!
      include "PPICLF"
      include "ELIBM.inc"
!
      integer*4 i, j
      real*8 yi(PPICLF_LRS), rpropi(PPICLF_LRP)
      real*8 yj(PPICLF_LRS), rpropj(PPICLF_LRP)
!
      real*8 f(3), rmi, rmj
      integer*4 k
!
      if (j .eq. 0) return
      if (nint(rpropi(PPICLF_R_JTYPE)) .ne. PPICLF_TYPE_EL) return
      if (nint(rpropj(PPICLF_R_JTYPE)) .ne. PPICLF_TYPE_EL) return

      rmi = 0.5235987755982988d0*rpropi(PPICLF_R_JDP)**3
     >     *rpropi(PPICLF_R_JRHOP)
      rmj = 0.5235987755982988d0*rpropj(PPICLF_R_JDP)**3
     >     *rpropj(PPICLF_R_JRHOP)
      call elibm_ContactForce(yi(PPICLF_JX),yi(PPICLF_JVX),
     >   0.5d0*rpropi(PPICLF_R_JDP),rmi,
     >   yj(PPICLF_JX),yj(PPICLF_JVX),0.5d0*rpropj(PPICLF_R_JDP),rmj,f)

      do k=1,3
         ppiclf_ydotc(PPICLF_JVX-1+k,i) =
     >      ppiclf_ydotc(PPICLF_JVX-1+k,i) + f(k)
      enddo

      return
      end
!-----------------------------------------------------------------------
