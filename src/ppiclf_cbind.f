!-----------------------------------------------------------------------
! C bindings for ppiclF. build_ppiclf.sh appends this file to the
! case's ppiclf_user.f, so it is compiled against the same
! PPICLF_USER.h sizes as the library.
!
! Upstream's -DPPICLC mode renames the Fortran routines themselves to
! their C labels, which breaks ppiclF's own internal calls (e.g.
! InitParticle -> AddParticles). We build without it and wrap instead.
! Routines taking character arguments (ReadParticleVTU, ReadWallVTK,
! InitTargetBins) are not wrapped.
!-----------------------------------------------------------------------
      subroutine ppiclc_comm_InitMPI(comm,id,np)
     > bind(C, name="ppiclc_comm_InitMPI")
      implicit none
      integer*4 comm
      integer*4 id
      integer*4 np
      call ppiclf_comm_InitMPI(comm,id,np)
      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclc_comm_InitOverlapMesh(ncell,lx1,ly1,lz1,xgrid,
     >   ygrid,zgrid)
     > bind(C, name="ppiclc_comm_InitOverlapMesh")
      implicit none
      integer*4 ncell
      integer*4 lx1
      integer*4 ly1
      integer*4 lz1
      real*8 xgrid(*)
      real*8 ygrid(*)
      real*8 zgrid(*)
      call ppiclf_comm_InitOverlapMesh(ncell,lx1,ly1,lz1,xgrid,ygrid,
     >   zgrid)
      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclc_solve_InitParticle(imethod,ndim,iendian,npart,y,
     >   rprop)
     > bind(C, name="ppiclc_solve_InitParticle")
      implicit none
      integer*4 imethod
      integer*4 ndim
      integer*4 iendian
      integer*4 npart
      real*8 y(*)
      real*8 rprop(*)
      call ppiclf_solve_InitParticle(imethod,ndim,iendian,npart,y,rprop)
      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclc_solve_AddParticles(npart,y,rprop)
     > bind(C, name="ppiclc_solve_AddParticles")
      implicit none
      integer*4 npart
      real*8 y(*)
      real*8 rprop(*)
      call ppiclf_solve_AddParticles(npart,y,rprop)
      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclc_solve_InitNeighborBin(width)
     > bind(C, name="ppiclc_solve_InitNeighborBin")
      implicit none
      real*8 width
      call ppiclf_solve_InitNeighborBin(width)
      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclc_solve_InitPeriodicX(xl,xr)
     > bind(C, name="ppiclc_solve_InitPeriodicX")
      implicit none
      real*8 xl
      real*8 xr
      call ppiclf_solve_InitPeriodicX(xl,xr)
      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclc_solve_InitPeriodicY(yl,yr)
     > bind(C, name="ppiclc_solve_InitPeriodicY")
      implicit none
      real*8 yl
      real*8 yr
      call ppiclf_solve_InitPeriodicY(yl,yr)
      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclc_solve_InitPeriodicZ(zl,zr)
     > bind(C, name="ppiclc_solve_InitPeriodicZ")
      implicit none
      real*8 zl
      real*8 zr
      call ppiclf_solve_InitPeriodicZ(zl,zr)
      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclc_solve_InitGaussianFilter(filt,alpha,iwallm)
     > bind(C, name="ppiclc_solve_InitGaussianFilter")
      implicit none
      real*8 filt
      real*8 alpha
      integer*4 iwallm
      call ppiclf_solve_InitGaussianFilter(filt,alpha,iwallm)
      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclc_solve_InitBoxFilter(filt,iwallm)
     > bind(C, name="ppiclc_solve_InitBoxFilter")
      implicit none
      real*8 filt
      integer*4 iwallm
      call ppiclf_solve_InitBoxFilter(filt,iwallm)
      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclc_solve_InitWall(xp1,xp2,xp3)
     > bind(C, name="ppiclc_solve_InitWall")
      implicit none
      real*8 xp1(*)
      real*8 xp2(*)
      real*8 xp3(*)
      call ppiclf_solve_InitWall(xp1,xp2,xp3)
      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclc_solve_IntegrateParticle(istep,iostep,dt,time)
     > bind(C, name="ppiclc_solve_IntegrateParticle")
      implicit none
      integer*4 istep
      integer*4 iostep
      real*8 dt
      real*8 time
      call ppiclf_solve_IntegrateParticle(istep,iostep,dt,time)
      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclc_solve_Project()
     > bind(C, name="ppiclc_solve_Project")
!
! Re-bin, exchange ghosts and project with the current particle state.
! Inside IntegrateParticle the last projection happens before the final
! RHS evaluation, so it lags the forces left in rprop by one RK stage.
!
      implicit none
!
      call ppiclf_solve_InitSolve

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclc_solve_InterpFieldUser(jp,infld)
     > bind(C, name="ppiclc_solve_InterpFieldUser")
!
      implicit none
!
      include "PPICLF"
!
      integer*4 jp
      real*8    infld(*)
!
      call ppiclf_solve_InterpFieldUser(jp,infld)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclc_solve_GetProFld(m,fld)
     > bind(C, name="ppiclc_solve_GetProFld")
!
! Copy projected field m (1-based) on the local fluid elements to fld,
! laid out as (LEX,LEY,LEZ,PPICLF_NEE).
!
      implicit none
!
      include "PPICLF"
!
      integer*4 m
      real*8    fld(*)
!
      integer*4 n
!
      n = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ*PPICLF_NEE
      call ppiclf_copy(fld,ppiclf_pro_fld(1,1,1,1,m),n)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclc_solve_GetNumParticles(npart)
     > bind(C, name="ppiclc_solve_GetNumParticles")
!
      implicit none
!
      include "PPICLF"
!
      integer*4 npart
!
      npart = ppiclf_npart

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclc_solve_GetParticles(y,rprop)
     > bind(C, name="ppiclc_solve_GetParticles")
!
! Copy the rank-local particles, laid out as (LRS,npart), (LRP,npart).
!
      implicit none
!
      include "PPICLF"
!
      real*8 y(*)
      real*8 rprop(*)
!
      call ppiclf_copy(y,ppiclf_y,PPICLF_LRS*ppiclf_npart)
      call ppiclf_copy(rprop,ppiclf_rprop,PPICLF_LRP*ppiclf_npart)

      return
      end
!-----------------------------------------------------------------------
