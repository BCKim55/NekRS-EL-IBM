!     Parameters and rigid-body state for src/ppiclf_user_elibm.f, set
!     from C++ (src/elibm.hpp) through the ppiclc_elibm_* bindings.
      real*8 elibm_mu, elibm_rhof, elibm_grav(3), elibm_dtf,
     >       elibm_kn, elibm_erest,
     >       elibm_sigma, elibm_theta, elibm_vbfrac
      integer*4 elibm_nbody, elibm_ifnb, elibm_istage, elibm_iattach,
     >          elibm_iplanar, elibm_ibridge_el, elibm_ibridge_ibm
      common /elibm_par_r/ elibm_mu, elibm_rhof, elibm_grav, elibm_dtf,
     >                     elibm_kn, elibm_erest,
     >                     elibm_sigma, elibm_theta, elibm_vbfrac
      common /elibm_par_i/ elibm_nbody, elibm_ifnb, elibm_istage,
     >                     elibm_iattach, elibm_iplanar,
     >                     elibm_ibridge_el, elibm_ibridge_ibm

      integer*4 elibm_nbmax
      parameter (elibm_nbmax = 64)
      real*8 elibm_bx(3,elibm_nbmax),
     >       elibm_bv(3,elibm_nbmax),
     >       elibm_bw(3,elibm_nbmax),
     >       elibm_brot(3,3,elibm_nbmax),
     >       elibm_brad(elibm_nbmax),
     >       elibm_bmass(elibm_nbmax)
      common /elibm_body_r/ elibm_bx, elibm_bv, elibm_bw, elibm_brot,
     >                      elibm_brad, elibm_bmass
      integer*4 elibm_bcyl(elibm_nbmax)
      common /elibm_body_i/ elibm_bcyl
