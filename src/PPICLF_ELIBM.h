/* Particle layout shared by src/ppiclf_user_elibm.f (Fortran) and
 * src/elibm.hpp (C++). A case's PPICLF_USER.h sets the sizes
 * (PPICLF_LPART, PPICLF_LEE, PPICLF_LEX/LEY/LEZ) and then includes this file.
 *
 * Two kinds of ppiclF particles:
 * EL point particles     (rprop JTYPE = PPICLF_TYPE_EL)
 * IBM surface markers    (rprop JTYPE = PPICLF_TYPE_MARKER), one set per
 * rigid body; the body state lives outside ppiclF. */

#define PPICLF_LRS 6
#define PPICLF_LRP 17
#define PPICLF_LRP_INT 3
#define PPICLF_LRP_PRO 7

#define PPICLF_TYPE_EL 0
#define PPICLF_TYPE_MARKER 1

/* state y */
#define PPICLF_JX 1
#define PPICLF_JY 2
#define PPICLF_JZ 3
#define PPICLF_JVX 4
#define PPICLF_JVY 5
#define PPICLF_JVZ 6

/* properties rprop */
#define PPICLF_R_JTYPE 1
#define PPICLF_R_JRHOP 2
#define PPICLF_R_JDP 3
/* EL: particle volume; marker: surface share * Gaussian thickness */
#define PPICLF_R_JVOLP 4
/* fluid velocity interpolated at the particle */
#define PPICLF_R_JUX 5
#define PPICLF_R_JUY 6
#define PPICLF_R_JUZ 7
/* force the particle exerts on the fluid */
#define PPICLF_R_JFX 8
#define PPICLF_R_JFY 9
#define PPICLF_R_JFZ 10
/* EL: global particle id; marker: body index (0-based) */
#define PPICLF_R_JID 11
/* marker: offset from the body centre in the body frame */
#define PPICLF_R_JRX 12
#define PPICLF_R_JRY 13
#define PPICLF_R_JRZ 14
/* EL: contact force exerted by IBM bodies on the particle */
#define PPICLF_R_JCX 15
#define PPICLF_R_JCY 16
#define PPICLF_R_JCZ 17

/* projected fields */
#define PPICLF_P_JFELX 1
#define PPICLF_P_JFELY 2
#define PPICLF_P_JFELZ 3
#define PPICLF_P_JFIBX 4
#define PPICLF_P_JFIBY 5
#define PPICLF_P_JFIBZ 6
#define PPICLF_P_JPHIP 7
