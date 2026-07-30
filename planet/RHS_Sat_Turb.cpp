/*
 * Saturn Atmosphere General Circulation Modell(ATSAT) applied to turbulent flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with additional transport equations for the condensable species and for the
 * turbulent kinetic energy k* and its dissipation dis* (epsilon* or omega*)
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 *
 * class to combine the right hand sides of the differential equations for the Runge-Kutta scheme
 *
 * Restructured to ATJUP's arrangement (RHS_Jup_Turb.cpp): the geometry of the cell arrives
 * precomputed in a CellGeometry, and the derivatives are named local variables rather than
 * entries in a parallel pair of std::vector<Array*> indexed by an enum. See the note below.
*/

#include "cSaturnModel.h"

using namespace std;


void cSaturnModel::RHSSat(int i, int j, int k, const CellGeometry& geo){

    // ===== Why this function no longer builds its own geometry or its own array tables =====
    //
    // What it used to do, per cell and per RK4 stage: call sin() and cos(), form eight products
    // and quotients of them, then build two std::vector<Array*> of 24 and 17 entries, allocate
    // six std::vector<double> to hold the derivatives, and walk them with an enum index. RHSSat
    // is called four times per cell per iteration, so all of that — six heap allocations and two
    // transcendental calls per call — happened four times for every cell of the grid, to produce
    // numbers that depend only on i and j.
    //
    // The geometry now arrives in the CellGeometry that RungeKuttaSat builds once per (i,j)
    // column, reciprocals already taken, and the derivatives are plain named doubles. That is
    // ATJUP's arrangement and the reason it is worth copying is not only speed: with the values
    // named, the metric factor in each term is visible at the point of use, which is what one
    // wants when the question is whether a 1/r belongs where it stands.
    //
    // NOT BIT-IDENTICAL, and it cannot be. The whole point of the precomputed struct is that
    // x/rm becomes x*inv_rm, and a division and a multiplication by the reciprocal do not round
    // the same way. Every term in the model moves in the last bits. The measurement is in the
    // commit message; nothing here changes any formula.

    // All geometric quantities come from the precomputed struct —
    // NO sin(), cos(), division, or reciprocal computation here.
    const double rm                   = geo.rm;
    const double sinthe               = geo.sinthe;
    const double costhe               = geo.costhe;
    const double cotanthe             = geo.cotanthe;
    const double inv_rm               = geo.inv_rm;
    const double inv_rm2              = geo.inv_rm2;
    const double inv_rmsinthe         = geo.inv_rmsinthe;
    const double inv_rm2sinthe        = geo.inv_rm2sinthe;
    const double inv_rm2sinthe2       = geo.inv_rm2sinthe2;
    const double costhe_inv_rm2sinthe = geo.costhe_inv_rm2sinthe;

    const double inv_2dr   = geo.inv_2dr;
    const double inv_2dthe = geo.inv_2dthe;
    const double inv_2dphi = geo.inv_2dphi;
    const double inv_dr2   = geo.inv_dr2;
    const double inv_dthe2 = geo.inv_dthe2;
    const double inv_dphi2 = geo.inv_dphi2;
    const double exp_rm    = geo.exp_rm;
    const double exp_2_rm  = geo.exp_2_rm;

    // Cache local cell values
    const double u_ijk = u.x[i][j][k];
    const double v_ijk = v.x[i][j][k];
    const double w_ijk = w.x[i][j][k];

    // ---- First-order derivative storage ----
    double dudr, dvdr, dwdr, dtdr, dpdr;
    double dh2odr, dh2ocdr, dh2oidr;
    double dh2sdr;
    double dnh3dr, dnh3cdr, dnh3idr, dnh4shdr;
    double dch4dr, dch4cdr, dch4idr;
    double dtkedr, ddisdr;

    double dudthe, dvdthe, dwdthe, dtdthe, dpdthe;
    double dh2odthe, dh2ocdthe, dh2oidthe;
    double dh2sdthe;
    double dnh3dthe, dnh3cdthe, dnh3idthe, dnh4shdthe;
    double dch4dthe, dch4cdthe, dch4idthe;
    double dtkedthe, ddisdthe;

    double dudphi, dvdphi, dwdphi, dtdphi, dpdphi;
    double dh2odphi, dh2ocdphi, dh2oidphi;
    double dh2sdphi;
    double dnh3dphi, dnh3cdphi, dnh3idphi, dnh4shdphi;
    double dch4dphi, dch4cdphi, dch4idphi;
    double dtkedphi, ddisdphi;

    // ---- Second-order derivative storage ----
    double d2udr2, d2vdr2, d2wdr2, d2tdr2;
    double d2h2odr2, d2h2ocdr2, d2h2oidr2;
    double d2h2sdr2;
    double d2nh3dr2, d2nh3cdr2, d2nh3idr2, d2nh4shdr2;
    double d2ch4dr2, d2ch4cdr2, d2ch4idr2;
    double d2tkedr2, d2disdr2;

    double d2udthe2, d2vdthe2, d2wdthe2, d2tdthe2;
    double d2h2odthe2, d2h2ocdthe2, d2h2oidthe2;
    double d2h2sdthe2;
    double d2nh3dthe2, d2nh3cdthe2, d2nh3idthe2, d2nh4shdthe2;
    double d2ch4dthe2, d2ch4cdthe2, d2ch4idthe2;
    double d2tkedthe2, d2disdthe2;

    double d2udphi2, d2vdphi2, d2wdphi2, d2tdphi2;
    double d2h2odphi2, d2h2ocdphi2, d2h2oidphi2;
    double d2h2sdphi2;
    double d2nh3dphi2, d2nh3cdphi2, d2nh3idphi2, d2nh4shdphi2;
    double d2ch4dphi2, d2ch4cdphi2, d2ch4idphi2;
    double d2tkedphi2, d2disdphi2;


    // ===== R-direction derivatives (central differences) =====
    // exp_rm and exp_2_rm are the coordinate-stretching factors. ATSAT does not stretch its
    // radial coordinate (coord_stretching is false, init_layer_heights is linear), so both are
    // exactly 1.0 and the multiplications below are exact. They are carried anyway, in ATJUP's
    // form, so that switching the stretching on is a change to one place and not to this file.
    #define COMPUTE_DR(FIELD, d1, d2) \
        d1 = (FIELD.x[i+1][j][k] - FIELD.x[i-1][j][k]) * inv_2dr * exp_rm; \
        d2 = (FIELD.x[i+1][j][k] - 2.0*FIELD.x[i][j][k] + FIELD.x[i-1][j][k]) * inv_dr2 * exp_2_rm;

    COMPUTE_DR(u,         dudr,     d2udr2)
    COMPUTE_DR(v,         dvdr,     d2vdr2)
    COMPUTE_DR(w,         dwdr,     d2wdr2)
    COMPUTE_DR(t,         dtdr,     d2tdr2)
    COMPUTE_DR(h2o,       dh2odr,   d2h2odr2)
    COMPUTE_DR(h2o_cloud, dh2ocdr,  d2h2ocdr2)
    COMPUTE_DR(h2o_ice,   dh2oidr,  d2h2oidr2)
    COMPUTE_DR(h2s,       dh2sdr,   d2h2sdr2)
    COMPUTE_DR(nh3,       dnh3dr,   d2nh3dr2)
    COMPUTE_DR(nh3_cloud, dnh3cdr,  d2nh3cdr2)
    COMPUTE_DR(nh3_ice,   dnh3idr,  d2nh3idr2)
    COMPUTE_DR(ch4,       dch4dr,   d2ch4dr2)
    COMPUTE_DR(ch4_cloud, dch4cdr,  d2ch4cdr2)
    COMPUTE_DR(ch4_ice,   dch4idr,  d2ch4idr2)
    COMPUTE_DR(nh4sh,     dnh4shdr, d2nh4shdr2)
    COMPUTE_DR(tke,       dtkedr,   d2tkedr2)
    COMPUTE_DR(dis,       ddisdr,   d2disdr2)
    dpdr = (p_dyn.x[i+1][j][k] - p_dyn.x[i-1][j][k]) * inv_2dr * exp_rm;
    #undef COMPUTE_DR


    // ===== Theta-direction derivatives (central differences) =====
    #define COMPUTE_DTHE(FIELD, d1, d2) \
        d1 = (FIELD.x[i][j+1][k] - FIELD.x[i][j-1][k]) * inv_2dthe; \
        d2 = (FIELD.x[i][j+1][k] - 2.0*FIELD.x[i][j][k] + FIELD.x[i][j-1][k]) * inv_dthe2;

    COMPUTE_DTHE(u,         dudthe,     d2udthe2)
    COMPUTE_DTHE(v,         dvdthe,     d2vdthe2)
    COMPUTE_DTHE(w,         dwdthe,     d2wdthe2)
    COMPUTE_DTHE(t,         dtdthe,     d2tdthe2)
    COMPUTE_DTHE(h2o,       dh2odthe,   d2h2odthe2)
    COMPUTE_DTHE(h2o_cloud, dh2ocdthe,  d2h2ocdthe2)
    COMPUTE_DTHE(h2o_ice,   dh2oidthe,  d2h2oidthe2)
    COMPUTE_DTHE(h2s,       dh2sdthe,   d2h2sdthe2)
    COMPUTE_DTHE(nh3,       dnh3dthe,   d2nh3dthe2)
    COMPUTE_DTHE(nh3_cloud, dnh3cdthe,  d2nh3cdthe2)
    COMPUTE_DTHE(nh3_ice,   dnh3idthe,  d2nh3idthe2)
    COMPUTE_DTHE(ch4,       dch4dthe,   d2ch4dthe2)
    COMPUTE_DTHE(ch4_cloud, dch4cdthe,  d2ch4cdthe2)
    COMPUTE_DTHE(ch4_ice,   dch4idthe,  d2ch4idthe2)
    COMPUTE_DTHE(nh4sh,     dnh4shdthe, d2nh4shdthe2)
    COMPUTE_DTHE(tke,       dtkedthe,   d2tkedthe2)
    COMPUTE_DTHE(dis,       ddisdthe,   d2disdthe2)
    dpdthe = (p_dyn.x[i][j+1][k] - p_dyn.x[i][j-1][k]) * inv_2dthe;
    #undef COMPUTE_DTHE


    // ===== Phi-direction derivatives (central differences) =====
    #define COMPUTE_DPHI(FIELD, d1, d2) \
        d1 = (FIELD.x[i][j][k+1] - FIELD.x[i][j][k-1]) * inv_2dphi; \
        d2 = (FIELD.x[i][j][k+1] - 2.0*FIELD.x[i][j][k] + FIELD.x[i][j][k-1]) * inv_dphi2;

    COMPUTE_DPHI(u,         dudphi,     d2udphi2)
    COMPUTE_DPHI(v,         dvdphi,     d2vdphi2)
    COMPUTE_DPHI(w,         dwdphi,     d2wdphi2)
    COMPUTE_DPHI(t,         dtdphi,     d2tdphi2)
    COMPUTE_DPHI(h2o,       dh2odphi,   d2h2odphi2)
    COMPUTE_DPHI(h2o_cloud, dh2ocdphi,  d2h2ocdphi2)
    COMPUTE_DPHI(h2o_ice,   dh2oidphi,  d2h2oidphi2)
    COMPUTE_DPHI(h2s,       dh2sdphi,   d2h2sdphi2)
    COMPUTE_DPHI(nh3,       dnh3dphi,   d2nh3dphi2)
    COMPUTE_DPHI(nh3_cloud, dnh3cdphi,  d2nh3cdphi2)
    COMPUTE_DPHI(nh3_ice,   dnh3idphi,  d2nh3idphi2)
    COMPUTE_DPHI(ch4,       dch4dphi,   d2ch4dphi2)
    COMPUTE_DPHI(ch4_cloud, dch4cdphi,  d2ch4cdphi2)
    COMPUTE_DPHI(ch4_ice,   dch4idphi,  d2ch4idphi2)
    COMPUTE_DPHI(nh4sh,     dnh4shdphi, d2nh4shdphi2)
    COMPUTE_DPHI(tke,       dtkedphi,   d2tkedphi2)
    COMPUTE_DPHI(dis,       ddisdphi,   d2disdphi2)
    dpdphi = (p_dyn.x[i][j][k+1] - p_dyn.x[i][j][k-1]) * inv_2dphi;
    #undef COMPUTE_DPHI

    // ATJUP follows the derivatives with two obstacle corrections — a mirrored pressure gradient
    // and a Neumann condition on k*/dis* — at the faces of its SeaMount. ATSAT has no solid body
    // (BC_seamount is never called, is_land returns false outright), so neither has anything to
    // act on and neither is carried over. If an obstacle is ever built, they are what this file is
    // missing.


// influence of the Coriolis force
    double Coriolis_rad = - 2.0 * omega * sinthe * w_ijk;
    double Coriolis_the = + 2.0 * omega * costhe * w_ijk;
    double Coriolis_phi = + 2.0 * omega * (- costhe * v_ijk
        + sinthe * u_ijk);


// ===== influence of the centrifugal force =====
//
// Centrifugal acceleration = Omega^2 * s * s_hat, with s = r*sin(theta) the distance from the
// ROTATION AXIS and s_hat = sin(theta)*e_r + cos(theta)*e_theta the unit vector pointing AWAY
// from it:
//     a_r     = +Omega^2 * r * sin^2(theta)
//     a_theta = +Omega^2 * r * sin(theta) * cos(theta)
//
// All three parts of this were wrong, and ATJUP found the same three in the same lines
// (RHS_Jup_Turb.cpp). ATSAT had Omega^2*r and Omega^2*r*|sin(theta)| entering rhs_u and rhs_v
// with a MINUS, so:
//   - the force pointed TOWARD the axis instead of away from it;
//   - the radial part had no sin^2 at all, i.e. full strength at the poles, where the
//     centrifugal force must vanish because there is no distance from the axis there;
//   - the meridional part had |sin| where sin*cos belongs, which is neither the right
//     magnitude nor equator-directed — and being an absolute value it could not change sign
//     at the equator, while the true term must, since it points toward the equator in BOTH
//     hemispheres.
//
// The true sine is reconstructed from the cosine rather than taken from geo.sinthe. In ATSAT
// the two are equal today, because RungeKuttaSat builds its sine table with no polar floor.
// ATJUP's table IS floored for the 1/sin^2 metric divisions, and such a floor has no business
// in a body force; writing it this way means that if ATSAT ever floors its metric sine the
// force does not silently inherit it. theta runs 0..pi so sin(theta) >= 0 and the positive
// root is the right one.
    const double sinthe_true = std::sqrt(std::max(0.0, 1.0 - costhe * costhe));
    double centrifugal_rad = omega * omega * rm * sinthe_true * sinthe_true;
    double centrifugal_the = omega * omega * rm * sinthe_true * costhe;

    double coeff_energy_p = u_0 * u_0/(cp_mix * t_ref); // coefficient for the source terms = 2.33e-4 (Eckert-number)


// transport terms in the Navier-Stokes-equations
    double pressure_t = coeff_energy_p * (u_ijk * dpdr
        + v_ijk * dpdthe * inv_rm
        + w_ijk * dpdphi * inv_rmsinthe);

    double transport_t = u_ijk * dtdr + v_ijk * dtdthe * inv_rm
        + w_ijk * dtdphi * inv_rmsinthe;

    double transport_u = u_ijk * dudr + v_ijk * dudthe * inv_rm
        + w_ijk * dudphi * inv_rmsinthe
        - (v_ijk * v_ijk + w_ijk * w_ijk) * inv_rm;
    double transport_v = u_ijk * dvdr + v_ijk * dvdthe * inv_rm
        + w_ijk * dvdphi * inv_rmsinthe
        + (u_ijk * v_ijk
        - w_ijk * w_ijk * cotanthe) * inv_rm;
    double transport_w = u_ijk * dwdr + v_ijk * dwdthe * inv_rm
        + w_ijk * dwdphi * inv_rmsinthe
        + (w_ijk * u_ijk
        + v_ijk * w_ijk * cotanthe) * inv_rm;

    double transport_ch4 = u_ijk * dch4dr + v_ijk * dch4dthe * inv_rm
        + w_ijk * dch4dphi * inv_rmsinthe;
    double transport_ch4_cloud = u_ijk * dch4cdr + v_ijk * dch4cdthe * inv_rm
        + w_ijk * dch4cdphi * inv_rmsinthe;
    double transport_ch4_ice = u_ijk * dch4idr + v_ijk * dch4idthe * inv_rm
        + w_ijk * dch4idphi * inv_rmsinthe;

    double transport_h2o = u_ijk * dh2odr + v_ijk * dh2odthe * inv_rm
        + w_ijk * dh2odphi * inv_rmsinthe;
    double transport_h2o_cloud = u_ijk * dh2ocdr + v_ijk * dh2ocdthe * inv_rm
        + w_ijk * dh2ocdphi * inv_rmsinthe;
    double transport_h2o_ice = u_ijk * dh2oidr + v_ijk * dh2oidthe * inv_rm
        + w_ijk * dh2oidphi * inv_rmsinthe;

    double transport_h2s = u_ijk * dh2sdr + v_ijk * dh2sdthe * inv_rm
        + w_ijk * dh2sdphi * inv_rmsinthe;

    double transport_nh3 = u_ijk * dnh3dr + v_ijk * dnh3dthe * inv_rm
        + w_ijk * dnh3dphi * inv_rmsinthe;
    double transport_nh3_cloud = u_ijk * dnh3cdr + v_ijk * dnh3cdthe * inv_rm
        + w_ijk * dnh3cdphi * inv_rmsinthe;
    double transport_nh3_ice = u_ijk * dnh3idr + v_ijk * dnh3idthe * inv_rm
        + w_ijk * dnh3idphi * inv_rmsinthe;

    double transport_nh4sh = u_ijk * dnh4shdr + v_ijk * dnh4shdthe * inv_rm
        + w_ijk * dnh4shdphi * inv_rmsinthe;


// diffusion terms in the Navier-Stokes-equations
    double two_inv_rm = 2.0 * inv_rm;
    double v_metric   = (1.0 + costhe / geo.sinthe2) * inv_rm2;

    double diffusion_t = (d2tdr2 + dtdr * two_inv_rm + d2tdthe2 * inv_rm2
        + dtdthe * costhe_inv_rm2sinthe + d2tdphi2 * inv_rm2sinthe2);

    double diffusion_u = d2udr2 + 2.0 * u_ijk * inv_rm2 + d2udthe2 * inv_rm2
        + 4.0 * dudr * inv_rm + dudthe * costhe_inv_rm2sinthe
        + d2udphi2 * inv_rm2sinthe2;
    double diffusion_v = d2vdr2 + dvdr * two_inv_rm + d2vdthe2 * inv_rm2
        + dvdthe * costhe_inv_rm2sinthe
        - v_metric * v_ijk
        + d2vdphi2 * inv_rm2sinthe2 + 2.0 * dudthe * inv_rm2
        - dwdphi * 2.0 * costhe * inv_rm2sinthe2;
    double diffusion_w = d2wdr2 + dwdr * two_inv_rm + d2wdthe2 * inv_rm2
        + dwdthe * costhe_inv_rm2sinthe
        - v_metric * w_ijk
        + d2wdphi2 * inv_rm2sinthe2 + 2.0 * dudphi * inv_rm2sinthe
        + dvdphi * 2.0 * costhe * inv_rm2sinthe2;

    double diffusion_ch4 = d2ch4dr2 + dch4dr * two_inv_rm + d2ch4dthe2 * inv_rm2
        + dch4dthe * costhe_inv_rm2sinthe
        + d2ch4dphi2 * inv_rm2sinthe2;
    double diffusion_ch4_cloud = d2ch4cdr2 + dch4cdr * two_inv_rm + d2ch4cdthe2 * inv_rm2
        + dch4cdthe * costhe_inv_rm2sinthe
        + d2ch4cdphi2 * inv_rm2sinthe2;
    double diffusion_ch4_ice = d2ch4idr2 + dch4idr * two_inv_rm + d2ch4idthe2 * inv_rm2
        + dch4idthe * costhe_inv_rm2sinthe
        + d2ch4idphi2 * inv_rm2sinthe2;

    double diffusion_h2o = d2h2odr2 + dh2odr * two_inv_rm + d2h2odthe2 * inv_rm2
        + dh2odthe * costhe_inv_rm2sinthe
        + d2h2odphi2 * inv_rm2sinthe2;
    double diffusion_h2o_cloud = d2h2ocdr2 + dh2ocdr * two_inv_rm + d2h2ocdthe2 * inv_rm2
        + dh2ocdthe * costhe_inv_rm2sinthe
        + d2h2ocdphi2 * inv_rm2sinthe2;
    double diffusion_h2o_ice = d2h2oidr2 + dh2oidr * two_inv_rm + d2h2oidthe2 * inv_rm2
        + dh2oidthe * costhe_inv_rm2sinthe
        + d2h2oidphi2 * inv_rm2sinthe2;

    double diffusion_h2s = d2h2sdr2 + dh2sdr * two_inv_rm + d2h2sdthe2 * inv_rm2
        + dh2sdthe * costhe_inv_rm2sinthe
        + d2h2sdphi2 * inv_rm2sinthe2;

    double diffusion_nh3 = d2nh3dr2 + dnh3dr * two_inv_rm + d2nh3dthe2 * inv_rm2
        + dnh3dthe * costhe_inv_rm2sinthe
        + d2nh3dphi2 * inv_rm2sinthe2;
    double diffusion_nh3_cloud = d2nh3cdr2 + dnh3cdr * two_inv_rm + d2nh3cdthe2 * inv_rm2
        + dnh3cdthe * costhe_inv_rm2sinthe
        + d2nh3cdphi2 * inv_rm2sinthe2;
    double diffusion_nh3_ice = d2nh3idr2 + dnh3idr * two_inv_rm + d2nh3idthe2 * inv_rm2
        + dnh3idthe * costhe_inv_rm2sinthe
        + d2nh3idphi2 * inv_rm2sinthe2;

    double diffusion_nh4sh = d2nh4shdr2 + dnh4shdr * two_inv_rm + d2nh4shdthe2 * inv_rm2
        + dnh4shdthe * costhe_inv_rm2sinthe
        + d2nh4shdphi2 * inv_rm2sinthe2;



// ===== Turbulence transport equations =====
//
// dk*/dt   = -v.grad k*   + div((1/re + nue*/sigma_k) grad k*)   + (P_k - Y_k)
// ddis*/dt = -v.grad dis* + div((1/re + nue*/sigma_w) grad dis*) + (P_w - Y_w + D_w)
//
// The source terms are what TurbulenceSat computed and left in tke_source / dis_source; this file
// only transports and diffuses them, exactly as it does for t and the species. The advection and
// the Laplacian are built in the same idiom as diffusion_t above, so the metric factors are the
// model's own and follow ATSAT_METRIC_RADIUS with everything else.
//
// With the closure off, nue, tke_source and dis_source are identically zero AND this block is
// skipped, so rhs_tke/rhs_dis stay zero and RungeKutta leaves k*/dis* untouched.
//
// sigma_k and sigma_w are the standard k-omega constants. TurbulenceSat carries its own sig_w2
// for the SST cross-diffusion; these two are the transport Prandtl numbers of the closure and are
// not derived from it.
    if(turb_active){
        constexpr double sig_k = 0.85, sig_w = 0.5;
        const double nue_here = nue.x[i][j][k];

        const double transport_tke = u_ijk * dtkedr + v_ijk * dtkedthe * inv_rm
                                   + w_ijk * dtkedphi * inv_rmsinthe;
        const double transport_dis = u_ijk * ddisdr + v_ijk * ddisdthe * inv_rm
                                   + w_ijk * ddisdphi * inv_rmsinthe;

        const double diffusion_tke = d2tkedr2 + dtkedr * two_inv_rm
                                   + d2tkedthe2 * inv_rm2
                                   + dtkedthe * costhe_inv_rm2sinthe
                                   + d2tkedphi2 * inv_rm2sinthe2;
        const double diffusion_dis = d2disdr2 + ddisdr * two_inv_rm
                                   + d2disdthe2 * inv_rm2
                                   + ddisdthe * costhe_inv_rm2sinthe
                                   + d2disdphi2 * inv_rm2sinthe2;

        rhs_tke.x[i][j][k] = - transport_tke
                           + diffusion_tke * (1.0/re + nue_here/sig_k)
                           + tke_source.x[i][j][k];
        rhs_dis.x[i][j][k] = - transport_dis
                           + diffusion_dis * (1.0/re + nue_here/sig_w)
                           + dis_source.x[i][j][k];
    } else {
        rhs_tke.x[i][j][k] = 0.0;
        rhs_dis.x[i][j][k] = 0.0;
    }

// ===== Turbulent (eddy) diffusion from the closure (opt-in) =====
//
// This is the step the closure was built for. Until it was wired, TurbulenceSat filled nue* and
// nothing read it: k* and dis* were transported, a viscosity was computed from them, and the
// momentum and scalar equations went on diffusing at their molecular rates. The closure could
// therefore be measured but not felt.
//
// nue is the DIMENSIONLESS eddy viscosity nue* = nue_phys/(u_0*L_atm), the same normalisation as
// 1/re (re = u_0*L_atm/nue_mol), so the two are directly additive and the effective momentum
// diffusivity is simply 1/re + nue*. Heat and the species get nue*/Pr_t added to their own
// 1/(sc*re): Pr_t = 0.9 is the standard turbulent Prandtl (Schmidt) number for shear-driven
// turbulence, meaning eddies mix a scalar slightly faster than they mix momentum. It is a closure
// constant, not a Saturn measurement — nothing in this model determines it.
//
// Each term below is written as the molecular term plus a separate eddy term rather than as one
// diffusion_x * (1/(sc*re) + nue*). Same physics, and it keeps the coupling-off path an exact
// addition of zero rather than a re-rounding of the molecular coefficient.
//
// Gated by ATSAT_TURB_COUPLING (default 0 = off). nue is nonzero only when the closure is active,
// so this knob alone does nothing.
//
// HOW BIG THIS IS, measured rather than assumed: with re = 1000 the molecular background is
// 1/re = 1.0e-3, and the closure's peak nue* over an 8-iteration run is 1.3e-5. The eddy
// viscosity is therefore about 77x SMALLER than the laminar one it is being added to, and at
// full strength this coupling moves the momentum diffusivity by ~1%. That is the opposite of
// the situation in ATJUP, whose comment at RHS_Jup_Turb.cpp:618 this port otherwise follows —
// the number belongs to re and to the state the closure is being run on, not to the scheme.
//
// So the knob being a double is not, here, a way to creep up on an unstable term. It is a way to
// scale nue* UP: at 1x this coupling cannot be expected to change the solution visibly, and a
// deliberate multiplier is the honest way to ask what an eddy viscosity of a believable size
// would do, rather than adjusting re or the closure constants until the closure looks important.
//
// ATSAT has no obstacle and therefore no wall-adjacent eddy viscosity: ATJUP adds a wall_nue term
// here that its computeWallViscosity() fills around the SeaMount. ATSAT's wall_nue array exists
// (it came with the port) and is identically zero, so it is not read.
    static const double turb_coupling = [](){
        const char* e = getenv("ATSAT_TURB_COUPLING"); return e ? atof(e) : 0.0; }();
    constexpr double Pr_t = 0.9;
    const double nue_t   = (turb_coupling != 0.0 && std::isfinite(nue.x[i][j][k]))
                         ? turb_coupling * std::max(0.0, nue.x[i][j][k]) : 0.0;
    const double nue_t_s = nue_t / Pr_t;      // scalar (heat / species) eddy diffusivity

// ===== Radiative heating source (opt-in) =====
//
// RadiationSat has been filling Q_rad since it was ported and nothing read it: Q_rad appeared
// only in the module that computes it, in the array allocation, and in printMinMax. The grey
// radiation scheme could be measured but not felt, exactly as the turbulence closure could not
// before ATSAT_TURB_COUPLING.
//
// Physically dT/dt = Q_rad/(rho*cp). Nondimensionalised by this model's energy scaling — lengths
// by L_rad, velocity by u_0, temperature by t_ref —
//
//     radiation_t = rad_coupling * Q_rad * L_rad / (rho * cp_mix * u_0 * t_ref).
//
// The density is the LOCAL one, read from rho_mix DIRECTLY rather than through rho_at(). That
// bypass is deliberate and is what ATJUP does here too: rho_at() is gated by ATSAT_LOCAL_RHO and
// returns the constant r_mix by default, which is right for the schemes that saturate on it and
// wrong for this one. Not a detail: the
// thin, cold upper atmosphere is where Q_rad > 0 does its heating, and it is the small rho there
// that lets those layers respond quickly and relax toward radiative equilibrium. Using r_mix
// would flatten exactly the part of the profile the scheme exists to set.
//
// rad_coupling = 1.0 IS THE PHYSICALLY CORRECT VALUE, not a tuning starting point — the
// expression is the exact nondimensional form of dT/dt = Q/(rho*cp) under this scaling. A value
// above 1 does not repair a scaling error, it is a deliberate ACCELERATION of the radiative
// timescale against the advective one, and it should be named as such when used. The temptation
// is real and worth stating: ATSAT's default timestep is 3.01e-5 nondimensional, which is 0.032 s
// of Saturn time, while radiative relaxation here is of order 1e7 s. At coupling = 1 equilibration
// would need ~1e9 iterations. Nothing in a run of practical length will show it.
    static const double rad_coupling = [](){
        const char* e = getenv("ATSAT_RAD_COUPLING"); return e ? atof(e) : 0.0; }();
    double radiation_t = 0.0;
    if(rad_coupling != 0.0){
        const double rho_f = rho_mix.x[i][j][k];               // [kg/m3], computeMixtureDensity()
        const double rho   = (rho_f > 0.0 && std::isfinite(rho_f)) ? rho_f : r_mix;
        const double L_rad = L_atm * 1.0e3;                    // shell thickness [m]
        if(rho > 0.0 && cp_mix > 0.0){
            radiation_t = rad_coupling * Q_rad.x[i][j][k] * L_rad
                        / (rho * cp_mix * u_0 * t_ref);
            // Explicit-scheme guard. The 1/rho factor can make this blow up in a very thin cell
            // at the top or in the polar corner. Test for non-finite FIRST — both comparisons
            // below are false for a NaN and would let it through.
            constexpr double rad_t_max = 0.5;
            if(!std::isfinite(radiation_t)) radiation_t = 0.0;
            else if(radiation_t >  rad_t_max) radiation_t =  rad_t_max;
            else if(radiation_t < -rad_t_max) radiation_t = -rad_t_max;
        }
    }

// ===== Latent-heat source from the precipitation microphysics (opt-in) =====
//
// Same conversion, and the same standing of Q_precip before this: filled by PrecipitationSat,
// read by nobody. Positive where riming and freezing release fusion heat, negative where melting
// or rain evaporation absorb it.
//
// THE DENSITY HERE IS r_mix, NOT the local one, and the difference is not a preference. Q_rad
// comes from real radiative fluxes, so dividing it by the local density is right. Q_precip is
// built from the condensate fields, and in this model those are mass concentrations defined as
// r_mix*ep*E/p — a mixing ratio scaled by the REFERENCE density — whose latent heat
// SaturationAdjustmentSat itself converts with /(cp_mix*r_mix) (Weather_Sat.cpp:163). The r_mix
// therefore cancels and leaves the true mixing-ratio tendency. Using rho_at() here instead would
// inflate the heating by r_mix/rho_local, which at Saturn's cloud decks is a large factor in the
// wrong direction. PrecipitationSat.h states the same convention at its own site.
    static const double precip_coupling = [](){
        const char* e = getenv("ATSAT_PRECIP_COUPLING"); return e ? atof(e) : 0.0; }();
    double precip_t = 0.0;
    if(precip_coupling != 0.0){
        const double L_rad = L_atm * 1.0e3;                    // shell thickness [m]
        if(r_mix > 0.0 && cp_mix > 0.0){
            precip_t = precip_coupling * Q_precip.x[i][j][k] * L_rad
                     / (r_mix * cp_mix * u_0 * t_ref);
            // r_mix is a constant, so there is no 1/rho blow-up to guard against here. The
            // limiter stays anyway: latent heating is stiff and locally concentrated — it
            // switches on hard at a phase boundary — and this is an explicit scheme.
            constexpr double precip_t_max = 0.5;
            if(!std::isfinite(precip_t)) precip_t = 0.0;
            else if(precip_t >  precip_t_max) precip_t =  precip_t_max;
            else if(precip_t < -precip_t_max) precip_t = -precip_t_max;
        }
    }

// right hand sides of the Navier-Stokes equations
    rhs_t.x[i][j][k] =
        + pressure_t
        - transport_t
        + diffusion_t/(re * pr) + diffusion_t * nue_t_s
        - chemical_reaction * thermalmassflux.x[i][j][k]
        + radiation_t
        + precip_t;

    rhs_u.x[i][j][k] =
        - dpdr
        + buoyancy * g * (p_stat.x[i][j][k] + p_dyn.x[i][j][k])  //  in kg/m³    (rho * g)
                        /(r_mix * R_mix * t.x[i][j][k] * t_ref)
//        + buoyancy * g * (1.0 - (t.x[i][j][k] - 1.0))          //  rho0 * g - rho0 * (t - t0)/t0 * g    for   del_rho << rho0
        - transport_u
        + diffusion_u/re + diffusion_u * nue_t
        - Coriolis * Coriolis_rad
        + centrifugal * centrifugal_rad;

    rhs_v.x[i][j][k] =
        - dpdthe * inv_rm
        - transport_v
        + diffusion_v/re + diffusion_v * nue_t
        - Coriolis * Coriolis_the
        + centrifugal * centrifugal_the;

    rhs_w.x[i][j][k] =
        - dpdphi * inv_rmsinthe
        - transport_w
        + diffusion_w/re + diffusion_w * nue_t
        - Coriolis * Coriolis_phi;

    rhs_ch4.x[i][j][k] =
        - transport_ch4
        + diffusion_ch4/(sc_ch4 * re) + diffusion_ch4 * nue_t_s;

    rhs_ch4_cloud.x[i][j][k] =
        - transport_ch4_cloud
        + diffusion_ch4_cloud/(sc_ch4 * re) + diffusion_ch4_cloud * nue_t_s;

    rhs_ch4_ice.x[i][j][k] =
        - transport_ch4_ice
        + diffusion_ch4_ice/(sc_ch4 * re) + diffusion_ch4_ice * nue_t_s;

    rhs_h2o.x[i][j][k] =
        - transport_h2o
        + diffusion_h2o/(sc_h2o * re) + diffusion_h2o * nue_t_s;

    rhs_h2o_cloud.x[i][j][k] =
        - transport_h2o_cloud
        + diffusion_h2o_cloud/(sc_h2o * re) + diffusion_h2o_cloud * nue_t_s;

    rhs_h2o_ice.x[i][j][k] =
        - transport_h2o_ice
        + diffusion_h2o_ice/(sc_h2o * re) + diffusion_h2o_ice * nue_t_s;

    rhs_h2s.x[i][j][k] =
        - transport_h2s
        + diffusion_h2s/(sc_h2s * re) + diffusion_h2s * nue_t_s
//        + chemical_reaction * w_h2s.x[i][j][k];
        + chemical_reaction * massflux_h2s.x[i][j][k];

    rhs_nh3.x[i][j][k] =
        - transport_nh3
        + diffusion_nh3/(sc_nh3 * re) + diffusion_nh3 * nue_t_s
//        + chemical_reaction * w_nh3.x[i][j][k];
        + chemical_reaction * massflux_nh3.x[i][j][k];

    rhs_nh3_cloud.x[i][j][k] =
        - transport_nh3_cloud
        + diffusion_nh3_cloud/(sc_nh3 * re) + diffusion_nh3_cloud * nue_t_s;

    rhs_nh3_ice.x[i][j][k] =
        - transport_nh3_ice
        + diffusion_nh3_ice/(sc_nh3 * re) + diffusion_nh3_ice * nue_t_s;

    rhs_nh4sh.x[i][j][k] =
        - transport_nh4sh
        + fluxlim_nh4sh.x[i][j][k]
        + diffusion_nh4sh/(sc_nh4sh * re) + diffusion_nh4sh * nue_t_s
//        + chemical_reaction * w_nh4sh.x[i][j][k];
        + chemical_reaction * massflux_nh4sh.x[i][j][k];


    aux_u.x[i][j][k] = rhs_u.x[i][j][k] + dpdr;
    aux_v.x[i][j][k] = rhs_v.x[i][j][k] + dpdthe * inv_rm;
    aux_w.x[i][j][k] = rhs_w.x[i][j][k] + dpdphi * inv_rmsinthe;

    return;
}
/*
*
*/
