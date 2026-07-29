#include "cSaturnModel.h"

using namespace std;

#define dxdr_a(X, dx) \
    ((X->x[i+1][j][k] - X->x[i-1][j][k])/(2.0 * dx))
#define d2xdr2_a(X, dx2) \
    ((X->x[i+1][j][k] - 2.0 * X->x[i][j][k] + X->x[i-1][j][k])/dx2)

#define dxdthe_a(X, dx) \
    ((X->x[i][j+1][k] - X->x[i][j-1][k])/(2.0 * dx))
#define d2xdthe2_a(X, dx2) \
    ((X->x[i][j+1][k] - 2.0 * X->x[i][j][k] + X->x[i][j-1][k])/dx2)

#define dxdphi_a(X, dx) \
    ((X->x[i][j][k+1] - X->x[i][j][k-1])/(2.0 * dx))
#define d2xdphi2_a(X, dx2) \
    ((X->x[i][j][k+1] - 2.0 * X->x[i][j][k] + X->x[i][j][k-1])/dx2)


void cSaturnModel::RHSSat(int i, int j, int k){
    double dr2, dthe2, dphi2, rm2;
    double rm;
    double sinthe, sinthe2;
    double costhe, rmsinthe, rm2sinthe, rm2sinthe2;
    double cotanthe;

    dr2 = dr * dr;
    dthe2 = dthe * dthe;
    dphi2 = dphi * dphi;
    rm = metricRadius(rad.z[i]);
    rm2 = rm * rm;
    sinthe = sin(the.z[j]);
    sinthe2 = sinthe * sinthe;
    costhe = cos(the.z[j]);
    cotanthe = costhe/sinthe;
    rmsinthe = rm * sinthe;
    rm2sinthe = rm2 * sinthe;
    rm2sinthe2 = rm2 * sinthe2;


    std::vector<Array*> arrays_1{&u, &v, &w, &t, &p_dyn,
        &ch4, &ch4_cloud, &ch4_ice,
        &h2o, &h2o_cloud, &h2o_ice,
        &nh3, &nh3_cloud, &nh3_ice,
        &h2s, &nh4sh,
        &j_nh4sh, &j_h2s, &j_nh3,
        &jT_nh4sh, &jT_h2s, &jT_nh3,
        &tke, &dis};

    std::vector<Array*> arrays_2{&u, &v, &w, &t,
        &ch4, &ch4_cloud, &ch4_ice,
        &h2o, &h2o_cloud, &h2o_ice,
        &nh3, &nh3_cloud, &nh3_ice,
        &h2s, &nh4sh,
        &tke, &dis};

    enum array_index_1{i_u_1, i_v_1, i_w_1, i_t_1, i_p_1,
        i_ch4_1, i_ch4_cloud_1, i_ch4_ice_1,
        i_h2o_1, i_h2o_cloud_1, i_h2o_ice_1,
        i_nh3_1, i_nh3_cloud_1, i_nh3_ice_1,
        i_h2s_1, i_nh4sh_1,
        i_j_nh4sh_1, i_j_h2s_1, i_j_nh3_1,
        i_jT_nh4sh_1, i_jT_h2s_1, i_jT_nh3_1,
        i_tke_1, i_dis_1,
        last_array_index_1};

    enum array_index_2{i_u_2, i_v_2, i_w_2, i_t_2,
        i_ch4_2, i_ch4_cloud_2, i_ch4_ice_2,
        i_h2o_2, i_h2o_cloud_2, i_h2o_ice_2,
        i_nh3_2, i_nh3_cloud_2, i_nh3_ice_2,
        i_h2s_2, i_nh4sh_2,
        i_tke_2, i_dis_2,
        last_array_index_2};

    std::vector<double> dxdr_vals(last_array_index_1), 
                        dxdthe_vals(last_array_index_1), 
                        dxdphi_vals(last_array_index_1),
                        d2xdr2_vals(last_array_index_2),
                        d2xdthe2_vals(last_array_index_2),
                        d2xdphi2_vals(last_array_index_2);

// field gradients
    for(int n=0; n<last_array_index_1; n++){
        dxdr_vals[n] = dxdr_a(arrays_1[n], dr); // 1. order accurate
        dxdthe_vals[n] = dxdthe_a(arrays_1[n], dthe); // 1. order accurate
        dxdphi_vals[n] = dxdphi_a(arrays_1[n], dphi); // 1. order accurate
    }
    for(int n=0; n<last_array_index_2; n++){
        d2xdr2_vals[n] = d2xdr2_a(arrays_2[n], dr2); // 2. order accurate
        d2xdthe2_vals[n] = d2xdthe2_a(arrays_2[n], dthe2); // 2. order accurate
        d2xdphi2_vals[n] = d2xdphi2_a(arrays_2[n], dphi2); // 2. order accurate
    }

// 1. order derivatives
    double dudr = dxdr_vals[i_u_1], dvdr = dxdr_vals[i_v_1],
           dwdr = dxdr_vals[i_w_1], dtdr = dxdr_vals[i_t_1],
           dpdr = dxdr_vals[i_p_1],
           dch4dr = dxdr_vals[i_ch4_1], dch4cdr = dxdr_vals[i_ch4_cloud_1],
           dch4idr = dxdr_vals[i_ch4_ice_1],
           dh2odr = dxdr_vals[i_h2o_1], dh2ocdr = dxdr_vals[i_h2o_cloud_1],
           dh2oidr = dxdr_vals[i_h2o_ice_1],
           dh2sdr = dxdr_vals[i_h2s_1],
           dnh3dr = dxdr_vals[i_nh3_1], dnh3cdr = dxdr_vals[i_nh3_cloud_1],
           dnh3idr = dxdr_vals[i_nh3_ice_1],
           dnh4shdr = dxdr_vals[i_nh4sh_1];

    double dudthe = dxdthe_vals[i_u_1], dvdthe = dxdthe_vals[i_v_1],
           dwdthe = dxdthe_vals[i_w_1], dtdthe = dxdthe_vals[i_t_1],
           dpdthe = dxdthe_vals[i_p_1],
           dch4dthe = dxdthe_vals[i_ch4_1], dch4cdthe = dxdthe_vals[i_ch4_cloud_1],
           dch4idthe = dxdthe_vals[i_ch4_ice_1],
           dh2odthe = dxdthe_vals[i_h2o_1], dh2ocdthe = dxdthe_vals[i_h2o_cloud_1],
           dh2oidthe = dxdthe_vals[i_h2o_ice_1],
           dh2sdthe = dxdthe_vals[i_h2s_1],
           dnh3dthe = dxdthe_vals[i_nh3_1], dnh3cdthe = dxdthe_vals[i_nh3_cloud_1],
           dnh4shdthe = dxdthe_vals[i_nh4sh_1],
           dnh3idthe = dxdthe_vals[i_nh3_ice_1];

    double dudphi = dxdphi_vals[i_u_1], dvdphi = dxdphi_vals[i_v_1],
           dwdphi = dxdphi_vals[i_w_1], dtdphi = dxdphi_vals[i_t_1],
           dpdphi = dxdphi_vals[i_p_1],
           dch4dphi = dxdphi_vals[i_ch4_1], dch4cdphi = dxdphi_vals[i_ch4_cloud_1],
           dch4idphi = dxdphi_vals[i_ch4_ice_1],
           dh2odphi = dxdphi_vals[i_h2o_1], dh2ocdphi = dxdphi_vals[i_h2o_cloud_1],
           dh2oidphi = dxdphi_vals[i_h2o_ice_1],
           dh2sdphi = dxdphi_vals[i_h2s_1],
           dnh3dphi = dxdphi_vals[i_nh3_1], dnh3cdphi = dxdphi_vals[i_nh3_cloud_1],
           dnh3idphi = dxdphi_vals[i_nh3_ice_1],
           dnh4shdphi = dxdphi_vals[i_nh4sh_1];
           

// 2. order derivatives
    double d2udr2 = d2xdr2_vals[i_u_2], d2vdr2 = d2xdr2_vals[i_v_2],
           d2wdr2 = d2xdr2_vals[i_w_2], d2tdr2 = d2xdr2_vals[i_t_2],
           d2ch4dr2 = d2xdr2_vals[i_ch4_2], d2ch4cdr2 = d2xdr2_vals[i_ch4_cloud_2],
           d2ch4idr2 = d2xdr2_vals[i_ch4_ice_2],
           d2h2odr2 = d2xdr2_vals[i_h2o_2], d2h2ocdr2 = d2xdr2_vals[i_h2o_cloud_2],
           d2h2oidr2 = d2xdr2_vals[i_h2o_ice_2],
           d2h2sdr2 = d2xdr2_vals[i_h2s_2],
           d2nh3dr2 = d2xdr2_vals[i_nh3_2], d2nh3cdr2 = d2xdr2_vals[i_nh3_cloud_2],
           d2nh3idr2 = d2xdr2_vals[i_nh3_ice_2],
           d2nh4shdr2 = d2xdr2_vals[i_nh4sh_2];

    double d2udthe2 = d2xdthe2_vals[i_u_2], d2vdthe2 = d2xdthe2_vals[i_v_2],
           d2wdthe2 = d2xdthe2_vals[i_w_2], d2tdthe2 = d2xdthe2_vals[i_t_2],
           d2ch4dthe2 = d2xdthe2_vals[i_ch4_2], d2ch4cdthe2 = d2xdthe2_vals[i_ch4_cloud_2],
           d2ch4idthe2 = d2xdthe2_vals[i_ch4_ice_2],
           d2h2odthe2 = d2xdthe2_vals[i_h2o_2], d2h2ocdthe2 = d2xdthe2_vals[i_h2o_cloud_2],
           d2h2oidthe2 = d2xdthe2_vals[i_h2o_ice_2],
           d2h2sdthe2 = d2xdthe2_vals[i_h2s_2],
           d2nh3dthe2 = d2xdthe2_vals[i_nh3_2], d2nh3cdthe2 = d2xdthe2_vals[i_nh3_cloud_2],
           d2nh3idthe2 = d2xdthe2_vals[i_nh3_ice_2],
           d2nh4shdthe2 = d2xdthe2_vals[i_nh4sh_2];

    double d2udphi2 = d2xdphi2_vals[i_u_2], d2vdphi2 = d2xdphi2_vals[i_v_2],
           d2wdphi2 = d2xdphi2_vals[i_w_2], d2tdphi2 = d2xdphi2_vals[i_t_2],
           d2ch4dphi2 = d2xdphi2_vals[i_ch4_2], d2ch4cdphi2 = d2xdphi2_vals[i_ch4_cloud_2],
           d2ch4idphi2 = d2xdphi2_vals[i_ch4_ice_2],
           d2h2odphi2 = d2xdphi2_vals[i_h2o_2], d2h2ocdphi2 = d2xdphi2_vals[i_h2o_cloud_2],
           d2h2oidphi2 = d2xdphi2_vals[i_h2o_ice_2],
           d2h2sdphi2 = d2xdphi2_vals[i_h2s_2],
           d2nh3dphi2 = d2xdphi2_vals[i_nh3_2], d2nh3cdphi2 = d2xdphi2_vals[i_nh3_cloud_2],
           d2nh3idphi2 = d2xdphi2_vals[i_nh3_ice_2],
           d2nh4shdphi2 = d2xdphi2_vals[i_nh4sh_2];


// influence of the Coriolis force
    double Coriolis_rad = - 2.0 * omega * sinthe * w.x[i][j][k];
    double Coriolis_the = + 2.0 * omega * costhe * w.x[i][j][k];
    double Coriolis_phi = + 2.0 * omega * (- costhe * v.x[i][j][k] 
        + sinthe * u.x[i][j][k]);


// influence of the centrifugal force
    double centrifugal_rad = omega * omega * rm;
    double centrifugal_the = omega * omega * rm * fabs(sinthe);

    double coeff_energy_p = u_0 * u_0/(cp_mix * t_ref); // coefficient for the source terms = 2.33e-4 (Eckert-number)


// transport terms in the Navier-Stokes-equations
    double pressure_t = coeff_energy_p * (u.x[i][j][k] * dpdr
        + v.x[i][j][k] * dpdthe/rm 
        + w.x[i][j][k] * dpdphi/rmsinthe);

    double transport_t = u.x[i][j][k] * dtdr + v.x[i][j][k] * dtdthe/rm
        + w.x[i][j][k] * dtdphi/rmsinthe; 

    double transport_u = u.x[i][j][k] * dudr + v.x[i][j][k] * dudthe/rm 
        + w.x[i][j][k] * dudphi/rmsinthe
        - (v.x[i][j][k] * v.x[i][j][k] + w.x[i][j][k] * w.x[i][j][k])/rm;
    double transport_v = u.x[i][j][k] * dvdr + v.x[i][j][k] * dvdthe/rm
        + w.x[i][j][k] * dvdphi/rmsinthe 
        + (u.x[i][j][k] * v.x[i][j][k] 
        - w.x[i][j][k] * w.x[i][j][k] * cotanthe)/rm;
    double transport_w = u.x[i][j][k] * dwdr + v.x[i][j][k] * dwdthe/rm
        + w.x[i][j][k] * dwdphi/rmsinthe 
        + (w.x[i][j][k] * u.x[i][j][k]
        + v.x[i][j][k] * w.x[i][j][k] * cotanthe)/rm;

    double transport_ch4 = u.x[i][j][k] * dch4dr + v.x[i][j][k] * dch4dthe/rm
        + w.x[i][j][k] * dch4dphi/rmsinthe;
    double transport_ch4_cloud = u.x[i][j][k] * dch4cdr + v.x[i][j][k] * dch4cdthe/rm
        + w.x[i][j][k] * dch4cdphi/rmsinthe;
    double transport_ch4_ice = u.x[i][j][k] * dch4idr + v.x[i][j][k] * dch4idthe/rm
        + w.x[i][j][k] * dch4idphi/rmsinthe;

    double transport_h2o = u.x[i][j][k] * dh2odr + v.x[i][j][k] * dh2odthe/rm
        + w.x[i][j][k] * dh2odphi/rmsinthe;
    double transport_h2o_cloud = u.x[i][j][k] * dh2ocdr + v.x[i][j][k] * dh2ocdthe/rm
        + w.x[i][j][k] * dh2ocdphi/rmsinthe;
    double transport_h2o_ice = u.x[i][j][k] * dh2oidr + v.x[i][j][k] * dh2oidthe/rm
        + w.x[i][j][k] * dh2oidphi/rmsinthe;

    double transport_h2s = u.x[i][j][k] * dh2sdr + v.x[i][j][k] * dh2sdthe/rm
        + w.x[i][j][k] * dh2sdphi/rmsinthe;

    double transport_nh3 = u.x[i][j][k] * dnh3dr + v.x[i][j][k] * dnh3dthe/rm
        + w.x[i][j][k] * dnh3dphi/rmsinthe;
    double transport_nh3_cloud = u.x[i][j][k] * dnh3cdr + v.x[i][j][k] * dnh3cdthe/rm
        + w.x[i][j][k] * dnh3cdphi/rmsinthe;
    double transport_nh3_ice = u.x[i][j][k] * dnh3idr + v.x[i][j][k] * dnh3idthe/rm
        + w.x[i][j][k] * dnh3idphi/rmsinthe;

    double transport_nh4sh = u.x[i][j][k] * dnh4shdr + v.x[i][j][k] * dnh4shdthe/rm
        + w.x[i][j][k] * dnh4shdphi/rmsinthe;


// diffusion terms in the Navier-Stokes-equations
    double diffusion_t = (d2tdr2 + dtdr * 2.0/rm + d2tdthe2/rm2
        + dtdthe * costhe/rm2sinthe + d2tdphi2/rm2sinthe2);

    double diffusion_u = d2udr2 + 2.0 * u.x[i][j][k]/rm2 + d2udthe2/rm2
        + 4.0 * dudr/rm + dudthe * costhe/rm2sinthe 
        + d2udphi2/rm2sinthe2;
    double diffusion_v = d2vdr2 + dvdr * 2.0/rm + d2vdthe2/rm2 
        + dvdthe/rm2sinthe * costhe
        - (1.0 + costhe/sinthe2)/rm2 * v.x[i][j][k] 
        + d2vdphi2/rm2sinthe2 + 2.0 * dudthe/rm2 
        - dwdphi * 2.0 * costhe/rm2sinthe2;
    double diffusion_w = d2wdr2 + dwdr * 2.0/rm + d2wdthe2/rm2
        + dwdthe/rm2sinthe * costhe 
        - (1.0 + costhe/sinthe2)/rm2 * w.x[i][j][k]
        + d2wdphi2/rm2sinthe2 + 2.0 * dudphi/rm2sinthe 
        + dvdphi * 2.0 * costhe/rm2sinthe2;

    double diffusion_ch4 = d2ch4dr2 + dch4dr * 2.0/rm + d2ch4dthe2/rm2
        + dch4dthe * costhe/rm2sinthe
        + d2ch4dphi2/rm2sinthe2;
    double diffusion_ch4_cloud = d2ch4cdr2 + dch4cdr * 2.0/rm + d2ch4cdthe2/rm2
        + dch4cdthe * costhe/rm2sinthe
        + d2ch4cdphi2/rm2sinthe2;
    double diffusion_ch4_ice = d2ch4idr2 + dch4idr * 2.0/rm + d2ch4idthe2/rm2
        + dch4idthe * costhe/rm2sinthe
        + d2ch4idphi2/rm2sinthe2;

    double diffusion_h2o = d2h2odr2 + dh2odr * 2.0/rm + d2h2odthe2/rm2
        + dh2odthe * costhe/rm2sinthe
        + d2h2odphi2/rm2sinthe2;
    double diffusion_h2o_cloud = d2h2ocdr2 + dh2ocdr * 2.0/rm + d2h2ocdthe2/rm2
        + dh2ocdthe * costhe/rm2sinthe 
        + d2h2ocdphi2/rm2sinthe2;
    double diffusion_h2o_ice = d2h2oidr2 + dh2oidr * 2.0/rm + d2h2oidthe2/rm2
        + dh2oidthe * costhe/rm2sinthe 
        + d2h2oidphi2/rm2sinthe2;

    double diffusion_h2s = d2h2sdr2 + dh2sdr * 2.0/rm + d2h2sdthe2/rm2
        + dh2sdthe * costhe/rm2sinthe 
        + d2h2sdphi2/rm2sinthe2;

    double diffusion_nh3 = d2nh3dr2 + dnh3dr * 2.0/rm + d2nh3dthe2/rm2
        + dnh3dthe * costhe/rm2sinthe 
        + d2nh3dphi2/rm2sinthe2;
    double diffusion_nh3_cloud = d2nh3cdr2 + dnh3cdr * 2.0/rm + d2nh3cdthe2/rm2
        + dnh3cdthe * costhe/rm2sinthe 
        + d2nh3cdphi2/rm2sinthe2;
    double diffusion_nh3_ice = d2nh3idr2 + dnh3idr * 2.0/rm + d2nh3idthe2/rm2
        + dnh3idthe * costhe/rm2sinthe 
        + d2nh3idphi2/rm2sinthe2;

    double diffusion_nh4sh = d2nh4shdr2 + dnh4shdr * 2.0/rm + d2nh4shdthe2/rm2
        + dnh4shdthe * costhe/rm2sinthe 
        + d2nh4shdphi2/rm2sinthe2;



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
// skipped, so rhs_tke/rhs_dis stay zero, RungeKutta leaves k*/dis* untouched, and the run is
// bit-identical to the model before the closure existed.
//
// sigma_k and sigma_w are the standard k-omega constants. TurbulenceSat carries its own sig_w2
// for the SST cross-diffusion; these two are the transport Prandtl numbers of the closure and are
// not derived from it.
    static const int turb_on = [](){
        const char* e = getenv("ATSAT_TURB"); return e ? atoi(e) : 0; }();
    if(turb_on != 0){
        constexpr double sig_k = 0.85, sig_w = 0.5;
        const double nue_here = nue.x[i][j][k];

        const double dtkedr   = dxdr_vals[i_tke_1],   ddisdr   = dxdr_vals[i_dis_1];
        const double dtkedthe = dxdthe_vals[i_tke_1], ddisdthe = dxdthe_vals[i_dis_1];
        const double dtkedphi = dxdphi_vals[i_tke_1], ddisdphi = dxdphi_vals[i_dis_1];

        const double transport_tke = u.x[i][j][k] * dtkedr + v.x[i][j][k] * dtkedthe/rm
                                   + w.x[i][j][k] * dtkedphi/rmsinthe;
        const double transport_dis = u.x[i][j][k] * ddisdr + v.x[i][j][k] * ddisdthe/rm
                                   + w.x[i][j][k] * ddisdphi/rmsinthe;

        const double diffusion_tke = d2xdr2_vals[i_tke_2] + dtkedr * 2.0/rm
                                   + d2xdthe2_vals[i_tke_2]/rm2
                                   + dtkedthe * costhe/(rm2 * sinthe)
                                   + d2xdphi2_vals[i_tke_2]/rm2sinthe2;
        const double diffusion_dis = d2xdr2_vals[i_dis_2] + ddisdr * 2.0/rm
                                   + d2xdthe2_vals[i_dis_2]/rm2
                                   + ddisdthe * costhe/(rm2 * sinthe)
                                   + d2xdphi2_vals[i_dis_2]/rm2sinthe2;

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

// right hand sides of the Navier-Stokes equations
    rhs_t.x[i][j][k] = 
        + pressure_t
        - transport_t 
//        + diffusion_t/(re * pr);
        + diffusion_t/(re * pr)
        - chemical_reaction * thermalmassflux.x[i][j][k];

    rhs_u.x[i][j][k] = 
        - dpdr
        + buoyancy * g * (p_stat.x[i][j][k] + p_dyn.x[i][j][k])  //  in kg/m³    (rho * g)
                        /(r_mix * R_mix * t.x[i][j][k] * t_ref)
//        + buoyancy * g * (1.0 - (t.x[i][j][k] - 1.0))          //  rho0 * g - rho0 * (t - t0)/t0 * g    for   del_rho << rho0
        - transport_u 
        + diffusion_u/re
        - Coriolis * Coriolis_rad
        - centrifugal * centrifugal_rad;

    rhs_v.x[i][j][k] = 
        - dpdthe/rm
        - transport_v
        + diffusion_v/re
        - Coriolis * Coriolis_the
        - centrifugal * centrifugal_the;

    rhs_w.x[i][j][k] = 
        - dpdphi/rmsinthe
        - transport_w
        + diffusion_w/re
        - Coriolis * Coriolis_phi;

    rhs_ch4.x[i][j][k] =
        - transport_ch4
        + diffusion_ch4/(sc_ch4 * re);

    rhs_ch4_cloud.x[i][j][k] =
        - transport_ch4_cloud
        + diffusion_ch4_cloud/(sc_ch4 * re);

    rhs_ch4_ice.x[i][j][k] =
        - transport_ch4_ice
        + diffusion_ch4_ice/(sc_ch4 * re);

    rhs_h2o.x[i][j][k] =
        - transport_h2o
        + diffusion_h2o/(sc_h2o * re);

    rhs_h2o_cloud.x[i][j][k] = 
        - transport_h2o_cloud
        + diffusion_h2o_cloud/(sc_h2o * re);

    rhs_h2o_ice.x[i][j][k] = 
        - transport_h2o_ice
        + diffusion_h2o_ice/(sc_h2o * re);

    rhs_h2s.x[i][j][k] = 
        - transport_h2s
        + diffusion_h2s/(sc_h2s * re)
//        + chemical_reaction * w_h2s.x[i][j][k];
        + chemical_reaction * massflux_h2s.x[i][j][k];

    rhs_nh3.x[i][j][k] = 
        - transport_nh3
        + diffusion_nh3/(sc_nh3 * re)
//        + chemical_reaction * w_nh3.x[i][j][k];
        + chemical_reaction * massflux_nh3.x[i][j][k];

    rhs_nh3_cloud.x[i][j][k] = 
        - transport_nh3_cloud
        + diffusion_nh3_cloud/(sc_nh3 * re);

    rhs_nh3_ice.x[i][j][k] = 
        - transport_nh3_ice
        + diffusion_nh3_ice/(sc_nh3 * re);

    rhs_nh4sh.x[i][j][k] =
        - transport_nh4sh
        + fluxlim_nh4sh.x[i][j][k]
        + diffusion_nh4sh/(sc_nh4sh * re)
//        + chemical_reaction * w_nh4sh.x[i][j][k];
        + chemical_reaction * massflux_nh4sh.x[i][j][k];


    aux_u.x[i][j][k] = rhs_u.x[i][j][k] + dpdr;
    aux_v.x[i][j][k] = rhs_v.x[i][j][k] + dpdthe/rm;
    aux_w.x[i][j][k] = rhs_w.x[i][j][k] + dpdphi/rmsinthe;

    return;
}
/*
*
*/

