/*
 * Saturn Atmosphere General Circulation Modell(ATSAT)
 * 4th order Runge-Kutta scheme for the prognostic fields of the model:
 * t, u, v, w, the condensable species, and — when the closure is active — the turbulent
 * kinetic energy k* and its dissipation dis*.
 *
 * Restructured to ATJUP's arrangement (RungeKutta_Jup_Turb.cpp): the trigonometry and the
 * reciprocals of the grid spacing are hoisted out of the cell loop into a CellGeometry built
 * once per (i,j) column and handed to RHSSat, and the fifty-six stage variables that used to be
 * declared at function scope and carried through an OpenMP private() clause are now declared
 * where they are used.
*/

#include "cSaturnModel.h"

using namespace std;

/*
 * Area-weighted horizontal mean, at each level, of the very expression the buoyancy takes the
 * anomaly of. Mirrors ATJUP's computeBuoyancyRefLevel().
 *
 * The mean has to be the mean OF the quantity whose anomaly is taken, or the anomaly no longer
 * has zero mean at that height — which is the whole point: with it subtracted, only horizontal
 * density contrasts drive vertical motion and hydrostatic balance is left to carry the mean.
 * sin(theta) is the spherical area weight; cells with no positive temperature are skipped, which
 * also catches NaN, and a non-finite contribution is dropped rather than poisoning the level.
 */
void cSaturnModel::computeBuoyancyRefLevel(){
    if((int)buoy_ref_level.size() != im) buoy_ref_level.assign(im, 0.0);

    #pragma omp parallel for schedule(static)
    for(int i = 0; i < im; i++){
        double sum = 0.0, wsum = 0.0;
        for(int j = 0; j < jm; j++){
            const double wgt = sin(the.z[j]);              // spherical area weight
            for(int k = 0; k < km; k++){
                if(!(t.x[i][j][k] > 0.0)) continue;        // also catches NaN
                const double b = g * p_stat.x[i][j][k]
                               / (r_mix * R_mix * t.x[i][j][k] * t_ref);
                if(!std::isfinite(b)) continue;
                sum  += wgt * b;
                wsum += wgt;
            }
        }
        buoy_ref_level[i] = (wsum > 0.0) ? sum / wsum : 0.0;
    }
}

void cSaturnModel::RungeKuttaSat(){
    cout << endl << "      ATSAT: RungeKuttaSat" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // The buoyancy reference level, refreshed once per RK4 step, before any stage reads it.
    computeBuoyancyRefLevel();

    // ===== Geometry hoisted out of the cell loop =====
    //
    // These depend on i and j alone, and RHSSat used to recompute all of them — two library
    // trigonometric calls and eight divisions — on every one of its four calls per cell. The
    // tables below are built once per invocation and the per-column struct once per (i,j).
    //
    // sinthe is taken straight from sin(the.z[j]), with no polar floor and no absolute value.
    // ATJUP clamps its equivalent away from zero because its own loop runs nearer the pole and
    // the 1/sin terms overflow there; ATSAT has never done so and this restructure is not the
    // place to start. If ATSAT's polar metric is ever revisited, this is the line that decides it.
    std::vector<double> sinthe_tbl(jm), costhe_tbl(jm);
    for(int j = 0; j < jm; j++){
        sinthe_tbl[j] = sin(the.z[j]);
        costhe_tbl[j] = cos(the.z[j]);
    }

    const double inv_2dr   = 1.0 / (2.0 * dr);
    const double inv_2dthe = 1.0 / (2.0 * dthe);
    const double inv_2dphi = 1.0 / (2.0 * dphi);
    const double inv_dr2   = 1.0 / (dr * dr);
    const double inv_dthe2 = 1.0 / (dthe * dthe);
    const double inv_dphi2 = 1.0 / (dphi * dphi);

    // Turbulence integration. k* and dis* become prognostic here: RHS_Sat_Turb.cpp assembles
    // rhs_tke/rhs_dis and the four RK4 stages below advance them alongside t, u, v, w. Gated by
    // cSaturnModel::turb_active, which is ATSAT_TURB and the configured turb_model resolved into
    // one flag. With the closure off the tendencies are identically zero, so this leaves k*/dis*
    // at their initial values.
    //
    // The two clamps are the ones ATJUP uses. k* cannot be negative and cannot exceed a generous
    // physical ceiling; dis* has a floor because it appears in denominators throughout the
    // closure (nue = k/dis among them) and a zero there is an infinity one step later.
    const bool turb_on_rk = turb_active;
    const double tke_max_nd = 1000.0 / (u_0 * u_0);   // 1000 m2/s2
    constexpr double dis_min_nd = 1.0e-10;            // matches TurbulenceSat::dis_min

    #pragma omp parallel for collapse(2) schedule(static)
    for(int i = 1; i < im-1; i++){
//        for(int j = 1; j < jm-1; j++){
        for(int j = 2; j < jm-2; j++){
//        for(int j = 3; j < jm-3; j++){

            // Build the geometry struct once per (i,j) column.
            CellGeometry geo;
            geo.rm       = metricRadius(rad.z[i]);
            geo.rm2      = geo.rm * geo.rm;
            geo.exp_rm   = coord_stretching ? 1.0 / (geo.rm + 1.0) : 1.0;
            geo.exp_2_rm = geo.exp_rm * geo.exp_rm;
            geo.sinthe   = sinthe_tbl[j];
            geo.sinthe2  = geo.sinthe * geo.sinthe;
            geo.costhe   = costhe_tbl[j];
            geo.cotanthe             = geo.costhe / geo.sinthe;
            geo.inv_rm               = 1.0 / geo.rm;
            geo.inv_rm2              = 1.0 / geo.rm2;
            geo.inv_rmsinthe         = 1.0 / (geo.rm * geo.sinthe);
            geo.inv_rm2sinthe        = geo.inv_rm2 / geo.sinthe;
            geo.inv_rm2sinthe2       = geo.inv_rm2 / geo.sinthe2;
            geo.costhe_inv_rm2sinthe = geo.costhe * geo.inv_rm2sinthe;
            geo.inv_2dr   = inv_2dr;
            geo.inv_2dthe = inv_2dthe;
            geo.inv_2dphi = inv_2dphi;
            geo.inv_dr2   = inv_dr2;
            geo.inv_dthe2 = inv_dthe2;
            geo.inv_dphi2 = inv_dphi2;

            for(int k = 1; k < km-1; k++){

                // Start-of-step values, read once. All four stages integrate from these.
                const double tn_ijk     = tn.x[i][j][k];
                const double un_ijk     = un.x[i][j][k];
                const double vn_ijk     = vn.x[i][j][k];
                const double wn_ijk     = wn.x[i][j][k];
                const double h2on_ijk   = h2on.x[i][j][k];
                const double h2ocn_ijk  = h2o_cloudn.x[i][j][k];
                const double h2oin_ijk  = h2o_icen.x[i][j][k];
                const double ch4n_ijk   = ch4n.x[i][j][k];
                const double ch4cn_ijk  = ch4_cloudn.x[i][j][k];
                const double ch4in_ijk  = ch4_icen.x[i][j][k];
                const double h2sn_ijk   = h2sn.x[i][j][k];
                const double nh3n_ijk   = nh3n.x[i][j][k];
                const double nh3cn_ijk  = nh3_cloudn.x[i][j][k];
                const double nh3in_ijk  = nh3_icen.x[i][j][k];
                const double nh4shn_ijk = nh4shn.x[i][j][k];
                const double tken_ijk   = tken.x[i][j][k];
                const double disn_ijk   = disn.x[i][j][k];

                // ----- RK stage 1 -----
                cSaturnModel::RHSSat(i, j, k, geo);

                const double kt1     = rhs_t.x[i][j][k];
                const double ku1     = rhs_u.x[i][j][k];
                const double kv1     = rhs_v.x[i][j][k];
                const double kw1     = rhs_w.x[i][j][k];
                const double kc1     = rhs_h2o.x[i][j][k];
                const double kcloud1 = rhs_h2o_cloud.x[i][j][k];
                const double kice1   = rhs_h2o_ice.x[i][j][k];
                const double kch41       = rhs_ch4.x[i][j][k];
                const double kch4_cloud1 = rhs_ch4_cloud.x[i][j][k];
                const double kch4_ice1   = rhs_ch4_ice.x[i][j][k];
                const double kh2s1       = rhs_h2s.x[i][j][k];
                const double knh31       = rhs_nh3.x[i][j][k];
                const double knh3_cloud1 = rhs_nh3_cloud.x[i][j][k];
                const double knh3_ice1   = rhs_nh3_ice.x[i][j][k];
                const double knh4sh1     = rhs_nh4sh.x[i][j][k];
                const double ktke1 = rhs_tke.x[i][j][k], kdis1 = rhs_dis.x[i][j][k];

                t.x[i][j][k] = tn_ijk + kt1 * 0.5 * dt;
                if(turb_on_rk){
                    tke.x[i][j][k] = std::min(std::max(tken_ijk + ktke1 * 0.5 * dt, 0.0), tke_max_nd);
                    dis.x[i][j][k] = std::max(disn_ijk + kdis1 * 0.5 * dt, dis_min_nd);
                }
                u.x[i][j][k] = un_ijk + ku1 * 0.5 * dt;
                v.x[i][j][k] = vn_ijk + kv1 * 0.5 * dt;
                w.x[i][j][k] = wn_ijk + kw1 * 0.5 * dt;

                h2o.x[i][j][k] = h2on_ijk + kc1 * 0.5 * dt;
                h2o_cloud.x[i][j][k] = h2ocn_ijk + kcloud1 * 0.5 * dt;
                h2o_ice.x[i][j][k] = h2oin_ijk + kice1 * 0.5 * dt;

                ch4.x[i][j][k] = ch4n_ijk + kch41 * 0.5 * dt;
                ch4_cloud.x[i][j][k] = ch4cn_ijk + kch4_cloud1 * 0.5 * dt;
                ch4_ice.x[i][j][k] = ch4in_ijk + kch4_ice1 * 0.5 * dt;

                h2s.x[i][j][k] = h2sn_ijk + kh2s1 * 0.5 * dt;

                nh3.x[i][j][k] = nh3n_ijk + knh31 * 0.5 * dt;
                nh3_cloud.x[i][j][k] = nh3cn_ijk + knh3_cloud1 * 0.5 * dt;
                nh3_ice.x[i][j][k] = nh3in_ijk + knh3_ice1 * 0.5 * dt;

                nh4sh.x[i][j][k] = nh4shn_ijk + knh4sh1 * 0.5 * dt;


                // ----- RK stage 2 -----
                cSaturnModel::RHSSat(i, j, k, geo);

                const double kt2     = rhs_t.x[i][j][k];
                const double ku2     = rhs_u.x[i][j][k];
                const double kv2     = rhs_v.x[i][j][k];
                const double kw2     = rhs_w.x[i][j][k];
                const double kc2     = rhs_h2o.x[i][j][k];
                const double kcloud2 = rhs_h2o_cloud.x[i][j][k];
                const double kice2   = rhs_h2o_ice.x[i][j][k];
                const double kch42       = rhs_ch4.x[i][j][k];
                const double kch4_cloud2 = rhs_ch4_cloud.x[i][j][k];
                const double kch4_ice2   = rhs_ch4_ice.x[i][j][k];
                const double kh2s2       = rhs_h2s.x[i][j][k];
                const double knh32       = rhs_nh3.x[i][j][k];
                const double knh3_cloud2 = rhs_nh3_cloud.x[i][j][k];
                const double knh3_ice2   = rhs_nh3_ice.x[i][j][k];
                const double knh4sh2     = rhs_nh4sh.x[i][j][k];
                const double ktke2 = rhs_tke.x[i][j][k], kdis2 = rhs_dis.x[i][j][k];

                t.x[i][j][k] = tn_ijk + kt2 * 0.5 * dt;
                if(turb_on_rk){
                    tke.x[i][j][k] = std::min(std::max(tken_ijk + ktke2 * 0.5 * dt, 0.0), tke_max_nd);
                    dis.x[i][j][k] = std::max(disn_ijk + kdis2 * 0.5 * dt, dis_min_nd);
                }
                u.x[i][j][k] = un_ijk + ku2 * 0.5 * dt;
                v.x[i][j][k] = vn_ijk + kv2 * 0.5 * dt;
                w.x[i][j][k] = wn_ijk + kw2 * 0.5 * dt;

                h2o.x[i][j][k] = h2on_ijk + kc2 * 0.5 * dt;
                h2o_cloud.x[i][j][k] = h2ocn_ijk + kcloud2 * 0.5 * dt;
                h2o_ice.x[i][j][k] = h2oin_ijk + kice2 * 0.5 * dt;

                ch4.x[i][j][k] = ch4n_ijk + kch42 * 0.5 * dt;
                ch4_cloud.x[i][j][k] = ch4cn_ijk + kch4_cloud2 * 0.5 * dt;
                ch4_ice.x[i][j][k] = ch4in_ijk + kch4_ice2 * 0.5 * dt;

                h2s.x[i][j][k] = h2sn_ijk + kh2s2 * 0.5 * dt;

                nh3.x[i][j][k] = nh3n_ijk + knh32 * 0.5 * dt;
                nh3_cloud.x[i][j][k] = nh3cn_ijk + knh3_cloud2 * 0.5 * dt;
                nh3_ice.x[i][j][k] = nh3in_ijk + knh3_ice2 * 0.5 * dt;

                nh4sh.x[i][j][k] = nh4shn_ijk + knh4sh2 * 0.5 * dt;


                // ----- RK stage 3 -----
                cSaturnModel::RHSSat(i, j, k, geo);

                const double kt3     = rhs_t.x[i][j][k];
                const double ku3     = rhs_u.x[i][j][k];
                const double kv3     = rhs_v.x[i][j][k];
                const double kw3     = rhs_w.x[i][j][k];
                const double kc3     = rhs_h2o.x[i][j][k];
                const double kcloud3 = rhs_h2o_cloud.x[i][j][k];
                const double kice3   = rhs_h2o_ice.x[i][j][k];
                const double kch43       = rhs_ch4.x[i][j][k];
                const double kch4_cloud3 = rhs_ch4_cloud.x[i][j][k];
                const double kch4_ice3   = rhs_ch4_ice.x[i][j][k];
                const double kh2s3       = rhs_h2s.x[i][j][k];
                const double knh33       = rhs_nh3.x[i][j][k];
                const double knh3_cloud3 = rhs_nh3_cloud.x[i][j][k];
                const double knh3_ice3   = rhs_nh3_ice.x[i][j][k];
                const double knh4sh3     = rhs_nh4sh.x[i][j][k];
                const double ktke3 = rhs_tke.x[i][j][k], kdis3 = rhs_dis.x[i][j][k];

                t.x[i][j][k] = tn_ijk + kt3 * dt;
                if(turb_on_rk){
                    tke.x[i][j][k] = std::min(std::max(tken_ijk + ktke3 * dt, 0.0), tke_max_nd);
                    dis.x[i][j][k] = std::max(disn_ijk + kdis3 * dt, dis_min_nd);
                }
                u.x[i][j][k] = un_ijk + ku3 * dt;
                v.x[i][j][k] = vn_ijk + kv3 * dt;
                w.x[i][j][k] = wn_ijk + kw3 * dt;

                h2o.x[i][j][k] = h2on_ijk + kc3 * dt;
                h2o_cloud.x[i][j][k] = h2ocn_ijk + kcloud3 * dt;
                h2o_ice.x[i][j][k] = h2oin_ijk + kice3 * dt;

                ch4.x[i][j][k] = ch4n_ijk + kch43 * dt;
                ch4_cloud.x[i][j][k] = ch4cn_ijk + kch4_cloud3 * dt;
                ch4_ice.x[i][j][k] = ch4in_ijk + kch4_ice3 * dt;

                h2s.x[i][j][k] = h2sn_ijk + kh2s3 * dt;

                nh3.x[i][j][k] = nh3n_ijk + knh33 * dt;
                nh3_cloud.x[i][j][k] = nh3cn_ijk + knh3_cloud3 * dt;
                nh3_ice.x[i][j][k] = nh3in_ijk + knh3_ice3 * dt;

                nh4sh.x[i][j][k] = nh4shn_ijk + knh4sh3 * dt;


                // ----- RK stage 4 and the weighted update -----
                cSaturnModel::RHSSat(i, j, k, geo);

                const double kt4     = rhs_t.x[i][j][k];
                const double ku4     = rhs_u.x[i][j][k];
                const double kv4     = rhs_v.x[i][j][k];
                const double kw4     = rhs_w.x[i][j][k];
                const double kc4     = rhs_h2o.x[i][j][k];
                const double kcloud4 = rhs_h2o_cloud.x[i][j][k];
                const double kice4   = rhs_h2o_ice.x[i][j][k];
                const double kch44       = rhs_ch4.x[i][j][k];
                const double kch4_cloud4 = rhs_ch4_cloud.x[i][j][k];
                const double kch4_ice4   = rhs_ch4_ice.x[i][j][k];
                const double kh2s4       = rhs_h2s.x[i][j][k];
                const double knh34       = rhs_nh3.x[i][j][k];
                const double knh3_cloud4 = rhs_nh3_cloud.x[i][j][k];
                const double knh3_ice4   = rhs_nh3_ice.x[i][j][k];
                const double knh4sh4     = rhs_nh4sh.x[i][j][k];

                t.x[i][j][k] = tn_ijk + dt * (kt1 + 2.0 * kt2
                    + 2.0 * kt3 + kt4)/6.0;
                if(turb_on_rk){
                    const double ktke4 = rhs_tke.x[i][j][k], kdis4 = rhs_dis.x[i][j][k];
                    tke.x[i][j][k] = std::min(std::max(tken_ijk
                        + dt * (ktke1 + 2.0*ktke2 + 2.0*ktke3 + ktke4)/6.0, 0.0), tke_max_nd);
                    dis.x[i][j][k] = std::max(disn_ijk
                        + dt * (kdis1 + 2.0*kdis2 + 2.0*kdis3 + kdis4)/6.0, dis_min_nd);
                }
                u.x[i][j][k] = un_ijk + dt * (ku1 + 2.0 * ku2
                    + 2.0 * ku3 + ku4)/6.0;
                v.x[i][j][k] = vn_ijk + dt * (kv1 + 2.0 * kv2
                    + 2.0 * kv3 + kv4)/6.0;
                w.x[i][j][k] = wn_ijk + dt * (kw1 + 2.0 * kw2
                    + 2.0 * kw3 + kw4)/6.0;

                h2o.x[i][j][k] = h2on_ijk
                    + dt * (kc1 + 2.0 * kc2
                    + 2.0 * kc3 + kc4)/6.0;
                h2o_cloud.x[i][j][k] = h2ocn_ijk
                    + dt * (kcloud1 + 2.0 * kcloud2
                    + 2.0 * kcloud3 + kcloud4)/6.0;
                h2o_ice.x[i][j][k] = h2oin_ijk
                    + dt * (kice1 + 2.0 * kice2
                    + 2.0 * kice3 + kice4)/6.0;

                ch4.x[i][j][k] = ch4n_ijk
                    + dt * (kch41 + 2.0 * kch42
                    + 2.0 * kch43 + kch44)/6.0;
                ch4_cloud.x[i][j][k] = ch4cn_ijk
                    + dt * (kch4_cloud1 + 2.0 * kch4_cloud2
                    + 2.0 * kch4_cloud3 + kch4_cloud4)/6.0;
                ch4_ice.x[i][j][k] = ch4in_ijk
                    + dt * (kch4_ice1 + 2.0 * kch4_ice2
                    + 2.0 * kch4_ice3 + kch4_ice4)/6.0;

                h2s.x[i][j][k] = h2sn_ijk
                    + dt * (kh2s1 + 2.0 * kh2s2
                    + 2.0 * kh2s3 + kh2s4)/6.0;

                nh3.x[i][j][k] = nh3n_ijk
                    + dt * (knh31 + 2.0 * knh32
                    + 2.0 * knh33 + knh34)/6.0;
                nh3_cloud.x[i][j][k] = nh3cn_ijk
                    + dt * (knh3_cloud1 + 2.0 * knh3_cloud2
                    + 2.0 * knh3_cloud3 + knh3_cloud4)/6.0;
                nh3_ice.x[i][j][k] = nh3in_ijk
                    + dt * (knh3_ice1 + 2.0 * knh3_ice2
                    + 2.0 * knh3_ice3 + knh3_ice4)/6.0;

                nh4sh.x[i][j][k] = nh4shn_ijk
                    + dt * (knh4sh1 + 2.0 * knh4sh2
                    + 2.0 * knh4sh3 + knh4sh4)/6.0;
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for solveRungeKutta\n", elapsed.count() * 1e-9);

    cout << "      ATSAT: RungeKuttaSat ended" << endl;
    return;
}
/*
*
*/
