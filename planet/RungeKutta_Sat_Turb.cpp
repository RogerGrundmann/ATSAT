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
    const double s_min = sinthe_min();          // 0.0 by default: no floor, as ATSAT always had
    for(int j = 0; j < jm; j++){
        sinthe_tbl[j] = sin(the.z[j]);
        if(sinthe_tbl[j] < s_min) sinthe_tbl[j] = s_min;
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
    // ===== Temperature limiter (ATSAT_T_LIMITER, default 0 = off) =====
    //
    // t is the only prognostic field in this model with no bound of any kind: the species have
    // FluxLimiterNH4SH and damp_wiggles on their mass fluxes, k* and dis* are clamped in this
    // very loop, and t has nothing. Advection cannot create a new extremum, so a temperature
    // that leaves the range spanned by its own neighbourhood is a scheme error rather than
    // physics — and the sharpest fronts in ATSAT are at the poles and the model top.
    //
    // What this does, after the RK4 update: clip t to the range spanned by the cell and its six
    // face neighbours IN THE OLD STATE (tn, untouched during the update). That is the clipping
    // step of an FCT scheme. It cannot create a new local extremum, and it leaves any update
    // that already stays inside the local bounds exactly as it was — so it is inert on a
    // well-behaved run and only bites where the scheme has already failed. It reports how often
    // it fires, because a limiter that acts silently hides the front it is standing in for.
    //
    // ATJUP additionally skips solid neighbours, whose tn is an extrapolated value rather than a
    // state. ATSAT has no solid body, so every neighbour takes part.
    static const int t_limiter_on = [](){
        const char* e = getenv("ATSAT_T_LIMITER"); return e ? atoi(e) : 0; }();
    long t_clip_hits = 0;

    const bool turb_on_rk = turb_active;
    const double tke_max_nd = 1000.0 / (u_0 * u_0);   // 1000 m2/s2
    constexpr double dis_min_nd = 1.0e-10;            // matches TurbulenceSat::dis_min

    // ===== RK4, WITH THE FOUR STAGES SEPARATED =====
    //
    // Every stage is now two passes over the whole grid with an implicit barrier between them:
    // one that fills rhs_* from the current stage input, and one that folds rhs_* into the
    // accumulator and forms the next stage input. Eight passes per step instead of one.
    //
    // WHY IT HAD TO CHANGE. All four stages used to run inside a single cell loop: RHSSat(i,j,k)
    // then an immediate overwrite of t, u, v, w and every species AT THAT CELL, four times over,
    // before the loop moved on. But RHSSat DIFFERENTIATES those same live fields — COMPUTE_DR,
    // COMPUTE_DTHE and COMPUTE_DPHI read them at i+-1, j+-1, k+-1 — so a cell computing its
    // stage 2 was reading neighbours that had already been advanced to their own stage 1, or had
    // not been touched yet, depending entirely on where the loop had got to.
    //
    // Under OpenMP that is a data race and it was THE source of ATSAT's run-to-run
    // irreproducibility: with every other parallel region in the model serialised and only this
    // loop left parallel, four runs produced four different states; serialised, the model is
    // bit-reproducible. It was found by bisection, one region at a time.
    //
    // Serially it was not a race but it was still not RK4: each stage differentiated neighbours
    // sitting at inconsistent stages, so the scheme was a pointwise four-substep update wearing
    // RK4's coefficients. Both defects have the same cause and the same fix, which is why this is
    // one commit and not two.
    //
    // The cost is memory: y_n lives in the *n arrays, the stage input in the live fields, and the
    // running sum k1 + 2k2 + 2k3 + k4 now needs a third place, the acc_* arrays.
    for(int stage = 0; stage < 4; stage++){

        // Offset of the NEXT stage's input from y_n, and this stage's weight in the sum.
        const double c_in = (stage == 0 || stage == 1) ? 0.5 * dt : (stage == 2 ? dt : 0.0);
        const double wgt  = (stage == 0 || stage == 3) ? 1.0 : 2.0;

        // ---- pass A: tendencies everywhere, from one consistent state ----
        // Reads the live fields with a stencil and writes rhs_*, which are different arrays, so
        // no cell can see a neighbour half-updated.
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 1; i < im-1; i++){
            for(int j = 2; j < jm-2; j++){
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
                for(int k = 1; k < km-1; k++)
                    cSaturnModel::RHSSat(i, j, k, geo);
            }
        }

        // ---- pass B: fold into the accumulator, then form the next stage input ----
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 1; i < im-1; i++){
            for(int j = 2; j < jm-2; j++){
                for(int k = 1; k < km-1; k++){
                    acc_t.x[i][j][k] = (stage == 0) ? wgt * rhs_t.x[i][j][k]
                                      : acc_t.x[i][j][k] + wgt * rhs_t.x[i][j][k];
                    acc_u.x[i][j][k] = (stage == 0) ? wgt * rhs_u.x[i][j][k]
                                      : acc_u.x[i][j][k] + wgt * rhs_u.x[i][j][k];
                    acc_v.x[i][j][k] = (stage == 0) ? wgt * rhs_v.x[i][j][k]
                                      : acc_v.x[i][j][k] + wgt * rhs_v.x[i][j][k];
                    acc_w.x[i][j][k] = (stage == 0) ? wgt * rhs_w.x[i][j][k]
                                      : acc_w.x[i][j][k] + wgt * rhs_w.x[i][j][k];
                    acc_h2o.x[i][j][k] = (stage == 0) ? wgt * rhs_h2o.x[i][j][k]
                                      : acc_h2o.x[i][j][k] + wgt * rhs_h2o.x[i][j][k];
                    acc_h2o_cloud.x[i][j][k] = (stage == 0) ? wgt * rhs_h2o_cloud.x[i][j][k]
                                      : acc_h2o_cloud.x[i][j][k] + wgt * rhs_h2o_cloud.x[i][j][k];
                    acc_h2o_ice.x[i][j][k] = (stage == 0) ? wgt * rhs_h2o_ice.x[i][j][k]
                                      : acc_h2o_ice.x[i][j][k] + wgt * rhs_h2o_ice.x[i][j][k];
                    acc_ch4.x[i][j][k] = (stage == 0) ? wgt * rhs_ch4.x[i][j][k]
                                      : acc_ch4.x[i][j][k] + wgt * rhs_ch4.x[i][j][k];
                    acc_ch4_cloud.x[i][j][k] = (stage == 0) ? wgt * rhs_ch4_cloud.x[i][j][k]
                                      : acc_ch4_cloud.x[i][j][k] + wgt * rhs_ch4_cloud.x[i][j][k];
                    acc_ch4_ice.x[i][j][k] = (stage == 0) ? wgt * rhs_ch4_ice.x[i][j][k]
                                      : acc_ch4_ice.x[i][j][k] + wgt * rhs_ch4_ice.x[i][j][k];
                    acc_h2s.x[i][j][k] = (stage == 0) ? wgt * rhs_h2s.x[i][j][k]
                                      : acc_h2s.x[i][j][k] + wgt * rhs_h2s.x[i][j][k];
                    acc_nh3.x[i][j][k] = (stage == 0) ? wgt * rhs_nh3.x[i][j][k]
                                      : acc_nh3.x[i][j][k] + wgt * rhs_nh3.x[i][j][k];
                    acc_nh3_cloud.x[i][j][k] = (stage == 0) ? wgt * rhs_nh3_cloud.x[i][j][k]
                                      : acc_nh3_cloud.x[i][j][k] + wgt * rhs_nh3_cloud.x[i][j][k];
                    acc_nh3_ice.x[i][j][k] = (stage == 0) ? wgt * rhs_nh3_ice.x[i][j][k]
                                      : acc_nh3_ice.x[i][j][k] + wgt * rhs_nh3_ice.x[i][j][k];
                    acc_nh4sh.x[i][j][k] = (stage == 0) ? wgt * rhs_nh4sh.x[i][j][k]
                                      : acc_nh4sh.x[i][j][k] + wgt * rhs_nh4sh.x[i][j][k];
                    if(turb_on_rk){
                        acc_tke.x[i][j][k] = (stage == 0) ? wgt * rhs_tke.x[i][j][k]
                                          : acc_tke.x[i][j][k] + wgt * rhs_tke.x[i][j][k];
                        acc_dis.x[i][j][k] = (stage == 0) ? wgt * rhs_dis.x[i][j][k]
                                          : acc_dis.x[i][j][k] + wgt * rhs_dis.x[i][j][k];
                    }

                    if(stage < 3){
                        t.x[i][j][k] = tn.x[i][j][k] + c_in * rhs_t.x[i][j][k];
                        u.x[i][j][k] = un.x[i][j][k] + c_in * rhs_u.x[i][j][k];
                        v.x[i][j][k] = vn.x[i][j][k] + c_in * rhs_v.x[i][j][k];
                        w.x[i][j][k] = wn.x[i][j][k] + c_in * rhs_w.x[i][j][k];
                        h2o.x[i][j][k] = h2on.x[i][j][k] + c_in * rhs_h2o.x[i][j][k];
                        h2o_cloud.x[i][j][k] = h2o_cloudn.x[i][j][k] + c_in * rhs_h2o_cloud.x[i][j][k];
                        h2o_ice.x[i][j][k] = h2o_icen.x[i][j][k] + c_in * rhs_h2o_ice.x[i][j][k];
                        ch4.x[i][j][k] = ch4n.x[i][j][k] + c_in * rhs_ch4.x[i][j][k];
                        ch4_cloud.x[i][j][k] = ch4_cloudn.x[i][j][k] + c_in * rhs_ch4_cloud.x[i][j][k];
                        ch4_ice.x[i][j][k] = ch4_icen.x[i][j][k] + c_in * rhs_ch4_ice.x[i][j][k];
                        h2s.x[i][j][k] = h2sn.x[i][j][k] + c_in * rhs_h2s.x[i][j][k];
                        nh3.x[i][j][k] = nh3n.x[i][j][k] + c_in * rhs_nh3.x[i][j][k];
                        nh3_cloud.x[i][j][k] = nh3_cloudn.x[i][j][k] + c_in * rhs_nh3_cloud.x[i][j][k];
                        nh3_ice.x[i][j][k] = nh3_icen.x[i][j][k] + c_in * rhs_nh3_ice.x[i][j][k];
                        nh4sh.x[i][j][k] = nh4shn.x[i][j][k] + c_in * rhs_nh4sh.x[i][j][k];
                        if(turb_on_rk){
                            tke.x[i][j][k] = std::min(std::max(tken.x[i][j][k]
                                + c_in * rhs_tke.x[i][j][k], 0.0), tke_max_nd);
                            dis.x[i][j][k] = std::max(disn.x[i][j][k]
                                + c_in * rhs_dis.x[i][j][k], dis_min_nd);
                        }
                    }
                }
            }
        }
    }

    // ===== Final assembly: y_{n+1} = y_n + dt/6 (k1 + 2k2 + 2k3 + k4) =====
    {
        const double sixth_dt = dt / 6.0;
        #pragma omp parallel for collapse(2) schedule(static) reduction(+:t_clip_hits)
        for(int i = 1; i < im-1; i++){
            for(int j = 2; j < jm-2; j++){
                for(int k = 1; k < km-1; k++){

                    const double tn_ijk = tn.x[i][j][k];
                    double t_new = tn_ijk + sixth_dt * acc_t.x[i][j][k];
                    if(t_limiter_on){
                        // Local bounds from the OLD state.
                        double lo = tn_ijk, hi = tn_ijk;
                        auto take = [&](int ii, int jj, int kk){
                            const double val = tn.x[ii][jj][kk];
                            if(val < lo) lo = val;
                            if(val > hi) hi = val;
                        };
                        take(i-1,j,k); take(i+1,j,k);
                        take(i,j-1,k); take(i,j+1,k);
                        take(i,j,k-1); take(i,j,k+1);
                        if(t_new < lo){ t_new = lo; ++t_clip_hits; }
                        else if(t_new > hi){ t_new = hi; ++t_clip_hits; }
                    }
                    t.x[i][j][k] = t_new;

                    if(turb_on_rk){
                        tke.x[i][j][k] = std::min(std::max(tken.x[i][j][k]
                            + sixth_dt * acc_tke.x[i][j][k], 0.0), tke_max_nd);
                        dis.x[i][j][k] = std::max(disn.x[i][j][k]
                            + sixth_dt * acc_dis.x[i][j][k], dis_min_nd);
                    }

                u.x[i][j][k] = un.x[i][j][k] + sixth_dt * acc_u.x[i][j][k];
                v.x[i][j][k] = vn.x[i][j][k] + sixth_dt * acc_v.x[i][j][k];
                w.x[i][j][k] = wn.x[i][j][k] + sixth_dt * acc_w.x[i][j][k];
                h2o.x[i][j][k] = h2on.x[i][j][k] + sixth_dt * acc_h2o.x[i][j][k];
                h2o_cloud.x[i][j][k] = h2o_cloudn.x[i][j][k] + sixth_dt * acc_h2o_cloud.x[i][j][k];
                h2o_ice.x[i][j][k] = h2o_icen.x[i][j][k] + sixth_dt * acc_h2o_ice.x[i][j][k];
                ch4.x[i][j][k] = ch4n.x[i][j][k] + sixth_dt * acc_ch4.x[i][j][k];
                ch4_cloud.x[i][j][k] = ch4_cloudn.x[i][j][k] + sixth_dt * acc_ch4_cloud.x[i][j][k];
                ch4_ice.x[i][j][k] = ch4_icen.x[i][j][k] + sixth_dt * acc_ch4_ice.x[i][j][k];
                h2s.x[i][j][k] = h2sn.x[i][j][k] + sixth_dt * acc_h2s.x[i][j][k];
                nh3.x[i][j][k] = nh3n.x[i][j][k] + sixth_dt * acc_nh3.x[i][j][k];
                nh3_cloud.x[i][j][k] = nh3_cloudn.x[i][j][k] + sixth_dt * acc_nh3_cloud.x[i][j][k];
                nh3_ice.x[i][j][k] = nh3_icen.x[i][j][k] + sixth_dt * acc_nh3_ice.x[i][j][k];
                nh4sh.x[i][j][k] = nh4shn.x[i][j][k] + sixth_dt * acc_nh4sh.x[i][j][k];
                }
            }
        }
    }

    if(t_clip_hits > 0)
        printf("      ATSAT: t limiter clipped %ld cells this iteration"
               " (ATSAT_T_LIMITER=0 to lift it)\n", t_clip_hits);

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for solveRungeKutta\n", elapsed.count() * 1e-9);

    cout << "      ATSAT: RungeKuttaSat ended" << endl;
    return;
}
/*
*
*/
