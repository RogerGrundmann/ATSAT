/*
 * Atmosphere General Circulation Modell(ATNEPT) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in aa spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and nh3 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to prepare the boundary and initial conditions for diverse variables
*/
#include "cSaturnModel.h"
#include "ATPhys.h"   // saturation_vapour_pressure
#include "SaturationAdjustmentSat.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;

// SaturationAdjustmentSat::run() now lives in SaturationAdjustmentSat.cpp, which carries the
// mirrored ATJUP algorithm and dispatches back to cSaturnModel::Saturation_Adjustment() below
// only with ATSAT_SATADJ=0. That routine is the inherited one and was the default until
// 2026-10-10; the header records the ten points on which the two differ.

//Tao, W.-K., Simpson, J., and McCumber, M.: 
//An Ice-Water Saturation Adjustment, American Meteorological Society, Notes and 
//Correspondence, Volume 1, 321–235, 1988, 1988.
//Ice_Water_SaturationAdjustment, distribution of cloud ice and cloud water 
//dependent on water vapour amount and temperature

void cSaturnModel::Saturation_Adjustment(std::string gas, 
    double &coeff_A, double &coeff_B, double &coeff_A_i, double &coeff_B_i, 
    double &t_0, double &t_00, 
    double &ep, double &lv, double &ls, double &cp, double &r, 
    double &C, double &L0, double &R, 
    double &del_alf, double &del_bet, double &m,
    Array &c, Array &cloud, Array &ice){

    cout << endl << "      SaturationAdjustment of "  << gas << " begin" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    cout.precision(9);

    bool sat_found = false;

    // ===== sat_found IS THE ONE SHARED DIAGNOSTIC HERE THAT WAS ACTUALLY RACING =====
    //
    // It is written `true` where a cell converges and `false` where it does not, from inside the
    // parallel region below, and it is NOT in that region's private() clause. Every thread wrote
    // it with no synchronisation, so the value that survived was whichever thread wrote last:
    // two runs of the same binary at the same thread count disagreed about whether saturation
    // was found at all, flipping the report between "saturation ... found" and "NO saturation
    // ... found". SaturationAdjustmentSat.h records this as point 8 of its own list, as "one
    // concrete, locatable contributor to ATSAT not reproducing itself run to run". It is.
    //
    // The OTHER shared diagnostics named there — i_sat, j_sat, k_sat, height_sat, saturation and
    // the *_sat pairs — were never racing under the OLD capture condition, because it fired only
    // at (j == 90 && k == 180), one collapsed (k,j) iteration and therefore one thread. That
    // condition is gone (see point 7 further down), so they are now captured per thread and
    // resolved the same way sat_found is.
    //
    // Serial semantics is "the last assignment in traversal order wins", so that is REPRODUCED
    // here rather than replaced. Each thread keeps the (key, value) of its own last assignment
    // and the largest key wins at the end. Within a thread the collapsed index k*jm + j ascends
    // and i ascends inside it, so keys are non-decreasing, and `>=` is what records the last
    // write for a given cell — which matters, because a cell may write `false` on several
    // convergence passes before writing `true` and breaking. Entries are strided so no two
    // threads share a cache line; this runs per convergence pass, not merely per cell.
    //
    // iter_prec_found IS ALSO FIXED, and it was a third bug rather than a variant of the other
    // two. It sits in the private() clause, so the copy printed after the region was never the
    // one the loop wrote: it always read 0, at any thread count, which is what the logs showed.
    // It is now carried out of the region in the capture record below and reports the converged
    // cell's actual pass count.
#ifdef _OPENMP
    const int sat_nthr = omp_get_max_threads();
#else
    const int sat_nthr = 1;
#endif
    const int sat_stride = 8;                       // 8 * sizeof(long long) = one cache line
    std::vector<long long> sat_tkey(sat_nthr * sat_stride, -1);
    std::vector<long long> sat_tval(sat_nthr * sat_stride,  0);

    // Per-thread capture of the last converged cell. Strided by the same amount so two threads
    // never share a cache line, and deliberately holding the RAW inputs rather than the derived
    // report: height, p_sat and the two latent differences are computed once after the region,
    // from the winner, so the arithmetic that builds the report exists in exactly one place.
    struct SatCapture {
        long long key;
        int    i, j, k, iter;
        double t_u, p_u, T, q_v_b, q_c_b, q_i_b, q_Rain;
    };
    std::vector<SatCapture> sat_trec(sat_nthr * sat_stride, SatCapture{-1,0,0,0,0,0,0,0,0,0,0,0});

    double saturation = 0.0;
    double height_sat = 0.0;
    double t_latent = 0.0;
    double p_latent = 0.0;

    int i_sat = 0;
    int j_sat = 0;
    int k_sat = 0;

    int iter_prec_end = 30;
    int iter_prec = 0;
    int iter_prec_found = 0;

    double q_v_hyp = 0.0;
    double p_u = 0.0;
    double t_u = 0.0;
    double T = 0.0;
    double d_t = 0.0;
    double E_Rain = 0.0;
    double E_Ice = 0.0;
    double q_Rain = 0.0;
    double q_Ice = 0.0;
    double q_v_b = 0.0;
    double q_c_b = 0.0;
    double q_i_b = 0.0;
    double q_v_b_sat = 0.0;
    double q_c_b_sat = 0.0;
    double q_i_b_sat = 0.0;
    double t_sat = 0.0;
    double p_sat = 0.0;
    double t_u_sat = 0.0;
    double p_u_sat = 0.0;
    double CND = 0.0;
    double DEP = 0.0;
    double d_q_v = 0.0;
    double d_q_c = 0.0;
    double d_q_i = 0.0;
    double q_diff = 0.0;
    double exp_pressure = g/(gam * R_ref);


// setting water vapour, cloud water and cloud ice into the proper thermodynamic ratio based on the local temperatures
// starting from a guessed parabolic temperature and water vapour distribution in north/south direction
    #pragma omp parallel for collapse(2) schedule(static) \
        private(iter_prec, iter_prec_found, q_v_hyp, p_u, t_u, T, d_t, \
                E_Rain, E_Ice, q_Rain, q_Ice, q_v_b, q_c_b, q_i_b, \
                CND, DEP, d_q_v, d_q_c, d_q_i, q_diff)
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
/*
    for(int k = 1; k < km-1; k++){
        for(int j = 1; j < jm-1; j++){
            for(int i = 1; i < im-1; i++){
*/
                t_u = t.x[i][j][k] * t_ref; // in K
                p_u = p_stat.x[i][j][k]; // in hPa

                if(t_u > t_0)  ice.x[i][j][k] = 0.0;

                if(c.x[i][j][k] < 0.0)  c.x[i][j][k] = 0.0;
                if(cloud.x[i][j][k] < 0.0)  cloud.x[i][j][k] = 0.0;
                if(ice.x[i][j][k] < 0.0)  ice.x[i][j][k] = 0.0;

                T = t_u; // in K
//                E_Rain = 1e3 * cNeptuneModel::Clausius_Clapeyron(T, coeff_A, coeff_B);  // saturation vapour pressure for the liquid phase at t > 0°C in bar
//                E_Ice = 1e3 * cNeptuneModel::Clausius_Clapeyron(T, coeff_A_i, coeff_B_i);  // saturation vapour pressure for the ice phase in bar

                E_Rain = ATPhys::saturation_vapour_pressure
                    (t_u, C, L0, R, del_alf, del_bet);
                E_Ice = ATPhys::saturation_vapour_pressure
                    (t_u, C, L0, R, del_alf, del_bet);

//                q_Rain = ep * E_Rain/(p_u - E_Rain); // relativ vapour contents in kg/m³
//                q_Ice = ep * E_Ice/(p_u - E_Ice); // relativ vapour contents in kg/m³

                q_Rain = r_mix * ep * E_Rain/p_u; // relativ vapour contents in kg/m³
                q_Ice = r_mix * ep * E_Ice/p_u; // relativ vapour contents in kg/m³
//                q_Rain = r * ep * E_Rain/p_u; // relativ vapour contents in kg/m³
//                q_Ice = r * ep * E_Ice/p_u; // relativ vapour contents in kg/m³

                q_v_b = c.x[i][j][k];  // in kg/m³
                q_c_b = cloud.x[i][j][k];
                q_i_b = ice.x[i][j][k];

                q_v_hyp = q_v_b;

                if(c.x[i][j][k] >= q_Rain){ // condition for cloud and ice formation, available vapor greater than at saturation

                    // §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§     iterations for mixed cloud phase     §§§§§§§§§§§§§§§§§§§§§

                    for(iter_prec = 1; iter_prec <= iter_prec_end; iter_prec++){ // iter_prec = 2 given by COSMO

                    // condensation ==> water vapor saturation for cloud water formation, deposition ==> ice crystal for cloud ice formation
                        CND = (T - t_00)/(t_0 - t_00);
                        DEP = (t_0 - T)/(t_0 - t_00);
                        if(T <= t_00){
                            CND = 0.0;
                            DEP = 1.0;
                        }
                        if(T >= t_0){
                            CND = 1.0;
                            DEP = 0.0;
                        }
                        d_q_v = q_v_hyp - q_v_b;  // changes in vapour causing cloud and ice
                        d_q_c = - d_q_v * CND;
                        d_q_i = - d_q_v * DEP;

//                        d_t = (lv * d_q_c + ls * d_q_i)/(cp * r_mix); // in K, temperature changes
                        d_t = (lv * d_q_c + ls * d_q_i)/(cp_mix * r_mix); // in K, temperature changes
//                        d_t = (lv * d_q_c + ls * d_q_i)/(cp * r); // in K, temperature changes
                        T = T + d_t; // in K

                        q_v_b = q_v_b + d_q_v;  // new values
                        q_c_b = q_c_b + d_q_c;
                        q_i_b = q_i_b + d_q_i;

                        if(q_c_b <= 0.0)  q_c_b = 0.0;
                        if(q_i_b <= 0.0)  q_i_b = 0.0;

//                        E_Rain = 1e3 * cNeptuneModel::Clausius_Clapeyron(T, coeff_A, coeff_B);  // saturation vapour pressure for the liquid phase at t > 0°C in bar
//                        E_Ice = 1e3 * cNeptuneModel::Clausius_Clapeyron(T, coeff_A_i, coeff_B_i);  // saturation vapour pressure for the ice phase in bar

                        E_Rain = ATPhys::saturation_vapour_pressure
                            (t_u, C, L0, R, del_alf, del_bet);
                        E_Ice = ATPhys::saturation_vapour_pressure
                            (t_u, C, L0, R, del_alf, del_bet);

//                        q_Rain = ep * E_Rain/(p_u - E_Rain); // relativ vapour contents in kg/m³
//                        q_Ice = ep * E_Ice/(p_u - E_Ice); // relativ vapour contents in kg/m³

                        q_Rain = r_mix * ep * E_Rain/p_u; // relativ vapour contents in kg/m³
                        q_Ice = r_mix * ep * E_Ice/p_u; // relativ vapour contents in kg/m³
//                        q_Rain = r * ep * E_Rain/p_u; // relativ vapour contents in kg/m³
//                        q_Ice = r * ep * E_Ice/p_u; // relativ vapour contents in kg/m³

                        if((q_c_b > 0.0)&&(q_i_b > 0.0))
                            q_v_hyp = (q_c_b * q_Rain + q_i_b * q_Ice) 
                                /(q_c_b + q_i_b);
                        if((q_c_b >= 0.0)&&(q_i_b == 0.0))  q_v_hyp = q_Rain;
                        if((q_c_b == 0.0)&&(q_i_b >= 0.0))  q_v_hyp = q_Ice;


                        if(T >= t_0) q_i_b = 0.0;

                        q_diff = fabs(q_v_b/q_v_hyp - 1.0);
                        if(q_diff <= 1.0e-4){
                            {   // last assignment in traversal order wins — see the note above
#ifdef _OPENMP
                                const int sat_x = omp_get_thread_num() * sat_stride;
#else
                                const int sat_x = 0;
#endif
                                const long long sat_k = ((long long)k * jm + j) * im + i;
                                if(sat_k >= sat_tkey[sat_x]){
                                    sat_tkey[sat_x] = sat_k;
                                    sat_tval[sat_x] = 1;
                                }
                            }
                            iter_prec_found = iter_prec;

                            // ===== POINT 7: CAPTURE THE LAST CONVERGED CELL, NOT A MEANINGLESS ONE =====
                            //
                            // The condition here used to be
                            //
                            //     if((j == 90) && (k == 180) && (i == iter_prec_found))
                            //
                            // which SaturationAdjustmentSat.h calls meaningless in point 7 of its
                            // own list, and it is: i is a GRID INDEX and iter_prec_found is a
                            // CONVERGENCE-PASS COUNTER, so which cell got reported depended on how
                            // many passes that cell happened to take. It was pinned to a single
                            // column besides (j = 90, k = 180), so the report described one point
                            // on one meridian and was printed as if it described the model.
                            //
                            // It now captures the LAST CONVERGED CELL in traversal order, which is
                            // what ATJUP and the shared SaturationAdjustment<Planet> report and what
                            // point 7 names as the behaviour to mirror. THIS CHANGES THE REPORTED
                            // NUMBERS, deliberately: no output file moves, but these log lines do,
                            // and they should — they were not reporting what they claimed to.
                            //
                            // No lock is needed. Within a thread the collapsed index k*jm + j
                            // ascends and i ascends inside it, so a thread's latest converged cell
                            // is simply its most recent one and it can overwrite its own slot
                            // unconditionally. The largest key across threads wins after the region.
                            {
#ifdef _OPENMP
                                const int sat_c = omp_get_thread_num() * sat_stride;
#else
                                const int sat_c = 0;
#endif
                                SatCapture &rec = sat_trec[sat_c];
                                rec.key    = ((long long)k * jm + j) * im + i;
                                rec.i      = i;
                                rec.j      = j;
                                rec.k      = k;
                                rec.t_u    = t_u;
                                rec.p_u    = p_u;
                                rec.T      = T;
                                rec.q_v_b  = q_v_b;
                                rec.q_c_b  = q_c_b;
                                rec.q_i_b  = q_i_b;
                                rec.q_Rain = q_Rain;
                                rec.iter   = iter_prec;
                            }
                            break;
                        }
                        else{
                            q_v_hyp = 0.5 * (q_v_hyp + q_v_b);  // has smoothing effect
                            {   // last assignment in traversal order wins — see the note above
#ifdef _OPENMP
                                const int sat_x = omp_get_thread_num() * sat_stride;
#else
                                const int sat_x = 0;
#endif
                                const long long sat_k = ((long long)k * jm + j) * im + i;
                                if(sat_k >= sat_tkey[sat_x]){
                                    sat_tkey[sat_x] = sat_k;
                                    sat_tval[sat_x] = 0;
                                }
                            }
                        }                            


/*
                        cout.precision(10);
                        cout.setf(ios::fixed);
                        if((j == 90)&&(k == 180))  cout << endl
                            << "  SaturationAdjustment of " << gas << endl 
                            << "  height_sat = "<< height_sat << endl
                            << "  saturation = "<< saturation * 1e3 << endl
                            << "  i = " << i << "  j = " << j << "  k = " << k << endl
                            << "  iprec = "<< iter_prec << endl
                            << "  CND = " << CND
                            << "  DEP = " << DEP << endl
                            << "  p_u = " << p_u
//                            << "  p_dyn = " << p_dyn.x[i][j][k] << endl
//                            << "  r_dry = " << r_dry.x[i][j][k] 
//                            << "  r_humid = " << r_humid.x[i][j][k] << endl
                            << "  dt = " << d_t 
                            << "  t_latent = " << t_latent
                            << "  T = " << T - t_0 
                            << "  t = " << t_u - t_0 << endl
                            << "  q_Rain = " << q_Rain * 1e3 
                            << "  q_Ice = " << q_Ice * 1e3 << endl
                            << "  q_diff = " << q_diff << endl
                            << "  d_q_v = " << d_q_v * 1e3 
                            << "  d_q_c = " << d_q_c * 1e3 
                            << "  d_q_i = " << d_q_i * 1e3 << endl
                            << "  q_v_hyp = " << q_v_hyp * 1e3
                            << "  q_v_b = " << q_v_b * 1e3 
                            << "  q_c_b = " << q_c_b * 1e3 
                            << "  q_i_b = " << q_i_b * 1e3 << endl
                            << "  c = " << c.x[i][j][k] * 1e3 
                            << "  cloud = " << q_v_b * 1e3 
                            << "  ice = " << q_c_b * 1e3 << endl;
*/

                    } // iter_prec end

                    c.x[i][j][k] = q_v_b;  // new values achieved after converged iterations
                    cloud.x[i][j][k] = q_c_b;
                    ice.x[i][j][k] = q_i_b;

                    t.x[i][j][k] = T/t_ref;

                    p_stat.x[i][j][k] = p_ref 
                        * pow((t.x[i][j][k]), exp_pressure); // in bar

                } // iterations for mixed cloud phase
//                sat_adjust = true;

            } // end i
        } // end j
    } // end k

    // Resolve sat_found: the largest key across threads is the last assignment in traversal
    // order, which is exactly the value a single thread would have been left holding. If no cell
    // assigned at all every key is still -1 and sat_found keeps its initial false, which is also
    // what the serial loop would have done.
    {
        long long sat_best = -1;
        for(int t = 0; t < sat_nthr; t++){
            const int sat_x = t * sat_stride;
            if(sat_tkey[sat_x] > sat_best){
                sat_best  = sat_tkey[sat_x];
                sat_found = (sat_tval[sat_x] != 0);
            }
        }
    }

    // And the capture: the largest key is the last converged cell in traversal order. The report
    // is derived here, once, from that cell's raw values.
    {
        long long cap_best = -1;
        const SatCapture *win = nullptr;
        for(int t = 0; t < sat_nthr; t++){
            const SatCapture &rec = sat_trec[t * sat_stride];
            if(rec.key > cap_best){ cap_best = rec.key; win = &rec; }
        }
        if(win != nullptr){
            i_sat = win->i;
            j_sat = win->j;
            k_sat = win->k;
            iter_prec_found = win->iter;
            height_sat = get_layer_height(i_sat);

            t_u_sat = win->t_u;
            p_u_sat = win->p_u;

            t_sat = win->T;
            p_sat = p_ref * pow(t_sat/t_ref, exp_pressure);

            t_latent = t_sat - t_u_sat;
            p_latent = p_sat - p_u_sat;

            q_v_b_sat = win->q_v_b;
            q_c_b_sat = win->q_c_b;
            q_i_b_sat = win->q_i_b;

            saturation = win->q_v_b - win->q_Rain;
        }
    }

    #pragma omp parallel for collapse(3) schedule(static)
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                if(c.x[i][j][k] <= 0.0) c.x[i][j][k] = 0.0;
                if(cloud.x[i][j][k] <= 0.0) cloud.x[i][j][k] = 0.0;
                if(ice.x[i][j][k] <= 0.0) ice.x[i][j][k] = 0.0;
            }
        }
    }

    if(sat_found == false)  
        cout << "      NO saturation of water vapour in SaturationAdjustment of " << gas << " found" 
        << endl
        << "      iter_prec_end = " << iter_prec_end
        << "   iter_prec = " << iter_prec << endl
        << endl;
    else
        cout << "      saturation of water vapour in SaturationAdjustment of " << gas << " found"
        << endl
        << "      iter_prec_found = " << iter_prec_found
        << "   iter_prec_end = " << iter_prec_end
        << "   iter_prec = " << iter_prec << endl

        << "      i_sat = " << i_sat
        << "   j_sat = " << j_sat
        << "   k_sat = " << k_sat
        << "   height_sat[km] = " << height_sat << endl

        << "      p_stat[bar] = " << p_sat
        << "   p_u[bar] = " << p_u_sat
        << "   p_latent[bar] = " <<  p_latent << endl

        << "      T[°C] = " << t_sat - t_ref
        << "   t_u[°C] = " << t_u_sat - t_ref
        << "   t_latent[°C] = " <<  t_latent << endl

        << "      saturation[g/m³] = " << saturation * 1e3 << endl

        << "      " << gas << " humid solution[g/m³] = " << q_v_b_sat * 1e3
        << "   " << gas << "-cloud[g/m³] = " << q_c_b_sat * 1e3
        << "   " << gas << "-ice[g/m³] = " << q_i_b_sat * 1e3 << endl;


    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for SaturationAdjustment\n", elapsed.count() * 1e-9);

    cout << "      SaturationAdjustment of " << gas << " ended" << endl;
    return;
}
/*
*
*/
// One-Category-Ice-Scheme, COSMO-module from the German Weather Forecast, 
// resulting the precipitation distribution formed of rain and snow
/*
void cNeptuneModel::OneCategoryIceScheme(){   
    cout << endl << "      OneCategoryIceScheme" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // constant coefficients for the transport of cloud water and cloud ice amount, 
    // rain and snow in the parameterization procedures
    bool rain = false;
    bool snow = false;
    double P_Rain = 0.0;
    double P_Snow = 0.0;
    double Rain = 0.0;
    double Snow = 0.0;
    int i_rain = 0;
    int i_snow = 0;
    int j_rain = 0;
    int j_snow = 0;
    int k_rain = 0;
    int k_snow = 0;
    double height_rain = 0.0;
    double height_snow = 0.0;
    double maxValue_rain = 0.0;
    double maxValue_snow = 0.0;
    double N_r_0 = 8.0e+6,  // in 1/m-4
           N_s_0 = 4.0e+5,  // in 1/m-4
           v_s_0 = 130.0,  // in m/s
           v_r_0 = 4.9,  // in m/s
           c_r_t = 12.63,  // in 
           c_s_t = 2.87,  // in 
           a_s_m = 0.038,  // in kg/m2
           c_ac = 0.24,  // m2/kg
           c_rim = 0.69,  // m2/kg
           b_ev = 5.98,  // m2*s/kg
           a_melt = 3.90e-6, // K/(kg/kg)
           b_melt = 10.50,  // m2*s/kg
           b_dep = 10.50,  // m2*s/kg
           a_if = 1.92e-6,
           a_cf = 3.97e-5,
           E_cf = 5.0e-3,
           N_cf_0_surf = 2.0e5,  // 1/m3
           N_cf_0_top = 1.0e4,  // 1/m3
//           tau_r = 1.0e4,  // s
           tau_r = 3.3e3,  // s          can be adjusted to fit best the average NASA-precipitation of 2.68 mm/d
           tau_s = 1.0e3,  // s
           a_mc = 0.08,  // kg/m2
           a_mv = 0.02;  // kg/m2
    double N_cf_0 = 0.0;
    double eps_t = 0.0;
    double a_m = 0.0;
    double a_ev = 0.0;
    double a_dep = 0.0;
    double N_cf = 0.0;
    double t_m1 = 0.5 * (t_0_h2o + t_000);
    double t_m2 = 0.5 * (t_0_h2o + t_00_h2o);
    double q_Rain = 0.0,
           E_Rain = 0.0;
    double q_Ice = 0.0,
           E_Ice = 0.0;
    double A_r = r_0_water * M_PI * N_r_0;
    double A_s = 2.0 * a_s_m * N_s_0;
    double B_r = r_0_water * M_PI * N_r_0 * v_r_0 * tgamma(4.5)/tgamma(4.0); 
    double B_s = N_s_0 * v_s_0 * a_s_m * tgamma(3.25);
    double S_nuc, S_frz, S_cf_frz, S_if_frz, S_dep, S_au, S_ac, S_rim, S_shed, S_ev, S_melt;
    double r_q_r, r_q_s, v_r_t, v_s_t;
    std::vector<double> step(im, 0.0);
    // rain and snow distribution based on parameterization schemes adopted from the COSMO code used by the German Weather Forecast
    // the choosen scheme is a Two Category Ice Scheme
    // besides the transport equation for the water vapour exist two equations for the cloud water and the cloud ice transport
    // since the diagnostic version of the code is applied the rain and snow mass transport is computed by column equilibrium integral equation

//    #pragma omp parallel for
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                if(h2o.x[i][j][k] < 0.0)  h2o.x[i][j][k] = 0.0;
                if(h2o_cloud.x[i][j][k] < 0.0)  h2o_cloud.x[i][j][k] = 0.0;
                if(h2o_ice.x[i][j][k] < 0.0)  h2o_ice.x[i][j][k] = 0.0;
                if(P_rain.x[i][j][k] < 0.0)  P_rain.x[i][j][k] = 0.0;
                if(P_snow.x[i][j][k] < 0.0)  P_snow.x[i][j][k] = 0.0;
                P_rain.x[i][j][k] = 0.0;
                P_snow.x[i][j][k] = 0.0;
            }
        }
    }
    for(int k = 1; k < km-1; k++){
        for(int j = 1; j < jm-1; j++){
            P_rain.x[im-1][j][k] = 0.0;
            P_snow.x[im-1][j][k] = 0.0;
            P_rain.x[im-2][j][k] = 0.0;
            P_snow.x[im-2][j][k] = 0.0;
            S_r.x[im-1][j][k] = 0.0;
            S_s.x[im-1][j][k] = 0.0;
            S_r.x[im-2][j][k] = 0.0;
            S_s.x[im-2][j][k] = 0.0;
            for(int i = im-2; i >= 0; i--){
                Rain = P_rain.x[i][j][k];
                Snow = P_snow.x[i][j][k];
                r_q_r = A_r * pow(B_r,-(8.0/9.0)) 
                    * pow(P_rain.x[i][j][k],(8.0/9.0));
                r_q_s = A_s * pow(B_s,-(12.0/13.0)) 
                    * pow(P_snow.x[i][j][k],(12.0/13.0));
                v_r_t = c_r_t * pow(r_q_r, (1.0/8.0));
                v_s_t = c_s_t * pow(r_q_s, (1.0/12.0));
                double t_u = t.x[i][j][k] * t_0_h2o;
                step[i] = get_layer_height(i) - get_layer_height(i-1);  // local atmospheric shell thickness
//                E_Rain = hp * exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
//                E_Ice = hp * exp_func(t_u, 21.8746, 7.66);

                E_Rain = cNeptuneModel::saturation_vapour_pressure(t_u, 
                    C_h2o, L0_h2o, R_h2o, del_alf_h2o, del_bet_h2o);
                E_Ice = cNeptuneModel::saturation_vapour_pressure(t_u, 
                    C_h2o, L0_h2o, R_h2o, del_alf_h2o, del_bet_h2o);
                q_Rain = ep * E_Rain/(p_stat.x[i][j][k] - E_Rain); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
                q_Ice = ep * E_Ice/(p_stat.x[i][j][k] - E_Ice); // relativ water vapour contents on ocean surface reduced by factor in kg/kg
    // mass size relation of circular plates
                if((t_u < t_0_h2o)&&(t_u >= t_m1))
                    a_m = a_mc - a_mv * (1.0 + cos(2.0 * M_PI 
                        * (t_u - t_m1)/(t_0_h2o - t_000)));
                else  a_m = a_mc;
    // autoconversion and nucleation process, epsilon(T)
                if((t_u > t_00_h2o)&&(t_u <= t_0_h2o))
                    eps_t = 0.5 * (1.0 + sin(M_PI*(t_m2 - t_u)/(t_0_h2o - t_00_h2o)));
                else
                    eps_t = 0.0;
//                if (h2o_cloud.x[i][j][k] >= 0.0){
//                    S_au = (1.0 - eps_t)/tau_r * h2o_cloud.x[i][j][k];
//                    S_nuc = eps_t/tau_s * h2o_cloud.x[i][j][k];
                    S_au = (1.0 - eps_t)/tau_r * max(0.0, h2o_cloud.x[i][j][k]);
                    S_nuc = eps_t/tau_s * max(0.0, h2o_cloud.x[i][j][k]);
//                }
    // collection mechanisms accretion, riming and shedding
    // accretion of cloud water by raindrops
                S_ac = (1.0 - eps_t) * c_ac * h2o_cloud.x[i][j][k]  // c_ac = 0.24, in m2/kg
                    * pow(Rain,(7.0/9.0));
                if(t_u < t_0_h2o)
                    S_rim = c_rim/a_m * h2o_cloud.x[i][j][k] * Snow;  // c_rim = 18.6,  m2/kg
                else  S_rim = 0.0;  // riming rate of snow mass due to collection of supercooled cloud droplets, < VIII >
                                   // by falling snow particles
                if(t_u >= t_0_h2o)  
                    S_shed = c_rim/a_m * h2o_cloud.x[i][j][k] * Snow;
                else  S_shed = 0.0;  // rate of water shed by melting wet snow particles, < IX >
                                    // collecting cloud droplets to produce rain
    // diffusional growth of rain and snow
    // evaporation of rain water
                a_ev = 2.76e-3 * exp(0.055 * (t_0_h2o - t_u));
                S_ev = a_ev * (1.0 + b_ev * pow(Rain,(1.0/6.0)))  // evaporation of rain due to water vapour diffusion, < XIII >
                    * (q_Rain - h2o.x[i][j][k]) 
                    * pow(Rain, (4.0/9.0));
    // deposition growth and sublimation of cloud ice
                a_dep = 1.13e-3 * exp(0.073 * (t_0_h2o - t_u));
                S_dep = a_dep/pow(a_m, - 0.5) * (1.0 + b_dep   // melting rate of snow to form rain, < XVI >
                    * pow(a_m, - 0.25) * pow(Snow, (0.9/4.3)))  // c_s_melt = 8.43e-5, (m2*s)/(K*kg)
                    * ((h2o.x[i][j][k] - q_Ice) * pow(Snow, (5.0/8.6)));
    // melting of snow to form cloud water
                S_melt = a_melt/pow(a_mc, - 0.5) * (1.0 + b_melt   // melting rate of snow to form rain, < XVI >
                    * pow(a_mc, - 0.25) * pow(Snow, (0.9/4.3)))  // c_s_melt = 8.43e-5, (m2*s)/(K*kg)
                    * ((t_u - t_0_h2o) * pow(Snow, (5.0/8.6)));
    // freezing of rain to form snow
                N_cf_0 = N_cf_0_surf + (N_cf_0_surf - N_cf_0_top)
                    /(p_stat.x[0][j][k] - 500.0) 
                    * (p_stat.x[i][j][k] - p_stat.x[0][j][k]);
                if(p_stat.x[i][j][k] <= 500.0)  N_cf_0 = N_cf_0_top;
                if(t_u < 270.16){
                    N_cf = N_cf_0 * pow((270.16 - t_u), 1.3);
                }
                else  N_cf = 0.0;
                if((t_u >= t_0_h2o)&&(t_u <= t_00_h2o)){
                    S_if_frz = a_if * (exp(a_if * (t_u - t_0_h2o)) - 1.0)  // immersion freezing
                        * pow(Rain, (14.0/9.0));
                    S_cf_frz = a_cf * E_cf * N_cf                       // contact freezing nucleation
                        * pow(Rain, (13.0/9.0));
                    S_frz = S_if_frz + S_cf_frz;                        // mixing ratio of snow
                }
                else  S_frz = 0.0;
    // sinks and sources
                S_v.x[i][j][k] = - S_c_c.x[i][j][k] + S_ev - S_dep;
                S_c.x[i][j][k] = S_c_c.x[i][j][k] - S_au - S_ac 
                    - S_nuc - S_rim - S_shed;
                S_r.x[i][j][k] = S_au + S_ac - S_ev + S_shed 
                    - S_frz + S_melt;
                S_s.x[i][j][k] = S_nuc + S_rim + S_dep + S_frz - S_melt;
    // rain and snow integration
                if(t_u >= t_0_h2o)
                    P_rain.x[i][j][k] = P_rain.x[i+1][j][k]
                        + r_humid.x[i+1][j][k] * S_r.x[i+1][j][k] 
                        * step[i];  // in kg/(m2 * s) == mm/s 
                else  P_rain.x[i][j][k] = 0.0;
                if(P_rain.x[i][j][k] >= 15.0/8.64e4)  P_rain.x[i][j][k] = 15.0/8.64e4;
                if(P_rain.x[i][j][k] < 0.0)  P_rain.x[i][j][k] = 0.0;
                if(P_rain.x[i][j][k] > 0.0){
                    rain = true;
                    if(P_rain.x[i][j][k] > maxValue_rain){
                        maxValue_rain = P_rain.x[i][j][k];
                        P_Rain = maxValue_rain * 8.64e4; // in mm/d
                        i_rain = i;
                        j_rain = j;
                        k_rain = k;
                        height_rain = get_layer_height(i);
                    }
                }
                if((t_u < t_0_h2o)&&(t_u >= t_00_h2o))
                    P_snow.x[i][j][k] = P_snow.x[i+1][j][k]
                        + r_humid.x[i+1][j][k] * S_s.x[i+1][j][k] 
                        * step[i];  // in kg/(m2 * s) == mm/s 
                else  P_snow.x[i][j][k] = 0.0;
                if(P_snow.x[i][j][k] < 0.0)  P_snow.x[i][j][k] = 0.0;
                if(P_snow.x[i][j][k] > 0.0){
                    snow = true;
                    if(P_snow.x[i][j][k] > maxValue_snow){
                        maxValue_snow = P_snow.x[i][j][k];
                        P_Snow = maxValue_snow * 8.64e4; // in mm/d
                        i_snow = i;
                        j_snow = j;
                        k_snow = k;
                        height_snow = get_layer_height(i);
                    }
                }
            }  // end i RainSnow
        }  // end j
    }  // end k

    #pragma omp parallel for
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            P_rain.x[im-1][j][k] = c43 * P_rain.x[im-2][j][k] 
                - c13 * P_rain.x[im-3][j][k];
            P_snow.x[im-1][j][k] = c43 * P_snow.x[im-2][j][k] 
                - c13 * P_snow.x[im-3][j][k];
                }
            }

    #pragma omp parallel for
    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
            P_rain.x[i][0][k] = c43 * P_rain.x[i][1][k] 
                - c13 * P_rain.x[i][2][k];
            P_rain.x[i][jm-1][k] = c43 * P_rain.x[i][jm-2][k] 
                - c13 * P_rain.x[i][jm-3][k];
            P_snow.x[i][0][k] = c43 * P_snow.x[i][1][k] 
                - c13 * P_snow.x[i][2][k];
            P_snow.x[i][jm-1][k] = c43 * P_snow.x[i][jm-2][k] 
                - c13 * P_snow.x[i][jm-3][k];
        }
    }

    #pragma omp parallel for
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            P_rain.x[i][j][0] = c43 * P_rain.x[i][j][1] 
                - c13 * P_rain.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            P_rain.x[i][j][km-1] = c43 * P_rain.x[i][j][km-2] 
                - c13 * P_rain.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
            P_rain.x[i][j][0] = P_rain.x[i][j][km-1] 
                = (P_rain.x[i][j][0] + P_rain.x[i][j][km-1])/2.0;
            P_snow.x[i][j][0] = c43 * P_snow.x[i][j][1] 
                - c13 * P_snow.x[i][j][2];  // von Neumann boundary condition dt/dphi = 0.0
            P_snow.x[i][j][km-1] = c43 * P_snow.x[i][j][km-2] 
                - c13 * P_snow.x[i][j][km-3];  // von Neumann boundary condition dt/dphi = 0.0
        }
    }

//    #pragma omp parallel for
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_mount = i_topography[j][k];
            for(int i = i_mount; i > 0; i--){
                if((is_land(SeaMount, i, j, k))&&(i < i_mount)){
                    P_rain.x[i][j][k] = 0.0;
                    P_snow.x[i][j][k] = 0.0;
                }
                P_rain.x[0][j][k] = P_rain.x[i_mount][j][k];
                P_snow.x[0][j][k] = P_snow.x[i_mount][j][k];
            }
        }
    }
    if(rain == false)  
        cout << "      no rain fall in OneCategoryIceScheme found" 
        << endl;
    else
        cout << "      rain fall in OneCategoryIceScheme found" 
        << endl
        << "      i_rain = " << i_rain
        << "      j_rain = " << j_rain
        << "   k_rain = " << k_rain
        << "   height_rain[m] = " << height_rain
        << "   P_Rain[mm/d] = " << P_Rain << endl;
    if(snow == false)  
        cout << "      no snow fall in OneCategoryIceScheme found" 
        << endl;
    else
        cout << "      snow fall in OneCategoryIceScheme found" 
        << endl
        << "      i_snow = " << i_snow
        << "      j_snow = " << j_snow
        << "   k_snow = " << k_snow
        << "   height_snow[m] = " << height_snow
        << "   P_Snow[mm/d] = " << P_Snow << endl;
    cout << "      OneCategoryIceScheme ended" << endl;

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for OneCategoryIceScheme\n", elapsed.count() * 1e-9);

    return;
}
*/
/*
*
*/
double cSaturnModel::Clausius_Clapeyron(double &T_K, double &A, double &B){
        return exp(A/T_K + B);  // temperature in °K   result pressure in bar
}
/*
*
*/
double cSaturnModel::Humility_critical(double &x, double Hu_cr_max, 
    double Hu_cr_mid){
    return (Hu_cr_max - Hu_cr_mid) * (x * x - 2.0 * x) 
        + Hu_cr_max;
    }
/*
*
*/
/*
*
*/
