/*
 * Saturn Atmosphere Circulation Model (ATSAT)
 * Dry convective adjustment — see the long note in ConvectiveAdjustmentSat.h.
 */

#include "cSaturnModel.h"
#include "ConvectiveAdjustmentSat.h"
#include "Utils.h"

using namespace std;

void ConvectiveAdjustmentSat::run(){
    cout << endl << "      ATSAT: ConvectiveAdjustmentSat (dry, Manabe-Strickler)" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // ---- The critical lapse rate, as a temperature drop per grid layer ----
    //
    // dT/dz = -g/cp is the dry adiabat; over one layer of thickness dz that is a drop of
    // (g/cp)*dz kelvin, and the model's temperature is stored as T/t_ref, so the number compared
    // against below is (g/cp)*dz/t_ref. cp is the mixture value the model computes for itself in
    // ChemistrySat::ThermalPropertiesSat, which runs once during initialisation and is therefore
    // filled before the first call here. Do not substitute a textbook cp: in ATJUP the difference
    // between the model's own cp_mix and a plausible H2/He value moved the criterion by 0.26 K/km,
    // which was larger than the entire superadiabatic excess that model developed.
    //
    // ATSAT_CONV_ADJ_LAPSE scales it, for asking what a different critical lapse rate would do —
    // 0 gives an isothermal criterion, values above 1 make the scheme stricter than the dry
    // adiabat, which is one crude way to stand in for a moist adiabat in a condensing region.
    const double dz_m       = (m.L_atm * 1.0e3) / double(m.im - 1);
    const double lapse_fac  = [](){
        const char* e = getenv("ATSAT_CONV_ADJ_LAPSE"); return e ? atof(e) : 1.0; }();
    const double dT_ad_nd   = lapse_fac * (m.g / m.cp_mix) * dz_m / m.t_ref;

    // A column is left alone unless it is superadiabatic by more than this, so that round-off
    // does not make the scheme fire on a column that is already neutral.
    constexpr double tol_nd = 1.0e-12;

    // Sweeps are repeated until the column is stable. A deep unstable layer needs one sweep per
    // layer in the worst case; the cap only exists so a pathological column cannot spin here.
    static const int max_pass = [](){
        const char* e = getenv("ATSAT_CONV_ADJ_PASSES");
        const int v = e ? atoi(e) : 64;
        return v > 0 ? v : 64; }();

    // ---- Diagnostics, reduced across the whole grid ----
    long long n_columns_adjusted = 0;    // columns that needed at least one sweep
    long long n_layers_mixed     = 0;    // individual layers put back on the adiabat
    int       max_passes_used    = 0;    // worst column, to show whether the cap was reached
    double    max_dT_K           = 0.0;  // largest temperature change applied, in kelvin
    double    max_rel_drift      = 0.0;  // worst enthalpy non-conservation, relative

    std::vector<double> w(m.im), t_col(m.im);

    #pragma omp parallel for collapse(2) schedule(dynamic, 8) firstprivate(w, t_col) \
        reduction(+:n_columns_adjusted, n_layers_mixed) \
        reduction(max:max_passes_used, max_dT_K, max_rel_drift)
    for(int j = 0; j < m.jm; j++){
        for(int k = 0; k < m.km; k++){

            // ---- The fluid part of this column ----
            // Solid cells hold boundary values, not a fluid state, so the column starts above the
            // topography and stops at the first solid cell above it (there should be none, but a
            // column that is walled in higher up must not be mixed across the wall).
            int i0 = 0;
            while(i0 < m.im && m.SeaMount.x[i0][j][k] == 1.0) i0++;
            int i1 = i0;
            while(i1 + 1 < m.im && m.SeaMount.x[i1 + 1][j][k] != 1.0) i1++;
            if(i1 - i0 < 1) continue;                       // nothing to mix

            // ---- Mass weights ----
            // The enthalpy per unit area of a layer is (cp/g) * dp, so the conserved quantity is
            // the dp-weighted temperature. dp is taken from the hydrostatic p_stat with faces at
            // the midpoints, and the two end layers get their half-cell. Using dp rather than a
            // density avoids a trap: rho = p/(R*T) makes rho*T identically p/R, so a density
            // weight would conserve nothing at all.
            bool weights_ok = true;
            for(int i = i0; i <= i1; i++){
                const double p_lo = (i > i0) ? 0.5 * (m.p_stat.x[i-1][j][k] + m.p_stat.x[i][j][k])
                                             : m.p_stat.x[i0][j][k];
                const double p_hi = (i < i1) ? 0.5 * (m.p_stat.x[i][j][k] + m.p_stat.x[i+1][j][k])
                                             : m.p_stat.x[i1][j][k];
                w[i] = p_lo - p_hi;
                if(!(w[i] > 0.0) || !std::isfinite(w[i])) weights_ok = false;
                t_col[i] = m.t.x[i][j][k];
                if(!std::isfinite(t_col[i]) || t_col[i] <= 0.0) weights_ok = false;
            }
            // A column with a non-monotonic p_stat or a bad temperature is left untouched: this
            // routine is not the place to repair either, and mixing across a NaN would spread it
            // through the whole column.
            if(!weights_ok) continue;

            double sum_before = 0.0;
            for(int i = i0; i <= i1; i++) sum_before += w[i] * t_col[i];

            // ---- Sweeps ----
            //
            // Whole unstable SEGMENTS are mixed at once, not adjacent pairs. Both converge to the
            // same profile, but pairwise mixing moves heat one layer per sweep, so a deep unstable
            // block needs as many sweeps as it has layers: in ATJUP a test with an artificially
            // strict criterion ran 137 million pair adjustments and still hit the 64-sweep cap.
            // Mixing the segment settles it in one step.
            //
            // A segment [a..b] is put on the adiabat, T_q = C - dT_ad*(q-a), with C chosen so the
            // dp-weighted temperature of the segment is unchanged:
            //     C = [ sum w_q T_q + dT_ad * sum w_q (q-a) ] / sum w_q
            // Mixing can destabilise the joint with the layer below or above, so the segment is
            // grown in whichever direction is still too steep and re-mixed. It can only grow, and
            // only within the column, so the inner loop terminates.
            int passes = 0;
            long long layers_mixed = 0;
            bool changed = true;
            while(changed && passes < max_pass){
                changed = false;
                passes++;
                int i = i0;
                while(i < i1){
                    if(t_col[i] - t_col[i+1] <= dT_ad_nd + tol_nd){ i++; continue; }

                    int a = i, b = i + 1;
                    for(;;){
                        double sw = 0.0, swt = 0.0, swk = 0.0;
                        for(int q = a; q <= b; q++){
                            sw  += w[q];
                            swt += w[q] * t_col[q];
                            swk += w[q] * double(q - a);
                        }
                        const double C = (swt + dT_ad_nd * swk) / sw;
                        for(int q = a; q <= b; q++){
                            const double t_new = C - dT_ad_nd * double(q - a);
                            const double dK = std::fabs(t_new - t_col[q]) * m.t_ref;
                            if(dK > max_dT_K) max_dT_K = dK;
                            t_col[q] = t_new;
                        }

                        bool extended = false;
                        if(a > i0 && t_col[a-1] - t_col[a] > dT_ad_nd + tol_nd){ a--; extended = true; }
                        if(b < i1 && t_col[b] - t_col[b+1] > dT_ad_nd + tol_nd){ b++; extended = true; }
                        if(!extended) break;
                    }

                    layers_mixed += (b - a + 1);
                    changed = true;
                    i = b;                       // carry on above the block just mixed
                }
            }

            if(layers_mixed == 0) continue;

            double sum_after = 0.0;
            for(int i = i0; i <= i1; i++) sum_after += w[i] * t_col[i];
            const double drift = (sum_before != 0.0)
                               ? std::fabs(sum_after - sum_before) / std::fabs(sum_before) : 0.0;
            if(drift > max_rel_drift) max_rel_drift = drift;

            for(int i = i0; i <= i1; i++) m.t.x[i][j][k] = t_col[i];

            n_columns_adjusted++;
            n_layers_mixed += layers_mixed;
            if(passes > max_passes_used) max_passes_used = passes;
        }
    }

    const double frac = 100.0 * double(n_columns_adjusted) / double(m.jm * m.km);
    printf("      ATSAT: convective adjustment — critical drop %.3f K per %.1f km layer;"
           " %lld of %d columns (%.2f %%), %lld layers, worst column %d sweeps of %d,"
           " max dT %.3f K, enthalpy drift %.2e\n",
           dT_ad_nd * m.t_ref, dz_m * 1.0e-3,
           n_columns_adjusted, m.jm * m.km, frac, n_layers_mixed,
           max_passes_used, max_pass, max_dT_K, max_rel_drift);
    if(max_passes_used >= max_pass)
        cout << "      ATSAT: WARNING - the sweep cap was reached, a column may still be "
                "superadiabatic (raise ATSAT_CONV_ADJ_PASSES)" << endl;

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for ConvectiveAdjustmentSat\n", elapsed.count() * 1e-9);
    cout << "      ATSAT: ConvectiveAdjustmentSat ended" << endl;
}
