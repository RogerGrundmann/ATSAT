#pragma once

#include "cSaturnModel.h"
#include "ATPhys.h"
#include "FluxLimiter.h"

#include <cmath>
#include <chrono>
#include <cstdio>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

class ChemistrySat {
public:
    explicit ChemistrySat(cSaturnModel& model) : m(model) {}

    // -----------------------------------------------------------------------

    // ATSAT_CHEM_ENTHALPY — what the second term of thermalmassflux is. Default 0 = unchanged.
    //
    // The term reads  t * (w_nh3*k_nh3 + w_h2s*k_h2s + w_nh4sh*k_nh4sh), and it is dimensionally
    // impossible: k_* is a THERMAL CONDUCTIVITY, W/(m K), and w_* a mass production rate,
    // kg/(m3 s), so the product is kg W/(m4 s K) while the term it is added to is W/m3. They
    // differ by kg/(m s). The stray factor of t is wrong too — no enthalpy source scales with
    // absolute temperature.
    //
    // What it was meant to be is written in the code's own comment: "source term by chemical
    // reaction (reaction enthalpy)". That is  w * dh  with dh a SPECIFIC ENTHALPY in J/kg, which
    // gives W/m3. =1 makes it that, using dh_nh4sh.
    //
    // ONE species, not three: w_nh3, w_h2s and w_nh4sh are the SAME reaction expressed in three
    // masses (the chemistry forms them as -m_nh3*R, -m_h2s*R, +m_nh4sh*R from one rate R), so
    // summing three enthalpy terms would count the same heat three times. The source is tied to
    // NH4SH's production rate alone.
    //
    // THIS ONE CHANGES THE ANSWER, unlike the other units work: thermalmassflux is read by
    // RHS_*, not only by the writers. Hence a knob, off by default, and a measurement.
    static int chem_enthalpy(){
        static const int v = [](){ const char* e = getenv("ATSAT_CHEM_ENTHALPY"); return e ? atoi(e) : 0; }();
        return v;
    }

    void ChemMassRateSat()
    {
        using namespace std;
        cout << endl << "      ATSAT: ChemMassRateSat" << endl;

        auto begin = chrono::high_resolution_clock::now();

        const int im = m.im, jm = m.jm, km = m.km;

        // reaction constants are the same for all cells.
        //
        // NOTE A and B DIFFER FROM ATJUP'S — 1500/-0.3 here against 15000/-0.5 there, an order of
        // magnitude in the prefactor, with the same dissociation temperature. Neither file records
        // which is Saturn's and which is inherited. That is a calibration question, it is NOT
        // touched by either knob below, and it has to be settled before the two chemistries can
        // share an implementation.
        const double keq = m.m_nh4sh / (m.m_nh3 * m.m_h2s);
        const double A   = 1500.0;
        const double B   = -0.3;
        const double T_d = 3020.0;

        // ====================================================================================
        // TWO DEFECTS PORTED FROM ATJUP AS KNOBS, both default 0 = exactly what ATSAT does today.
        // ATJUP ships (1) ON, on Jupiter evidence; whether that evidence transfers to Saturn is
        // what these are here to measure. ATJUP's own notes are at ChemistryJup.h:45.
        //
        // (1) ATSAT_CHEM_GATE_ZERO — the temperature gate has no else branch. w_nh3, w_h2s and
        //     w_nh4sh are only assigned inside `if(t_00 <= t_u <= t_0)`, and they are persistent
        //     Arrays, so a cell that LEAVES the formation window keeps the last rate it ever had
        //     and goes on applying it to the tendency every iteration for the rest of the run.
        //     On ATJUP every NH4SH hot spot sat OUTSIDE the window, i.e. was made entirely of
        //     frozen rates, and enabling this dropped them by a factor ~20. Set to 1 to zero the
        //     three rates outside the window instead of freezing them.
        //
        // (2) ATSAT_CHEM_MOLAR_CONC — which concentration the quadratic rate law uses. ATSAT and
        //     ATJUP do not merely differ here, they are three different things, so this knob has
        //     three values rather than ATJUP's two:
        //
        //       0  (default)  c_x = w_x / m_x                  — ATSAT today: NO density factor
        //       1             c_x = rho_mix * w_x / m_x        — the plain molar concentration
        //                                                        [kmol/m3]; ATJUP's =1
        //       2             c_x = (r_mix / sum_c) * w_x/m_x  — ATJUP's DEFAULT normalisation,
        //                                                        sum_c over nh3, h2s, nh4sh
        //
        //     Value 2 is provided to make the two models comparable, NOT because it is right:
        //     ATJUP's note records that it forces c_nh3 + c_h2s + c_nh4sh = r_mix, i.e. rescales
        //     three trace species to carry the entire mixture density, and that NH4SH then holds
        //     99.99 % of sum_c so the back reaction saturates and the rate law stops seeing its
        //     own product. It is also a unit mismatch. Value 1 is the dimensionally correct one.
        //
        //     rho_mix is filled by computeMixtureDensity(), which runs AFTER this routine on the
        //     first call, so fall back to the scalar r_mix then — as ATJUP does for the same
        //     reason.
        //
        // A THIRD DIFFERENCE, DELIBERATELY NOT MADE A KNOB HERE because it is not ATJUP's: ATSAT
        // divides both rates by m.dt (see react_rate_ordinary). ATJUP does not. That makes the
        // "rate" timestep-dependent — the chemistry applies a fixed INCREMENT per step rather
        // than a rate — so halving dt does not halve the chemical change per unit time. It is
        // flagged, not settled, and it is the reason a straight A/B against ATJUP's numbers would
        // not mean what it looks like.
        // ====================================================================================
        // (3) ATSAT_CHEM_DIFFLUX_PLUS — the sign of the diffusive term in massflux_*. DEFAULT 1 since 2026-11-03
        //     with ATSAT_CHEM_DIFFLUX_LAPLACIAN=1 (the correct pair, see DiffMassFluxSat); =0 restores the minus.
        //
        //     WHAT ATSAT COMPUTES, which is NOT what ATJUP computes. ATJUP forms difflux_* in one
        //     step as D_x * laplacian_spherical(c_x), the diffusive tendency directly. ATSAT does
        //     it in two passes (DiffMassFluxSat, this file, lines ~298 and ~325):
        //         pass 1   j_x       = (c_mix/r_mix) * (m_x * D_x * grad_x - jT_x)
        //         pass 2   difflux_x = dj_x/dr + |dj_x/dthe|/rm + dj_x/dphi/(rm*sinthe)
        //     so difflux_* here is the DIVERGENCE OF A FLUX, not a Laplacian. Do not read ATJUP's
        //     description of its own difflux_* onto this one.
        //
        //     WHY THE SIGN IS STILL WRONG IN r AND phi. Fick's law is j = -D grad(c) and the
        //     tendency is -div(j) = +D lap(c). Pass 1 builds j as +D grad(c) — Fick's minus is
        //     missing — and pass 2 takes +div of it, so the two compose to +D lap(c), the same
        //     quantity ATJUP's difflux_* holds. massflux_* = w_* MINUS difflux_* therefore
        //     subtracts a tendency that should be added, and the multicomponent diffusion acts as
        //     an ANTI-diffusion, sharpening every species gradient instead of smoothing it.
        //
        //     WHAT THE KNOB CANNOT FIX, and this is larger than it first looked. difflux_* is not
        //     a divergence at all: pass 1 SUMS the three components of grad(c) into a scalar, so
        //     j_* is a scalar rather than a flux vector, and pass 2 differentiates that scalar and
        //     sums again. A spherical divergence carries (1/r^2) d(r^2 F_r)/dr and
        //     (1/(r sinthe)) d(F_the sinthe)/dthe; none of those metric factors are present. The
        //     std::abs() on the theta component is deliberate and documented at DiffMassFluxSat —
        //     it symmetrises an otherwise pole-antisymmetric term — not an oversight.
        //
        //     So =1 flips the sign of a quantity that is not the diffusive tendency in the first
        //     place. It is worth having as a diagnostic, and the growth measured below is real,
        //     but it cannot be read as "the diffusion now has the right sign".
        //
        //     THE ACTUAL FIX is ATSAT_CHEM_DIFFLUX_LAPLACIAN=1, which replaces the two-pass scalar
        //     with a real D * grad^2(c) — ATJUP's formulation. The physically correct combination
        //     is that knob AND this one:
        //         ATSAT_CHEM_DIFFLUX_LAPLACIAN=1 ATSAT_CHEM_DIFFLUX_PLUS=1
        //     An earlier version of this comment said the fix was "Fick's minus plus dropping the
        //     abs()". That was wrong: neither change turns a summed-components scalar into a
        //     Laplacian.
        //
        //     Why it survived unnoticed in both models: the coefficient is tiny. D_x is built as
        //     mue_x/(rg_x*sc_x) with rg_h2s the density of the CONDENSED phase rather than the
        //     gas, giving D ~ 1e-8 m2/s, and difflux_* prints as 0.000000 against massflux_* of
        //     order 1e-4. Tiny is not the same as harmless: anti-diffusion is SELF-AMPLIFYING —
        //     it sharpens a gradient, the sharper gradient raises the Laplacian, which sharpens
        //     it further — so the right test is whether the difference GROWS, not how big it is
        //     at one moment. ATJUP measured 1.27 % on nh3 and 1.56 % on nh4sh at 100 iterations.
        //     The measurement for Saturn is in the commit message.
        //
        //     The questionable D_x is left alone: it belongs to a separate question. So does the
        //     fact that `chemical_reaction` in the RHS multiplies the WHOLE of massflux_*, so
        //     setting that switch to 0 to disable the chemistry silently disables this diffusion
        //     as well.
        static const int gate_zero  = [](){
            const char* e = getenv("ATSAT_CHEM_GATE_ZERO");  return e ? atoi(e) : 0; }();
        static const int molar_conc = [](){
            const char* e = getenv("ATSAT_CHEM_MOLAR_CONC"); return e ? atoi(e) : 0; }();
        static const double dsign = [](){
            const char* e = getenv("ATSAT_CHEM_DIFFLUX_PLUS");      // DEFAULT 1 since 2026-11-03, see below
            return (e ? atoi(e) != 0 : true) ? 1.0 : -1.0; }();
        // ATSAT_NH4SH_EQ=1: the rate law contributes nothing (w_* = 0); the diffusion term stays.
        const bool eq_on = nh4sh_eq();

        #pragma omp parallel for collapse(3) schedule(static)
        for(int k = ATPhys::phi_first(ATPhys::phi_periodic<cSaturnModel>()); k < km-1; k++){   // ATSAT_BC_PHI_PERIODIC
            for(int j = 1; j < jm-1; j++){
                for(int i = 1; i < im-1; i++){
                    const double t_u = m.t.x[i][j][k] * m.t_ref;
                    double kf = react_rate_const(t_u, T_d, A, B);
                    double kb = kf / keq;

                    if(eq_on){
                        // ATSAT_NH4SH_EQ: EquilibriumNH4SH forms / decomposes the solid; no rate law.
                        m.w_nh3.x[i][j][k]   = 0.0;
                        m.w_h2s.x[i][j][k]   = 0.0;
                        m.w_nh4sh.x[i][j][k] = 0.0;
                    }
                    else if((t_u <= m.t_0_nh4sh) && (t_u >= m.t_00_nh4sh))
                        react_rate_ordinary(i, j, k, kf, kb, molar_conc,
                            m.nh3, m.h2s, m.nh4sh, m.w_nh3, m.w_h2s, m.w_nh4sh);
                    else if(gate_zero){
                        // Outside the formation window there is no reaction, so the rates are
                        // zero — not "whatever they were the last time this cell was inside it".
                        m.w_nh3.x[i][j][k]   = 0.0;
                        m.w_h2s.x[i][j][k]   = 0.0;
                        m.w_nh4sh.x[i][j][k] = 0.0;
                    }

                    // (3) ATSAT_CHEM_DIFFLUX_PLUS — the sign of the multicomponent diffusion.
                    // See the note above the loop. dsign is +1.0 or -1.0, so the default path is
                    // the subtraction it always was, exactly.
                    m.massflux_h2s.x[i][j][k]   = m.w_h2s.x[i][j][k]   + dsign * m.difflux_h2s.x[i][j][k];
                    m.massflux_nh3.x[i][j][k]   = m.w_nh3.x[i][j][k]   + dsign * m.difflux_nh3.x[i][j][k];
                    m.massflux_nh4sh.x[i][j][k] = m.w_nh4sh.x[i][j][k] + dsign * m.difflux_nh4sh.x[i][j][k];
                }
            }
        }

        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for ChemMassRateSat\n", elapsed.count() * 1e-9);
        cout << "      ATSAT: ChemMassRateSat ended" << endl;
    }

    // The TVD limiter is the SHARED FluxLimiter<Planet> template — 97 % identical between
    // the two models, the whole difference being metricRadius(). See FluxLimiter.h for why
    // this one routine could be shared while the rest of the chemistry cannot.
    void FluxLimiterNH4SH(){ FluxLimiter<cSaturnModel>(m).nh4sh(); }

    // ===================================================================================================
    // ATSAT_NH4SH_EQ (2026-11-03, row SAT-CHEM-DIVERGE; DEFAULT 1 since 2026-11-03, =0 restores the rate law):
    // NH4SH as an EQUILIBRIUM ADJUSTMENT, replacing the rate law of ChemMassRateSat (whose w_* are then
    // zero). Ported from ATJUP (ATJUP_NH4SH_EQ, ChemistryJup.h, default 1 there since 2026-10-26).
    //
    // WHY. The rate law has no source for any of its numbers: kb = kf / keq with keq = m_nh4sh /
    // (m_nh3 m_h2s), a ratio of molar masses, not an equilibrium constant; A / B / T_d = 1500 / -0.3 / 3020
    // here against 15000 / -0.5 / 3020 on ATJUP; the window 140-230 K here against 200-230 K; both rates
    // divided by m.dt (an increment per step, not a rate); and no heat is released. The equilibrium makes
    // all of those moot: K_p is a property of the reaction, not of the planet.
    //
    // WHAT. NH3(g) + H2S(g) <-> NH4SH(s) is in equilibrium where p_NH3 p_H2S = K_p(T),
    //     log10 K_p [atm^2] = 14.82 - 4705 / T          (Lewis 1969; Weidenschilling & Lewis 1973),
    // with p_x = (c_x / m_x) R* T [Pa] and c_x the species density [kg/m3] (the convention of the shared
    // saturation adjustment, target rho ep E / p). Per fluid cell x [kmol/m3] of NH4SH forms (x > 0,
    // product above K_p) or decomposes (x < 0, product below K_p and solid present, at most all of it),
    // heating the cell by L x m_nh4sh / (cp_mix rho): L = 4705 ln(10) R* / m_nh4sh = 1.762e6 J/kg, the
    // van't Hoff enthalpy of the same K_p (90.1 kJ/mol; dh_nh4sh in cSaturnModel.h is 1.7687e6), so the
    // heat and the equilibrium agree. T and K_p are solved together (bisection on x; f(x) = (a-x)(b-x)
    // (R* T(x))^2 - K_p(T(x)) falls monotonically). m_nh3 + m_h2s = m_nh4sh, so N, S and the mass are
    // conserved by construction. No window. Called after the saturation adjustments (set-up and loop);
    // H2S and NH4SH join ATSAT_RK_SCALAR_SYNC; the MCA carries the enthalpy in NH4SH's tracer L.
    // MEASURED (satchk/chemeq_2026-11-03, 1000 it from scratch, K0 = rate law vs KE = equilibrium, diffusion pair on in
    // both): NH4SH max settles at 21.9-22.2 g/m3 from it 300 (rate law 3.67 -> 1.92 g/m3, still falling); NH3 ice max
    // 2.9-3.1 (rate law 4.24 -> 5.39, growing); per call ~1.42 M cells form / ~1.19 M decompose, constant 100 -> 1000;
    // T, u, OLR 1.375, residuum identical to the printed digit; the MCA fires in 77-94 % of the columns (rate law
    // 4-39 %), max dT per call <= 0.44 K.
    // ===================================================================================================
    static double nh4sh_formation_heat(){ return 4705.0 * std::log(10.0) * 8314.46 / 51.1114; }

    static bool nh4sh_eq(){
        static const bool v = [](){ const char* e = getenv("ATSAT_NH4SH_EQ"); return e ? atoi(e) != 0 : true; }();
        return v;
    }

    // One cell to equilibrium; returns x [kmol/m3] formed (> 0) or decomposed (< 0) and applies it.
    double equilibrateNH4SHCell(int i, int j, int k, bool& all_decomposed, double& sat_entry, double& dT)
    {
        constexpr double R_univ = 8314.46, atm = 101325.0;
        static const double L_f = nh4sh_formation_heat();
        auto Kp = [&](double T){ return std::pow(10.0, 14.82 - 4705.0 / T) * atm * atm; };   // Pa^2
        all_decomposed = false; sat_entry = 0.0; dT = 0.0;
        const double T0 = m.t.x[i][j][k] * m.t_ref;
        if (!(T0 > 0.0)) return 0.0;
        double rho = m.rho_mix.x[i][j][k];
        if (!(rho > 0.0) || !std::isfinite(rho)) rho = m.r_mix;
        const double a = std::max(0.0, m.nh3.x[i][j][k])   / m.m_nh3;          // kmol/m3
        const double b = std::max(0.0, m.h2s.x[i][j][k])   / m.m_h2s;
        const double s = std::max(0.0, m.nh4sh.x[i][j][k]) / m.m_nh4sh;
        const double heat = L_f * m.m_nh4sh / (m.cp_mix * rho);                // K per kmol/m3 formed
        auto f = [&](double x){
            const double T = T0 + heat * x, RT = R_univ * T;
            return (a - x) * (b - x) * RT * RT - Kp(T); };
        const double f0 = f(0.0);
        double lo, hi;
        if (f0 > 0.0)         { lo = 0.0; hi = std::min(a, b); sat_entry = f0 / Kp(T0); }
        else if (s > 0.0) {
            if (f(-s) <= 0.0) { lo = hi = -s; all_decomposed = true; }
            else              { lo = -s; hi = 0.0; }
        } else return 0.0;
        if (hi > lo) for (int it = 0; it < 80; it++) {
            const double xm = 0.5 * (lo + hi);
            if (f(xm) > 0.0) lo = xm; else hi = xm;
            if (hi - lo <= 1.0e-15 * std::max(a, std::max(b, s))) break;
        }
        const double x = (f0 > 0.0) ? lo : hi;                                 // never past the root
        if (x == 0.0) return 0.0;
        m.nh3.x[i][j][k]   = std::max(0.0, m.nh3.x[i][j][k]   - m.m_nh3   * x);
        m.h2s.x[i][j][k]   = std::max(0.0, m.h2s.x[i][j][k]   - m.m_h2s   * x);
        m.nh4sh.x[i][j][k] = std::max(0.0, m.nh4sh.x[i][j][k] + m.m_nh4sh * x);
        dT = heat * x;
        m.t.x[i][j][k] = (T0 + dT) / m.t_ref;
        return x;
    }

    // where == "init": the set-up's column construction (Weidenschilling & Lewis 1973): each column is
    // swept bottom-up with a rising parcel; a cell takes NH3 and H2S as the parcel's remaining mass
    // fractions (capped by its own set-up value), is equilibrated, and passes what is left upward. Not
    // conservative by design (it builds a state). Any other where: every fluid cell in place, conservative.
    // Each column is independent, so the result does not depend on the thread count.
    void EquilibriumNH4SH(const char* where)
    {
        using namespace std;
        auto begin = chrono::high_resolution_clock::now();
        const bool sweep = std::string(where) == "init";
        long n_form = 0, n_decomp = 0, n_all = 0;
        double dm_form = 0.0, dm_decomp = 0.0, dT_max = 0.0, s_max = 0.0;
        #pragma omp parallel for collapse(2) schedule(dynamic, 8) \
            reduction(+:n_form, n_decomp, n_all, dm_form, dm_decomp) reduction(max:dT_max, s_max)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                double q_nh3 = 1.0e300, q_h2s = 1.0e300;                            // the parcel, mass fractions
                for (int i = 0; i < m.im; i++) {
                    if (m.SeaMount.x[i][j][k] == 1.0) continue;
                    double rho = m.rho_mix.x[i][j][k];
                    if (!(rho > 0.0) || !isfinite(rho)) rho = m.r_mix;
                    if (sweep) {
                        m.nh3.x[i][j][k] = min(m.nh3.x[i][j][k], q_nh3 * rho);
                        m.h2s.x[i][j][k] = min(m.h2s.x[i][j][k], q_h2s * rho);
                    }
                    bool all_ = false; double sat_ = 0.0, dT_ = 0.0;
                    const double x = equilibrateNH4SHCell(i, j, k, all_, sat_, dT_);
                    if (sweep) { q_nh3 = m.nh3.x[i][j][k] / rho; q_h2s = m.h2s.x[i][j][k] / rho; }
                    s_max = max(s_max, sat_); dT_max = max(dT_max, fabs(dT_));
                    if (all_) n_all++;
                    if (x > 0.0)      { n_form++;   dm_form   += m.m_nh4sh * x; }
                    else if (x < 0.0) { n_decomp++; dm_decomp -= m.m_nh4sh * x; }
                }
            }
        }
        auto end = chrono::high_resolution_clock::now();
        printf("      ATSAT: EquilibriumNH4SH (%s%s, Lewis K_p, L %.4e J/kg): formed in %ld cells (sum %.4e kg/m3),"
               " decomposed in %ld (sum %.4e kg/m3, %ld completely), max p_NH3 p_H2S / K_p - 1 at entry %.3e,"
               " max dT %.4f K; %.3f s\n", where, sweep ? ", parcel sweep" : "",
               nh4sh_formation_heat(), n_form, dm_form, n_decomp, dm_decomp, n_all, s_max, dT_max,
               chrono::duration_cast<chrono::nanoseconds>(end - begin).count() * 1e-9);
    }

    // -----------------------------------------------------------------------
    // ATSAT_CHEM_DIFFLUX_LAPLACIAN — what difflux_* actually is. DEFAULT 1 since 2026-11-03 together with
    // ATSAT_CHEM_DIFFLUX_PLUS=1 (row SAT-CHEM-DIVERGE item 4); =0 + PLUS=0 restores the old two-pass scalar.
    // 200 it from scratch (satchk/chemeq_2026-11-03 L1 vs D0): the spurious massflux_nh4sh at the lid
    // 652 -> 7 ug/(m3 s), NH4SH / NH3 / H2S maxima within 0.01 %, u / v / w / residuum / OLR unchanged
    // (1.382 -> 1.380 W/m2), top T min -0.03 K.
    //
    // WHAT THE DEFAULT COMPUTES, and it is worth being exact because it is easy to misread as a
    // divergence. Pass 1 builds a SCALAR per species,
    //     dc = dc/dr + |dc/dthe|/rm + dc/dphi/(rm*sinthe)
    // by SUMMING the three components of grad(c) as if they were scalars — the |.| is deliberate
    // and documented below, to symmetrise an otherwise pole-antisymmetric term. j_* is then that
    // scalar, not a flux vector. Pass 2 differentiates j_* and sums the components the same way.
    //
    // So difflux_* is neither grad nor div: a spherical divergence carries (1/r^2) d(r^2 F_r)/dr
    // and (1/(r sinthe)) d(F_the sinthe)/dthe, and none of those metric factors appear. Adding
    // Fick's missing minus sign or dropping the |.| would replace one ill-defined scalar with a
    // different ill-defined scalar; neither turns this into a Laplacian.
    //
    // WHAT =1 COMPUTES: difflux_x = D_x * laplacian_spherical(c_x), a real spherical Laplacian
    // with all five terms including 2/rm * dc/dr and cotthe/rm^2 * dc/dthe. This is ATJUP's
    // formulation, ported unchanged apart from asking metricRadius() for rm. It is the diffusive
    // TENDENCY, so it belongs with ATSAT_CHEM_DIFFLUX_PLUS=1: the pair
    //     ATSAT_CHEM_DIFFLUX_LAPLACIAN=1 ATSAT_CHEM_DIFFLUX_PLUS=1
    // is the physically correct combination, and either alone is not. They are kept as separate
    // knobs rather than one because they answer separate questions — what the term IS, and which
    // way it enters — and collapsing them would make the measurement unattributable.
    //
    // Pass 1 and thermalmassflux are untouched either way: j_* still feeds the thermal term.
    //
    // This is also what blocks sharing DiffMassFlux with ATJUP, which forms D*laplacian directly
    // and carries Lewis-number groups. The two routines overlap 25 % line-for-line today; with
    // this knob on they would agree about the quantity, which is the first of several steps.
    void DiffMassFluxSat()
    {
        using namespace std;
        cout << endl << "      ATSAT: DiffMassFluxSat" << endl;

        static const int difflux_laplacian = [](){
            const char* e = getenv("ATSAT_CHEM_DIFFLUX_LAPLACIAN"); return e ? atoi(e) : 1; }();   // DEFAULT 1 since 2026-11-03

        auto begin = chrono::high_resolution_clock::now();

        const int im = m.im, jm = m.jm, km = m.km;

        // Diffusion coefficients are spatially constant — hoist out of the loop.
        m.D_nh3    = m.mue_nh3   / (m.rg_nh3   * m.sc_nh3);
        m.D_h2s    = m.mue_h2s   / (m.rg_h2s   * m.sc_h2s);
        m.D_nh4sh  = m.mue_nh4sh / (m.rg_nh4sh * m.sc_nh4sh);
        m.DT_nh3   = m.D_nh3;
        m.DT_h2s   = m.D_h2s;
        m.DT_nh4sh = m.D_nh4sh;

        // Pass 1: compute jT_* and j_* (ordinary + thermal diffusion fluxes).
        // Writes only at [i][j][k]; reads neighbouring cells of t/h2s/nh3/nh4sh.
        #pragma omp parallel for collapse(3) schedule(static)
        for(int k = ATPhys::phi_first(ATPhys::phi_periodic<cSaturnModel>()); k < km-1; k++){   // ATSAT_BC_PHI_PERIODIC
            for(int j = 1; j < jm-1; j++){
                for(int i = 1; i < im-1; i++){
                    const double rm       = m.metricRadius(m.rad.z[i]);
                    const double sinthe   = sin(m.the.z[j]);
                    const double rmsinthe = rm * sinthe;

                    double dtdr = 0.0,     dtdthe = 0.0,     dtdphi = 0.0;
                    double dh2sdr = 0.0,   dh2sdthe = 0.0,   dh2sdphi = 0.0;
                    double dnh3dr = 0.0,   dnh3dthe = 0.0,   dnh3dphi = 0.0;
                    double dnh4shdr = 0.0, dnh4shdthe = 0.0, dnh4shdphi = 0.0;

                    derivative_1_order(i, j, k, dtdr,     dtdthe,     dtdphi,     m.t);
                    derivative_1_order(i, j, k, dh2sdr,   dh2sdthe,   dh2sdphi,   m.h2s);
                    derivative_1_order(i, j, k, dnh3dr,   dnh3dthe,   dnh3dphi,   m.nh3);
                    derivative_1_order(i, j, k, dnh4shdr, dnh4shdthe, dnh4shdphi, m.nh4sh);

                    m.jT_nh3.x[i][j][k]   = m.DT_nh3   / m.t.x[i][j][k] * dtdr;
                    m.jT_h2s.x[i][j][k]   = m.DT_h2s   / m.t.x[i][j][k] * dtdr;
                    m.jT_nh4sh.x[i][j][k] = m.DT_nh4sh / m.t.x[i][j][k] * dtdr;

                    // Pole-symmetric sum of gradient components: dXdr and dXdphi are
                    // pole-symmetric for symmetric inputs, dXdthe is pole-antisymmetric.
                    // std::abs() on the theta term symmetrizes the resulting scalar.
                    const double dnh3   = dnh3dr   + std::abs(dnh3dthe)/rm   + dnh3dphi/rmsinthe;
                    const double dh2s   = dh2sdr   + std::abs(dh2sdthe)/rm   + dh2sdphi/rmsinthe;
                    const double dnh4sh = dnh4shdr + std::abs(dnh4shdthe)/rm + dnh4shdphi/rmsinthe;

                    m.j_nh3.x[i][j][k]   = m.c_mix/m.r_mix * (m.m_nh3   * m.D_nh3   * dnh3   - m.jT_nh3.x[i][j][k]);
                    m.j_h2s.x[i][j][k]   = m.c_mix/m.r_mix * (m.m_h2s   * m.D_h2s   * dh2s   - m.jT_h2s.x[i][j][k]);
                    m.j_nh4sh.x[i][j][k] = m.c_mix/m.r_mix * (m.m_nh4sh * m.D_nh4sh * dnh4sh - m.jT_nh4sh.x[i][j][k]);
                }
            }
        }

        // Pass 2: compute difflux_* (divergence of j_*) and thermalmassflux.
        // Reads j_* from neighbouring cells — all written by pass 1.
        #pragma omp parallel for collapse(3) schedule(static)
        for(int k = ATPhys::phi_first(ATPhys::phi_periodic<cSaturnModel>()); k < km-1; k++){   // ATSAT_BC_PHI_PERIODIC
            for(int j = 1; j < jm-1; j++){
                for(int i = 1; i < im-1; i++){
                    const double rm       = m.metricRadius(m.rad.z[i]);
                    const double sinthe   = sin(m.the.z[j]);
                    const double rmsinthe = rm * sinthe;

                    double dtdr = 0.0,       dtdthe = 0.0,       dtdphi = 0.0;
                    double dj_h2sdr = 0.0,   dj_h2sdthe = 0.0,   dj_h2sdphi = 0.0;
                    double dj_nh3dr = 0.0,   dj_nh3dthe = 0.0,   dj_nh3dphi = 0.0;
                    double dj_nh4shdr = 0.0, dj_nh4shdthe = 0.0, dj_nh4shdphi = 0.0;

                    derivative_1_order(i, j, k, dtdr,       dtdthe,       dtdphi,       m.t);
                    derivative_1_order(i, j, k, dj_h2sdr,   dj_h2sdthe,   dj_h2sdphi,   m.j_h2s);
                    derivative_1_order(i, j, k, dj_nh3dr,   dj_nh3dthe,   dj_nh3dphi,   m.j_nh3);
                    derivative_1_order(i, j, k, dj_nh4shdr, dj_nh4shdthe, dj_nh4shdphi, m.j_nh4sh);

                    if(difflux_laplacian){
                        // A genuine D * grad^2(c). See the note above this routine.
                        m.difflux_h2s.x[i][j][k]   = m.D_h2s   * laplacian_spherical(i, j, k, m.h2s);
                        m.difflux_nh3.x[i][j][k]   = m.D_nh3   * laplacian_spherical(i, j, k, m.nh3);
                        m.difflux_nh4sh.x[i][j][k] = m.D_nh4sh * laplacian_spherical(i, j, k, m.nh4sh);
                    }else{
                        m.difflux_h2s.x[i][j][k]   = dj_h2sdr   + std::abs(dj_h2sdthe)/rm   + dj_h2sdphi/rmsinthe;
                        m.difflux_nh3.x[i][j][k]   = dj_nh3dr   + std::abs(dj_nh3dthe)/rm   + dj_nh3dphi/rmsinthe;
                        m.difflux_nh4sh.x[i][j][k] = dj_nh4shdr + std::abs(dj_nh4shdthe)/rm + dj_nh4shdphi/rmsinthe;
                    }

                    m.thermalmassflux.x[i][j][k] =
                         (m.j_nh3.x[i][j][k]   * m.cp_nh3
                        + m.j_h2s.x[i][j][k]   * m.cp_h2s
                        + m.j_nh4sh.x[i][j][k] * m.cp_nh4sh)
                        * (dtdr + std::abs(dtdthe)/rm + dtdphi/rmsinthe)
                        + (chem_enthalpy()
                             ? m.w_nh4sh.x[i][j][k] * m.dh_nh4sh
                             : m.t.x[i][j][k] * (m.w_nh3.x[i][j][k]   * m.k_nh3
                                          + m.w_h2s.x[i][j][k]   * m.k_h2s
                                          + m.w_nh4sh.x[i][j][k] * m.k_nh4sh));
                }
            }
        }

        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for DiffMassFluxSat\n", elapsed.count() * 1e-9);
        cout << "      ATSAT: DiffMassFluxSat ended" << endl;
    }

    // -----------------------------------------------------------------------
    void ThermalPropertiesSat()
    {
        using namespace std;
        cout << endl << "      ATSAT: ThermalPropertiesSat" << endl;

        auto begin = chrono::high_resolution_clock::now();

        double r_mixture = 1.326;

        m.rg_mix = m.rg_h2 + m.rg_he + m.rg_h2s + m.rg_nh3 + m.rg_nh4sh + m.rg_h2o + m.rg_ch4;
        m.r_mix  = m.r_h2  + m.r_he  + m.r_h2s  + m.r_nh3  + m.r_nh4sh  + m.r_h2o  + m.r_ch4;
        m.c_mix  = m.c_h2  + m.c_he  + m.c_h2s  + m.c_nh3  + m.c_nh4sh  + m.c_h2o  + m.c_ch4;

        double M_mix = m.r_mix / m.c_mix;

        // MASS-FRACTION WEIGHTING, which is ATJUP's form (ChemistryJup.h:441) and was not ATSAT's.
        //
        // These three are per-unit-MASS properties, so a mixture value is sum(mass fraction * x),
        // i.e. sum(r_x * X_x) / r_mix. What stood here divided each species by its OWN MOLAR MASS
        // instead of by the mixture density: sum(r_x/m_x * X_x). Those weights are mole densities
        // and they sum to c_mix = 0.4833 kmol/m3, not to 1, so the result was neither a mass
        // average nor a molar one — and every species heavier than H2 was suppressed by the ratio
        // m_x/r_mix, which for CH4 is a factor of 11.6.
        //
        // THIS IS A DEFECT AND NOT A MODELLING CHOICE, on the routine's own evidence. R_mix below
        // is already normalised by the density sum, so the two lines disagreed with each other;
        // and for an ideal gas cp - cv = R, which makes the old cp_mix say
        //     gamma = cp/(cp - R_mix) = 6357.57/(6357.57 - 2909.25) = 1.8428
        // for an atmosphere that is 62.5 % H2 by mass. No ideal gas can exceed the monatomic 5/3;
        // a H2-dominated one should be near 1.4. The corrected value gives 1.4066.
        //
        // Measured, ATSAT, from ThermalPropertiesSat's own printout:
        //     cp_mix   6357.5720524226 -> 10058.5352... J/(kg K)     +58.2 %
        //     mue_mix  0.0001419220    -> see commit message
        //     k_mix    0.0606394882    -> see commit message
        // ATSAT's gam is a config constant (0.7 K/km), so init_temperature is NOT affected here.
        // What moves is every consumer of cp_mix: the ConvectiveAdjustment reference adiabat
        // g/cp_mix, 1.5729 -> 0.9942 K/km, plus Q_Sensible, the RHS energy coefficient and the
        // saturation adjustment. ATNEPT and ATURAN carry the same defect and are fixed alongside;
        // on those two gam IS g*1e3/cp_mix, so their initial profiles move as well.
        m.cp_mix  = (m.r_h2 * m.cp_h2 + m.r_h2s * m.cp_h2s
                  +  m.r_he * m.cp_he + m.r_nh3 * m.cp_nh3
                  +  m.r_nh4sh * m.cp_nh4sh + m.r_h2o * m.cp_h2o
                  +  m.r_ch4 * m.cp_ch4) / m.r_mix;
        m.mue_mix = (m.r_h2 * m.mue_h2 + m.r_h2s * m.mue_h2s
                  +  m.r_he * m.mue_he + m.r_nh3 * m.mue_nh3
                  +  m.r_nh4sh * m.mue_nh4sh + m.r_h2o * m.mue_h2o
                  +  m.r_ch4 * m.mue_ch4) / m.r_mix;
        m.k_mix   = (m.r_h2 * m.k_h2 + m.r_h2s * m.k_h2s
                  +  m.r_he * m.k_he + m.r_nh3 * m.k_nh3
                  +  m.r_nh4sh * m.k_nh4sh + m.r_h2o * m.k_h2o
                  +  m.r_ch4 * m.k_ch4) / m.r_mix;
        m.R_mix   = (m.r_h2 * m.R_h2 + m.r_he * m.R_he + m.r_h2o * m.R_h2o
                  + m.r_h2s * m.R_h2s + m.r_nh3 * m.R_nh3 + m.r_nh4sh * m.R_nh4sh
                  + m.r_ch4 * m.R_ch4)
                  / (m.r_h2 + m.r_he + m.r_h2o + m.r_h2s + m.r_nh3 + m.r_nh4sh + m.r_ch4);

        cout.precision(10);
        cout.setf(ios::fixed);
        cout << endl
            << "     r_mixture[kg/m³] = " << r_mixture << endl << endl
            << "     rg_mix[kg/m³] = " << m.rg_mix << endl << endl
            << "     r_mix[kg/m³] = " << m.r_mix << endl << endl
            << "     c_mix[kmol/m³] = " << m.c_mix << endl << endl
            << "     M_mix[kg/Kmol] = " << M_mix << endl << endl
            << "     cp_mix[J/(kg*K)] = " << m.cp_mix << endl
            << "     mue_mix[Ns/m²] = " << m.mue_mix << endl
            << "     k_mix[W/(m*K)] = " << m.k_mix << endl
            << "     R_mix[J/(kg*K)] = " << m.R_mix << endl << endl;

        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for ThermalPropertiesSat\n", elapsed.count() * 1e-9);
        cout << "      ATSAT: ThermalPropertiesSat ended" << endl;
    }

private:
    cSaturnModel& m;

    // -----------------------------------------------------------------------
    void derivative_1_order(int i, int j, int k,
        double &dcdr, double &dcdthe, double &dcdphi, Array &c)
    {
        // The phi neighbours: k-1, k+1, k-2, k+2, taken round the ring at the seam when
        // ATSAT_BC_PHI_PERIODIC=1 (ATPhys.h).
        const bool phi_per = ATPhys::phi_periodic<cSaturnModel>();
        const int kw  = ATPhys::phi_wrap(phi_per, k-1, m.km), ke  = ATPhys::phi_wrap(phi_per, k+1, m.km);
        const int kw2 = ATPhys::phi_wrap(phi_per, k-2, m.km), ke2 = ATPhys::phi_wrap(phi_per, k+2, m.km);

        // radial gradients — SeaMount-aware
        if((m.SeaMount.x[i][j][k]==1.0) && (m.SeaMount.x[i+1][j][k]==0.0) && (m.SeaMount.x[i+2][j][k]==0.0)){
            c.x[i][j][k] = c.x[i+3][j][k] - 3.0*c.x[i+2][j][k] + 3.0*c.x[i+1][j][k];
            dcdr = (-3.0*c.x[i][j][k] + 4.0*c.x[i+1][j][k] - c.x[i+2][j][k]) / (2.0*m.dr);
        }
        if((m.SeaMount.x[i-1][j][k]==1.0) && (m.SeaMount.x[i][j][k]==0.0) && (m.SeaMount.x[i+1][j][k]==0.0))
            dcdr = (c.x[i+1][j][k] - c.x[i-1][j][k]) / (2.0*m.dr);
        if((m.SeaMount.x[i-1][j][k]==0.0) && (m.SeaMount.x[i][j][k]==0.0) && (m.SeaMount.x[i+1][j][k]==0.0))
            dcdr = (c.x[i+1][j][k] - c.x[i-1][j][k]) / (2.0*m.dr);

        // south (j+1) gradients
        if((m.SeaMount.x[i][j][k]==1.0) && (m.SeaMount.x[i][j+1][k]==0.0) && (m.SeaMount.x[i][j+2][k]==0.0)){
            c.x[i][j][k] = c43*c.x[i][j+1][k] - c13*c.x[i][j+2][k];
            dcdthe = (-3.0*c.x[i][j][k] + 4.0*c.x[i][j+1][k] - c.x[i][j+2][k]) / (2.0*m.dthe);
        }
        if((m.SeaMount.x[i][j-1][k]==1.0) && (m.SeaMount.x[i][j][k]==0.0) && (m.SeaMount.x[i][j+1][k]==0.0))
            dcdthe = (c.x[i][j+1][k] - c.x[i][j-1][k]) / (2.0*m.dthe);

        // north (j-1) gradients
        if((m.SeaMount.x[i][j][k]==1.0) && (m.SeaMount.x[i][j-1][k]==0.0) && (m.SeaMount.x[i][j-2][k]==0.0)){
            c.x[i][j][k] = c43*c.x[i][j-1][k] - c13*c.x[i][j-2][k];
            dcdthe = (-3.0*c.x[i][j][k] + 4.0*c.x[i][j-1][k] - c.x[i][j-2][k]) / (2.0*m.dthe);
        }
        if((m.SeaMount.x[i][j-1][k]==1.0) && (m.SeaMount.x[i][j][k]==0.0) && (m.SeaMount.x[i][j+1][k]==0.0))
            dcdthe = (c.x[i][j+1][k] - c.x[i][j-1][k]) / (2.0*m.dthe);
        if((m.SeaMount.x[i][j-1][k]==0.0) && (m.SeaMount.x[i][j][k]==0.0) && (m.SeaMount.x[i][j+1][k]==0.0))
            dcdthe = (c.x[i][j+1][k] - c.x[i][j-1][k]) / (2.0*m.dthe);

        // east (k+1) gradients
        if((m.SeaMount.x[i][j][k]==1.0) && (m.SeaMount.x[i][j][ke]==0.0) && (m.SeaMount.x[i][j][ke2]==0.0)){
            c.x[i][j][k] = c43*c.x[i][j][ke] - c13*c.x[i][j][ke2];
            dcdphi = (-3.0*c.x[i][j][k] + 4.0*c.x[i][j][ke] - c.x[i][j][ke2]) / (2.0*m.dphi);
        }
        if((m.SeaMount.x[i][j][kw]==1.0) && (m.SeaMount.x[i][j][k]==0.0) && (m.SeaMount.x[i][j][ke]==0.0))
            dcdphi = (c.x[i][j][ke] - c.x[i][j][kw]) / (2.0*m.dphi);

        // west (k-1) gradients
        if((m.SeaMount.x[i][j][k]==1.0) && (m.SeaMount.x[i][j][kw]==0.0) && (m.SeaMount.x[i][j][kw2]==0.0)){
            c.x[i][j][k] = c43*c.x[i][j][kw] - c13*c.x[i][j][kw2];
            dcdphi = (-3.0*c.x[i][j][k] + 4.0*c.x[i][j][kw] - c.x[i][j][kw2]) / (2.0*m.dphi);
        }
        if((m.SeaMount.x[i][j][kw]==1.0) && (m.SeaMount.x[i][j][k]==0.0) && (m.SeaMount.x[i][j][ke]==0.0))
            dcdphi = (c.x[i][j][ke] - c.x[i][j][kw]) / (2.0*m.dphi);
        if((m.SeaMount.x[i][j][kw]==0.0) && (m.SeaMount.x[i][j][k]==0.0) && (m.SeaMount.x[i][j][ke]==0.0))
            dcdphi = (c.x[i][j][ke] - c.x[i][j][kw]) / (2.0*m.dphi);
    }

    // -----------------------------------------------------------------------
    void derivative_1_order_boundary(int i, int j, int k,
        double &dcdr, double &dcdthe, double &dcdphi, Array &c)
    {
        if(i == 0)       dcdr   = (-3.0*c.x[0][j][k]    + 4.0*c.x[1][j][k]    - c.x[2][j][k])    / (2.0*m.dr);
        if(i == m.im-1)  dcdr   = (-3.0*c.x[m.im-1][j][k] + 4.0*c.x[m.im-2][j][k] - c.x[m.im-3][j][k]) / (2.0*m.dr);
        if(j == 0)       dcdthe = (-3.0*c.x[i][0][k]    + 4.0*c.x[i][1][k]    - c.x[i][2][k])    / (2.0*m.dthe);
        if(j == m.jm-1)  dcdthe = (-3.0*c.x[i][m.jm-1][k] + 4.0*c.x[i][m.jm-2][k] - c.x[i][m.jm-3][k]) / (2.0*m.dthe);
        if(k == 0)       dcdphi = (-3.0*c.x[i][j][0]    + 4.0*c.x[i][j][1]    - c.x[i][j][2])    / (2.0*m.dphi);
        if(k == m.km-1)  dcdphi = (-3.0*c.x[i][j][m.km-1] + 4.0*c.x[i][j][m.km-2] - c.x[i][j][m.km-3]) / (2.0*m.dphi);
    }

    // -----------------------------------------------------------------------
    // Spherical Laplacian, ported from ATJUP's ChemistryJup::laplacian_spherical. Used only when
    // ATSAT_CHEM_DIFFLUX_LAPLACIAN is set; see the note above DiffMassFluxSat.
    //
    // Two ATSAT adaptations, both deliberate:
    //   rm comes from metricRadius(), so the term follows ATSAT_METRIC_RADIUS like everything else
    //   here. coord_stretching is false in ATSAT, so exp_rm folds to 1 — the branch is kept so the
    //   file stays a faithful port and would still be right if that ever changed.
    //
    //   sinthe is floored at 0.4, as ATJUP floors it, because d2c/dphi2 / (rm*sinthe)^2 diverges
    //   at the pole otherwise. NOTE this is a LOCAL floor and NOT cSaturnModel::sinthe_min(),
    //   which defaults to 0.0 — ATSAT deliberately has no global polar floor. Using the global
    //   one here would make this term blow up at the default, so the local 0.4 matches both ATJUP
    //   and FluxLimiter.h. That the two floors disagree is a real question for ATSAT's polar
    //   treatment, and it is flagged rather than settled.
    double laplacian_spherical(int i, int j, int k, Array& c)
    {
        const double rm       = m.metricRadius(m.rad.z[i]);
        const double exp_rm   = m.coord_stretching ? 1.0 / (rm + 1.0) : 1.0;
        const double exp_2_rm = exp_rm * exp_rm;
        const double sinthe   = sin(m.the.z[j]);
        const double costhe   = cos(m.the.z[j]);
        const double rmsinthe = rm * std::max(sinthe,
                                    ATPhys::polar_divisor_floor<cSaturnModel>());

        const double d2cdr2   = (c.x[i+1][j][k] - 2.0*c.x[i][j][k] + c.x[i-1][j][k]) / (m.dr   * m.dr) * exp_2_rm;
        const double d2cdthe2 = (c.x[i][j+1][k] - 2.0*c.x[i][j][k] + c.x[i][j-1][k]) / (m.dthe * m.dthe);
        const bool phi_per = ATPhys::phi_periodic<cSaturnModel>();          // ATSAT_BC_PHI_PERIODIC
        const int kw = ATPhys::phi_wrap(phi_per, k-1, m.km), ke = ATPhys::phi_wrap(phi_per, k+1, m.km);
        const double d2cdphi2 = (c.x[i][j][ke] - 2.0*c.x[i][j][k] + c.x[i][j][kw]) / (m.dphi * m.dphi);

        const double dcdr   = (c.x[i+1][j][k] - c.x[i-1][j][k]) / (2.0 * m.dr) * exp_rm;
        const double dcdthe = (c.x[i][j+1][k] - c.x[i][j-1][k]) / (2.0 * m.dthe);

        return d2cdr2
             + 2.0 / rm * dcdr
             + d2cdthe2 / (rm * rm)
             + costhe / (rm * rmsinthe) * dcdthe
             + d2cdphi2 / (rmsinthe * rmsinthe);
    }

    void react_rate_ordinary(int i, int j, int k,
        double &f, double &b, int molar_conc,
        Array &molecule_nh3, Array &molecule_h2s, Array &molecule_nh4sh,
        Array &mass_nh3,     Array &mass_h2s,     Array &mass_nh4sh)
    {
        // ATSAT_CHEM_MOLAR_CONC — see the note in ChemMassRateSat(). denom = 1.0 is ATSAT's own
        // form and multiplying by it is exact, so the default path is bit-identical.
        double denom = 1.0;
        if(molar_conc == 1){
            double rho = m.rho_mix.x[i][j][k];
            if(!(rho > 0.0)) rho = m.r_mix;   // first call: computeMixtureDensity has not run
            denom = rho;
        }else if(molar_conc == 2){
            const double sum_c = molecule_nh3.x[i][j][k]   / m.m_nh3
                               + molecule_h2s.x[i][j][k]   / m.m_h2s
                               + molecule_nh4sh.x[i][j][k] / m.m_nh4sh;
            denom = (sum_c == 0.0) ? 0.0 : m.r_mix / sum_c;
        }

        double c_nh3   = molecule_nh3.x[i][j][k]   / m.m_nh3   * denom;
        double c_h2s   = molecule_h2s.x[i][j][k]   / m.m_h2s   * denom;
        double c_nh4sh = molecule_nh4sh.x[i][j][k] / m.m_nh4sh * denom;

        double Rf = f * (c_nh3 * c_h2s) / m.dt;
        double Rb = b * c_nh4sh         / m.dt;
        double R_diff = Rf - Rb;

        mass_nh3.x[i][j][k]   = -m.m_nh3   * R_diff;
        mass_h2s.x[i][j][k]   = -m.m_h2s   * R_diff;
        mass_nh4sh.x[i][j][k] =  m.m_nh4sh * R_diff;
    }

    // -----------------------------------------------------------------------
    static double react_rate_const(double T_K, const double T_d,
        const double A, const double B)
    {
        return A * std::pow(T_K, B) * std::exp(-T_d / T_K);
    }

    // cSaturnModel::c43 / c13 shorthands
    static constexpr double c43 = 4.0/3.0;
    static constexpr double c13 = 1.0/3.0;

    // Superbee: most compressive TVD limiter — best for sharp cloud fronts.

    // Van Leer: smooth, differentiable — good general-purpose alternative.
};
