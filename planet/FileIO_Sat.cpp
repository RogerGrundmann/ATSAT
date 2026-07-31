/*
 * Atmosphere General Circulation Modell(ATSAT) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in aa spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and nh3 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to prepare the boundary and initial conditions for diverse variables
*/
#include "cSaturnModel.h"
#include "Utils.h"

#include <vector>
#include <cstdint>

using namespace std;
using namespace AtomUtils;

void cSaturnModel::writeData(){
    cout << endl << "      ATSAT: writeData" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

//    int i_radial = 0;
//    int i_radial = 3;
    int i_radial = 20;
//    int i_radial = 36;
//    int i_radial = 39;
//    int i_radial = 40;
    paraview_vtk_radial(iter_n, i_radial);
    int j_longal = 112;
    paraview_vtk_longal(iter_n, j_longal);
    int k_zonal = 180;
//    int k_zonal = 0;
    paraview_vtk_zonal(iter_n, k_zonal);

    // The full 3D panorama .vts, which is the only output carrying every field over the whole
    // volume. It used to fire on panorama_print alone, which is a different cadence from the
    // checkpoint that writes the .vtk slices and from the 100-iteration restart stride — so a
    // run could produce slices at 100 and 200 and a panorama at neither, which is what the
    // 2026-07-31 200-iteration run had to work around by switching panorama_print to 999.
    // The 100-stride is added so the three output kinds land on the same iterations: with
    // checkpoint = 100 a run now writes .vtk, .vts and .bin together at 100 and at 200.
    // panorama_print is kept as an additional trigger so existing configurations still get
    // what they asked for.
    if(paraview_flag && (iter_n % panorama_print == 0 || iter_n % 100 == 0)){
        paraview_panorama_vts(iter_n);
//        paraview_sphere_vts(iter_n);
    }

    SaturnPlotData();

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for writeData\n", elapsed.count() * 1e-9);

    cout << "      ATSAT: writeData ended" << endl;

    return;
}
/*
*
*/
void cSaturnModel::writeResults(){
    cout << endl << "      ATSAT: writeResults" << endl;
//    double coeff_mmWS = r_mix/r_h2o;  // coeff_mmWS = 1.2041/0.0094 [kg/m³/kg/m³] = 128,0827 [/]
//    double coeff_lv = lv_h2o /(cp_h2o * t_0);  // coefficient for the specific latent Evaporation heat(Condensation heat), coeff_lv = 9.1069 in [/]
//    double coeff_ls = ls /(cp_h2o * t_0);  // coefficient for the specific latent Evaporation heat(Condensation heat), coeff_ls = 10.3091 in [/]
//    double a, e;

    auto begin = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            precipitable_water.y[j][k] = 0.;  // precipitable water
        }
    }
    #pragma omp parallel for
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            Q_Latent.x[0][j][k] = c43 * Q_Latent.x[1][j][k] - c13 * Q_Latent.x[2][j][k];
            Q_Latent.x[im-1][j][k] = c43 * Q_Latent.x[im-2][j][k] - c13 * Q_Latent.x[im-3][j][k];
            Q_Sensible.x[0][j][k] = c43 * Q_Sensible.x[1][j][k] - c13 * Q_Sensible.x[2][j][k];
            Q_Sensible.x[im-1][j][k] = c43 * Q_Sensible.x[im-2][j][k] - c13 * Q_Sensible.x[im-3][j][k];
            BuoyancyForce.x[0][j][k] = c43 * BuoyancyForce.x[1][j][k] - c13 * BuoyancyForce.x[2][j][k];
            BuoyancyForce.x[im-1][j][k] = c43 * BuoyancyForce.x[im-2][j][k] - c13 * BuoyancyForce.x[im-3][j][k];
/*
            nh3.x[0][j][k] = c43 * nh3.x[1][j][k] - c13 * nh3.x[2][j][k];
            nh3.x[im-1][j][k] = c43 * nh3.x[im-2][j][k] - c13 * nh3.x[im-3][j][k];
            nh3_cloud.x[0][j][k] = c43 * nh3_cloud.x[1][j][k] - c13 * nh3_cloud.x[2][j][k];
            nh3_cloud.x[im-1][j][k] = c43 * nh3_cloud.x[im-2][j][k] - c13 * nh3_cloud.x[im-3][j][k];
            nh3_ice.x[0][j][k] = c43 * nh3_ice.x[1][j][k] - c13 * nh3_ice.x[2][j][k];
            nh3_ice.x[im-1][j][k] = c43 * nh3_ice.x[im-2][j][k] - c13 * nh3_ice.x[im-3][j][k];
*/
        }
    }
    #pragma omp parallel for
    for(int k = 0; k < km; k++){
        for(int i = 0; i < im; i++){
            t.x[i][0][k] = c43 * t.x[i][1][k] - c13 * t.x[i][2][k];
            t.x[i][jm-1][k] = c43 * t.x[i][jm-2][k] - c13 * t.x[i][jm-3][k];
            Q_Latent.x[i][0][k] = c43 * Q_Latent.x[i][1][k] - c13 * Q_Latent.x[i][2][k];
            Q_Latent.x[i][jm-1][k] = c43 * Q_Latent.x[i][jm-2][k] - c13 * Q_Latent.x[i][jm-3][k];
            Q_Sensible.x[i][0][k] = c43 * Q_Sensible.x[i][1][k] - c13 * Q_Sensible.x[i][2][k];
            Q_Sensible.x[i][jm-1][k] = c43 * Q_Sensible.x[i][jm-2][k] - c13 * Q_Sensible.x[i][jm-3][k];
            BuoyancyForce.x[i][0][k] = c43 * BuoyancyForce.x[i][1][k] - c13 * BuoyancyForce.x[i][2][k];
            BuoyancyForce.x[i][jm-1][k] = c43 * BuoyancyForce.x[i][jm-2][k] - c13 * BuoyancyForce.x[i][jm-3][k];
/*
            nh3.x[i][0][k] = c43 * nh3.x[i][1][k] - c13 * nh3.x[i][2][k];
            nh3.x[i][jm-1][k] = c43 * nh3.x[i][jm-2][k] - c13 * nh3.x[i][jm-3][k];
            nh3_cloud.x[i][0][k] = c43 * nh3_cloud.x[i][1][k] - c13 * nh3_cloud.x[i][2][k];
            nh3_cloud.x[i][jm-1][k] = c43 * nh3_cloud.x[i][jm-2][k] - c13 * nh3_cloud.x[i][jm-3][k];
            nh3_ice.x[i][0][k] = c43 * nh3_ice.x[i][1][k] - c13 * nh3_ice.x[i][2][k];
            nh3_ice.x[i][jm-1][k] = c43 * nh3_ice.x[i][jm-2][k] - c13 * nh3_ice.x[i][jm-3][k];
            h2o.x[i][0][k] = c43 * h2o.x[i][1][k] - c13 * h2o.x[i][2][k];
            h2o.x[i][jm-1][k] = c43 * h2o.x[i][jm-2][k] - c13 * h2o.x[i][jm-3][k];
            h2o_cloud.x[i][0][k] = c43 * h2o_cloud.x[i][1][k] - c13 * h2o_cloud.x[i][2][k];
            h2o_cloud.x[i][jm-1][k] = c43 * h2o_cloud.x[i][jm-2][k] - c13 * h2o_cloud.x[i][jm-3][k];
            h2o_ice.x[i][0][k] = c43 * h2o_ice.x[i][1][k] - c13 * h2o_ice.x[i][2][k];
            h2o_ice.x[i][jm-1][k] = c43 * h2o_ice.x[i][jm-2][k] - c13 * h2o_ice.x[i][jm-3][k];
*/
/*
            P_rain.x[i][0][k] = c43 * P_rain.x[i][1][k] - c13 * P_rain.x[i][2][k];
            P_rain.x[i][jm-1][k] = c43 * P_rain.x[i][jm-2][k] - c13 * P_rain.x[i][jm-3][k];
            P_snow.x[i][0][k] = c43 * P_snow.x[i][1][k] - c13 * P_snow.x[i][2][k];
            P_snow.x[i][jm-1][k] = c43 * P_snow.x[i][jm-2][k] - c13 * P_snow.x[i][jm-3][k];
*/
        }
    }
    #pragma omp parallel for
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            t.x[i][j][0] = c43 * t.x[i][j][1] - c13 * t.x[i][j][2];
            t.x[i][j][km-1] = c43 * t.x[i][j][km-2] - c13 * t.x[i][j][km-3];
            t.x[i][j][0] = t.x[i][j][km-1] =(t.x[i][j][0] + t.x[i][j][km-1])/ 2.;
            Q_Latent.x[i][j][0] = c43 * Q_Latent.x[i][j][1] - c13 * Q_Latent.x[i][j][2];
            Q_Latent.x[i][j][km-1] = c43 * Q_Latent.x[i][j][km-2] - c13 * Q_Latent.x[i][j][km-3];
            Q_Latent.x[i][j][0] = Q_Latent.x[i][j][km-1] =(Q_Latent.x[i][j][0] + Q_Latent.x[i][j][km-1])/ 2.;
            Q_Sensible.x[i][j][0] = c43 * Q_Sensible.x[i][j][1] - c13 * Q_Sensible.x[i][j][2];
            Q_Sensible.x[i][j][km-1] = c43 * Q_Sensible.x[i][j][km-2] - c13 * Q_Sensible.x[i][j][km-3];
            Q_Sensible.x[i][j][0] = Q_Sensible.x[i][j][km-1] =(Q_Sensible.x[i][j][0] + Q_Sensible.x[i][j][km-1])/ 2.;
            BuoyancyForce.x[i][j][0] = c43 * BuoyancyForce.x[i][j][1] - c13 * BuoyancyForce.x[i][j][2];
            BuoyancyForce.x[i][j][km-1] = c43 * BuoyancyForce.x[i][j][km-2] - c13 * BuoyancyForce.x[i][j][km-3];
            BuoyancyForce.x[i][j][0] = BuoyancyForce.x[i][j][km-1] =(BuoyancyForce.x[i][j][0] + BuoyancyForce.x[i][j][km-1])/ 2.;
/*
            nh3.x[i][j][0] = c43 * nh3.x[i][j][1] - c13 * nh3.x[i][j][2];
            nh3.x[i][j][km-1] = c43 * nh3.x[i][j][km-2] - c13 * nh3.x[i][j][km-3];
            nh3.x[i][j][0] = nh3.x[i][j][km-1] =(nh3.x[i][j][0] + nh3.x[i][j][km-1])/ 2.;
            nh3_cloud.x[i][j][0] = c43 * nh3_cloud.x[i][j][1] - c13 * nh3_cloud.x[i][j][2];
            nh3_cloud.x[i][j][km-1] = c43 * nh3_cloud.x[i][j][km-2] - c13 * nh3_cloud.x[i][j][km-3];
            nh3_cloud.x[i][j][0] = nh3_cloud.x[i][j][km-1] =(nh3_cloud.x[i][j][0] + nh3_cloud.x[i][j][km-1])/ 2.;
            nh3_ice.x[i][j][0] = c43 * nh3_ice.x[i][j][1] - c13 * nh3_ice.x[i][j][2];
            nh3_ice.x[i][j][km-1] = c43 * nh3_ice.x[i][j][km-2] - c13 * nh3_ice.x[i][j][km-3];
            nh3_ice.x[i][j][0] = nh3_ice.x[i][j][km-1] =(nh3_ice.x[i][j][0] + nh3_ice.x[i][j][km-1])/ 2.;
            h2o_cloud.x[i][j][0] = c43 * h2o_cloud.x[i][j][1] - c13 * h2o_cloud.x[i][j][2];
            h2o_cloud.x[i][j][km-1] = c43 * h2o_cloud.x[i][j][km-2] - c13 * h2o_cloud.x[i][j][km-3];
            h2o_cloud.x[i][j][0] = h2o_cloud.x[i][j][km-1] =(h2o_cloud.x[i][j][0] + h2o_cloud.x[i][j][km-1])/ 2.;
            h2o_ice.x[i][j][0] = c43 * h2o_ice.x[i][j][1] - c13 * h2o_ice.x[i][j][2];
            h2o_ice.x[i][j][km-1] = c43 * h2o_ice.x[i][j][km-2] - c13 * h2o_ice.x[i][j][km-3];
            h2o_ice.x[i][j][0] = h2o_ice.x[i][j][km-1] =(h2o_ice.x[i][j][0] + h2o_ice.x[i][j][km-1])/ 2.;
*/
/*
            P_rain.x[i][j][0] = c43 * P_rain.x[i][j][1] - c13 * P_rain.x[i][j][2];
            P_rain.x[i][j][km-1] = c43 * P_rain.x[i][j][km-2] - c13 * P_rain.x[i][j][km-3];
            P_rain.x[i][j][0] = P_rain.x[i][j][km-1] =(P_rain.x[i][j][0] + P_rain.x[i][j][km-1])/ 2.;
            P_snow.x[i][j][0] = c43 * P_snow.x[i][j][1] - c13 * P_snow.x[i][j][2];
            P_snow.x[i][j][km-1] = c43 * P_snow.x[i][j][km-2] - c13 * P_snow.x[i][j][km-3];
            P_snow.x[i][j][0] = P_snow.x[i][j][km-1] =(P_snow.x[i][j][0] + P_snow.x[i][j][km-1])/ 2.;
*/
        }
    }
/*
    precipitable_water.y[0][0] = 0.;
    precipitable_water.y[0][0] = 0.;
    nh3_total.y[0][0] = 0.;
    nh3_cloud_total.y[0][0] = 0.;
    nh3_ice_total.y[0][0] = 0.;
    double precipitablewater_average = 0.;
    double precipitation_average = 0.;
    double h2_vegetation_average = 0.;
    double he_vegetation_average = 0.;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            nh3_total.y[j][k] = nh3.x[0][j][k];
            nh3_cloud_total.y[j][k] = nh3_cloud.x[0][j][k];
            nh3_ice_total.y[j][k] = nh3_ice.x[0][j][k];
            for(int i = 0; i < im; i++){
                e = h2o.x[i][j][k] * p_stat.x[i][j][k]/ep_h2o;  // water vapour pressure in hPa
                a = 216.6 * e /(t.x[i][j][k] * t_0);  // absolute humidity in kg/m3
                precipitable_water.y[j][k] += a * L_atm /(double)(im - 1);  // kg/m³ * m
            }
        }
    }
    double coeff_prec = 86400.;  // dimensions see below
// surface values of precipitation and precipitable water
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            Precipitation.y[j][k] = coeff_prec * (P_rain.x[0][j][k] + P_snow.x[0][j][k]);
            // 60 s * 60 min * 24 h = 86400 s == 1 d
            // Precipitation, P_rain and P_snow in kg/ ( m² * s ) = mm/s
            // Precipitation in 86400. * kg/ ( m² * d ) = 86400 mm/d
            // kg/ ( m² * s ) == mm/s ( Kraus, p. 94 )
            if(Precipitation.y[j][k] >= 25.)  Precipitation.y[j][k] = 25.;
            if(Precipitation.y[j][k] <= 0)  Precipitation.y[j][k] = 0.;
        }
    }
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            precipitablewater_average += precipitable_water.y[j][k];
            precipitation_average += Precipitation.y[j][k];
            h2_vegetation_average += nh3_total.y[j][k];
            he_vegetation_average += nh3_cloud_total.y[j][k];
        }
    }
    h2_vegetation_average = h2_vegetation_average /(double)((jm-1)*(km-1));
    he_vegetation_average = he_vegetation_average /(double)((jm-1)*(km-1));
    precipitablewater_average = precipitablewater_average /(double)((jm-1)*(km-1));
    precipitation_average = 365. * precipitation_average /(double)((jm-1)*(km-1));
    cout.precision(2);
    string level = "m";
    string deg_north = "°N";
    string deg_south = "°S";
    string deg_west = "°W";
    string deg_east = "°E";
    string name_Value_7 = " precipitable water average ";
    string name_Value_8 = " precipitation average per year ";
    string name_Value_9 = " precipitation average per day ";
    string name_Value_22 = " H2_average ";
    string name_Value_25 = " He_average ";
    string name_unit_mmd = " mm/d";
    string name_unit_mm = " mm";
    string name_unit_mma = " mm/a";
    string name_unit_ppm = " kg/kg";
    cout << endl;
    double Value_7 = precipitablewater_average;
    double Value_8 = precipitation_average;
    cout << setw(6)<< setiosflags(ios::left)<< setw(40)<< setfill('.')<< name_Value_7 << " = " << resetiosflags(ios::left)<< setw(7)<< fixed << setfill(' ')<< Value_7 << setw(6)<< name_unit_mm << "   " << setiosflags(ios::left)<< setw(40)<< setfill('.')<< name_Value_8 << " = " << resetiosflags(ios::left)<< setw(7)<< fixed << setfill(' ')<< Value_8 << setw(6)<< name_unit_mma << "   " << setiosflags(ios::left)<< setw(40)<< setfill('.')<< name_Value_9 << " = " << resetiosflags(ios::left)<< setw(7)<< fixed << setfill(' ')<< Value_8/365. << setw(6)<< name_unit_mmd << endl;
    double Value_9 = h2_vegetation_average;
    double Value_25 = he_vegetation_average;
    cout << setw(6)<< setiosflags(ios::left)<< setw(40)<< setfill('.')<< name_Value_22 << " = " << resetiosflags(ios::left)<< setw(7)<< fixed << setfill(' ')<< Value_9 << setw(6)<< name_unit_ppm << "   " << setiosflags(ios::left)<< setw(40)<< setfill('.')<< name_Value_25 << " = " << resetiosflags(ios::left)<< setw(7)<< fixed << setfill(' ')<< Value_25 << setw(6)<< name_unit_ppm << endl << endl << endl;
    return;
*/

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for writeResults\n", elapsed.count() * 1e-9);

    cout << "      ATSAT: writeResults ended" << endl;

}
/*
*
*/

/*
 * ATSAT_NANCHECK — a CENSUS of non-finite cells, not just the first one.
 *
 * Ported from cJupiterModel::nan_watch. The reason it counts rather than reporting the first
 * hit: the scan runs i outermost, so a NaN born in the interior gets reported at the i=0
 * boundary cell that merely inherited it through the bcRadius extrapolation. Counting per field
 * and recording the index extent shows at a glance whether this is one interior cell, a whole
 * boundary plane, or already everywhere — which is the difference between a local physics
 * problem and a global one.
 *
 * Worth having at all because NaN is INVISIBLE to printMinMax: searchMinMax_3D compares with a
 * bare > and every comparison against a NaN is false, so a non-finite cell is silently skipped
 * and a thoroughly broken field can print a perfectly reasonable maximum.
 *
 * The bit test is exact and branch-free: exponent all ones means inf or NaN, whatever the
 * payload, and it does not depend on the compiler's floating-point flags the way isnan() can.
 *
 * ATJUP scans the arrays its restart serialises. ATSAT has no restart (that is item 7 of the
 * ATJUP/ATSAT gap list), so the list is written out here — the prognostic fields plus the two
 * pressures, which is what determines whether the next iteration is meaningful.
 */
bool cSaturnModel::nan_watch(int iter){
    struct Named { const char* name; Array* a; };
    Named fields[] = {
        {"t",         &t},         {"u",         &u},         {"v",     &v},
        {"w",         &w},         {"p_dyn",     &p_dyn},     {"p_stat", &p_stat},
        {"h2o",       &h2o},       {"h2o_cloud", &h2o_cloud}, {"h2o_ice", &h2o_ice},
        {"ch4",       &ch4},       {"ch4_cloud", &ch4_cloud}, {"ch4_ice", &ch4_ice},
        {"nh3",       &nh3},       {"nh3_cloud", &nh3_cloud}, {"nh3_ice", &nh3_ice},
        {"h2s",       &h2s},       {"nh4sh",     &nh4sh},
        {"tke",       &tke},       {"dis",       &dis},       {"nue",     &nue},
        {"rho_mix",   &rho_mix},
    };
    const int nf = (int)(sizeof(fields)/sizeof(fields[0]));

    long total = 0;
    bool any = false;

    for(int a = 0; a < nf; a++){
        long n = 0;
        int i_lo = im, i_hi = -1, j_lo = jm, j_hi = -1, k_lo = km, k_hi = -1;
        int fi = -1, fj = -1, fk = -1;
        for(int i = 0; i < im; i++)
            for(int j = 0; j < jm; j++)
                for(int k = 0; k < km; k++){
                    std::uint64_t bits;
                    std::memcpy(&bits, &fields[a].a->x[i][j][k], sizeof(bits));
                    if((bits & 0x7FF0000000000000ULL) == 0x7FF0000000000000ULL){
                        n++;
                        if(i < i_lo) i_lo = i;
                        if(i > i_hi) i_hi = i;
                        if(j < j_lo) j_lo = j;
                        if(j > j_hi) j_hi = j;
                        if(k < k_lo) k_lo = k;
                        if(k > k_hi) k_hi = k;
                        if(fi < 0){ fi = i; fj = j; fk = k; }
                    }
                }
        if(n > 0){
            if(!any){
                printf("      ATSAT: ===== NAN WATCH: state went non-finite at iteration %d =====\n",
                       iter);
                any = true;
            }
            printf("        %-12s %8ld cells   i[%d..%d] j[%d..%d] k[%d..%d]   first (%d,%d,%d)"
                   "  t=%g p_stat=%g\n",
                   fields[a].name, n, i_lo, i_hi, j_lo, j_hi, k_lo, k_hi, fi, fj, fk,
                   t.x[fi][fj][fk], p_stat.x[fi][fj][fk]);
            total += n;
        }
    }
    if(any) printf("        total %ld non-finite cells\n", total);
    return !any;
}
/*
*
*/

/*
 * ===== Restart checkpoint (item 7 of the ATJUP/ATSAT gap list) =====
 *
 * Ported from cJupiterModel::save_state / load_state / restart_state_is_clean. ATSAT has had no
 * restart of any kind: every run began at iteration 0 and any result past the first few hundred
 * iterations had to be paid for again from scratch. The file is a raw dump — a 5-int header
 * followed by the arrays, each written as im*jm contiguous rows of km doubles.
 *
 * Only the genuinely PROGNOSTIC arrays are stored, with the same reasoning as ATJUP: the
 * reaction rates, the diffusive and thermal mass fluxes, the forces, the radiation and
 * precipitation fluxes and the latent and sensible heat fields are all recomputed from these at
 * the top of the next iteration, so storing them would add bulk and one more way for the file to
 * disagree with itself.
 *
 * p_stat is stored although it is quasi-static: the buoyancy term and the whole saturation chain
 * read it, and it has to match the temperature field it was built with.
 *
 * The magic is "SAT1", not ATJUP's "JUP1", so a Jupiter restart handed to ATSAT is rejected by
 * the header check instead of being read as 41x181x361 doubles of nonsense — the two models
 * have the same grid dimensions, so the dimension check alone would not catch it.
 *
 * p_dynn IS stored, but only since the commit that revived steadyQuery. It used to be declared
 * and never allocated — no p_dynn.initArray() anywhere — so its data pointer was the NULL that
 * Array's default constructor leaves, and putting it in this list segfaulted the first
 * 100-iteration checkpoint. It is now allocated and refreshed by restoreVar with the other
 * n-copies, which is what makes it worth serialising: a run resumed from a restart has to know
 * the previous iteration's pressure or its first steady-state report is meaningless.
 */
std::vector<Array*> cSaturnModel::restart_arrays(){
    return { &t,   &u,   &v,   &w,
             &tn,  &un,  &vn,  &wn,
             &h2o,  &h2o_cloud,  &h2o_ice,
             &h2on, &h2o_cloudn, &h2o_icen,
             &h2s,  &h2sn,
             &nh3,  &nh3_cloud,  &nh3_ice,
             &nh3n, &nh3_cloudn, &nh3_icen,
             &ch4,  &ch4_cloud,  &ch4_ice,
             &ch4n, &ch4_cloudn, &ch4_icen,
             &nh4sh, &nh4shn,
             &p_dyn, &p_dynn, &p_stat,
             &tke, &dis, &tken, &disn, &nue };
}

// Every array in restart_arrays() must actually own storage. An Array that was declared but
// never initArray'd has x == NULL, and the three routines below walk these pointers without
// looking — which is how p_dynn took down the first checkpoint. Checked once, on the first call,
// and reported by index so the offending entry is identifiable without a debugger.
static bool restart_arrays_allocated(const std::vector<Array*>& arrs, const char* who){
    for(size_t a = 0; a < arrs.size(); a++){
        if(arrs[a]->x == NULL){
            printf("      ATSAT: %s ABORTED - restart array #%zu was never initArray'd"
                   " (declared but not allocated); fix restart_arrays()\n", who, a);
            return false;
        }
    }
    return true;
}

void cSaturnModel::save_state(int iter){
    const string fn = output_path + "/sat_restart_" + std::to_string(iter) + ".bin";
    std::ofstream f(fn, std::ios::binary);
    if(!f){
        cout << "      ATSAT: save_state FAILED to open " << fn << endl;
        return;
    }
    // Header: magic, grid dimensions, and the iteration this state belongs to. The grid is
    // checked on load, so a restart written at another resolution is rejected rather than read
    // as garbage.
    const int32_t hdr[5] = { 0x53415431 /*"SAT1"*/, im, jm, km, iter };
    f.write(reinterpret_cast<const char*>(hdr), sizeof(hdr));

    std::vector<Array*> arrs = restart_arrays();
    if(!restart_arrays_allocated(arrs, "save_state")) return;
    for(size_t a = 0; a < arrs.size(); a++)
        for(int i = 0; i < im; i++)
            for(int j = 0; j < jm; j++)
                f.write(reinterpret_cast<const char*>(arrs[a]->x[i][j]), km * sizeof(double));

    if(!f){
        cout << "      ATSAT: save_state FAILED while writing " << fn
             << " (disk full?)" << endl;
        return;
    }
    const double mb = (double)(sizeof(hdr) + arrs.size() * (size_t)im * jm * km * sizeof(double))
                    / (1024.0 * 1024.0);
    printf("      ATSAT: save_state wrote %zu arrays (%.1f MB) to %s\n",
           arrs.size(), mb, fn.c_str());
}

bool cSaturnModel::load_state(int iter){
    const string fn = output_path + "/sat_restart_" + std::to_string(iter) + ".bin";
    std::ifstream f(fn, std::ios::binary);
    if(!f){
        cout << "      ATSAT: load_state: no file " << fn
             << " - running from scratch" << endl;
        return false;
    }
    int32_t hdr[5];
    f.read(reinterpret_cast<char*>(hdr), sizeof(hdr));
    if(!f || hdr[0] != 0x53415431 || hdr[1] != im || hdr[2] != jm || hdr[3] != km){
        cout << "      ATSAT: load_state: bad header / grid mismatch in " << fn
             << " - running from scratch" << endl;
        return false;
    }

    std::vector<Array*> arrs = restart_arrays();
    if(!restart_arrays_allocated(arrs, "load_state")) return false;
    for(size_t a = 0; a < arrs.size(); a++)
        for(int i = 0; i < im; i++)
            for(int j = 0; j < jm; j++){
                f.read(reinterpret_cast<char*>(arrs[a]->x[i][j]), km * sizeof(double));
                if(!f){
                    cout << "      ATSAT: load_state: truncated file " << fn
                         << " - running from scratch" << endl;
                    return false;
                }
            }

    cout << "      ATSAT: load_state restored " << arrs.size() << " arrays from "
         << fn << " (resuming after iteration " << hdr[4] << ")" << endl;
    return true;
}

// True when every serialised prognostic field is finite everywhere. This guards the periodic
// checkpoint: a diverged state must never overwrite a good restart point, because being able to
// resume from it is the file's whole value. The bit test is the one nan_watch uses above —
// exponent all ones means inf or NaN whatever the payload, and unlike isfinite() it does not
// depend on the compiler's floating-point flags. Which matters here: the ATSAT_DT=0.001 run of
// 2026-07-31 was non-finite in 512797 cells at iteration 1 and would otherwise have written a
// 1.2 GB file of NaN over its predecessor.
bool cSaturnModel::restart_state_is_clean(){
    std::vector<Array*> arrs = restart_arrays();
    if(!restart_arrays_allocated(arrs, "restart_state_is_clean")) return false;
    bool clean = true;
    for(size_t a = 0; a < arrs.size() && clean; a++)
        for(int i = 0; i < im && clean; i++)
            for(int j = 0; j < jm && clean; j++)
                for(int k = 0; k < km; k++){
                    std::uint64_t bits;
                    std::memcpy(&bits, &arrs[a]->x[i][j][k], sizeof(bits));
                    if((bits & 0x7FF0000000000000ULL) == 0x7FF0000000000000ULL){
                        clean = false;
                        break;
                    }
                }
    return clean;
}
/*
*
*/

/*
 * ===== Zero floor for the condensable species, with accounting =====
 *
 * Ported from cJupiterModel::clampNegativeSpecies. WHY IT IS NEEDED HERE, measured rather than
 * assumed: the 200-iteration run of 2026-07-31 (config_run200.xml, default timestep) reports
 *
 *     min h2o    -0.000002 g/m3      min nh3   -0.0517 -> -0.0318 g/m3
 *     min nh4sh  -0.0316 g/m3 at iteration 100, -0.0626 at iteration 200
 *
 * i.e. NH4SH's negative minimum DOUBLES over 100 iterations. A mass density below zero has no
 * meaning, the saturation formulas take logarithms of these fields, and the chemistry multiplies
 * two of them together — a negative concentration there produces a reaction rate with the wrong
 * SIGN, which is a source where there should be a sink.
 *
 * The cause is the same one ATJUP identified: the centred differences of the transport terms
 * undershoot wherever a species has a sharp edge, and the (4/3,-1/3) extrapolation at the radial
 * boundary planes undershoots a field that is already essentially zero there. ATJUP measured the
 * clipped amount as small enough that a plain floor is the right answer rather than a
 * mass-conserving filler; the measurement for ATSAT is in the commit that adds this.
 *
 * A floor is in principle a mass SOURCE, so the clipped amount is accumulated per field and
 * reported next to printMinMax rather than left invisible. Read that report as GROSS clipping,
 * not net mass gained: the saturation adjustment re-partitions vapour and condensate on the next
 * pass and hands most of it straight back. What the counter is FOR is the day that stops being
 * true. ATSAT_NO_CLAMP=1 turns the floor off and restores the old behaviour exactly.
 */
void cSaturnModel::clampNegativeSpecies(){
    static const int off = [](){ const char* e = getenv("ATSAT_NO_CLAMP"); return e ? atoi(e) : 0; }();
    if(off) return;

    Array* fields[] = {
        &h2o, &h2o_cloud, &h2o_ice,
        &h2s,
        &nh3, &nh3_cloud, &nh3_ice,
        &ch4, &ch4_cloud, &ch4_ice,
        &nh4sh };
    const int nf = (int)(sizeof(fields) / sizeof(fields[0]));

    if((int)clamp_added.size() != nf){
        clamp_added.assign(nf, 0.0);
        clamp_cells.assign(nf, 0);
    }

    for(int f = 0; f < nf; f++){
        Array& F = *fields[f];
        double added = 0.0;
        long   cells = 0;
        #pragma omp parallel for collapse(2) schedule(static) reduction(+:added,cells)
        for(int i = 0; i < im; i++){
            for(int j = 0; j < jm; j++){
                for(int k = 0; k < km; k++){
                    const double v = F.x[i][j][k];
                    // Written as !(v >= 0.0) so a NaN is caught here too rather than carried on.
                    if(!(v >= 0.0)){
                        if(std::isfinite(v)){ added -= v; cells++; }
                        F.x[i][j][k] = 0.0;
                    }
                }
            }
        }
        clamp_added[f] += added;
        clamp_cells[f] += cells;
    }
}

// Companion report, called from printMinMax so it shares the checkpoint cadence.
void cSaturnModel::reportClampBudget(){
    static const char* const names[] = {
        "h2o","h2o_cloud","h2o_ice","h2s","nh3","nh3_cloud","nh3_ice",
        "ch4","ch4_cloud","ch4_ice","nh4sh" };
    const int nf = (int)(sizeof(names)/sizeof(names[0]));
    if((int)clamp_added.size() != nf) return;

    Array* fields[] = {
        &h2o, &h2o_cloud, &h2o_ice, &h2s, &nh3, &nh3_cloud, &nh3_ice,
        &ch4, &ch4_cloud, &ch4_ice, &nh4sh };

    bool any = false;
    for(int f = 0; f < nf; f++) if(clamp_cells[f] > 0) any = true;
    if(!any) return;

    printf("\n      ATSAT: negative-value clamp, cumulative since start\n");
    for(int f = 0; f < nf; f++){
        if(clamp_cells[f] == 0) continue;
        double pos = 0.0;
        #pragma omp parallel for collapse(2) schedule(static) reduction(+:pos)
        for(int i = 0; i < im; i++)
            for(int j = 0; j < jm; j++)
                for(int k = 0; k < km; k++)
                    if(fields[f]->x[i][j][k] > 0.0) pos += fields[f]->x[i][j][k];
        printf("        %-12s gross %.4e over %10ld clippings = %8.4f %% of the current"
               " field mass (gross, not net — see the note in FileIO_Sat.cpp)\n",
               names[f], clamp_added[f], clamp_cells[f],
               (pos > 0.0) ? 100.0 * clamp_added[f] / pos : 0.0);
    }
}
/*
*
*/
