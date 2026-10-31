/*
 * Saturn General Circulation Modell(ATSAT)applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and nh3 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to write sequel, transfer and paraview files
*/

#include "cSaturnModel.h"
#include "ParaViewWriter.h"

using namespace std;
//using namespace AtomUtils;

// The five dumpers moved to the SHARED ParaViewWriter.h, as namespace ParaViewIO.
// They were already byte-identical in both models; see that header.

// ===== The fields ATJUP's ParaView writes and ATSAT's did not =====
//
// Mirrored from ParaView_Jup.cpp: the mixture density, the radiation pair, the eight
// precipitation fluxes and the six turbulence fields. Everything ATSAT already wrote is kept —
// the NH3 condensate pair and the w_/j_/jT_/massflux_/difflux_ quintets have no counterpart in
// ATJUP's writer, so this is the union of the two field sets and not a replacement.
//
// ONE list, not four. ATJUP repeats its field list once per slice writer, which is how its own
// radial and zonal outputs came to disagree (the radial one has NH3Cloud commented out and the
// zonal one does not). dump_radial, dump_zonal and dump_longal share a signature, so the list
// lives here once and each writer passes its own dumper and slice index. Adding a field is one
// line and it appears in every slice.
//
// UNITS. These streams run precision(4) with ios::fixed, so a raw SI value below 5e-5 prints as
// 0.0000 and the field is lost. Three groups are therefore scaled, with the unit in the comment:
//   - precipitation fluxes to mm/day (x86400 from kg/m2/s; 1 kg/m2 of water == 1 mm of depth),
//   - Q_rad and Q_precip to mW/m3,
//   - tke to m2/s2 (x u_0^2) and nue_t to m2/s (x u_0 * L_atm in metres).
// dis stays DIMENSIONLESS on purpose: its conversion depends on which closure ran (eps* uses
// u_0^3/L, omega* uses u_0/L) and the writer cannot know that. prod and the two source terms
// are left dimensionless for the same reason; their range prints without loss.
//
// All of these are ZERO unless the module that fills them is switched on — ATSAT_PRECIP for the
// P_* and Q_precip, ATSAT_TURB for the turbulence six, ATSAT_RADIATION for the radiation pair.
// They are written regardless, so that a run with the knob off shows a field of zeros rather
// than a missing array, which in ParaView looks the same as a broken writer.
#define DUMP_EXTRA_FIELDS_3D(DUMP, IDX, F)                                    \
    DUMP("rho_mix",       rho_mix,    1.0,     IDX, F);                       \
    DUMP("Q_rad",         Q_rad,      1.0e3,   IDX, F);   /* mW/m3  */        \
    DUMP("Radiation",     radiation,  1.0,     IDX, F);   /* W/m2   */        \
    DUMP("Emissivity",    epsilon,    1.0,     IDX, F);   /* 0..1   */        \
    DUMP("P_rain",        P_rain,     86400.0, IDX, F);   /* mm/day */        \
    DUMP("P_snow",        P_snow,     86400.0, IDX, F);                       \
    DUMP("P_graupel",     P_graupel,  86400.0, IDX, F);                       \
    DUMP("P_nh3_rain",    P_nh3_rain, 86400.0, IDX, F);                       \
    DUMP("P_nh3_snow",    P_nh3_snow, 86400.0, IDX, F);                       \
    DUMP("P_ch4_rain",    P_ch4_rain, 86400.0, IDX, F);                       \
    DUMP("P_ch4_snow",    P_ch4_snow, 86400.0, IDX, F);                       \
    DUMP("P_nh4sh",       P_nh4sh,    86400.0, IDX, F);                       \
    DUMP("Q_precip",      Q_precip,   1.0e3,   IDX, F);   /* mW/m3  */        \
    DUMP("tke",           tke,        u_0 * u_0, IDX, F); /* m2/s2  */        \
    DUMP("disd",          dis,        1.0,     IDX, F);   /* nondimensional */\
    DUMP("nue_t",         nue,        u_0 * L_atm * 1.0e3, IDX, F); /* m2/s */\
    DUMP("prod",          prod,       1.0,     IDX, F);                       \
    DUMP("tke_source",    tke_source, 1.0,     IDX, F);                       \
    DUMP("dis_source",    dis_source, 1.0,     IDX, F);

// The same list for the panorama .vts, whose dumper takes no slice index. The names carry their
// unit as a suffix here because the .vts header has to name every scalar in one attribute
// string, where a bare "P_rain" gives the reader no way to know it is not kg/m2/s. ATJUP's
// panorama uses exactly these names; keeping them identical is what lets one ParaView state
// file open a Jupiter and a Saturn panorama.
#define DUMP_EXTRA_FIELDS_VTS(F)                                              \
    dump_array("rho_mix",           rho_mix,       1.0,     F);               \
    dump_array("Q_rad_mW_m3",       Q_rad,         1.0e3,   F);               \
    dump_array("Radiation",         radiation,     1.0,     F);               \
    dump_array("Emissivity",        epsilon,       1.0,     F);               \
    dump_array("P_rain_mmd",        P_rain,        86400.0, F);               \
    dump_array("P_snow_mmd",        P_snow,        86400.0, F);               \
    dump_array("P_graupel_mmd",     P_graupel,     86400.0, F);               \
    dump_array("P_nh3_rain_mmd",    P_nh3_rain,    86400.0, F);               \
    dump_array("P_nh3_snow_mmd",    P_nh3_snow,    86400.0, F);               \
    dump_array("P_ch4_rain_mmd",    P_ch4_rain,    86400.0, F);               \
    dump_array("P_ch4_snow_mmd",    P_ch4_snow,    86400.0, F);               \
    dump_array("P_nh4sh_mmd",       P_nh4sh,       86400.0, F);               \
    dump_array("Q_precip_mW_m3",    Q_precip,      1.0e3,   F);               \
    dump_array("tke_m2s2",          tke,           u_0 * u_0, F);             \
    dump_array("dis_nd",            dis,           1.0,     F);               \
    dump_array("nue_t_m2s",         nue,           u_0 * L_atm * 1.0e3, F);   \
    dump_array("prod_nd",           prod,          1.0,     F);               \
    dump_array("tke_source_nd",     tke_source,    1.0,     F);               \
    dump_array("dis_source_nd",     dis_source,    1.0,     F);

// The scalar names the .vts header must announce, in one place so the header and the writer
// cannot drift apart. A name listed here but not written (or the reverse) is not an error
// ParaView reports — it simply shows an empty array.
#define PANORAMA_EXTRA_SCALARS                                                \
    "rho_mix Q_rad_mW_m3 Radiation Emissivity P_rain_mmd P_snow_mmd "         \
    "P_graupel_mmd "                                                          \
    "P_nh3_rain_mmd P_nh3_snow_mmd "                        \
    "P_ch4_rain_mmd P_ch4_snow_mmd P_nh4sh_mmd "            \
    "Q_precip_mW_m3 tke_m2s2 dis_nd nue_t_m2s prod_nd tke_source_nd "         \
    "dis_source_nd "
/*
 * 
*/
void cSaturnModel::paraview_panorama_vts(int n){
    using namespace ParaViewIO;
    // Header, coordinates, the Velocity array and the Temperature array are the
    // SHARED ParaViewWriter.h. What stays here is the field list below and the
    // scalars string that has to agree with it.
    ParaViewWriter<cSaturnModel> pv(*this);
    // ATSAT_CHROM (ChromSat.h): the chromophore arrays only when the tracer is on, so a default panorama is unchanged.
    // The list names exactly the arrays written below (2026-10-31: it announced PressureDynamic, PressureStatic and
    // Q_Sensible, which were never written, and left out seven that were).
    const std::string pano_scalars = std::string(
        "Temperature u-component v-component w-component PressureDyn CoriolisForce BuoyancyForce "
        "CH4 CH4Cloud CH4Ice H2O H2OCloud H2OIce H2S w_h2s NH3 NH3Cloud NH3Ice w_nh3 NH4SH w_nh4sh ")
        + (chrom_on() ? "Chromophore_mgm3 ChromUV " : "") + "Q_Latent " PANORAMA_EXTRA_SCALARS;
    ofstream Saturn_panorama_vts_File = pv.open_panorama(n, pano_scalars.c_str());
    pv.panorama_velocity(Saturn_panorama_vts_File);
    pv.panorama_temperature(Saturn_panorama_vts_File);
    dump_array("u-component", u, u_0, Saturn_panorama_vts_File);
    dump_array("v-component", v, u_0, Saturn_panorama_vts_File);
    dump_array("w-component", w, u_0, Saturn_panorama_vts_File);

    dump_array("PressureDyn", p_dyn, 1.0, Saturn_panorama_vts_File);

    dump_array("CoriolisForce", CoriolisForce, 1.0, Saturn_panorama_vts_File);
    dump_array("BuoyancyForce", BuoyancyForce, 1.0, Saturn_panorama_vts_File);

    dump_array("CH4", ch4, 1e3, Saturn_panorama_vts_File);
    dump_array("CH4Cloud", ch4_cloud, 1e3, Saturn_panorama_vts_File);
    dump_array("CH4Ice", ch4_ice, 1e3, Saturn_panorama_vts_File);

    dump_array("H2O", h2o, 1e3, Saturn_panorama_vts_File);
    dump_array("H2OCloud", h2o_cloud, 1e3, Saturn_panorama_vts_File);
    dump_array("H2OIce", h2o_ice, 1e3, Saturn_panorama_vts_File);

    dump_array("H2S", h2s, 1e3, Saturn_panorama_vts_File);
    dump_array("w_h2s", w_h2s, 1e3, Saturn_panorama_vts_File);

    dump_array("NH3", nh3, 1e3, Saturn_panorama_vts_File);
    dump_array("NH3Cloud", nh3_cloud, 1e3, Saturn_panorama_vts_File);
    dump_array("NH3Ice", nh3_ice, 1e3, Saturn_panorama_vts_File);
    dump_array("w_nh3", w_nh3, 1e3, Saturn_panorama_vts_File);

    dump_array("NH4SH", nh4sh, 1e9, Saturn_panorama_vts_File);
    dump_array("w_nh4sh", w_nh4sh, 1e9, Saturn_panorama_vts_File);
    if(chrom_on()){ dump_array("Chromophore_mgm3", chrom, 1e6, Saturn_panorama_vts_File);
                    dump_array("ChromUV", chrom_uv, 1.0, Saturn_panorama_vts_File); }

    dump_array("Q_Latent", Q_Latent, 1.0, Saturn_panorama_vts_File);

    DUMP_EXTRA_FIELDS_VTS(Saturn_panorama_vts_File)

    pv.close_panorama(Saturn_panorama_vts_File, n);
    return;
}
/*
 * Saturn_cloudview_<n>.vtk (2026-10-30, ATSAT_CLOUDVIEW_STRIDE; ported from ATJUP's Jupiter_cloudview, d18961a /
 * 5593beb / e930d03): what a camera above the planet would see of the two upper cloud decks. Every radial
 * layer i is a horizontal (lat, lon) sheet holding the column mass ABOVE it, summed from the top down to and
 * including layer i:
 *
 *     X_above(i, j, k) = sum_{i' = i .. im-2} X(i', j, k) * layer_thickness_m(i')       [g/m2]
 *
 * so layer 0 is the whole column (the photo's view), and stepping i upward strips the deeper cloud away.
 * Layer im-1, the lid, has no thickness and holds 0. NH3Ice is the upper white deck; NH4SH the deeper layer
 * that shows through where the NH3 deck is thin. Saturn's deck structure is Jupiter's, so the arrays are the
 * same. Grid: (lat, lon) layout (x = 0.1 per latitude row, y = 0.1 per longitude), the layers stacked in z
 * at 0.1 each and FLIPPED -- layer i at z = 0.1 (im-1-i), so layer 0 is the TOP face and a view from above
 * shows the photo at once. With ATSAT_CHROM the chromophore column Chrom_above_mgm2 (mg/m2) is added.
 */
void cSaturnModel::paraview_cloudview(int n){
    const std::string name = output_path + "/Saturn_cloudview_" + std::to_string(n) + ".vtk";
    std::ofstream f;
    f.precision(4);
    f.setf(std::ios::fixed);
    f.open(name);
    if(!f.is_open()){
        std::cerr << "ERROR: could not open cloudview file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    const size_t npt = (size_t)im * jm * km;
    f << "# vtk DataFile Version 3.0" << std::endl;
    f << "Cloudview_Data_Sat_Circulation\n";
    f << "ASCII" << std::endl;
    f << "DATASET STRUCTURED_GRID" << std::endl;
    f << "DIMENSIONS " << km << " " << jm << " " << im << std::endl;
    f << "POINTS " << npt << " float" << std::endl;
    for(int i = 0; i < im; i++)
        for(int j = 0; j < jm; j++)
            for(int k = 0; k < km; k++)
                f << j * 0.1 << " " << k * 0.1 << " " << (im - 1 - i) * 0.1 << std::endl;   // layer 0 on top
    f << "POINT_DATA " << npt << std::endl;
    auto above = [&](const char* desc, Array& a, double unit = 1.0e3){
        std::vector<double> c((size_t)im * jm * km, 0.0);          // c[i][j][k], layer im-1 stays 0
        for(int j = 0; j < jm; j++) for(int k = 0; k < km; k++){
            double s = 0.0;
            for(int i = im - 2; i >= 0; i--){
                s += a.x[i][j][k] * layer_thickness_m(i) * unit;     // kg/m3 * m -> g/m2 (unit 1e3)
                c[((size_t)i * jm + j) * km + k] = s;
            }
        }
        f << "SCALARS " << desc << " float " << 1 << std::endl;
        f << "LOOKUP_TABLE default" << std::endl;
        for(size_t q = 0; q < c.size(); q++) f << c[q] << std::endl;
        return c;
    };
    const std::vector<double> c_ice = above("NH3Ice_above_gm2", nh3_ice);
    const std::vector<double> c_nhs = above("NH4SH_above_gm2", nh4sh);
    std::vector<double> c_chr;
    if(chrom_on()) c_chr = above("Chrom_above_mgm2", chrom, 1.0e6);     // ATSAT_CHROM, ChromSat.h: the colour

    // PhotoRGB: one colour per point, red / green / blue in 0..1 (ParaView: colour by PhotoRGB, untick
    // "Map Scalars"). A RECIPE, not radiative transfer -- ATJUP's, with Saturn's paler tints:
    //   b  brightness of the NH3 ice column             -> pale cream deck, B = 0.15 + 0.85 b
    //   n  NH4SH column where the NH3 deck is thin, n (1 - b) -> mixed toward tan (0.78, 0.66, 0.46)
    //   c  chromophore column                            -> mixed toward butterscotch (0.86, 0.66, 0.34)
    // Weaker mixing than on Jupiter (0.5 / 0.6 against 0.7 / 0.8): Saturn's haze mutes the contrasts.
    // Each column scaled between its 2nd and 98th percentile of LAYER 0 over |lat| <= 80 deg (the
    // chromophore between its 2nd percentile and its maximum), so the same scale holds on every layer.
    // The scale is per file: two files are not on one scale.
    {
        auto range = [&](const std::vector<double>& c, bool to_max){
            std::vector<double> v;
            for(int j = 10; j < jm - 10; j++) for(int k = 0; k < km - 1; k++) v.push_back(c[(size_t)j * km + k]);
            std::sort(v.begin(), v.end());
            const double lo = v[v.size() * 2 / 100], hi = to_max ? v.back() : v[v.size() * 98 / 100];
            return std::make_pair(lo, hi > lo ? hi : lo + 1.0);
        };
        auto unit01 = [](double x, const std::pair<double, double>& r){
            return std::min(1.0, std::max(0.0, (x - r.first) / (r.second - r.first))); };
        const auto r_ice = range(c_ice, false), r_nhs = range(c_nhs, false);
        const auto r_chr = c_chr.empty() ? std::make_pair(0.0, 1.0) : range(c_chr, true);
        f << "COLOR_SCALARS PhotoRGB 3" << std::endl;
        for(size_t q = 0; q < c_ice.size(); q++){
            const double b  = unit01(c_ice[q], r_ice);
            const double nb = unit01(c_nhs[q], r_nhs) * (1.0 - b);
            const double cc = c_chr.empty() ? 0.0 : unit01(c_chr[q], r_chr);
            double rgb[3] = { 1.0, 0.97, 0.90 };
            const double tan_[3] = { 0.78, 0.66, 0.46 }, butter[3] = { 0.86, 0.66, 0.34 };
            for(int m_ = 0; m_ < 3; m_++){
                rgb[m_] += 0.5 * nb * (tan_[m_] - rgb[m_]);
                rgb[m_] += 0.6 * cc * (butter[m_] - rgb[m_]);
                rgb[m_] *= 0.15 + 0.85 * b;
            }
            f << rgb[0] << " " << rgb[1] << " " << rgb[2] << std::endl;
        }
    }
    f.close();
    std::cout << "   File:  Saturn_cloudview_" << n << ".vtk"
              << "  has been written to Directory:  " << output_path << std::endl;
}
/*
 * 
*/
void cSaturnModel::paraview_vtk_radial(int n, int i_radial){
    using namespace ParaViewIO;
    ofstream Saturn_vtk_radial_File = ParaViewWriter<cSaturnModel>(*this)
        .open_slice("radial", "Radial", i_radial, n, km, jm, 0.1, false);
    // Third component of the in-plane vector below. A radial slice is a (j,k) plane, so the
    // out-of-plane component is flat zero; it used to be the geometry loop's leftover z.
    const double z = 0.0;
    dump_radial("u-Component", u, u_0, i_radial, Saturn_vtk_radial_File);
    dump_radial("v-Component", v, u_0, i_radial, Saturn_vtk_radial_File);
    dump_radial("w-Component", w, u_0, i_radial, Saturn_vtk_radial_File);
    Saturn_vtk_radial_File <<  "SCALARS Temperature float " << 1 << endl;
    Saturn_vtk_radial_File <<  "LOOKUP_TABLE default"  <<endl;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Saturn_vtk_radial_File << t.x[i_radial][j][k] * t_ref/10.0 << endl;
        }
    }

    dump_radial("thermalmassflux", thermalmassflux, 1.0, i_radial, Saturn_vtk_radial_File);
    dump_radial("CH4", ch4, 1e3, i_radial, Saturn_vtk_radial_File);
    dump_radial("CH4Cloud", ch4_cloud, 1e3, i_radial, Saturn_vtk_radial_File);
    dump_radial("CH4Ice", ch4_ice, 1e3, i_radial, Saturn_vtk_radial_File);

    dump_radial("H2O", h2o, 1e3, i_radial, Saturn_vtk_radial_File);
    dump_radial("H2OCloud", h2o_cloud, 1e3, i_radial, Saturn_vtk_radial_File);
    dump_radial("H2OIce", h2o_ice, 1e3, i_radial, Saturn_vtk_radial_File);

    dump_radial("H2S", h2s, 1e3, i_radial, Saturn_vtk_radial_File);
    dump_radial("w_h2s", w_h2s, 1e3, i_radial, Saturn_vtk_radial_File);
    dump_radial("j_h2s", j_h2s, 1e3, i_radial, Saturn_vtk_radial_File);
    dump_radial("jT_h2s", jT_h2s, 1e3, i_radial, Saturn_vtk_radial_File);
    dump_radial("massflux_h2s", massflux_h2s, 1e3, i_radial, Saturn_vtk_radial_File);
    dump_radial("difflux_h2s", difflux_h2s, 1e3, i_radial, Saturn_vtk_radial_File);

    dump_radial("NH3", nh3, 1e3, i_radial, Saturn_vtk_radial_File);
    dump_radial("NH3Cloud", nh3_cloud, 1e3, i_radial, Saturn_vtk_radial_File);
    dump_radial("NH3Ice", nh3_ice, 1e3, i_radial, Saturn_vtk_radial_File);
    dump_radial("w_nh3", w_nh3, 1e3, i_radial, Saturn_vtk_radial_File);
    dump_radial("j_nh3", j_nh3, 1e3, i_radial, Saturn_vtk_radial_File);
    dump_radial("jT_nh3", jT_nh3, 1e3, i_radial, Saturn_vtk_radial_File);
    dump_radial("massflux_nh3", massflux_nh3, 1e3, i_radial, Saturn_vtk_radial_File);
    dump_radial("difflux_nh3", difflux_nh3, 1e3, i_radial, Saturn_vtk_radial_File);

    dump_radial("NH4SH", nh4sh, 1e9, i_radial, Saturn_vtk_radial_File);
    if(chrom_on()){ dump_radial("Chromophore_mgm3", chrom, 1e6, i_radial, Saturn_vtk_radial_File);
                    dump_radial("ChromUV", chrom_uv, 1.0, i_radial, Saturn_vtk_radial_File); }
    dump_radial("w_nh4sh", w_nh4sh, 1e9, i_radial, Saturn_vtk_radial_File);
    dump_radial("j_nh4sh", j_nh4sh, 1e9, i_radial, Saturn_vtk_radial_File);
    dump_radial("jT_nh4sh", jT_nh4sh, 1e9, i_radial, Saturn_vtk_radial_File);
    dump_radial("massflux_nh4sh", massflux_nh4sh, 1e9, i_radial, Saturn_vtk_radial_File);
    dump_radial("difflux_nh4sh", difflux_nh4sh, 1e9, i_radial, Saturn_vtk_radial_File);


    dump_radial("PressureDyn", p_dyn, 1.0, i_radial, Saturn_vtk_radial_File);
    dump_radial("PressureStat", p_stat, 1.0, i_radial, Saturn_vtk_radial_File);
    dump_radial("rho_mix", rho_mix, 1.0, i_radial, Saturn_vtk_radial_File);

    dump_radial("CoriolisForce", CoriolisForce, 1.0, i_radial, Saturn_vtk_radial_File);
    dump_radial("CentrifugalForce", CentrifugalForce, 1e3, i_radial, Saturn_vtk_radial_File);
    dump_radial("BuoyancyForce", BuoyancyForce, 1.0, i_radial, Saturn_vtk_radial_File);
    dump_radial("PresGradForce", PresGradForce, 1.0, i_radial, Saturn_vtk_radial_File);

    dump_radial("Q_Latent", Q_Latent, 1.0, i_radial, Saturn_vtk_radial_File);
    dump_radial("Q_Sensible", Q_Sensible, 1.0, i_radial, Saturn_vtk_radial_File);
    DUMP_EXTRA_FIELDS_3D(dump_radial, i_radial, Saturn_vtk_radial_File)

    // All-species SURFACE precipitation map, mm/day. The fields above are sampled at this
    // file's one altitude i_radial and so miss any deck that does not sit there; these are 2D
    // and always taken at the base of the column, so they are the condensate mass flux actually
    // arriving at the bottom, with its per-species breakdown. A radial slice is a (j,k) plane,
    // which is the shape of these arrays — that is why they appear here and not in the zonal or
    // longitudinal writers.
    dump_radial_2d("Precip_total", precip_srf_total, 86400.0, Saturn_vtk_radial_File);
    dump_radial_2d("Precip_h2o",   precip_srf_h2o,   86400.0, Saturn_vtk_radial_File);
    dump_radial_2d("Precip_nh3",   precip_srf_nh3,   86400.0, Saturn_vtk_radial_File);
    dump_radial_2d("Precip_ch4",   precip_srf_ch4,   86400.0, Saturn_vtk_radial_File);
    dump_radial_2d("Precip_nh4sh", precip_srf_nh4sh, 86400.0, Saturn_vtk_radial_File);

    // Per-column friction velocity u_tau, from TurbulenceSat::compute_vel_star.
    dump_radial_2d("vel_star_ms", vel_star, 1.0, Saturn_vtk_radial_File);

    Saturn_vtk_radial_File <<  "VECTORS v-w-Cell float " << endl;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            Saturn_vtk_radial_File << v.x[i_radial][j][k] << " " << w.x[i_radial][j][k] << " " << z << endl;
        }
    }
    ParaViewWriter<cSaturnModel>(*this).close_slice(Saturn_vtk_radial_File, "radial", i_radial, n);
    return;
}
/*
 * 
*/
void cSaturnModel::paraview_vtk_zonal(int n, int k_zonal){
    using namespace ParaViewIO;
    ofstream Saturn_vtk_zonal_File = ParaViewWriter<cSaturnModel>(*this)
        .open_slice("zonal", "Zonal", k_zonal, n, jm, im, 0.05, false);
    const double z = 0.0;   // out-of-plane component, as in the radial writer
    dump_zonal("u-Component", u, u_0, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("v-Component", v, u_0, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("w-Component", w, u_0, k_zonal, Saturn_vtk_zonal_File);
    Saturn_vtk_zonal_File <<  "SCALARS Temperature float " << 1 << endl;
    Saturn_vtk_zonal_File <<  "LOOKUP_TABLE default"  <<endl;

    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            Saturn_vtk_zonal_File << t.x[i][j][k_zonal] * t_ref/10.0 << endl;
            aux.x[i][j][k_zonal] = get_layer_height(i);
        }
    }
    dump_zonal("thermalmassflux", thermalmassflux, 1.0, k_zonal, Saturn_vtk_zonal_File);

    dump_zonal("height", aux, 1.0, k_zonal, Saturn_vtk_zonal_File);

    dump_zonal("CH4", ch4, 1e3, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("CH4Cloud", ch4_cloud, 1e3, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("CH4Ice", ch4_ice, 1e3, k_zonal, Saturn_vtk_zonal_File);

    dump_zonal("H2O", h2o, 1e3, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("H2OCloud", h2o_cloud, 1e3, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("H2OIce", h2o_ice, 1e3, k_zonal, Saturn_vtk_zonal_File);

    dump_zonal("H2S", h2s, 1e3, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("w_h2s", w_h2s, 1e3, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("j_h2s", j_h2s, 1e3, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("jT_h2s", jT_h2s, 1e3, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("massflux_h2s", massflux_h2s, 1e3, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("difflux_h2s", difflux_h2s, 1e3, k_zonal, Saturn_vtk_zonal_File);

    dump_zonal("NH3", nh3, 1e3, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("NH3Cloud", nh3_cloud, 1e3, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("NH3Ice", nh3_ice, 1e3, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("w_nh3", w_nh3, 1e3, k_zonal, Saturn_vtk_zonal_File);
   dump_zonal("j_nh3", j_nh3, 1e3, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("jT_nh3", jT_nh3, 1e3, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("massflux_nh3", massflux_nh3, 1e3, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("difflux_nh3", difflux_nh3, 1e3, k_zonal, Saturn_vtk_zonal_File);

    dump_zonal("NH4SH", nh4sh, 1e9, k_zonal, Saturn_vtk_zonal_File);
    if(chrom_on()){ dump_zonal("Chromophore_mgm3", chrom, 1e6, k_zonal, Saturn_vtk_zonal_File);
                    dump_zonal("ChromUV", chrom_uv, 1.0, k_zonal, Saturn_vtk_zonal_File); }
    dump_zonal("w_nh4sh", w_nh4sh, 1e9, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("j_nh4sh", j_nh4sh, 1e9, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("jT_nh4sh", jT_nh4sh, 1e9, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("massflux_nh4sh", massflux_nh4sh, 1e9, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("difflux_nh4sh", difflux_nh4sh, 1e9, k_zonal, Saturn_vtk_zonal_File);

    dump_zonal("PressureDyn", p_dyn, 1.0, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("PressureStat", p_stat, 1.0, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("rho_mix", rho_mix, 1.0, k_zonal, Saturn_vtk_zonal_File);

    dump_zonal("CoriolisForce", CoriolisForce, 1.0, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("CentrifugalForce", CentrifugalForce, 1e3, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("BuoyancyForce", BuoyancyForce, 1.0, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("PresGradForce", PresGradForce, 1.0, k_zonal, Saturn_vtk_zonal_File);

    dump_zonal("Q_Latent", Q_Latent, 1.0, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("Q_Sensible", Q_Sensible, 1.0, k_zonal, Saturn_vtk_zonal_File);
    DUMP_EXTRA_FIELDS_3D(dump_zonal, k_zonal, Saturn_vtk_zonal_File)

    Saturn_vtk_zonal_File <<  "VECTORS u-v-Cell float" << endl;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            Saturn_vtk_zonal_File << u.x[i][j][k_zonal] << " " << v.x[i][j][k_zonal] << " " << z << endl;
        }
    }
    ParaViewWriter<cSaturnModel>(*this).close_slice(Saturn_vtk_zonal_File, "zonal", k_zonal, n);
    return;
}
/*
 * 
*/
void cSaturnModel::paraview_vtk_longal(int n, int j_longal){
    using namespace ParaViewIO;
    ofstream Saturn_vtk_longal_File = ParaViewWriter<cSaturnModel>(*this)
        .open_slice("longal", "Longitudinal", j_longal, n, km, im, 0.025, true);
    const double y = 0.0;   // out-of-plane component; longal advances z, so y stayed 0
    dump_longal("u-Component", u, u_0, j_longal, Saturn_vtk_longal_File);
    dump_longal("v-Component", v, u_0, j_longal, Saturn_vtk_longal_File);
    dump_longal("w-Component", w, u_0, j_longal, Saturn_vtk_longal_File);
    Saturn_vtk_longal_File <<  "SCALARS Temperature float " << 1 << endl;
    Saturn_vtk_longal_File <<  "LOOKUP_TABLE default"  <<endl;

    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            Saturn_vtk_longal_File << t.x[i][j_longal][k] * t_ref/10.0 << endl;
            aux.x[i][j_longal][k] = get_layer_height(i);
        }
    }

    dump_longal("thermalmassflux", thermalmassflux, 1.0, j_longal, Saturn_vtk_longal_File);

    dump_longal("height", aux, 1.0, j_longal, Saturn_vtk_longal_File);
    dump_longal("CH4", ch4, 1e3, j_longal, Saturn_vtk_longal_File);
    dump_longal("CH4Cloud", ch4_cloud, 1e3, j_longal, Saturn_vtk_longal_File);
    dump_longal("CH4Ice", ch4_ice, 1e3, j_longal, Saturn_vtk_longal_File);

    dump_longal("H2O", h2o, 1e3, j_longal, Saturn_vtk_longal_File);
    dump_longal("H2OCloud", h2o_cloud, 1e3, j_longal, Saturn_vtk_longal_File);
    dump_longal("H2OIce", h2o_ice, 1e3, j_longal, Saturn_vtk_longal_File);

    dump_longal("H2S", h2s, 1e3, j_longal, Saturn_vtk_longal_File);
    dump_longal("w_h2s", w_h2s, 1e3, j_longal, Saturn_vtk_longal_File);
    dump_longal("j_h2s", j_h2s, 1e3, j_longal, Saturn_vtk_longal_File);
    dump_longal("jT_h2s", jT_h2s, 1e3, j_longal, Saturn_vtk_longal_File);
    dump_longal("massflux_h2s", massflux_h2s, 1e3, j_longal, Saturn_vtk_longal_File);
    dump_longal("difflux_h2s", difflux_h2s, 1e3, j_longal, Saturn_vtk_longal_File);

    dump_longal("NH3", nh3, 1e3, j_longal, Saturn_vtk_longal_File);
    dump_longal("NH3Cloud", nh3_cloud, 1e3, j_longal, Saturn_vtk_longal_File);
    dump_longal("NH3Ice", nh3_ice, 1e3, j_longal, Saturn_vtk_longal_File);
    dump_longal("w_nh3", w_nh3, 1e3, j_longal, Saturn_vtk_longal_File);
    dump_longal("j_nh3", j_nh3, 1e3, j_longal, Saturn_vtk_longal_File);
    dump_longal("jT_nh3", jT_nh3, 1e3, j_longal, Saturn_vtk_longal_File);
    dump_longal("massflux_nh3", massflux_nh3, 1e3, j_longal, Saturn_vtk_longal_File);
    dump_longal("difflux_nh3", difflux_nh3, 1e3, j_longal, Saturn_vtk_longal_File);

    dump_longal("NH4SH", nh4sh, 1e9, j_longal, Saturn_vtk_longal_File);
    if(chrom_on()){ dump_longal("Chromophore_mgm3", chrom, 1e6, j_longal, Saturn_vtk_longal_File);
                    dump_longal("ChromUV", chrom_uv, 1.0, j_longal, Saturn_vtk_longal_File); }
    dump_longal("w_nh4sh", w_nh4sh, 1e9, j_longal, Saturn_vtk_longal_File);
    dump_longal("j_nh4sh", j_nh4sh, 1e9, j_longal, Saturn_vtk_longal_File);
    dump_longal("jT_nh4sh", jT_nh4sh, 1e9, j_longal, Saturn_vtk_longal_File);
    dump_longal("massflux_nh4sh", massflux_nh4sh, 1e9, j_longal, Saturn_vtk_longal_File);
    dump_longal("difflux_nh4sh", difflux_nh4sh, 1e9, j_longal, Saturn_vtk_longal_File);

    dump_longal("PressureDyn", p_dyn, 1.0, j_longal, Saturn_vtk_longal_File);
    dump_longal("PressureStat", p_stat, 1.0, j_longal, Saturn_vtk_longal_File);

    dump_longal("CoriolisForce", CoriolisForce, 1.0, j_longal, Saturn_vtk_longal_File);
    dump_longal("CentrifugalForce", CentrifugalForce, 1e3, j_longal, Saturn_vtk_longal_File);
    dump_longal("BuoyancyForce", BuoyancyForce, 1.0, j_longal, Saturn_vtk_longal_File);
    dump_longal("PresGradForce", PresGradForce, 1.0, j_longal, Saturn_vtk_longal_File);

    dump_longal("Q_Latent", Q_Latent, 1.0, j_longal, Saturn_vtk_longal_File);
    dump_longal("Q_Sensible", Q_Sensible, 1.0, j_longal, Saturn_vtk_longal_File);
    DUMP_EXTRA_FIELDS_3D(dump_longal, j_longal, Saturn_vtk_longal_File)

    Saturn_vtk_longal_File <<  "VECTORS u-w-Cell float" << endl;
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            Saturn_vtk_longal_File << u.x[i][j_longal][k] << " " 
                << y << " " << w.x[i][j_longal][k] << endl;
        }
    }
    ParaViewWriter<cSaturnModel>(*this).close_slice(Saturn_vtk_longal_File, "longal", j_longal, n);
    return;
}
/*
 * 
*/
void cSaturnModel::paraview_sphere_vts(int n){
    using namespace ParaViewIO;
    double x, y, z, sinthe, sinphi, costhe, cosphi;
    string Saturn_sphere_vts_File_Name = output_path + "/Saturn_sphere_" 
        + std::to_string(n) + ".vts";
    string file_name = "Saturn_sphere_" + std::to_string(n) + ".vts";
    ofstream Saturn_sphere_vts_File;
    Saturn_sphere_vts_File.precision(4);
    Saturn_sphere_vts_File.setf(ios::fixed);
    Saturn_sphere_vts_File.open(Saturn_sphere_vts_File_Name);
    if (!Saturn_sphere_vts_File.is_open()){
        cerr << "ERROR: could not open paraview_vts file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    Saturn_sphere_vts_File <<  "<?xml version=\"1.0\"?>\n"  << endl;
    Saturn_sphere_vts_File <<  "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"  << endl;
    Saturn_sphere_vts_File <<  " <StructuredGrid WholeExtent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Saturn_sphere_vts_File <<  "  <Piece Extent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Saturn_sphere_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"Temperature PressureDyn H2O H2S NH3 NH4SH H2OCloud NH3Cloud H2OIce NH3Ice CH4 CH4Cloud CH4Ice u-Component v-Component w-Component\">\n"  << endl;
    Saturn_sphere_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        sinphi = sin( phi.z[k]);
        cosphi = cos( phi.z[k]);
        for(int j = 0; j < jm; j++){
            sinthe = sin(the.z[j]);
            costhe = cos(the.z[j]);
            for(int i = 0; i < im; i++){
                // LOCALS, not aux_u/aux_v/aux_w. These are the spherical-to-Cartesian velocity
                // components and they are wanted only for the line printed immediately below —
                // but aux_* are not scratch. RHS_Sat_Turb.cpp:760 fills aux_u with the
                // intermediate velocity rhs_u + dpdr and PressureSolver differentiates it, so
                // writing them here left the projection reading rendering coordinates.
                //
                // MEASURED, and this is why the writer's call site was commented out rather than
                // the writer fixed: with the call enabled and aux_* used as below, three
                // iterations gave sat_restart_3.bin = a585481f... against b64cb36f... for the
                // same run without it. A control with the panorama writing every iteration and
                // only the sphere disabled reproduced b64cb36f... exactly, so the panorama is
                // innocent and this loop was the whole difference. An output routine that
                // changes the answer is worse than no output routine.
                const double cart_u = sinthe * cosphi * u.x[i][j][k] + costhe * cosphi * v.x[i][j][k] - sinphi * w.x[i][j][k];
                const double cart_v = sinthe * sinphi * u.x[i][j][k] + sinphi * costhe * v.x[i][j][k] + cosphi * w.x[i][j][k];
                const double cart_w = costhe * u.x[i][j][k] - sinthe * v.x[i][j][k];
                Saturn_sphere_vts_File << cart_u << " " << cart_v << " " << cart_w  << endl;
            }
            Saturn_sphere_vts_File <<  "\n"  << endl;
        }
        Saturn_sphere_vts_File <<  "\n"  << endl;
    }
    Saturn_sphere_vts_File <<  "\n"  << endl;
    Saturn_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Saturn_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
//                Saturn_sphere_vts_File << t.x[i][j][k] * t_ref << endl;
                Saturn_sphere_vts_File << t.x[i][j][k] * t_ref/10.0 << endl;
            }
            Saturn_sphere_vts_File <<  "\n"  << endl;
        }
        Saturn_sphere_vts_File <<  "\n"  << endl;
    }
    Saturn_sphere_vts_File <<  "\n"  << endl;
    Saturn_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Saturn_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"PressureDyn\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Saturn_sphere_vts_File << p_dyn.x[i][j][k] << endl;
            }
            Saturn_sphere_vts_File <<  "\n"  << endl;
        }
        Saturn_sphere_vts_File <<  "\n"  << endl;
    }
/*
    Saturn_sphere_vts_File <<  "\n"  << endl;
    Saturn_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Saturn_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"PressureStat\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Saturn_sphere_vts_File << 1e-3 * p_stat.x[i][j][k] << endl;
            }
            Saturn_sphere_vts_File <<  "\n"  << endl;
        }
        Saturn_sphere_vts_File <<  "\n"  << endl;
    }
*/
    Saturn_sphere_vts_File <<  "\n"  << endl;
    Saturn_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Saturn_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"H2O\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Saturn_sphere_vts_File << 1e3 * h2o.x[i][j][k] << endl;
            }
            Saturn_sphere_vts_File <<  "\n"  << endl;
        }
        Saturn_sphere_vts_File <<  "\n"  << endl;
    }
    Saturn_sphere_vts_File <<  "\n"  << endl;
    Saturn_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Saturn_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"H2S\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Saturn_sphere_vts_File << 1e3 * h2s.x[i][j][k] << endl;
            }
            Saturn_sphere_vts_File <<  "\n"  << endl;
        }
        Saturn_sphere_vts_File <<  "\n"  << endl;
    }
    Saturn_sphere_vts_File <<  "\n"  << endl;
    Saturn_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Saturn_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"NH3\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Saturn_sphere_vts_File << 1e3 * nh3.x[i][j][k] << endl;
            }
            Saturn_sphere_vts_File <<  "\n"  << endl;
        }
        Saturn_sphere_vts_File <<  "\n"  << endl;
    }
    Saturn_sphere_vts_File <<  "\n"  << endl;
    Saturn_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Saturn_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"NH4SH\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Saturn_sphere_vts_File << 1e3 * nh4sh.x[i][j][k] << endl;
            }
            Saturn_sphere_vts_File <<  "\n"  << endl;
        }
        Saturn_sphere_vts_File <<  "\n"  << endl;
    }
    Saturn_sphere_vts_File <<  "\n"  << endl;
    Saturn_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Saturn_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"H2OCloud\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Saturn_sphere_vts_File << 1e3 * h2o_cloud.x[i][j][k] << endl;
            }
            Saturn_sphere_vts_File <<  "\n"  << endl;
        }
        Saturn_sphere_vts_File <<  "\n"  << endl;
    }

    Saturn_sphere_vts_File <<  "\n"  << endl;
    Saturn_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Saturn_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"NH3Cloud\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Saturn_sphere_vts_File << 1e3 * nh3_cloud.x[i][j][k] << endl;
            }
            Saturn_sphere_vts_File <<  "\n"  << endl;
        }
        Saturn_sphere_vts_File <<  "\n"  << endl;
    }

    Saturn_sphere_vts_File <<  "\n"  << endl;

    Saturn_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Saturn_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"H2OIce\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Saturn_sphere_vts_File << 1e3 * h2o_ice.x[i][j][k] << endl;
            }
            Saturn_sphere_vts_File <<  "\n"  << endl;
        }
        Saturn_sphere_vts_File <<  "\n"  << endl;
    }

    Saturn_sphere_vts_File <<  "\n"  << endl;
    Saturn_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Saturn_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"NH3Ice\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Saturn_sphere_vts_File << 1e3 * nh3_ice.x[i][j][k] << endl;
            }
            Saturn_sphere_vts_File <<  "\n"  << endl;
        }
        Saturn_sphere_vts_File <<  "\n"  << endl;
    }
    Saturn_sphere_vts_File <<  "\n"  << endl;

    Saturn_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Saturn_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"CH4\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Saturn_sphere_vts_File << 1e3 * ch4.x[i][j][k] << endl;
            }
            Saturn_sphere_vts_File <<  "\n"  << endl;
        }
        Saturn_sphere_vts_File <<  "\n"  << endl;
    }
    Saturn_sphere_vts_File <<  "\n"  << endl;
    Saturn_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Saturn_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"CH4Cloud\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Saturn_sphere_vts_File << 1e3 * ch4_cloud.x[i][j][k] << endl;
            }
            Saturn_sphere_vts_File <<  "\n"  << endl;
        }
        Saturn_sphere_vts_File <<  "\n"  << endl;
    }
    Saturn_sphere_vts_File <<  "\n"  << endl;
    Saturn_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Saturn_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"CH4Ice\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Saturn_sphere_vts_File << 1e3 * ch4_ice.x[i][j][k] << endl;
            }
            Saturn_sphere_vts_File <<  "\n"  << endl;
        }
        Saturn_sphere_vts_File <<  "\n"  << endl;
    }
    Saturn_sphere_vts_File <<  "\n"  << endl;
    Saturn_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Saturn_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"u-Component\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Saturn_sphere_vts_File << u_0 * aux_u.x[i][j][k] << endl;
            }
            Saturn_sphere_vts_File <<  "\n"  << endl;
        }
        Saturn_sphere_vts_File <<  "\n"  << endl;
    }
    Saturn_sphere_vts_File <<  "\n"  << endl;
    Saturn_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Saturn_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"v-Component\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Saturn_sphere_vts_File << u_0 * aux_v.x[i][j][k] << endl;
            }
            Saturn_sphere_vts_File <<  "\n"  << endl;
        }
        Saturn_sphere_vts_File <<  "\n"  << endl;
    }
    Saturn_sphere_vts_File <<  "\n"  << endl;
    Saturn_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Saturn_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"w-Component\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Saturn_sphere_vts_File << u_0 * aux_w.x[i][j][k] << endl;
            }
            Saturn_sphere_vts_File <<  "\n"  << endl;
        }
        Saturn_sphere_vts_File <<  "\n"  << endl;
    }
/*
    Saturn_sphere_vts_File <<  "\n"  << endl;
    Saturn_sphere_vts_File <<  "    </DataArray>\n" << endl;
    Saturn_sphere_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Seamount\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Saturn_sphere_vts_File << SeaMount.x[i][j][k] << endl;
            }
            Saturn_sphere_vts_File <<  "\n"  << endl;
        }
        Saturn_sphere_vts_File <<  "\n"  << endl;
    }
*/
    Saturn_sphere_vts_File <<  "\n"  << endl;
    Saturn_sphere_vts_File <<  "    </DataArray>\n"  << endl;
    Saturn_sphere_vts_File <<  "   </PointData>\n" << endl;
    Saturn_sphere_vts_File <<  "   <Points>\n"  << endl;
    Saturn_sphere_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                x = rad.z[i] * sin( the.z[j])* cos(phi.z[k]);
                y = rad.z[i] * sin( the.z[j])* sin(phi.z[k]);
                z = rad.z[i] * cos( the.z[j]);
                Saturn_sphere_vts_File << x << " " << y << " " << z  << endl;
            }
            Saturn_sphere_vts_File <<  "\n"  << endl;
        }
        Saturn_sphere_vts_File <<  "\n"  << endl;
    }
    Saturn_sphere_vts_File <<  "    </DataArray>\n"  << endl;
    Saturn_sphere_vts_File <<  "   </Points>\n"  << endl;
    Saturn_sphere_vts_File <<  "  </Piece>\n"  << endl;
    Saturn_sphere_vts_File <<  " </StructuredGrid>\n"  << endl;
    Saturn_sphere_vts_File <<  "</VTKFile>\n"  << endl;
    Saturn_sphere_vts_File.close();
    cout << "   File:  " << "[" << file_name << "]_Sat_sphere_" 
        << n << ".vts" << "  has been written to Directory:  " 
        << output_path << endl;
}
/*
 * 
*/
void cSaturnModel::SaturnPlotData(){
    // The whole body moved to the SHARED ParaViewWriter.h — it was 100 % identical
    // between the two models apart from the planet name in the file name. The
    // header-row/column mismatch it carries is described there.
    ParaViewWriter<cSaturnModel>(*this).plot_data();
}


