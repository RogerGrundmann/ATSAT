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

using namespace std;
//using namespace AtomUtils;

namespace ParaViewSaturn{
    void dump_array(const string &name, Array &a, double multiplier, ofstream &f) {
        f <<  "    <DataArray type=\"Float32\" Name=\"" << name << "\" format=\"ascii\">\n";
        for (int k = 0; k < a.km; k++){
            for (int j = 0; j < a.jm; j++){
                for (int i = 0; i < a.im; i++){
                    f << (a.x[i][j][k] * multiplier) << endl;
                }
                f << "\n";
            }
            f << "\n";
        }
        f << "\n";
        f << "    </DataArray>\n";
    }
/*
 * 
*/
    void dump_radial(const string &desc, Array &a, double multiplier, int i, ofstream &f){
        f << "SCALARS " << desc << " float " << 1 << endl;
        f << "LOOKUP_TABLE default" << endl;
        for (int j = 0; j < a.jm; j++){
            for (int k = 0; k < a.km; k++){
                f << (a.x[i][j][k] * multiplier) << endl;
            }
        }
    }
/*
 * 
*/
    void dump_radial_2d(const string &desc, Array_2D &a, double multiplier, ofstream &f){
        f << "SCALARS " << desc << " float " << 1 << endl;
        f << "LOOKUP_TABLE default" << endl;
        for (int j = 0; j < a.jm; j++){
            for (int k = 0; k < a.km; k++){
                f << (a.y[j][k] * multiplier) << endl;
            }
        }
    }
/*
 * 
*/
    void dump_zonal(const string &desc, Array &a, double multiplier, int k, ofstream &f){
        f <<  "SCALARS " << desc << " float " << 1 << endl;
        f <<  "LOOKUP_TABLE default" << endl;
        for(int i = 0; i < a.im; i++){
            for(int j = 0; j < a.jm; j++){
                f << (a.x[i][j][k] * multiplier) << endl;
            }
        }
    }
/*
 * 
*/
    void dump_longal(const string &desc, Array &a, double multiplier, int j, ofstream &f){
        f << "SCALARS " << desc << " float " << 1 << endl;
        f << "LOOKUP_TABLE default" << endl;
        for (int i = 0; i < a.im; i++){
            for (int k = 0; k < a.km; k++){
                f << (a.x[i][j][k] * multiplier) << endl;
            }
        }
    }
}

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
    DUMP("P_rain",        P_rain,     86400.0, IDX, F);   /* mm/day */        \
    DUMP("P_snow",        P_snow,     86400.0, IDX, F);                       \
    DUMP("P_graupel",     P_graupel,  86400.0, IDX, F);                       \
    DUMP("P_nh3_rain",    P_nh3_rain, 86400.0, IDX, F);                       \
    DUMP("P_nh3_snow",    P_nh3_snow, 86400.0, IDX, F);                       \
    DUMP("P_nh3_graupel", P_nh3_graupel, 86400.0, IDX, F);                    \
    DUMP("P_ch4_rain",    P_ch4_rain, 86400.0, IDX, F);                       \
    DUMP("P_ch4_snow",    P_ch4_snow, 86400.0, IDX, F);                       \
    DUMP("P_ch4_graupel", P_ch4_graupel, 86400.0, IDX, F);                    \
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
    dump_array("P_rain_mmd",        P_rain,        86400.0, F);               \
    dump_array("P_snow_mmd",        P_snow,        86400.0, F);               \
    dump_array("P_graupel_mmd",     P_graupel,     86400.0, F);               \
    dump_array("P_nh3_rain_mmd",    P_nh3_rain,    86400.0, F);               \
    dump_array("P_nh3_snow_mmd",    P_nh3_snow,    86400.0, F);               \
    dump_array("P_nh3_graupel_mmd", P_nh3_graupel, 86400.0, F);               \
    dump_array("P_ch4_rain_mmd",    P_ch4_rain,    86400.0, F);               \
    dump_array("P_ch4_snow_mmd",    P_ch4_snow,    86400.0, F);               \
    dump_array("P_ch4_graupel_mmd", P_ch4_graupel, 86400.0, F);               \
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
    "rho_mix Q_rad_mW_m3 Radiation P_rain_mmd P_snow_mmd P_graupel_mmd "      \
    "P_nh3_rain_mmd P_nh3_snow_mmd P_nh3_graupel_mmd "                        \
    "P_ch4_rain_mmd P_ch4_snow_mmd P_ch4_graupel_mmd P_nh4sh_mmd "            \
    "Q_precip_mW_m3 tke_m2s2 dis_nd nue_t_m2s prod_nd tke_source_nd "         \
    "dis_source_nd "
/*
 * 
*/
void cSaturnModel::paraview_panorama_vts(int n){
    using namespace ParaViewSaturn;
    double x, y, z, dx, dy, dz;
    string Saturn_panorama_vts_File_Name = output_path + "/Saturn_panorama_" 
        + std::to_string(n) + ".vts";
    string file_name = "Saturn_sphere_" 
        + std::to_string(n) + ".vts";
    ofstream Saturn_panorama_vts_File;
    Saturn_panorama_vts_File.precision(4);
    Saturn_panorama_vts_File.setf(ios::fixed);
    Saturn_panorama_vts_File.open(Saturn_panorama_vts_File_Name);
    if(!Saturn_panorama_vts_File.is_open()){
        cerr << "ERROR: could not open shpere_vts file " << __FILE__ 
            << " at line " << __LINE__ << "\n";
        abort();
    }
    Saturn_panorama_vts_File <<  "<?xml version=\"1.0\"?>\n"  << endl;
    Saturn_panorama_vts_File <<  "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"  << endl;
    Saturn_panorama_vts_File <<  " <StructuredGrid WholeExtent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Saturn_panorama_vts_File <<  "  <Piece Extent=\"" << 1 << " "<< im << " "<< 1 << " " << jm << " "<< 1 << " " << km << "\">\n"  << endl;
    Saturn_panorama_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"Temperature PressureDynamic PressureStatic CH4 CH4Cloud CH4Ice NH3 NH3Cloud NH3Ice H2O H2OCloud H2OIce Q_Latent Q_Sensible "
        PANORAMA_EXTRA_SCALARS "BuoyancyForce \">\n"  << endl;

    Saturn_panorama_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                Saturn_panorama_vts_File << u.x[i][j][k] << " " << v.x[i][j][k] << " " << w.x[i][j][k] << endl;
            }
            Saturn_panorama_vts_File <<  "\n"  << endl;
        }
        Saturn_panorama_vts_File <<  "\n"  << endl;
    }
    Saturn_panorama_vts_File <<  "\n"  << endl;
    Saturn_panorama_vts_File <<  "    </DataArray>\n" << endl;
    Saturn_panorama_vts_File <<  "    <DataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
//                Saturn_panorama_vts_File << t.x[i][j][k] * t_ref - t_ref << endl;
//                Saturn_panorama_vts_File << t.x[i][j][k] * t_ref << endl;
                Saturn_panorama_vts_File << t.x[i][j][k] * t_ref/10.0 << endl;
            }
            Saturn_panorama_vts_File <<  "\n"  << endl;
        }
        Saturn_panorama_vts_File <<  "\n"  << endl;
    }
    Saturn_panorama_vts_File <<  "\n"  << endl;
    Saturn_panorama_vts_File <<  "    </DataArray>\n" << endl;
//    dump_array("Seamount", SeaMount, 1.0, Saturn_panorama_vts_File);
    dump_array("u-component", u, u_0, Saturn_panorama_vts_File);
    dump_array("v-component", v, u_0, Saturn_panorama_vts_File);
    dump_array("w-component", w, u_0, Saturn_panorama_vts_File);

    dump_array("PressureDyn", p_dyn, 1.0, Saturn_panorama_vts_File);
//    dump_array("PressureStat", p_stat, 1.0, Saturn_panorama_vts_File);
//    dump_array("Density", rho, 1.0, Saturn_panorama_vts_File);

    dump_array("CoriolisForce", CoriolisForce, 1.0, Saturn_panorama_vts_File);
//    dump_array("CentrifugalForce", CentrifugalForce, 1e3, Saturn_panorama_vts_File);
    dump_array("BuoyancyForce", BuoyancyForce, 1.0, Saturn_panorama_vts_File);
//    dump_array("PresGradForce", PresGradForce, 1.0, Saturn_panorama_vts_File);

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

    dump_array("NH4SH", nh4sh, 1e3, Saturn_panorama_vts_File);
    dump_array("w_nh4sh", w_nh4sh, 1e3, Saturn_panorama_vts_File);
//    dump_array("j_nh4sh", j_nh4sh, 1e3, Saturn_panorama_vts_File);
//    dump_array("jT_nh4sh", jT_nh4sh, 1e3, Saturn_panorama_vts_File);

    dump_array("Q_Latent", Q_Latent, 1.0, Saturn_panorama_vts_File);
//    dump_array("Q_Sensible", Q_Sensible, 1.0, Saturn_panorama_vts_File);

    DUMP_EXTRA_FIELDS_VTS(Saturn_panorama_vts_File)

    Saturn_panorama_vts_File <<  "   </PointData>\n" << endl;
    Saturn_panorama_vts_File <<  "   <Points>\n"  << endl;
    Saturn_panorama_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n"  << endl;
    x = 0.0;
    y = 0.0;
    z = 0.0;
    dx = 0.1;
    dy = 0.1;
    dz = 0.1;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                if(k == 0 || j == 0) x = 0.0;
                else x = x + dx;
                Saturn_panorama_vts_File << x << " " << y << " " << z  << endl;
            }
            x = 0;
            y = y + dy;
            Saturn_panorama_vts_File <<  "\n"  << endl;
        }
        y = 0.0;
        z = z + dz;
        Saturn_panorama_vts_File <<  "\n"  << endl;
    }
    Saturn_panorama_vts_File <<  "    </DataArray>\n"  << endl;
    Saturn_panorama_vts_File <<  "   </Points>\n"  << endl;
    Saturn_panorama_vts_File <<  "  </Piece>\n"  << endl;
    Saturn_panorama_vts_File <<  " </StructuredGrid>\n"  << endl;
    Saturn_panorama_vts_File <<  "</VTKFile>\n"  << endl;
    Saturn_panorama_vts_File.close();
    cout << "   File:  " << "[" << file_name << "]_Sat_panorama_" 
        << n << ".vts" << "  has been written to Directory:  " 
        << output_path << endl;
    return;
}
/*
 * 
*/
void cSaturnModel::paraview_vtk_radial(int n, int i_radial){
    using namespace ParaViewSaturn;
    double x, y, z, dx, dy;
    string Saturn_radial_File_Name = output_path + "/Saturn_radial_" 
        + std::to_string(i_radial) + "_" + std::to_string(n) + ".vtk";
    string file_name = "Saturn_radial_" 
        + std::to_string(i_radial) + "_" + std::to_string(n) + ".vtk";
    ofstream Saturn_vtk_radial_File;
    Saturn_vtk_radial_File.precision (4);
    Saturn_vtk_radial_File.setf(ios::fixed);
    Saturn_vtk_radial_File.open(Saturn_radial_File_Name);
    if(!Saturn_vtk_radial_File.is_open()){
        cerr << "ERROR: could not open paraview_vtk file " << __FILE__ 
            << " at line " << __LINE__ << "\n";
        abort();
    }
    Saturn_vtk_radial_File <<  "# vtk DataFile Version 3.0" << endl;
    Saturn_vtk_radial_File <<  "Radial_Data_Sat_Circulation\n";
    Saturn_vtk_radial_File <<  "ASCII" << endl;
    Saturn_vtk_radial_File <<  "DATASET STRUCTURED_GRID" << endl;
    Saturn_vtk_radial_File <<  "DIMENSIONS " << km << " "<< jm << " " << 1 << endl;
    Saturn_vtk_radial_File <<  "POINTS " << jm * km << " float" << endl;
    x = 0.0;
    y = 0.0;
    z = 0.0;
    dx = 0.1;
    dy = 0.1;
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            if(k == 0) y = 0.0;
            else y = y + dy;
            Saturn_vtk_radial_File << x << " " << y << " "<< z << endl;
        }
        y = 0.0;
        x = x + dx;
    }
    Saturn_vtk_radial_File <<  "POINT_DATA " << jm * km << endl;
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

    dump_radial("thermalmassflux", thermalmassflux, 1e3, i_radial, Saturn_vtk_radial_File);
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

    dump_radial("NH4SH", nh4sh, 1e3, i_radial, Saturn_vtk_radial_File);
    dump_radial("w_nh4sh", w_nh4sh, 1e3, i_radial, Saturn_vtk_radial_File);
    dump_radial("j_nh4sh", j_nh4sh, 1e3, i_radial, Saturn_vtk_radial_File);
    dump_radial("jT_nh4sh", jT_nh4sh, 1e3, i_radial, Saturn_vtk_radial_File);
    dump_radial("massflux_nh4sh", massflux_nh4sh, 1e3, i_radial, Saturn_vtk_radial_File);
    dump_radial("difflux_nh4sh", difflux_nh4sh, 1e3, i_radial, Saturn_vtk_radial_File);


    dump_radial("PressureDyn", p_dyn, 1.0, i_radial, Saturn_vtk_radial_File);
    dump_radial("PressureStat", p_stat, 1.0, i_radial, Saturn_vtk_radial_File);
//    dump_radial("Density", rho, 1.0, i_radial, Saturn_vtk_radial_File);
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
    Saturn_vtk_radial_File.close();
    cout << "   File:  " << "[" <<  file_name << "]_Sat_radial_" 
        << i_radial << "_" << n << ".vtk" 
        << "  has been written to Directory:  " << output_path << endl;
    return;
}
/*
 * 
*/
void cSaturnModel::paraview_vtk_zonal(int n, int k_zonal){
    using namespace ParaViewSaturn;
    double x, y, z, dx, dy;
    string Saturn_zonal_File_Name = output_path + "/Saturn_zonal_" 
        + std::to_string(k_zonal) + "_" + std::to_string(n) + ".vtk";
    string file_name = "Saturn_zonal_" 
        + std::to_string(k_zonal) + "_" + std::to_string(n) + ".vtk";
    ofstream Saturn_vtk_zonal_File;
    Saturn_vtk_zonal_File.precision(4);
    Saturn_vtk_zonal_File.setf(ios::fixed);
    Saturn_vtk_zonal_File.open(Saturn_zonal_File_Name);
    if(!Saturn_vtk_zonal_File.is_open()){
        cerr << "ERROR: could not open vtk_zonal file " << __FILE__ 
            << " at line " << __LINE__ << "\n";
        abort();
    }
    Saturn_vtk_zonal_File <<  "# vtk DataFile Version 3.0" << endl;
    Saturn_vtk_zonal_File <<  "Zonal_Data_Sat_Circulation\n";
    Saturn_vtk_zonal_File <<  "ASCII" << endl;
    Saturn_vtk_zonal_File <<  "DATASET STRUCTURED_GRID" << endl;
    Saturn_vtk_zonal_File <<  "DIMENSIONS " << jm << " "<< im << " " << 1 << endl;
    Saturn_vtk_zonal_File <<  "POINTS " << im * jm << " float" << endl;
    x = 0.0;
    y = 0.0;
    z = 0.0;
    dx = 0.1;
    dy = 0.05;
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            if(j == 0) y = 0.0;
            else y = y + dy;
            Saturn_vtk_zonal_File << x << " " << y << " "<< z << endl;
        }
        y = 0.0;
        x = x + dx;
    }
    Saturn_vtk_zonal_File <<  "POINT_DATA " << im * jm << endl;
//    dump_zonal("Seamount", SeaMount, 1.0, k_zonal, Saturn_vtk_zonal_File);
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
    dump_zonal("thermalmassflux", thermalmassflux, 1e3, k_zonal, Saturn_vtk_zonal_File);

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

    dump_zonal("NH4SH", nh4sh, 1e3, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("w_nh4sh", w_nh4sh, 1e3, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("j_nh4sh", j_nh4sh, 1e3, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("jT_nh4sh", jT_nh4sh, 1e3, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("massflux_nh4sh", massflux_nh4sh, 1e3, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("difflux_nh4sh", difflux_nh4sh, 1e3, k_zonal, Saturn_vtk_zonal_File);

    dump_zonal("PressureDyn", p_dyn, 1.0, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("PressureStat", p_stat, 1.0, k_zonal, Saturn_vtk_zonal_File);
    dump_zonal("Density", rho, 1.0, k_zonal, Saturn_vtk_zonal_File);

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
    Saturn_vtk_zonal_File.close();
    cout << "   File:  " << "[" << file_name << "]_Sat_zonal_" 
        << k_zonal << "_" << n << ".vtk" 
        << "  has been written to Directory:  " << output_path << endl;
    return;
}
/*
 * 
*/
void cSaturnModel::paraview_vtk_longal(int n, int j_longal){
    using namespace ParaViewSaturn;
    double x, y, z, dx, dz;
    string Saturn_longal_File_Name = output_path + "/Saturn_longal_" 
        + std::to_string(j_longal) + "_" + std::to_string(n) + ".vtk";
    string file_name = "Saturn_longal_" 
        + std::to_string(j_longal) + "_" + std::to_string(n) + ".vtk";
    ofstream Saturn_vtk_longal_File;
    Saturn_vtk_longal_File.precision(4);
    Saturn_vtk_longal_File.setf(ios::fixed);
    Saturn_vtk_longal_File.open(Saturn_longal_File_Name);
    if(!Saturn_vtk_longal_File.is_open()){
        cerr << "ERROR: could not open vtk_longal file " 
            << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    Saturn_vtk_longal_File <<  "# vtk DataFile Version 3.0" << endl;
    Saturn_vtk_longal_File <<  "Longitudinal_Data_Sat_Circulation\n";
    Saturn_vtk_longal_File <<  "ASCII" << endl;
    Saturn_vtk_longal_File <<  "DATASET STRUCTURED_GRID" << endl;
    Saturn_vtk_longal_File <<  "DIMENSIONS " << km << " "<< im << " " << 1 << endl;
    Saturn_vtk_longal_File <<  "POINTS " << im * km << " float" << endl;
    x = 0.0;
    y = 0.0;
    z = 0.0;
    dx = 0.1;
    dz = 0.025;
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            if(k == 0){
                z = 0.0;
            }else{
                z = z + dz;
            }
            Saturn_vtk_longal_File << x << " " << y << " "<< z << endl;
        }
        z = 0.0;
        x = x + dx;
    }
    Saturn_vtk_longal_File <<  "POINT_DATA " << im * km << endl;
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

    dump_longal("thermalmassflux", thermalmassflux, 1e3, j_longal, Saturn_vtk_longal_File);

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

    dump_longal("NH4SH", nh4sh, 1e3, j_longal, Saturn_vtk_longal_File);
    dump_longal("w_nh4sh", w_nh4sh, 1e3, j_longal, Saturn_vtk_longal_File);
    dump_longal("j_nh4sh", j_nh4sh, 1e3, j_longal, Saturn_vtk_longal_File);
    dump_longal("jT_nh4sh", jT_nh4sh, 1e3, j_longal, Saturn_vtk_longal_File);
    dump_longal("massflux_nh4sh", massflux_nh4sh, 1e3, j_longal, Saturn_vtk_longal_File);
    dump_longal("difflux_nh4sh", difflux_nh4sh, 1e3, j_longal, Saturn_vtk_longal_File);

    dump_longal("PressureDyn", p_dyn, 1.0, j_longal, Saturn_vtk_longal_File);
    dump_longal("PressureStat", p_stat, 1.0, j_longal, Saturn_vtk_longal_File);
//    dump_longal("Density", rho, 1.0, j_longal, Saturn_vtk_longal_File);

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
    Saturn_vtk_longal_File.close();
    cout << "   File:  " << "[" << file_name << "]_Sat_longal_" 
        << j_longal << "_" << n << ".vtk" 
        << "  has been written to Directory:  " << output_path << endl;
    return;
}
/*
 * 
*/
void cSaturnModel::paraview_sphere_vts(int n){
    using namespace ParaViewSaturn;
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
    Saturn_sphere_vts_File <<  "   <PointData Vectors=\"Velocity\" Scalars=\"Temperature PressureDyn PressureStat  H2S H2SCloud H2SIce CH4 CH4Cloud CH4Ice NH3 NH3Cloud NH3Ice H2O H2OCloud H2OIce \">\n"  << endl;
    Saturn_sphere_vts_File <<  "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n"  << endl;
    for(int k = 0; k < km; k++){
        sinphi = sin( phi.z[k]);
        cosphi = cos( phi.z[k]);
        for(int j = 0; j < jm; j++){
            sinthe = sin(the.z[j]);
            costhe = cos(the.z[j]);
            for(int i = 0; i < im; i++){
                aux_u.x[i][j][k] = sinthe * cosphi * u.x[i][j][k] + costhe * cosphi * v.x[i][j][k] - sinphi * w.x[i][j][k];
                aux_v.x[i][j][k] = sinthe * sinphi * u.x[i][j][k] + sinphi * costhe * v.x[i][j][k] + cosphi * w.x[i][j][k];
                aux_w.x[i][j][k] = costhe * u.x[i][j][k] - sinthe * v.x[i][j][k];
                Saturn_sphere_vts_File << aux_u.x[i][j][k] << " " << aux_v.x[i][j][k] << " " << aux_w.x[i][j][k]  << endl;
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
    string Name_PlotData_File = output_path + "/PlotData_Saturn.xyz";
    ofstream PlotData_File;
    PlotData_File.precision(4);
    PlotData_File.setf(ios::fixed);
    PlotData_File.open(Name_PlotData_File);
    if(!PlotData_File.is_open()){
        cerr << "ERROR: could not open PlotData file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }
    PlotData_File << "lons(deg)" << ", " << "lats(deg)" << ", " 
        << "topography" << ", " << "v-velocity(m/s)" << ", " 
        << "w-velocity(m/s)" << ", " << "velocity-mag(m/s)" << ", " 
        << "temperature(Celsius)" << ", " << "water_vapour(g/kg)" 
        << ", " << "precipitation(mm)" << ", " 
        <<  "precipitable water(mm)" << endl;
    double vel_mag;
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            vel_mag = sqrt(pow(v.x[0][j][k] * u_0, 2) + pow(w.x[0][j][k] * u_0, 2));
            PlotData_File << k << " " << j << " " << SeaMount.x[0][j][k] << " " 
                << v.x[0][j][k] * u_0 << " " << w.x[0][j][k] * u_0 << " " 
                << vel_mag << " " << t.x[0][j][k] * t_ref - t_ref << " " 
                << h2o.x[0][j][k] << " "<< nh3.x[0][j][k] <<  endl;
        }
    }
    PlotData_File.close();
    return;
}


