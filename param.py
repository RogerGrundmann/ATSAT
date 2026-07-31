# Given a parameter definition, generates necessary C++, Python and XML bindings
# coding=utf-8

def main():
    # read the input definition
    # name, description, datatype, default
    PARAMS = {
        'common': [
            ('bathymetry_path', '', 'string', '../data/Paleotopography_bathymetry/Golonka_rev210'),
            ('BathymetrySuffix', '', 'string', 'Ma_Golonka.xyz'),
            ('verbose', '', 'bool', False),
#            ('verbose', '', 'bool', True),
            ('output_path', 'directory where model outputs should be placed (must end in /)', 'string', 'output'),
#            ('paraview_flag','flag to control if create paraview panorama', 'bool', False),
            ('paraview_flag','flag to control if create paraview panorama', 'bool', True),
       ],








        'saturn': [
            ('nm', 'the maximum number of iterations', 'int', 224),
            ('checkpoint', "control when to write output files", 'int', 8),
            ('panorama_print', "control when to write panorama files", 'int', 32),
            ('checkpoint_save_iter', 'dump the full 3D prognostic state to output_path/sat_restart_<iter>.bin when iter_n reaches this; -1 disables', 'int', -1),
            ('restart_from_iter', 'load output_path/sat_restart_<iter>.bin and resume the iteration loop at iter+1, skipping the spin-up; -1 disables', 'int', -1),


            ('Coriolis', 'Coriolis force', 'double', 1),
            ('centrifugal', 'centrifugal force', 'double', 1),
            ('buoyancy', 'buoyancy force', 'double', 1),
            ('chemical_reaction', 'chemical reactions included', 'double', 1),
            ('epsres', 'accuracy of relative and absolute errors', 'double', 0.00001),

# Turbulence closure selection, as in ATJUP and ATOM. TurbulenceSat reads this; the module as a
# whole is still gated by ATSAT_TURB, and ATSAT_TURB_MODEL overrides this value at runtime.
# Uncomment one of the alternatives below to change the default.
            ('turb_model', 'turbulence model: none, k_epsilon, k_omega, k_omega_SST', 'string', 'k_omega_SST'),
#            ('turb_model', 'turbulence model: none, k_epsilon, k_omega, k_omega_SST', 'string', 'k_omega'),
#            ('turb_model', 'turbulence model: none, k_epsilon, k_omega, k_omega_SST', 'string', 'k_epsilon'),
#            ('turb_model', 'turbulence model: none, k_epsilon, k_omega, k_omega_SST', 'string', 'none'),

            ('L_atm', 'extension of the atmosphere shell in km, 120km/40 steps = 3.0km', 'double', 500.0),

            ('tropopause_pole', 'extension of the troposphere at the poles in km', 'double', 115.0),
            ('tropopause_equator', 'extension of the troposphere at the equator in km', 'double', 125.0),

            ('re', 'Reynolds number: ratio viscous to inertia forces, Re = u * L/nue', 'double', 1000.0),
            ('ec', 'Eckert number: ratio kinetic energy to enthalpy, Ec = u²/cp T', 'double', 0.00044),

            ('ep', 'ratio of the gas constants of dry air to water vapour [/]', 'double', 0.623),
            ('hp', 'water vapour pressure at T = 0°C: E = 6.1 hPa', 'double', 6.1078),
            ('lv', 'specific latent evaporation heat(condensation heat) in J/kg', 'double', 2.52e6),
            ('ls', 'specific latent vaporisation heat(sublimation heat) in J/kg', 'double', 2.83e6),
            ('cp_l', 'specific heat capacity of dry air at constant pressure and 20°C in J/(kg K)', 'double', 1005.0),

            ('pr', 'Prandtl number of h2oe for laminar flows', 'double', 0.69),
            ('g', 'gravitational acceleration of Saturn in m/s²', 'double', 10.0),
            ('omega', 'rotation number of Saturn in 1/s', 'double', 1.63e-4),

            ('gam', 'temperature lapse rate   gam = 0.7 K/km', 'double', 0.7),

            ('u_0', 'maximum value of velocity in m/s', 'double', 470.0),
            ('r_0_water', 'reference density of fresh water in kg/m3', 'double', 997.0),

            ('ua', 'initial velocity component in r-direction', 'double', 0.0),
            ('va', 'initial velocity component in theta-direction', 'double', 0.0),
            ('wa', 'initial velocity component in phi-direction', 'double', 0.0),
            ('pa', 'initial value for the pressure field', 'double', 0.0),
            ('ca', 'value 0.04 stands for the maximum value of 40 g/kg water vapour', 'double', 0.0),
            ('ta', 'initial value for the temperature field, 1.0 compares to 0° C compares to 273.15 K', 'double', 1.0),

            ('t_ref', 'temperature in K compare to 0°C', 'double', 134.0),

            ('t_equator', 'temperature at the equator 56.85°C compares to 330K', 'double', 400.0),
            ('t_pole', 'temperature at the poles 26.85°C compares to 300 K', 'double', 380.0),

            ('p_ref', 'pressure in bar', 'double', 1.0),
            ('R_ref', 'average gas constant in J/(g*K) after Sanchez et. al.', 'double', 3.89),

            ('t_tropopause', 'temperature in the tropopause -145°C compares to 110K', 'double', 0.4064),
            ('t_react_onset', 'temperature when reaction starts -33°C compares to 230K', 'double', 0.8420),

            ('p_tropopause', 'static pressure in the tropopause 0.14 bar', 'double', 0.14),

            ('h2_tropopause', 'minimum water vapour at tropopause h2_tropopause = 0.001 compares to 0.001 kg/kg', 'double', 0.0),
            ('he_tropopause', 'minimum water vapour at tropopause he_tropopause = 0.001 compares to 0.001 kg/kg', 'double', 0.0),
            ('h2o_tropopause', 'minimum water vapour at tropopause h2o_tropopause = 0.001 compares to 0.001 kg/kg', 'double', 0.0),
            ('h2s_tropopause', 'minimum water vapour at tropopause h2s_tropopause = 0.001 compares to 0.001 kg/kg', 'double', 0.0),
            ('nh3_tropopause', 'minimum rate nh3 at tropopause nh3_tropopause = 0.001 compares to 0.001 kg/kg', 'double', 0.0),
            ('nh4sh_tropopause', 'minimum rate nh4sh at tropopause nh4sh_tropopause = 0.001 compares to 0.001 kg/kg', 'double', 0.0),
            ('ch4_tropopause', 'minimum ch4 at tropopause kg/kg', 'double', 0.246),

       ],
 
    }

    XML_READ_FUNCS = {
        "string": "FillStringWithElement",
        "double": "FillDoubleWithElement",
        "int": "FillIntWithElement",
        "bool": "FillBoolWithElement"
    }

    def write_cpp_defaults(filename, classname, sections):
        with open(filename, 'w') as f:
            f.write("// header files\n")
            f.write("// THIS FILE IS AUTOMATICALLY GENERATED BY param.py\n")
            f.write("// ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME\n")
            f.write("\n")
            f.write("void %s::SetDefaultConfig() {\n" % classname)
            for section in sections:
                f.write('\n  // %s section\n' % section)
                for slug, desc, ctype, default in PARAMS[section]:
                    rhs = default
                    if ctype == 'string':
                        rhs = '"%s"' % default
                    elif ctype == 'bool':
                        if default:
                            rhs = 'true'
                        else:
                            rhs = 'false'
                    f.write('  %s = %s;\n' % (slug, rhs))
            f.write("}")

    def write_cpp_load_config(filename, classname, sections):
        with open(filename, 'w') as f:
            f.write("// config files\n")
            f.write("// THIS FILE IS AUTOMATICALLY GENERATED BY param.py\n")
            f.write("// ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME\n")
            f.write("\n")
            for section in sections:
                f.write('\n  // %s section\n' % section)
                element_var_name = 'elem_%s' % section
                f.write('\n  if (%s) {\n' % (element_var_name))
                for slug, desc, ctype, default in PARAMS [section]:
                    func_name = XML_READ_FUNCS [ctype]
                    f.write('    Config::%s(%s, "%s", %s);\n' % (func_name, element_var_name, slug, slug))
                f.write("  }\n")

    def write_cpp_ueaders(filename, sections, is_extern=False):
        with open(filename, 'w') as f:
            f.write("// header files\n")
            f.write("// THIS FILE IS AUTOMATICALLY GENERATED BY param.py\n")
            f.write("// ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME\n")
            f.write("\n")
            if is_extern:
                f.write("#include<string>\n\n")
                f.write("using namespace std;\n")
            for section in sections:
                f.write('\n// %s section\n' % section)
                for slug, desc, ctype, default in PARAMS [section]:
                    if is_extern:
                        f.write('   extern %s %s;\n' %(ctype, slug))
                    else:
                        f.write('%s %s;\n' %(ctype, slug))
           
            if is_extern:
                f.write("}\n")

    def write_pxi(input_filename, output_filename, substitutions):
        data = open(input_filename, 'r').read()
        indent = '    '
        for key, classname, sections in substitutions:
            rep = ''
            for section in sections:
                rep += '%s# %s section\n' % (indent, section)
                for slug, desc, ctype, default in PARAMS[section]:
                    rep += '%sproperty %s:\n' % (indent, slug)
                    rep += '%s    def __get__(%s self):\n' % (indent, classname)
                    rep += '%s        self._check_alive()\n' % indent
                    rep += '%s        return self._thisptr.%s\n' % (indent, slug)
                    rep += '%s\n' % indent
                    rep += '%s    def __set__(%s self, value):\n' % (indent, classname)
                    rep += '%s        self._check_alive()\n' % indent
                    rep += '%s        self._thisptr.%s = <%s> value\n' % (indent, slug, ctype)
                    rep += '%s\n' % indent
            data = data.replace('{{ %s }}' % key, rep)
        with open(output_filename, 'w') as f:
            f.write("""# pxi files\n""")
            f.write("# THIS FILE IS AUTOMATICALLY GENERATED BY param.py\n")
            f.write("# ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME\n")
            f.write(data)

    def write_pxd(filename, model, sections):
        with open(filename, 'w') as f:
            # Sadly, Cython docs are incorrect on usage of 'include', so we must include a whole lot of boilerplate
            f.write("""# pxd files\n""")
            f.write("""# THIS FILE IS AUTOMATICALLY GENERATED BY param.py
# ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME
from libcpp.vector cimport vector
cdef extern from "c%sModel.h":
    cppclass c%sModel:
        c%sModel() except +  # NB! std::bad_alloc will be converted to MemoryError
        void LoadConfig(const char *filename)
        void Run()
""" % (model, model, model))
            for section in sections:
                f.write('        # %s section\n' % section)
                for slug, desc, ctype, default in PARAMS [section]:
                    f.write('        %s %s\n' % (ctype, slug))

    def write_config_xml(filename, sections):
        with open(filename, 'w') as f:
            f.write("""<!-- THIS FILE IS GENERATED AUTOMATICALLY BY param.py. DO NOT EDIT. -->""")
            f.write('<atsat>')
            for section in sections:
                f.write('    <%s>\n' % section)
                for slug, desc, ctype, default in PARAMS [section]:
                    if ctype == 'bool':
                        default = str(default).lower()  # Python uses True/False, C++, uses true/false
                    f.write('        <%s>%s</%s>  <!-- %s (%s) -->\n' % (slug, default, slug, desc, ctype))
                f.write('    </%s>\n' % section)
            f.write('</atsat>')

    saturn_atmosphere_sections = ['common', 'saturn']

    for filename, classname, sections in [
       ('planet/cSaturnDefaults.cpp.inc', 'cSaturnModel', saturn_atmosphere_sections),
   ]:
        write_cpp_defaults(filename, classname, sections)

    for filename, classname, sections in [
        ('planet/SaturnLoadConfig.cpp.inc', 'cSaturnModel', saturn_atmosphere_sections),
   ]:
        write_cpp_load_config(filename, classname, sections)

    for filename, sections in [
        ('planet/SaturnParams.h.inc', saturn_atmosphere_sections),
   ]:
        write_cpp_ueaders(filename, sections)

    write_pxi ('python/pyatsat.pyx.template', 'python/pyatsat.pyx', [
        ('saturn_params', 'Saturn', saturn_atmosphere_sections)]
    )

    for filename, model, sections in [
        ('python/saturn_pxd.pxi', 'Saturn', saturn_atmosphere_sections),
   ]:
        write_pxd(filename, model, sections)

    for  filename, sections in [
        ('python/config_atsat.xml', saturn_atmosphere_sections),
   ]:
        write_config_xml(filename, sections)

    for  filename, sections in [
        ('cli/config_atsat.xml', saturn_atmosphere_sections),
   ]:
        write_config_xml(filename, sections)

    for  filename, sections in [
        ('saturn/config_atsat.xml', saturn_atmosphere_sections),
   ]:
        write_config_xml(filename, sections)

if __name__ == '__main__':
    main()
