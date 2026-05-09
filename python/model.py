#!/usr/bin/env python

from pyatsat import Saturn

class Model(object):
    """
    ATSAT Saturn's cell structure in the atmosphere
    """
    def __init__(self, cfg_xml):
        self.sat = Saturn()
        self.config_xml = cfg_xml

    def print_config_sat(self, config_xml):
        print("\n\n\n   Saturn's cell structure in the atmosphere")
        print("\n\n\n   Saturn atmosphere configuration file name = ", self.config_xml)
        print("   output path is           ", self.sat.output_path.decode('utf-8'))

    def run_Model_sat(self):
        print("\n   run_Model for the Saturn-Atmosphere code prepared")
        self.sat.run()
        print("\n    successfully terminated Saturn-Atmosphere code")
        print("\n")

sat = Model("config_atsat.xml")
sat.print_config_sat("config_atsat.xml")
sat.run_Model_sat()
