# This will contain the baseclass for different types of cells
"""This module defines a Hydrogen Electrolyzer Cell."""

from typing import Union

import numpy as np
import scipy
import pandas as pd
import rainflow
from attrs import field, define
from scipy.signal import tf2ss, cont2discrete
from scipy.constants import R, physical_constants
from electrolyzer.tools.type_dec import NDArrayFloat, FromDictMixin, array_converter
from electrolyzer.tools.validators import contains


F, _, _ = physical_constants["Faraday constant"] 

@define
class Cell(FromDictMixin):
   
    # n: int = 2

    # Constants
    #TODO: change to z_c
    z: int = 2 # number of electrons transferred in reaction
    F: float = 96485.34  # Faraday's Constant (C/mol) or [As/mol]
    R: float = 8.314  # Ideal Gas Constant (J/mol/K)

    M_H: float = 1.00784  # molecular weight of Hydrogen [g/mol]
    M_H2: float = 2.016 #[g/mol]
    M_O: float = 15.999  # molecular weight of Oxygen [g/mol]
    M_O2: float = 31.999 #[g/mol]
    M_K: float = 39.0983  # molecular weight of Potassium [g/mol]

    lhv: float = 33.33  # lower heating value of H2 [kWh/kg]
    hhv: float = 39.41  # higher heating value of H2 [kWh/kg]
    gibbs: float = 237.24e3  # Gibbs Energy of global reaction (J/mol)

    def calc_current_density(self,current):
        j = current/self.cell_area
        return j

    def calc_cell_voltage(self,current,temperature):
        E_cell = self.calc_open_circuit_voltage()
        V_act_a, V_act_c = self.calc_activation_overpotential(current,temperature)
        V_ohm = self.calc_ohmic_overpotential(current,temperature)
        V_cell = E_cell + V_act_a + V_act_c + V_ohm
        return V_cell

    def calc_ohmic_overpotential(self,current,temperature):
        R_tot = self.calc_total_resistance(current,temperature)
        V_ohm = R_tot*current 
        return V_ohm

    def calc_activation_overpotential(self,current,temperature):
        pass

    def calc_open_circuit_voltage(self,temperature):
        E_rev_0 = self.calc_reversible_voltage()
        return E_rev_0

    def calc_reversible_voltage(self):
        return self.gibbs / (self.z * F)

    def calc_total_resistance(self,current,temperature):
        R_mem = self.calc_membrane_resistance()
        R_c,R_a = self.calc_electrode_resistance()
        R_elec = self.calc_electrolyte_resistance()

        R_tot = R_mem + R_c + R_a + R_elec
        return R_tot

    def calc_membrane_resistance(self):
        pass
    
    def calc_electrode_resistance(self):
        pass

    def calc_electrolyte_resistance(self):
        pass
    
    def calc_thermoneutral_voltage(self):
        pass

    def calc_temperature(self):
        pass

    def calc_reverse_faradays(self,hydrogen_demand_per_cell):
        # TODO: replace 2 with self.z
        I_reqd_BOL_noFaradaicLoss=(hydrogen_demand_per_cell*1000*2*F)/(self.dt*self.M_H2)
        n_f = self.calc_faradaic_efficiency(I_reqd_BOL_noFaradaicLoss)
        I_reqd = (hydrogen_demand_per_cell*1000*2*F)/(n_f*self.dt*self.M_H2)
        return I_reqd

    def calc_faradaic_efficiency(self,current):
        j = self.calc_current_density(current)
        j *= 1000
        eta_F = self.f_2 * (j**2) / (self.f_1 + j**2)
        return eta_F

    def calc_H2_mass_flow_rate(self,current):
        # TODO: replace 2 with self.z
        eta_F = self.calc_faradaic_efficiency(current)
        h2_prod_mol = eta_F * current / (2 * F) #[mol/sec]
        mfr = self.M_H2 * h2_prod_mol  # [g/sec]
        mfr_H2 = self.dt*mfr / 1e3  # [kg/cell-dt]
        return mfr_H2

    def calc_O2_mass_flow_rate(self,current):
        # TODO: replace 4 with ....
        eta_F = self.calc_faradaic_efficiency(current)
        o2_prod_mol = eta_F * current / (4 * F) #[mol/sec]
        mfr = self.M_O2 * o2_prod_mol  # [g/sec]
        mfr_O2 = self.dt*mfr / 1e3  # [kg/cell-dt]
        return mfr_O2

