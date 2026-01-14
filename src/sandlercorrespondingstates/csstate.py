from __future__ import annotations
import numpy as np
from copy import deepcopy
from dataclasses import dataclass, field
from .charts import get_charts
from sandlermisc import GasConstant, DeltaH_IG, DeltaS_IG, StateReporter
from sandlerprops import Compound, get_database

import logging
logger = logging.getLogger(__name__)

@dataclass
class CSState:
    """
    Represents the thermodynamic state of a pure component using the corresponding states principle.
    """
    _cache: dict = field(default_factory=dict, init=False, repr=False)
    """ Cache for computed properties """
    _input_state: dict | None = field(default=None, init=False, repr=False)
    """ Snapshot of input field values for cache validation """

    _INPUT_FIELDS = frozenset(['T', 'P', 'Tc', 'Pc', 'x'])
    """ Fields that affect calculations and cache validity """

    _STATE_VAR_ORDERED_FIELDS = ['T', 'P', 'v', 's', 'h', 'u']

    _STATE_VAR_FIELDS = frozenset(_STATE_VAR_ORDERED_FIELDS)
    """ State variables for reporting """

    _UNIT_FIELDS = frozenset(['mass_unit', 'temperature_unit', 'pressure_unit', 'volume_unit'])
    """ Fields that define the units for caching purposes """

    def __setattr__(self, name, value):
        """ Clear cache on input field changes """
        if name in self._INPUT_FIELDS and hasattr(self, '_cache'):
            self._cache.clear()
        super().__setattr__(name, value)

    def _get_current_input_state(self):
        """ Snapshot of current input field values. """
        return {field: getattr(self, field) for field in self._INPUT_FIELDS}

    T: float  = None # Temperature in K
    P: float  = None # Pressure in MPa
    x: float  = None # Vapor fraction
    
    Liquid: CSState = None
    Vapor: CSState  = None

    Tc: float  = None # Critical temperature in K
    Pc: float  = None # Critical pressure in MPa
    
    z_read: float = None
    h_dep_read: float = None
    s_dep_read: float = None

    temperature_unit: str = "K"
    pressure_unit: str = "MPa"
    volume_unit: str = "m3"
    energy_unit: str = "J"

    Cp: float | list[float] | dict [str, float] | None = None
    """ heat capacity data for ideal-gas contributions """
    
    Tref: float = 298.15
    """ reference temperature for 'absolute' internal energy, enthalpy, and entropy calculations """
    Pref_MPa: float = 0.1
    """ reference pressure for 'absolute' entropy calculations, in MPa """

    round: int = 2

    def __post_init__(self):
        self._validate_inputs()
        self._input_state = self._get_current_input_state()

    def _validate_inputs(self):
        if self.T is None and self.P is None:
            raise ValueError("At least one of temperature (T) or pressure (P) must be provided.")
        if self.x is not None and self.T is not None and self.P is not None:
            raise ValueError("May not specify all three of T, P, and x.")
        if self.x is not None:
            if self.P is not None:
                if self.Pr > 1.0:
                    raise ValueError("Vapor fraction (x) can only be specified for subcritical pressures (Pr < 1).")
                else:
                    charts = get_charts()
                    Tsats = {}
                    for chart in ['Zchart', 'Hchart', 'Schart']:
                        Tsats[chart] = charts[chart].sat_indep_to_saturation_Tr(self.Pr) * self.Tc
                    self.T = np.mean(list(Tsats.values()))
                    logger.debug(f"Computed saturation temperature T={self.T} K for P={self.P} MPa and x={self.x} ({list(Tsats.values())})")
            if self.T is not None:
                if self.Tr > 1.0:
                    raise ValueError("Vapor fraction (x) can only be specified for subcritical temperatures (Tr < 1).")
                else:
                    charts = get_charts()
                    Psats = {}
                    for chart in ['Zchart', 'Hchart', 'Schart']:
                        Psats[chart] = charts[chart].sat_Tr_to_saturation_indep(self.Tr) * self.Pc
                    self.P = np.mean(list(Psats.values()))
                    logger.debug(f"Computed saturation pressure P={self.P} MPa for T={self.T} K and x={self.x} ({list(Psats.values())})")

    @property
    def R(self) -> GasConstant:
        return GasConstant("pa", "m3")

    @property
    def R_pv(self) -> GasConstant:
        return GasConstant(self.pressure_unit, self.volume_unit)

    @property
    def Tr(self) -> float:
        """Reduced temperature."""
        return self.T / self.Tc

    @property
    def Pr(self) -> float:
        """Reduced pressure."""
        return self.P / self.Pc

    @property
    def Z(self) -> np.ndarray:
        """Compressibility factor Z."""
        if 'Z' not in self._cache:
            charts = get_charts()
            self._cache['Z'] = charts['Zchart'].get_depvar(self.Pr, self.Tr, round=self.round)
        phase = self._cache['Z']['phase']
        logger.debug(f"Computed Z: {self._cache['Z']}")
        if phase == 'two_phase':
            if 'Liquid' not in self._cache:
                self.Liquid = self.clone(x=0.0)
                self._cache['Liquid'] = self.Liquid
            if 'Vapor' not in self._cache:
                self.Vapor = self.clone(x=1.0)
                self._cache['Vapor'] = self.Vapor
            return np.array([self._cache['Z']['depvar_vapor'], self._cache['Z']['depvar_liquid']])
        return np.array([self._cache['Z']['depvar']])
    
    @property
    def h_dep(self) -> np.ndarray:
        """Enthalpy departure in J/mol."""
        if 'h_dep_read' not in self._cache:
            charts = get_charts()
            self._cache['h_dep_read'] = charts['Hchart'].get_depvar(self.Pr, self.Tr, round=self.round)
        phase = self._cache['h_dep_read']['phase']
        logger.debug(f"Computed h_dep: {self._cache['h_dep_read']}")
        if phase == 'two_phase':
            if 'Liquid' not in self._cache:
                self.Liquid = self.clone(x=0.0)
                self._cache['Liquid'] = self.Liquid
            if 'Vapor' not in self._cache:
                self.Vapor = self.clone(x=1.0)
                self._cache['Vapor'] = self.Vapor
            return np.array([self._cache['h_dep_read']['depvar_vapor'] * -4.184 * self.Tc,
                             self._cache['h_dep_read']['depvar_liquid'] * -4.184 * self.Tc])
        return np.array([self._cache['h_dep_read']['depvar'] * -4.184 * self.Tc])

    @property
    def s_dep(self) -> np.ndarray:
        """Entropy departure in J/mol-K."""
        if 's_dep_read' not in self._cache:
            charts = get_charts()
            self._cache['s_dep_read'] = charts['Schart'].get_depvar(self.Pr, self.Tr, round=self.round)
        phase = self._cache['s_dep_read']['phase']
        logger.debug(f"Computed s_dep: {self._cache['s_dep_read']}")
        if phase == 'two_phase':
            if 'Liquid' not in self._cache:
                self.Liquid = self.clone(x=0.0)
                self._cache['Liquid'] = self.Liquid
            if 'Vapor' not in self._cache:
                self.Vapor = self.clone(x=1.0)
                self._cache['Vapor'] = self.Vapor
            return np.array([self._cache['s_dep_read']['depvar_vapor'] * -4.184,
                             self._cache['s_dep_read']['depvar_liquid'] * -4.184])
        return np.array([self._cache['s_dep_read']['depvar'] * -4.184])
    
    @property
    def v(self) -> np.ndarray:
        """Molar volume in m3/mol."""
        return self.Z * self.R_pv * self.T / self.P

    @property
    def Pv(self) -> np.ndarray:
        """Molar pressure-volume product in MPa-m3/mol."""
        return self.P * self.v

    @property
    def h(self) -> np.ndarray:
        """
        Computes absolute enthalpy at state T and P
        
        Returns
        -------
        np.ndarray
            Absolute enthalpy/enthalpies at state T and P
        """
        if self.Cp is None:
            raise ValueError("Cp data required for absolute enthalpy calculation.")
        dH_ideal = DeltaH_IG(self.Tref, self.T, self.Cp)
        logger.debug(f'dh_ideal: {dH_ideal} for Tref {self.Tref} to T {self.T} with Cp {self.Cp}')
        return self.h_dep + dH_ideal

    @property
    def u(self) -> np.ndarray:
        """
        Computes absolute internal energy at state T and P
        
        Returns
        -------
        np.ndarray
            Absolute internal energy/energies at state T and P
        """
        if self.Cp is None:
            raise ValueError("Cp data required for absolute internal energy calculation.")
        # u = h - pv
        return self.h - self.Pv * self.R / self.R_pv

    @property
    def Pref_local(self) -> float:
        """ 
        Returns reference pressure in local units

        Returns
        -------
        float
            Reference pressure in local units
        """
        if self.pressure_unit == 'mpa' or self.pressure_unit == 'MPa':
            return self.Pref_MPa
        elif self.pressure_unit == 'bar':
            return self.Pref_MPa * 10.0
        elif self.pressure_unit == 'pa' or self.pressure_unit == 'Pa':
            return self.Pref_MPa * 1.e6
        elif self.pressure_unit == 'kpa' or self.pressure_unit == 'kPa':
            return self.Pref_MPa * 1.e3
        elif self.pressure_unit == 'atm':
            return self.Pref_MPa * 9.86923
        else:
            raise ValueError(f"Unsupported pressure unit: {self.pressure_unit}")

    @property
    def s(self) -> np.ndarray:
        """
        Computes absolute entropy at state T and P
        
        Returns
        -------
        np.ndarray
            Absolute entropy/entropies at state T and P
        """
        if self.Cp is None:
            raise ValueError("Cp data required for absolute entropy calculation.")
        # make sure Pref_local is in correct units (same as self.P)
        dS_ideal = DeltaS_IG(self.Tref, self.Pref_local, self.T, self.P, self.Cp, self.R)
        return self.s_dep + dS_ideal

    def delta_h(self, other: CSState) -> np.ndarray:
        """
        Computes and returns enthalpy change from self to other state
        """
        self.unit_consistency(other)
        return other.h - self.h

    def delta_s(self, other: CSState) -> np.ndarray:
        """
        Computes and returns entropy change from self to other state
        """
        self.unit_consistency(other)
        return other.s - self.s
    
    def delta_pv(self, other: CSState) -> np.ndarray:
        """
        Returns Delta(PV) in thermal (not PV) units 
        """
        self.unit_consistency(other)
        return (other.Pv - self.Pv) * self.R / self.R_pv
    
    def delta_u(self, other: CSState) -> np.ndarray:
        """
        Returns Delta(U) (internal energy)
        """
        return self.delta_h(other) - self.delta_pv(other)

    def set_compound(self, compound: str | Compound):
        """
        Set critical properties and Cp data from a compound name.

        Parameters
        ----------
        compound: str | Compound
            Name of the compound to retrieve properties for
        """
        db = get_database()
        compound = db.get_compound(compound) if isinstance(compound, str) else compound
        compound_name = compound.Name if compound is not None else str(compound)
        if compound is None:
            raise ValueError(f"Compound '{compound_name}' not found in database.")
        return self.transfer_crits_from_compound(compound)

    def transfer_crits_from_compound(self, compound: Compound = None):
        """
        Set critical properties and Cp data from a Compound object.

        Parameters
        ----------
        compound: Compound
            Compound object containing critical properties and Cp data
        """
        if compound is not None:
            self.Tc = compound.Tc
            Pc_bar = compound.Pc
            if self.pressure_unit == 'MPa' or self.pressure_unit == 'mpa':
                self.Pc = Pc_bar / 10.0
            elif self.pressure_unit == 'bar':
                self.Pc = Pc_bar
            elif self.pressure_unit == 'kPa' or self.pressure_unit == 'kpa':
                self.Pc = Pc_bar * 100.0
            elif self.pressure_unit == 'Pa' or self.pressure_unit == 'pa':
                self.Pc = Pc_bar * 1.e5
            elif self.pressure_unit == 'atm':
                self.Pc = Pc_bar / 1.01325
            else:
                raise ValueError(f"Unsupported pressure unit: {self.pressure_unit}")
            self.Cp = deepcopy(compound.Cp)
        return self

    def report(self):
        reporter = StateReporter()
        for p in self._STATE_VAR_ORDERED_FIELDS + ['Pv']:
            if getattr(self, p) is not None:
                reporter.add_property(p, self._from_table_units(getattr(self, p), self.get_unit(p)), self.get_unit(p), self.get_formatter(p))
        if self.x is not None:
            reporter.add_property('x', self.x, f'{self.mass_unit} vapor/{self.mass_unit} total')
            for phase, state in [('L', self.Liquid), ('V', self.Vapor)]:
                for p in self._STATE_VAR_ORDERED_FIELDS + ['Pv']:
                    if not p in 'TP':
                        if getattr(state, p) is not None:
                            reporter.add_property(f'{p}{phase}', self._from_table_units(getattr(state, p), self.get_unit(p)), self.get_unit(p), self.get_formatter(p))
        return reporter.report()

    def __repr__(self):
        return f'CSState(T={self.T}, P={self.P}, v={self.v}, u={self.u}, h={self.h}, s={self.s}, x={self.x})'

    def clone(self, **kwargs) -> CSState:
        """ Create a copy of this State instance """
        new_state = CSState()
        for f in self._INPUT_FIELDS.union(self._UNIT_FIELDS):
            setattr(new_state, f, getattr(self, f))
        for k, v in kwargs.items():
            if k in self._INPUT_FIELDS.union(self._UNIT_FIELDS):
                setattr(new_state, k, v)
        return new_state