from __future__ import annotations
import numpy as np
from copy import deepcopy
from dataclasses import dataclass, field
from .charts import get_charts
from sandlermisc import R, DeltaH_IG, DeltaS_IG, ureg, cached_property, ThermodynamicState
from sandlerprops import Compound, get_database
import pint

import logging
logger = logging.getLogger(__name__)

@dataclass
class CSState(ThermodynamicState):
    """
    Represents the thermodynamic state of a pure component using the corresponding states principle.
    """
    phase: str = 'unspecified' 
    """ strict phase requirement flag; 'vapor', 'liquid', 'unspecified' (multiple roots) """

    name: str = 'cs'

    description: str = 'Corresponding states equation of state'

    Pc: pint.Quantity = None
    """ critical pressure """
    Tc: pint.Quantity = None
    """ critical temperature """

    Cp: float | list[float] | dict [str, float] | None = None
    """ heat capacity data for ideal-gas contributions """

    _PARAMETER_ORDERED_FIELDS = ['Tc', 'Pc', 'Cp']
    _PARAMETER_FIELDS = frozenset(_PARAMETER_ORDERED_FIELDS)
    """ Fields that define parameters for the EOS; to be defined in subclasses """
        
    Tref: pint.Quantity = 298.15 * ureg.kelvin
    """ reference temperature for 'absolute' internal energy, enthalpy, and entropy calculations """
    Pref: pint.Quantity = 0.1 * ureg.megapascal
    """ reference pressure for 'absolute' entropy calculations, in MPa """

    round: int = 2

    def get_default_unit(self, field_name: str) -> pint.Unit:
        """
        Get the default unit for a given field
        """
        _default_unit_map = {
            'P': ureg.megapascal,
            'T': ureg.kelvin,
            'v': ureg.meter**3 / ureg.mol,
            'u': ureg.joule / ureg.mol,
            'h': ureg.joule / ureg.mol,
            's': ureg.joule / (ureg.mol * ureg.kelvin),
        }
        return _default_unit_map.get(field_name, ureg.dimensionless)

    def resolve(self) -> bool:
        """
        Resolve all state variables from the two input variables.   This method
        is called automatically when any property assignment results
        in a fully specified state.  There should never be a situation in which
        this method is called explicitly; it is only called by __setattr__ or __post_init__.
        """
        hasT, hasP, hasx = self._is_specified_input_var('T'), self._is_specified_input_var('P'), self._is_specified_input_var('x')
        if hasT and hasP:
            assert not hasx, f'Cannot specify T, P, and x simultaneously.'
            resolved = self._resolve_TP()
        elif hasT:
            if hasx: # explicitly saturated state
                resolved = self._resolve_saturated_Tx()
            else:
                raise ValueError(f'Cannot resolve state with T specified but not P or x.')
        elif hasP:
            if hasx: # explicitly saturated state
                resolved = self._resolve_saturated_Px()
            else:
                raise ValueError(f'Cannot resolve state with P specified but not T or x.')
        if not resolved:
            raise ValueError(f'Unable to resolve state with inputs: {self.get_input_varnames}')
        return resolved

    def _resolve_TP(self) -> bool:
        """Resolve state variables from T and P inputs."""
        self.Z = self._solve_for_Z()
        return self._calculate_vhus()

    def _resolve_saturated_Tx(self) -> bool:
        if self.T > self.Tc:
            raise ValueError(f'Cannot resolve saturated state with T={self.T} above critical temperature Tc={self.Tc}.')
        self.P = self.Pvap
        if 0 < self.x < 1:
            self.Liquid = self.__class__(x=0.0, T=self.T, Pc=self.Pc, Tc=self.Tc, Cp=self.Cp)
            self.Vapor = self.__class__(x=1.0, T=self.T, Pc=self.Pc, Tc=self.Tc, Cp=self.Cp)
            for op in self._STATE_VAR_FIELDS - {'T', 'P', 'x'}:
                setattr(self, op, self.x * getattr(self.Vapor, op) + (1 - self.x) * getattr(self.Liquid, op))
            return True
        else:
            return self._resolve_TP()
    
    def _resolve_saturated_Px(self) -> bool:
        if self.P > self.Pc:
            raise ValueError(f'Cannot resolve saturated state with P={self.P} above critical pressure Pc={self.Pc}.')
        self.T = self.Tsat
        if 0 < self.x < 1:
            self.Liquid = self.__class__(x=0.0, P=self.P, Pc=self.Pc, Tc=self.Tc, Cp=self.Cp)
            self.Vapor = self.__class__(x=1.0, P=self.P, Pc=self.Pc, Tc=self.Tc, Cp=self.Cp)
            for op in self._STATE_VAR_FIELDS - {'T', 'P', 'x'}:
                setattr(self, op, self.x * getattr(self.Vapor, op) + (1 - self.x) * getattr(self.Liquid, op))
            return True
        else:
            return self._resolve_TP()

    @cached_property
    def Tr(self) -> float:
        """Reduced temperature"""
        return self._calc_Tr()

    def _calc_Tr(self) -> float:
        return self.T / self.Tc

    @cached_property
    def Pr(self) -> float:
        """Reduced pressure"""
        return self._calc_Pr()

    def _calc_Pr(self) -> float:
        return self.P / self.Pc

    @cached_property
    def Pvap(self) -> pint.Quantity:
        return self._calc_Pvap()

    def _calc_Pvap(self) -> pint.Quantity:
        if self.T >= self.Tc:
            raise ValueError(f'Cannot calculate Pvap at T={self.T} above critical temperature Tc={self.Tc}.')
        charts = get_charts()
        Psats = {}
        for chart in ['Zchart', 'Hchart', 'Schart']:
            Psats[chart] = charts[chart].sat_Tr_to_saturation_indep(self.Tr)
        Pvap = np.mean(list(Psats.values())) * ureg.megapascal
        logger.debug(f"Computed saturation pressure Pvap={Pvap} MPa for T={self.T} K and Tr={self.Tr} ({list(Psats.values())})")
        return Pvap

    @cached_property
    def Tsat(self) -> pint.Quantity:
        return self._calc_Tsat()    

    def _calc_Tsat(self) -> pint.Quantity:
        if self.P >= self.Pc:
            raise ValueError(f'Cannot calculate Tsat at P={self.P} above critical pressure Pc={self.Pc}.')
        charts = get_charts()
        Trsats = {}
        for chart in ['Zchart', 'Hchart', 'Schart']:
            Trsats[chart] = charts[chart].sat_indep_to_saturation_Tr(self.Pr)
        Trsat = np.mean(list(Trsats.values()))
        Tsat = Trsat * self.Tc
        logger.debug(f"Computed saturation temperature Tsat={Tsat} K for P={self.P} MPa and Pr={self.Pr} ({list(Trsats.values())})")
        return Tsat

    def _spawn_helper(self) -> CSState:
        """
        Create a helper CSState instance for internal calculations
        """
        # logger.debug(f'_spawn_helper: Creating helper state for {self.name}')
        # logger.debug(f'_spawn_helper: Creating helper state for {self.name} with T={self.T:.3f}, P={self.P:.3f}, Tc={self.Tc:.3f}, Pc={self.Pc:.3f}, omega={self.omega:.3f}, Cp={self.Cp}')
        helper_state = self.__class__.simple(
            name=f'{self.name}_helper' if self.name else 'CSState_helper',
            T=self.T,
            P=self.P,
            Tc=self.Tc,
            Pc=self.Pc,
            Cp=self.Cp,
        )
        # logger.debug(f'_spawn_helper: Created helper state {helper_state.name} with T={helper_state.T:.3f}, P={helper_state.P:.3f}, Tc={helper_state.Tc:.3f}, Pc={helper_state.Pc:.3f}, omega={helper_state.omega:.3f}, Cp={helper_state.Cp}')
        return helper_state

    def _calc_Hvap(self) -> pint.Quantity:
        charts = get_charts()
        h_chart = charts['Hchart']
        liquid_hdep = h_chart.sat_indep_to_liquid_depvar(self.Tr)
        vapor_hdep = h_chart.sat_indep_to_vapor_depvar(self.Tr)
        logger.debug(f"_calc_Hvap: liquid_hdep={liquid_hdep}, vapor_hdep={vapor_hdep} for T={self.T} K and Tr={self.Tr}")
        Hvap = -(vapor_hdep - liquid_hdep) * self.Tc * ureg.cal / ureg.mol / ureg.kelvin
        Hvap = Hvap.to(ureg.joule / ureg.mol)
        logger.debug(f"Computed Hvap={Hvap} J/mol for T={self.T} K and Tr={self.Tr} (liquid h_dep={liquid_hdep}, vapor h_dep={vapor_hdep})")
        return Hvap

    @cached_property
    def Hvap(self) -> pint.Quantity:
        return self._calc_Hvap()
   
    def _calc_Svap(self) -> pint.Quantity:
        charts = get_charts()
        s_chart = charts['Schart']
        liquid_sdep = s_chart.sat_indep_to_liquid_depvar(self.Tr)
        vapor_sdep = s_chart.sat_indep_to_vapor_depvar(self.Tr)
        logger.debug(f"_calc_Svap: liquid_sdep={liquid_sdep}, vapor_sdep={vapor_sdep} for T={self.T} K and Tr={self.Tr}")
        Svap = -(vapor_sdep - liquid_sdep) * ureg.cal / ureg.mol / ureg.kelvin
        Svap = Svap.to(ureg.joule / (ureg.mol * ureg.kelvin))
        logger.debug(f"Computed Svap={Svap} J/mol-K for T={self.T} K and Tr={self.Tr} (liquid s_dep={liquid_sdep}, vapor s_dep={vapor_sdep})")
        return Svap
    
    @cached_property
    def Svap(self) -> pint.Quantity:
        return self._calc_Svap()

    def _solve_for_Z(self) -> float:
        """Compressibility factor Z."""
        charts = get_charts()
        z_var = charts['Zchart'].get_depvar(self.Pr, self.Tr, round=self.round)
        return z_var['depvar']

    @cached_property
    def h_departure(self) -> pint.Quantity:
        return self._calc_h_departure()

    def _calc_h_departure(self) -> pint.Quantity:
        """Enthalpy departure"""
        charts = get_charts()
        h_vardict = charts['Hchart'].get_depvar(self.Pr, self.Tr, round=self.round)
        self.neg_h_dep_over_Tc_read = h_vardict['depvar'] * ureg.cal / ureg.mol / ureg.kelvin
        h_read = -1 * h_vardict['depvar'] * ureg.cal / ureg.mol / ureg.kelvin * self.Tc
        logger.debug(f"Computed h_dep: {h_read} for Pr={self.Pr}, Tr={self.Tr}")
        return h_read.to(ureg.joule / ureg.mol)

    @cached_property
    def s_departure(self) -> pint.Quantity:
        return self._calc_s_departure()

    def _calc_s_departure(self) -> pint.Quantity:
        """Entropy departure"""
        charts = get_charts()
        s_vardict = charts['Schart'].get_depvar(self.Pr, self.Tr, round=self.round)
        self.neg_s_dep_read = s_vardict['depvar'] * ureg.cal / ureg.mol / ureg.kelvin
        s_read = -1 * s_vardict['depvar'] * ureg.cal / ureg.mol / ureg.kelvin
        logger.debug(f"Computed s_dep: {s_read} for Pr={self.Pr}, Tr={self.Tr}")
        return s_read.to(ureg.joule / (ureg.mol * ureg.kelvin))

    def _calculate_op(self, op: str) -> bool:
        """
        Calculate other properties (v, h, u, s) from Z, T, P
        """
        vunit = self.get_default_unit('v')
        punit = self.get_default_unit('P')
        tunit = self.get_default_unit('T')
        eunit = self.get_default_unit('h')
        sunit = self.get_default_unit('s')
        if 'v' in op:
            v = self.Z * R.to(punit * vunit / tunit) * self.T / self.P
            setattr(self, 'v', v)
        if 'h' in op or 'u' in op or 's' in op:
            if self.Cp is None:
                raise ValueError("Cp data required for absolute enthalpy/entropy calculation.")
        if 'h' in op:
            dH_ideal = DeltaH_IG(self.Tref.to('K'), self.T.to('K'), self.Cp).to(eunit)
            h = self.h_departure + dH_ideal
            setattr(self, 'h', h)
            if 'u' in op:
                Pv = (self.P * self.v).to(eunit)
                u = h - Pv
                setattr(self, 'u', u)
                setattr(self, 'Pv', Pv)
        if 'u' in op and not 'h' in op:
            dH_ideal = DeltaH_IG(self.Tref.to('K'), self.T.to('K'), self.Cp).to(eunit)
            Pv = (self.P * self.v).to(eunit)
            u = h - Pv
            setattr(self, 'u', u)
            setattr(self, 'Pv', Pv)
        if 's' in op:
            dS_ideal = DeltaS_IG(self.Tref.to('K'), self.Pref.to(punit), 
                                self.T.to('K'), self.P.to(punit), 
                                self.Cp, R.to(eunit/tunit)).to(sunit)
            s = self.s_departure + dS_ideal
            setattr(self, 's', s)
        return True

    def _calculate_hus(self) -> bool:
        return self._calculate_op('hus')

    def _calculate_vhus(self) -> bool: # T, P
        return self._calculate_op('vhus')

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
            self.Pc = compound.Pc.to('MPa')
            self.omega = compound.Omega
            self.Cp = deepcopy(compound.Cp)
        return self

