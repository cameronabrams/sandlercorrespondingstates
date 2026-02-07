from __future__ import annotations
import numpy as np
import json
from scipy.interpolate import interp1d, Rbf
import warnings

import logging

from abc import ABC, abstractmethod
from importlib.resources import files
from pathlib import Path

from sandlermisc import StateReporter, R, ureg

import pint

logger = logging.getLogger(__name__)

class CorrStsChart(ABC):
    """
    general abstract class for all corresponding states charts
    """
    data_path = files('sandlercorrespondingstates') / 'data'

    indepvar = 'Pr' # reduced pressure

    indicials_json = None
    @property
    @abstractmethod
    def isotherms_json(self):
        pass

    @property
    @abstractmethod
    def depvar(self) -> str:
        pass

    @property
    @abstractmethod
    def config(self):
        pass

    def __init__(self):
        """Initialize with digitized data in JSON format from WebPlotDigitizer (https://automeris.io/) """
        self.load_data(self.isotherms_json)
        if self.indicials_json:
            self.warp_data(self.indicials_json)
        self._identify_phase_transitions(min_dindep=self.config['min_dindep'], min_jump=self.config['min_depjump'], liquid_sense=self.config['liquid_sense'])
        self._create_interpolators()

    def load_data(self, isotherms_file: Path | str):

        with open(isotherms_file, 'r') as f:
            isotherms_data = json.load(f)
        
        self.isotherms = {}
        
        for dataset in isotherms_data['datasetColl']:
            name: str = dataset['name']
            
            # Only process isotherm datasets  
            if not name.startswith(f'{self.depvar}_T_'):
                continue
            
            tr_str = name.replace(f'{self.depvar}_T_', '').replace('p', '.')
            Tr = float(tr_str)
            
            pixel_x = []
            pixel_y = []
            Indep = []
            Dep = []

            for point in dataset['data']:
                value = point['value']
                Indep.append(value[0])
                Dep.append(value[1])
                pixel_x.append(point['x'])
                pixel_y.append(point['y'])
            
            Indep = np.array(Indep)
            Dep = np.array(Dep)
            pixel_x = np.array(pixel_x)
            pixel_y = np.array(pixel_y)

            sorted_indices = np.argsort(Indep)
            Indep = Indep[sorted_indices]
            Dep = Dep[sorted_indices]
            pixel_x = pixel_x[sorted_indices]
            pixel_y = pixel_y[sorted_indices]
            self.isotherms[Tr] = {
                self.indepvar: np.array(Indep),
                self.depvar: np.array(Dep),
                'x': pixel_x,
                'y': pixel_y,
                'n_points': len(Indep)
            }
            
        self.Tr_values = np.array(sorted(self.isotherms.keys()))
        
        logger.debug(f"\nLoaded {len(self.isotherms)} isotherms")
        logger.debug(f"Tr range: {self.Tr_values.min():.2f} to {self.Tr_values.max():.2f}")

    def warp_data(self, indicials_json: Path | str = None):
        if not indicials_json:
            return
        
        with open(indicials_json, 'r') as f:
            indicials_data = json.load(f)

        # Extract indicials
        pixel_x = []
        pixel_y = []
        data_indep = []
        data_dep = []
        
        for point in indicials_data['data']:
            pixel_x.append(point['x'])
            pixel_y.append(point['y'])
            data_indep.append(point['value'][0])
            data_dep.append(point['value'][1])
        
        pixel_x = np.array(pixel_x)
        pixel_y = np.array(pixel_y)
        data_indep = np.array(data_indep)
        data_dep = np.array(data_dep)
        
        logger.debug(f"Loaded {len(pixel_x)} indicials")
        
        if self.config['indepvar_scale'] == 'log':
            log_data_indep = np.log10(data_indep)
            self.warp_to_indep = Rbf(pixel_x, pixel_y, log_data_indep, 
                                  function='thin_plate', smooth=0)
        else:
            self.warp_to_indep = Rbf(pixel_x, pixel_y, data_indep, 
                                  function='thin_plate', smooth=0)
            
        self.warp_to_dep = Rbf(pixel_x, pixel_y, data_dep, 
                            function='thin_plate', smooth=0)
        
        
        for Tr, isotherm in self.isotherms.items():
            # Apply coordinate transformation
            if self.config['indepvar_scale'] == 'log':
                warped_log_indep = self.warp_to_indep(isotherm['x'], isotherm['y'])
                warped_indep = 10**warped_log_indep
            else:
                warped_indep = self.warp_to_indep(isotherm['x'], isotherm['y'])
            warped_dep = self.warp_to_dep(isotherm['x'], isotherm['y'])
            isotherm[self.indepvar] = warped_indep
            isotherm[self.depvar] = warped_dep

    def _identify_phase_transitions(self, min_jump = 10.0, min_dindep = 0.1, liquid_sense: str ='lower'):
        """
        Identify phase transitions in subcritical isotherms.
        
        The data may not be sorted by Indepvar - it follows the curve visually,
        tracing the liquid branch, then jumping to vapor branch.
        We need to identify the discontinuity and split accordingly.
        """
        self.phase_transitions = {}
        
        for Tr in self.Tr_values:
            if Tr >= 1.0:
                # Supercritical - no phase transition
                self.phase_transitions[Tr] = None
                continue
            
            indep = self.isotherms[Tr][self.indepvar]
            depvar = self.isotherms[Tr][self.depvar]
            
            search_dep_idx = np.where(indep < 1.0)
            search_dep = depvar[search_dep_idx]

            # Look for large jumps in Z between consecutive points
            # This indicates the discontinuity
            ddep = np.diff(search_dep)
            large_jumps = np.where((np.abs(ddep) > min_jump))[0]
            logger.debug(f'Tr {Tr}: large jumps: {large_jumps}: {ddep[large_jumps]}')
            if len(large_jumps) == 0:
                # No clear phase transition found
                logger.debug(f'No transitions detected by jumps in Hdep at Tr {Tr}')
                self.phase_transitions[Tr] = None
                continue
            elif len(large_jumps) == 3:
                # this is a boomerang -- the first two points should be swapped
                swap_idx = large_jumps[1:]
                logger.debug(f'   swapping indices {swap_idx}')
                i, j = swap_idx
                sav = indep[i]
                indep[i] = indep[j]
                indep[j] = sav
                self.isotherms[Tr][self.indepvar] = indep
                sav = depvar[i]
                depvar[i] = depvar[j]
                depvar[j] = sav
                self.isotherms[Tr][self.depvar] = depvar
                large_jumps = np.array([large_jumps[1]])
            # Find the largest jump
            jump_idx = large_jumps[np.argmax(np.abs(depvar[large_jumps]))]
            
            # The discontinuity is between jump_idx and jump_idx+1
            # Before jump: one phase
            # After jump: other phase
            
            dep_before = depvar[jump_idx]
            dep_after = depvar[jump_idx + 1]
            
            indep_before = indep[jump_idx]
            indep_after = indep[jump_idx + 1]
            indep_sat = 0.5 * (indep_before + indep_after)
            # regularize
            indep[jump_idx] = indep_sat
            indep[jump_idx+1] = indep_sat

            start_idx = jump_idx - 3
            if start_idx < 0:
                start_idx = 0
            for i in range(start_idx, jump_idx + 3):
                logger.debug(f'{i:>2d} {indep[i]:.3f} {depvar[i]:.3f}')
                if i == jump_idx:
                    logger.debug("-"*50)

            if np.abs(indep_after - indep_before) > min_dindep:
                logger.debug(f'Tr {Tr}: discretization-jump {dep_before}->{dep_after} indep_before {indep_before} indep_after {indep_after}; no transition')
                self.phase_transitions[Tr] = None
                continue

            if liquid_sense == 'lower':
                # lower value of depvar corresponds to liquid
                if dep_before < dep_after:
                    # liquid branch is before
                    liquid_indices = np.arange(0, jump_idx + 1)
                    vapor_indices = np.arange(jump_idx + 1, len(indep))
                    dep_L = dep_before
                    dep_V = dep_after
                else:
                    # liquid branch is after
                    vapor_indices = np.arange(0, jump_idx + 1)
                    liquid_indices = np.arange(jump_idx + 1, len(indep)) 
                    dep_L = dep_after
                    dep_V = dep_before
            elif liquid_sense == 'higher':
                if dep_before > dep_after:
                    # liquid branch is before
                    liquid_indices = np.arange(0, jump_idx + 1)
                    vapor_indices = np.arange(jump_idx + 1, len(indep))
                    dep_L = dep_before
                    dep_V = dep_after
                else:
                    # liquid branch is after
                    vapor_indices = np.arange(0, jump_idx + 1)
                    liquid_indices = np.arange(jump_idx + 1, len(indep)) 
                    dep_L = dep_after
                    dep_V = dep_before

            liq_indep = indep[liquid_indices]
            liq_depvar = depvar[liquid_indices]
            liq_sort = np.argsort(liq_indep)
            
            vap_indep = indep[vapor_indices]
            vap_depvar = depvar[vapor_indices]
            vap_sort = np.argsort(vap_indep)
            
            self.phase_transitions[Tr] = {
                f'{self.indepvar}_sat': indep_sat,
                f'{self.depvar}_L': dep_L,
                f'{self.depvar}_V': dep_V,
                'liquid': {
                    self.indepvar: liq_indep[liq_sort],
                    self.depvar : liq_depvar[liq_sort]
                },
                'vapor': {
                    self.indepvar: vap_indep[vap_sort],
                    self.depvar: vap_depvar[vap_sort]
                }
            }
            
            logger.debug(f"Tr={Tr:.2f}: Phase transition at {self.indepvar}_sat≈{indep_sat:.3f}")
            logger.debug(f"  Discontinuity between indices {jump_idx} ({indep[jump_idx]:.3f}) and {jump_idx+1} ({indep[jump_idx+1]:.3f})")
            logger.debug(f"  Liquid: {len(liquid_indices)} points, Hdep: [{liq_depvar.min():.3f}, {liq_depvar.max():.3f}], Pr: [{liq_indep.min():.3f}, {liq_indep.max():.3f}]")
            logger.debug(f"  Vapor:  {len(vapor_indices)} points, Hdep: [{vap_depvar.min():.3f}, {vap_depvar.max():.3f}], Pr: [{vap_indep.min():.3f}, {vap_indep.max():.3f}]")

        # make interpolators for the liquid and vapor branch of the two-phase envelope
        both_branch_indep = [c[f'{self.indepvar}_sat'] for c in self.phase_transitions.values() if c is not None]
        both_branch_Tr = [Tr for Tr in self.Tr_values if self.phase_transitions[Tr] is not None]
        vapor_branch_depvar = [c[f'{self.depvar}_V'] for c in self.phase_transitions.values() if c is not None]
        liquid_branch_depvar = [c[f'{self.depvar}_L'] for c in self.phase_transitions.values() if c is not None]

        both_branch_indep.append(1.0) # critical point
        both_branch_Tr.append(1.0)
        # get the critical depvar value from Tr=1.0 isotherm and indepvar value 1.0
        crit_depvar = self.isotherms[1.0][self.depvar][np.argmin(np.abs(self.isotherms[1.0][self.indepvar]-1.0))]
        vapor_branch_depvar.append(crit_depvar)
        liquid_branch_depvar.append(crit_depvar)
        logger.debug(f'Critical value of {self.depvar} at Tr=1.0 is {crit_depvar:.3f}')

        logger.debug(f"Creating saturation curve interpolators for {self.depvar}")
        for i, v, l in zip(both_branch_indep, vapor_branch_depvar, liquid_branch_depvar):
            logger.debug(f'  indep: {i:.3f} vapor {self.depvar}: {v:.3f} liquid {self.depvar}: {l:.3f}')
        for t, i in zip(both_branch_Tr, both_branch_indep):
            logger.debug(f'  Tr: {t:.3f} indep: {i:.3f}')
        self.sat_indep_to_vapor_depvar = interp1d(both_branch_indep, vapor_branch_depvar, kind='linear', fill_value='extrapolate')
        self.sat_indep_to_liquid_depvar = interp1d(both_branch_indep, liquid_branch_depvar, kind='linear', fill_value='extrapolate')
        self.sat_indep_to_saturation_Tr = interp1d(both_branch_indep, both_branch_Tr, kind='linear', fill_value='extrapolate')
        self.sat_Tr_to_saturation_indep = interp1d(both_branch_Tr, both_branch_indep, kind='linear', fill_value='extrapolate')

    def _create_interpolators(self):
        """Create interpolators for each isotherm, handling phase transitions."""
        self.interpolators = {}
        
        for Tr in self.Tr_values:
            if self.phase_transitions[Tr] is None:
                # Single-phase region - simple interpolator
                indep = self.isotherms[Tr][self.indepvar]
                depvar = self.isotherms[Tr][self.depvar]
                self.interpolators[Tr] = {
                    'type': 'single_phase',
                    'phase': None,
                    'interp': interp1d(indep, depvar, kind='linear',
                                      bounds_error=False, fill_value='extrapolate')
                }
            else:
                # Two-phase region - separate interpolators
                trans = self.phase_transitions[Tr]
                
                self.interpolators[Tr] = {
                    'type': 'two_phase',
                    f'{self.indepvar}_sat': trans[f'{self.indepvar}_sat'],
                    'liquid': interp1d(trans['liquid'][self.indepvar], trans['liquid'][self.depvar],
                                      kind='linear', bounds_error=False, fill_value='extrapolate'),
                    'vapor': interp1d(trans['vapor'][self.indepvar], trans['vapor'][self.depvar],
                                     kind='linear', bounds_error=False, fill_value='extrapolate')
                }
        
        logger.debug(f"{self.depvar}: Created interpolators for {len(self.interpolators)} isotherms")
        
        # Count two-phase isotherms
        n_two_phase = sum(1 for Tr in self.Tr_values 
                         if self.phase_transitions[Tr] is not None)
        logger.debug(f"  Single-phase: {len(self.interpolators) - n_two_phase}")
        logger.debug(f"  Two-phase: {n_two_phase}")
    
    def get_depvar(self, indep: float, Tr: float, round: int = None) -> dict:
        """
        Determine dependent variable value(s) from independent variable and reduced temperature.
        
        Parameters
        ----------
        indep : float
            independent variable
        Tr : float
            Reduced temperature
        round : int, optional
            Number of decimal places to round the result to.

        Returns
        -------
        dict:
            Dependent variable value(s) with keys 'subcooled', 'satd-liquid', 'satd-vapor', 'superheated', 'supercritical' as applicable.
        
        """
        
        idx = np.searchsorted(self.Tr_values, Tr)
            
        if idx == 0 or idx >= len(self.Tr_values):
            # Extrapolation - use nearest isotherm
            raise ValueError(f'Tr={Tr} is out of bounds ({self.Tr_values[0]} to {self.Tr_values[-1]}) for interpolation.')
            
        if indep > 1.0:
            depvar = self._supercrit(indep, Tr, round)
        else:
            depvar = self._subcrit(indep, Tr, round)

        return depvar

    def _supercrit(self, indep, Tr, round: int = 2):
        """Handle supercritical region (indep >= 1.0)."""

        idx = np.searchsorted(self.Tr_values, Tr)
        if self.Tr_values[idx-1] == Tr:
            # Exactly on an isotherm
            depvar_value = self.interpolators[Tr]['interp'](indep)            
        else:
            # Between isotherms - linear interpolation
            Tr_low = self.Tr_values[idx-1]
            Tr_high = self.Tr_values[idx]
            logger.debug(f"Supercritical interpolation at Pr = {indep}: Tr={Tr:.3f} between {Tr_low:.3f} and {Tr_high:.3f}")
            logger.debug(f'Trlow {self.interpolators[Tr_low]}')
            logger.debug(f'Trhigh {self.interpolators[Tr_high]}')

            interp = 'interp' if 'interp' in self.interpolators[Tr_low] else 'liquid'
            depvar_low = self.interpolators[Tr_low].get(interp)(indep)
            interp = 'interp' if 'interp' in self.interpolators[Tr_high] else 'liquid'
            depvar_high = self.interpolators[Tr_high].get(interp)(indep)

            alpha = (Tr - Tr_low) / (Tr_high - Tr_low)
            depvar_value = depvar_low + alpha * (depvar_high - depvar_low)
        
        return dict(state='supercritical', phase='supercritical', depvar=np.round(depvar_value, round))

    def _subcrit(self, indep, Tr, round):
        """ handle subcritical region (indep < 1.0)."""
        depvar_sat_L = self.sat_indep_to_liquid_depvar(indep)
        depvar_sat_V = self.sat_indep_to_vapor_depvar(indep)
        Tr_sat = self.sat_indep_to_saturation_Tr(indep)
        logger.debug(f"Subcritical: indep={indep:.3f}, Tr={Tr:.3f}, Tr_sat({self.depvar})={Tr_sat:.3f}, depvar_sat_L={depvar_sat_L:.3f}, depvar_sat_V={depvar_sat_V:.3f}")

        if Tr < Tr_sat:
            # Subcooled liquid
            depvar_value = self._subcrit_choose(indep, Tr, 'liquid', round)
            return dict(
                state='subcooled',
                phase='liquid',
                depvar=np.round(depvar_value, round)
            )
        else:
            # Superheated vapor
            depvar_value = self._subcrit_choose(indep, Tr, 'vapor', round)
            return dict(
                state='superheated',
                phase='vapor',
                depvar=np.round(depvar_value, round)
            )
        
    def _subcrit_choose(self, indep, Tr, phase, round):
        """ Interpolate in subcooled liquid region. """
        idx = np.searchsorted(self.Tr_values, Tr)
        
        if self.Tr_values[idx-1] == Tr:
            # Exactly on an isotherm
            interp = self.interpolators[Tr]
            if interp['type'] == 'single_phase':
                depvar_value = interp['interp'](indep)
            else:
                depvar_value = interp[phase](indep)
        else:
            # Between isotherms - linear interpolation
            Tr_low = self.Tr_values[idx-1]
            Tr_high = self.Tr_values[idx]
            
            interp_low = self.interpolators[Tr_low]
            if interp_low['type'] == 'single_phase':
                depvar_low = interp_low['interp'](indep)
            else:
                depvar_low = interp_low[phase](indep)
                
            interp_high = self.interpolators[Tr_high]
            if interp_high['type'] == 'single_phase':
                depvar_high = interp_high['interp'](indep)
            else:
                depvar_high = interp_high[phase](indep)
            
            alpha = (Tr - Tr_low) / (Tr_high - Tr_low)
            depvar_value = depvar_low + alpha * (depvar_high - depvar_low)
        
        return depvar_value
    
    def plot_chart(self, Tr_curves=None, figsize=(12, 8), show_phases=True, show_points=[]):
        """Plot the compressibility chart with phase transitions marked."""
        import matplotlib.pyplot as plt
        from matplotlib import colormaps as cm

        cmap = cm['viridis']
        config = self.config

        if Tr_curves is None:
            Tr_curves = self.Tr_values
        
        fig, ax = plt.subplots(figsize=figsize)
        
        n_iso = len(Tr_curves)

        for idx, Tr in enumerate(Tr_curves):
            if Tr not in self.isotherms:
                continue
            
            if self.phase_transitions[Tr] is None:
                # Single phase
                data = self.isotherms[Tr]
                shortcode = 'o-'
                if 'phase' in data:
                    if data['phase'] == 'liquid':
                        shortcode = 'o--'

                ax.plot(data[self.indepvar], data[self.depvar], shortcode, 
                       markersize=2, linewidth=1.5, alpha=0.7,
                       label=f'Tr = {Tr:.2f}', color=cmap(idx/n_iso))
            else:
                # Two phase
                trans = self.phase_transitions[Tr]
                
                # Plot liquid branch
                liq = trans['liquid']
                ax.plot(liq[self.indepvar], liq[self.depvar], 'o-',
                       markersize=2, linewidth=1.5, alpha=0.7, color=cmap(idx/n_iso),
                       label=f'Tr = {Tr:.2f}')
                
                # Plot vapor branch  
                vap = trans['vapor']
                ax.plot(vap[self.indepvar], vap[self.depvar], 'o--',
                       markersize=2, linewidth=1.5, alpha=0.7, color=cmap(idx/n_iso))
                
                # Mark transition
                if show_phases:
                    ax.plot([trans[f'{self.indepvar}_sat'], trans[f'{self.indepvar}_sat']], [trans[f'{self.depvar}_L'], trans[f'{self.depvar}_V']], color=cmap(idx/n_iso), 
                              linestyle='--', linewidth=1.5, alpha=0.5)
        
        for point in show_points:
            indep_val, dep_val = point
            ax.plot(indep_val, dep_val, 'rx', markersize=8, markeredgewidth=2, color='red')

        ax.set_xlabel('Reduced Pressure, Pr', fontsize=12)
        ax.set_ylabel(config['depvar_axis_label'], fontsize=12)
        ax.set_title(f'Corresponding States Chart for {config["depvar_name"]} (Zc = 0.27)\n',
                    fontsize=14, fontweight='bold')
        ax.set_xscale(self.config['indepvar_scale'])
        ax.grid(True, which='both', alpha=0.3)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, ncol=1)
        ax.set_xlim(*(config['indepvar_lims']))
        ax.set_ylim(*(config['depvar_lims']))
        
        plt.tight_layout()
        return fig, ax

class CompressibilityChart(CorrStsChart):
    
    @property
    def isotherms_json(self):
        return CorrStsChart.data_path / 'Z_vs_Pr_isotherms.json'
    
    indicials_json = CorrStsChart.data_path / 'Z_vs_Pr_indicials.json'

    @property
    def depvar(self):
        return 'Z'
    
    @property
    def config(self):
        return dict(
            liquid_sense = 'lower',
            min_dindep = 0.1,
            min_depjump = 0.3,
            depvar_name = 'Compressibility Factor',
            depvar_axis_label = 'Z',
            indepvar_lims = [0.1, 50.0],
            indepvar_scale = 'log',
            depvar_lims = [0, 1.5],
            depvar_units = ureg('dimensionless')
        )

class EnthalpyDepartureChart(CorrStsChart):

    @property
    def isotherms_json(self):
        return CorrStsChart.data_path / 'Hdep_vs_Pr_isotherms.json'
    
    @property
    def depvar(self):
        return 'Hdep'
    
    @property
    def config(self):
        return dict(
            liquid_sense = 'higher',
            min_dindep = 0.1,
            min_depjump = 3.0,
            depvar_name = 'Enthaply Departure',
            depvar_axis_label = r'$(h^{\rm IG} -h)/T_c\ \ \frac{\rm cal}{\rm mol K}$ (1 cal = 4.184 J)',
            indepvar_lims = [0.1, 50.0],
            indepvar_scale = 'log',
            depvar_lims = [0, 15.0],
            depvar_units = ureg('cal/(mol K)')
        )

class EntropyDepartureChart(CorrStsChart):

    @property
    def isotherms_json(self):
        return CorrStsChart.data_path / 'Sdep_vs_Pr_isotherms.json'
    
    @property
    def depvar(self):
        return 'Sdep'
    
    @property
    def config(self):
        return dict(
            liquid_sense = 'higher',
            min_dindep = 0.1,
            min_depjump = 3.0,
            depvar_name = 'Entropy Departure',
            depvar_axis_label = r'$s^{\rm IG} -s\ \ \frac{\rm cal}{\rm mol K}$ (1 cal = 4.184 J)',
            indepvar_lims = [0.1, 50.0],
            indepvar_scale = 'log',
            depvar_lims = [0, 21.5],
            depvar_units = ureg('cal/(mol K)')
        )

_instance = None
def get_charts():
    global _instance
    if _instance is None:
        _instance = {
            'Zchart': CompressibilityChart(),
            'Hchart': EnthalpyDepartureChart(),
            'Schart': EntropyDepartureChart()
        }
    return _instance

def demo():
    import matplotlib.pyplot as plt
    
    charts = get_charts()
    Zchart = charts['Zchart']
    Hchart = charts['Hchart']
    Schart = charts['Schart']
    
    # Validation
    logger.debug("\n" + "="*70)
    logger.debug("VALIDATION TESTS")
    logger.debug("="*70)
    
    tests = [
        (1.0, 1.0, 0.27, "auto", "Critical point"),
        (1.0, 1.5, 0.92, "auto", "Reference point"),
        (1.0, 0.90, None, "vapor", "Tr=0.90 vapor (above Pr_sat)"),
        (0.3, 0.90, None, "liquid", "Tr=0.90 liquid (below Pr_sat)"),
        (5.0, 2.0, None, "auto", "Supercritical"),
    ]
    
    logger.debug(f"\n{'Pr':>6} {'Tr':>6} {'Phase':>10} {'Z':>8} {'Expected':>10} {'Description':<30}")
    logger.debug("-"*85)
    
    for pr, tr, expected, phase, desc in tests:
        z = Zchart.get_depvar(pr, tr, phase=phase)
        if np.isnan(z):
            logger.debug(f"{pr:>6.1f} {tr:>6.2f} {phase:>10} {'NaN':>8} {'--':>10} {desc:<30}")
        elif expected:
            error = abs(z - expected) / expected * 100
            logger.debug(f"{pr:>6.1f} {tr:>6.2f} {phase:>10} {z:>8.4f} {expected:>10.2f} ({error:>4.1f}%) {desc:<30}")
        else:
            logger.debug(f"{pr:>6.1f} {tr:>6.2f} {phase:>10} {z:>8.4f} {'--':>10} {desc:<30}")
    
    # Plot
    logger.debug("\n" + "="*70)
    logger.debug("GENERATING PLOTS")
    logger.debug("="*70)
    
    fig, ax = Zchart.plot_chart(figsize=(8,8))
    plt.savefig('Zchart_with_phases.png', dpi=150, bbox_inches='tight')
    logger.debug("Saved: Zchart_with_phases.png")
    plt.close()
    
    fig, ax = Hchart.plot_chart(figsize=(8,12))
    plt.savefig('Hchart_with_phases.png', dpi=150, bbox_inches='tight')
    logger.debug("Saved: Hchart_with_phases.png")
    plt.close()

    fig, ax = Schart.plot_chart(figsize=(8,12))
    plt.savefig('Schart_with_phases.png', dpi=150, bbox_inches='tight')
    logger.debug("Saved: Schart_with_phases.png")
    plt.close()

    logger.debug("\n" + "="*70)
    logger.debug("DEMO COMPLETE")
    logger.debug("="*70)

if __name__ == "__main__":
    demo()
