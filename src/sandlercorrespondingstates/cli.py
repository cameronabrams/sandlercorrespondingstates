from .charts import *
from sandlerprops.properties import PropertiesDatabase
from sandlermisc.gas_constant import GasConstant
from sandlermisc.thermals import DeltaH_IG, DeltaS_IG
from sandlermisc.statereporter import StateReporter
import argparse as ap
import shutil
import logging
import os

banner = """
   _____                 ____                                                  
  / ___/____ _____  ____/ / /__  _____                                         
  \__ \/ __ `/ __ \/ __  / / _ \/ ___/                                         
 ___/ / /_/ / / / / /_/ / /  __/ /                                             
/____/\__,_/_/ /_/\__,_/_/\___/_/                               ___            
          _________  _____________  _________  ____  ____  ____/ (_)___  ____ _
         / ___/ __ \/ ___/ ___/ _ \/ ___/ __ \/ __ \/ __ \/ __  / / __ \/ __ `/
        / /__/ /_/ / /  / /  /  __(__  ) /_/ / /_/ / / / / /_/ / / / / / /_/ / 
        \___/\____/_/  /_/   \___/____/ .___/\____/_/ /_/\__,_/_/_/ /_/\__, /  
                  _____/ /_____ _/ /_/_/  _____                       /____/   
                 / ___/ __/ __ `/ __/ _ \/ ___/                                
                (__  ) /_/ /_/ / /_/  __(__  )                                 
               /____/\__/\__,_/\__/\___/____/                                  
                                        
(c) 2025, Cameron F. Abrams <cfa22@drexel.edu>
"""

logger = logging.getLogger(__name__)

def setup_logging(args):    
    loglevel_numeric = getattr(logging, args.logging_level.upper())
    if args.log:
        if os.path.exists(args.log):
            shutil.copyfile(args.log, args.log+'.bak')
        logging.basicConfig(filename=args.log,
                            filemode='w',
                            format='%(asctime)s %(name)s %(message)s',
                            level=loglevel_numeric
        )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(levelname)s> %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

def state(args):
    db = PropertiesDatabase()
    component = db.get_compound(args.n)
    if component is None:
        print(f"Component '{args.n}' not found in database.")
        return
    cs = CorrespondingStatesChartReader()
    Rpv = GasConstant("mpa", "m3")
    result = cs.dimensionalized_lookup(
        T = args.T,
        P = args.P,
        Tc = component.Tc,
        Pc = component.Pc/10,
        R_pv = Rpv
    )
    if result is not None:
        print(result.report())
    else:
        print("Could not find corresponding states properties for the given inputs.")

def delta(args):
    db = PropertiesDatabase()
    component = db.get_compound(args.n) # pressures are in bars!
    if component is None:
        print(f"Component '{args.n}' not found in database.")
        return
    Cp = [component.CpA, component.CpB, component.CpC, component.CpD]
    cs = CorrespondingStatesChartReader()
    Rpv = GasConstant("mpa", "m3")
    state1 = cs.dimensionalized_lookup(
        T = args.T1,
        P = args.P1,
        Tc = component.Tc,
        Pc = component.Pc/10,
        R_pv = Rpv
    )
    state2 = cs.dimensionalized_lookup(
        T = args.T2,
        P = args.P2,
        Tc = component.Tc,
        Pc = component.Pc/10,
        R_pv = Rpv
    )
    if state1 is not None and state2 is not None:
        delta_State = StateReporter({})
        prop_State = StateReporter({})
        prop_State.add_property('Tc', component.Tc, 'K', fstring="{:.2f}")
        prop_State.add_property('Pc', component.Pc/10, 'MPa', fstring="{:.2f}")
        prop_State.pack_Cp(Cp, fmts=["{:.2f}", "{:.3e}", "{:.3e}", "{:.3e}"])
        hdep2 = state2.get_value('Hdep')
        sdep2 = state2.get_value('Sdep')
        hdep1 = state1.get_value('Hdep')
        sdep1 = state1.get_value('Sdep')
        delta_h = hdep2 + DeltaH_IG(args.T1, args.T2, Cp) - hdep1
        delta_s = sdep2 + DeltaS_IG(args.T1, args.P1, args.T2, args.P2, Cp, GasConstant("pa", "m3")) - sdep1
        delta_pv = args.P2 * state2.get_value('v') - args.P1 * state1.get_value('v') # mpa m3/mol
        delta_u = delta_h - delta_pv*GasConstant("pa", "m3")/Rpv
        delta_State.add_property('Delta H', delta_h, 'J/mol', fstring="{:.2f}")
        delta_State.add_property('Delta S', delta_s, 'J/mol-K', fstring="{:.2f}")
        delta_State.add_property('Delta U', delta_u, 'J/mol', fstring="{:.2f}")
        if args.show_states:
            print("State 1:")
            print(state1.report())
            print("\nState 2:")
            print(state2.report())
            print("\nProperty differences:")
        print(delta_State.report())
        print("\nConstants used for calculations:")
        print(prop_State.report())
        
def cli():
    subcommands = {
        'state': dict(
            func = state,
            help = 'work with corresponding states for a single state'
        ),
        'delta': dict(
            func = delta,
            help = 'work with property differences between two states'
        ),
    }
    parser = ap.ArgumentParser(
        prog='sandlercorrespondingstates',
        description='Interact with corresponding states in Sandler\'s textbook'
    )
    parser.add_argument(
        '-b',
        '--banner',
        default=False,
        action=ap.BooleanOptionalAction,
        help='toggle banner message'
    )
    parser.add_argument(
        '--logging-level',
        type=str,
        default='debug',
        choices=[None, 'info', 'debug', 'warning'],
        help='Logging level for messages written to diagnostic log'
    )
    parser.add_argument(
        '-l',
        '--log',
        type=str,
        default='',
        help='File to which diagnostic log messages are written'
    )
    subparsers = parser.add_subparsers(
        title="subcommands",
        dest="command",
        metavar="<command>",
        required=True,
    )
    command_parsers={}
    for k, specs in subcommands.items():
        command_parsers[k] = subparsers.add_parser(
            k,
            help=specs['help'],
            formatter_class=ap.RawDescriptionHelpFormatter
        )
        command_parsers[k].set_defaults(func=specs['func'])

    state_args = [
        ('P', 'pressure', 'pressure in MPa', float, True),
        ('T', 'temperature', 'temperature in K', float, True),
        ('Pc', 'critical_pressure', 'critical pressure in MPa (if component not specified)', float, False),
        ('Tc', 'critical_temperature', 'critical temperature in K (if component not specified)', float, False),
        ('n', 'component', 'component name (e.g., methane, ethane, etc.)', str, False)
    ]
    for prop, long_arg, explanation, arg_type, required in state_args:
        command_parsers['state'].add_argument(
            f'-{prop}',
            f'--{long_arg}',
            dest=prop,
            type=arg_type,
            required=required,
            help=explanation
        )
    
    delta_args = [
        ('P1', 'pressure1', 'pressure of state 1 in MPa', float, True),
        ('T1', 'temperature1', 'temperature of state 1 in K', float, True),
        ('P2', 'pressure2', 'pressure of state 2 in MPa', float, True),
        ('T2', 'temperature2', 'temperature of state 2 in K', float, True),
        ('Pc', 'critical_pressure', 'critical pressure in MPa (if component not specified)', float, False),
        ('Tc', 'critical_temperature', 'critical temperature in K (if component not specified)', float, False),
        ('n', 'component', 'component name (e.g., methane, ethane, etc.)', str, False)
    ]
    for prop, long_arg, explanation, arg_type, required in delta_args:
        command_parsers['delta'].add_argument(
            f'-{prop}',
            f'--{long_arg}',
            dest=prop,
            type=arg_type,
            required=required,
            help=explanation
        )
    command_parsers['delta'].add_argument(
        '--show-states',
        default=False,
        action=ap.BooleanOptionalAction,
        help='also show the full states for state 1 and state 2'
    )
    args = parser.parse_args()
    setup_logging(args)
    if args.banner:
        print(banner)
    if hasattr(args, 'func'):
        args.func(args)
    else:
        my_list = ', '.join(list(subcommands.keys()))
        print(f'No subcommand found. Expected one of {my_list}')
    if args.banner:
        print('Thanks for using sandlercorrespondingstates!')