from .csstate import CSState
from sandlermisc import GasConstant, StateReporter
import argparse as ap
import shutil
import logging
import os
from importlib.metadata import version

banner = r"""
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
                                        
(c) 2026, Cameron F. Abrams <cfa22@drexel.edu>
"""

logger = logging.getLogger(__name__)

def setup_logging(args):    
    loglevel_numeric = getattr(logging, args.log_level.upper())
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
    cs = CSState(T=args.T, P=args.P).set_compound(args.n)
    if cs is not None:
        print(cs.report())
    else:
        print("Could not find corresponding states properties for the given inputs.")

def delta(args):
    state1 = CSState(T=args.T1, P=args.P1).set_compound(args.n)
    state2 = CSState(T=args.T2, P=args.P2).set_compound(args.n)
    delta_h = state2.h - state1.h
    delta_s = state2.s - state1.s
    delta_u = state2.u - state1.u
    delta_State = StateReporter({})
    prop_State = StateReporter({})
    prop_State.add_property('Tc', state1.Tc, 'K', fstring="{:.2f}")
    prop_State.add_property('Pc', state1.Pc/10, 'MPa', fstring="{:.2f}")
    prop_State.pack_Cp(state1.Cp, fmts=["{:.2f}", "{:.3e}", "{:.3e}", "{:.3e}"])
    delta_State.add_property('Δh', delta_h, 'J/mol', fstring="{: 6g}")
    delta_State.add_property('Δs', delta_s, 'J/mol-K', fstring="{: 6g}")
    delta_State.add_property('Δu', delta_u, 'J/mol', fstring="{: 6g}")
    print(f"State-change calculations for {args.n} using corresponding states:")

    if args.show_states:
        print()
        two_states = ["State 1:                                      State 2:"]
        for line1, line2 in zip(state1.report().splitlines(), state2.report().splitlines()):
            two_states.append(f"{line1:<41s}     {line2}")
        print("\n".join(two_states))
        print()
        print("Property changes:")

    print(delta_State.report())
    print("\nConstants used for calculations:")
    print(prop_State.report())
        
def cli():
    subcommands = {
        'state': dict(
            func = state,
            help = 'Look up corresponding states properties for a single state'
        ),
        'delta': dict(
            func = delta,
            help = 'Compute corresponding states property differences between two states'
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
        '--log-level',
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
    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version=f'sandlercubics version {version("sandlercubics")}',
        help='show program version and exit'
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