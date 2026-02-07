from .csstate import CSState
from sandlermisc import R, StateReporter, ureg
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

def state_subcommand(args):
    """
    Calculate and report the state for a single condition using the corresponding states method.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    """
    csstate = CSState()
    if args.n is not None:
        csstate.set_compound(args.n)
    if args.Tc is not None:
        csstate.Tc = args.Tc
    if args.Pc is not None:
        csstate.Pc = args.Pc
    if args.Cp is not None:
        csstate.Cp = args.Cp
    for p in 'TPx':
        v = getattr(args, p, None)
        if v is not None:
            setattr(csstate, p, v)
    additional_vars = ['Tc', 'Pc', 'Tr', 'Pr', 'Z', 'h_departure', 's_departure', 'neg_h_dep_over_Tc_read', 'neg_s_dep_read']
    if csstate.T < csstate.Tc:
        additional_vars.extend(['Pvap', 'Hvap', 'Svap'])
    if csstate.P < csstate.Pc:
        additional_vars.extend(['Tsat'])
    property_notes = {
        'Pvap': f'at {csstate.T.to(ureg.kelvin):g}',
        'Hvap': f'at {csstate.T.to(ureg.kelvin):g}',
        'Svap': f'at {csstate.T.to(ureg.kelvin):g}',
        'Tsat': f'at {csstate.P.to(ureg.megapascal):g}',
    }
    
    print(csstate.report(additional_vars=additional_vars, show_parameters=args.show_props, 
                    property_notes=property_notes))

def delta_subcommand(args):
    csstate1 = CSState(T=args.T1, P=args.P1).set_compound(args.n)
    csstate2 = CSState(T=args.T2, P=args.P2).set_compound(args.n)
    if args.n is not None:
        csstate1.set_compound(args.n)
        csstate2.set_compound(args.n)
    if args.Tc is not None:
        csstate1.Tc = args.Tc
        csstate2.Tc = args.Tc
    if args.Pc is not None:
        csstate1.Pc = args.Pc
        csstate2.Pc = args.Pc
    if args.Cp is not None:
        csstate1.Cp = args.Cp
        csstate2.Cp = args.Cp

    for p in 'TPx':
        v = getattr(args, f'{p}1', None)
        if v is not None:
            setattr(csstate1, p, v)
        v = getattr(args, f'{p}2', None)
        if v is not None:
            setattr(csstate2, p, v)
    additional_vars = ['Tc', 'Pc', 'Tr', 'Pr', 'Z', 'h_departure', 's_departure', 'neg_h_dep_over_Tc_read', 'neg_s_dep_read']
    state_1 = csstate1.report(additional_vars=additional_vars)
    state_2 = csstate2.report(additional_vars=additional_vars)
    delta = csstate1.delta(csstate2, additional_vars=['Pv', 'Z'])
    print(f"State-change calculations for {args.n} using {csstate1.description}:")
    if args.show_props or args.show_states:
        print()
        two_states = ["State 1:                                                     State 2:"]
        for line1, line2 in zip(state_1.splitlines(), state_2.splitlines()):
            two_states.append(f"{line1:<56s}     {line2}")
        print("\n".join(two_states))
        print()
        print("Property changes:")
    for p in ['T', 'P', 'h', 's', 'u', 'v', 'Pv', 'Z']:
        if p in delta:
            val = delta[p]
            eq = ' =' if p == 'Pv' else '  ='
            print(f'Δ{p}{eq} {val: 6g}')

def cli():
    subcommands = {
        'state': dict(
            func = state_subcommand,
            help = 'Look up corresponding states properties for a single state'
        ),
        'delta': dict(
            func = delta_subcommand,
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

    crit_args = [
        ('n', 'component', 'component name (e.g., methane, ethane, etc.)', str, False),
        ('Pc', 'critical_pressure', 'critical pressure (if component not specified)', float, False),
        ('Tc', 'critical_temperature', 'critical temperature in K (if component not specified)', float, False),
    ]

    state_args = [
        ('P', 'pressure', 'pressure in MPa', float, False),
        ('T', 'temperature', 'temperature in K (always in K)', float, False),
        ('x', 'vapor_fraction', 'vapor fraction (dimensionless)', float, False),
    ]
    for prop, long_arg, explanation, arg_type, required in state_args + crit_args:
        command_parsers['state'].add_argument(
            f'-{prop}',
            f'--{long_arg}',
            dest=prop,
            type=arg_type,
            required=required,
            help=explanation
        )
    command_parsers['state'].add_argument(
        '--Cp',
        nargs=4,
        type=float,
        metavar=('CpA', 'CpB', 'CpC', 'CpD'),
        help='heat capacity polynomial coefficients A, B, C, D (J/mol-K, J/mol-K^2, J/mol-K^3, J/mol-K^4) (if component not specified)',
        default=None
    )
    command_parsers['state'].add_argument(
        '--show-props',
        default=False,
        action=ap.BooleanOptionalAction,
        help='also show all critical properties and Cp coefficients used'
    )
    
    delta_args = [
        ('P1', 'pressure1', 'pressure of state 1 in MPa', float, True),
        ('T1', 'temperature1', 'temperature of state 1 in K', float, True),
        ('P2', 'pressure2', 'pressure of state 2 in MPa', float, True),
        ('T2', 'temperature2', 'temperature of state 2 in K', float, True),
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
    for prop, long_arg, explanation, arg_type, required in crit_args:
        command_parsers['delta'].add_argument(
            f'-{prop}',
            f'--{long_arg}',
            dest=prop,
            type=arg_type,
            required=required,
            help=explanation
        )
    command_parsers['delta'].add_argument(
        '--Cp',
        nargs=4,
        type=float,
        metavar=('CpA', 'CpB', 'CpC', 'CpD'),
        help='heat capacity polynomial coefficients A, B, C, D (J/mol-K, J/mol-K^2, J/mol-K^3, J/mol-K^4) (if component not specified)',
        default=None
    )
    command_parsers['delta'].add_argument(
        '--show-props',
        default=False,
        action=ap.BooleanOptionalAction,
        help='show properties for both states in addition to the differences'
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