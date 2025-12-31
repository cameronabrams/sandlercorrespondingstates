from .charts import *
from sandlerprops.properties import PropertiesDatabase
from sandlermisc.gas_constant import GasConstant
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
    component = db.get_compound(args.component)
    if component is None:
        print(f"Component '{args.component}' not found in database.")
        return
    cs = CorrespondingStatesChartReader()
    Rpv = GasConstant("bar", "m3")
    result = cs.dimensionalized_lookup(
        T = args.temperature,
        P = args.pressure,
        Tc = component.Tc,
        Pc = component.Pc/10,
        R_pv = Rpv
    )
    if result is not None:
        for prop, value in result.items():
            print(f"{prop}: {value}")
    else:
        print("Could not find corresponding states properties for the given inputs.")

def cli():
    subcommands = {
        'state': dict(
            func = state,
            help = 'work with corresponding states for a single state'
        ),
        'delta': dict(
            func = None,
            help = 'work with property differences between two states (not implemented yet)'
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
    command_parsers['state'].add_argument(
        '-T',
        '--temperature',
        type=float,
        help='temperature in K'
    )
    command_parsers['state'].add_argument(
        '-P',
        '--pressure',
        type=float,
        help='pressure in MPa'
    )
    command_parsers['state'].add_argument(
        '-Tc',
        '--critical-temperature',
        type=float,
        help='critical temperature in K (if component not specified)'
    )
    command_parsers['state'].add_argument(
        '-Pc',
        '--critical-pressure',
        type=float,
        help='critical pressure in MPa (if component not specified)'
    )
    command_parsers['state'].add_argument(
        '-n',
        '--component',
        type=str,
        help='component name (e.g., methane, ethane, etc.)'
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