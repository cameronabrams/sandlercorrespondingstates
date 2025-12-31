# Sandlercorrespondingstates

> Corresponding states utilities from Sandler's 5th ed.

Sandlercorrespondingstates implements a python interface to a corresponding states calculations using charts from  _Chemical, Biochemical, and Engineering Thermodynamics_ (5th edition) by Stan Sandler (Wiley, USA). It should be used for educational purposes only.


## Installation 

Sandlercorrespondingstates is available via `pip`:

```sh
pip install sandlercorrespondingstates
```

## Usage

### Command-line

```sh
$ sandlercorrespondingstates state -n methane -P 7.5 -T 400 
T  =  400.00 K
P  =  7.50 mpa
Tc =  190.40 K
Pc =  4.60 mpa
Tr =  2.10
Pr =  1.63
v  =  0.000435 m3/mol
Z  =  0.98
Hdep = -438.15 J/mol = -0.55 cal/mol-K
Sdep = -1.13 J/mol-K = -0.27 cal/mol-K
```

### API

```python
>>> from sandlercorrespondingstates.charts import CorrespondingStatesChartReader
>>> from sandlermisc.gas_constant import GasConstant
>>> from sandlerprops.properties import PropertiesDatabase
>>> db = PropertiesDatabase()
>>> component = db.get_compound('methane')
>>> cs = CorrespondingStatesChartReader()
>>> Rpv = GasConstant("bar", "m3")
>>> result = cs.dimensionalized_lookup(T=400, P=0.5, Tc=component.Tc, Pc=component.Pc/10, R_pv=Rpv)
>>> print(result.report())
T  =  400.00 K
P  =  7.50 mpa
Tc =  190.40 K
Pc =  4.60 mpa
Tr =  2.10
Pr =  1.63
v  =  0.000435 m3/mol
Z  =  0.98
Hdep = -438.15 J/mol = -0.55 cal/mol-K
Sdep = -1.13 J/mol-K = -0.27 cal/mol-K
```

## Release History

* 0.3.0
    * `delta` subcommand added
* 0.2.0
    * `StateReporter` used
* 0.1.2
    * fixed messaging errors
* 0.1.0
    * Initial release

## Meta

Cameron F. Abrams – cfa22@drexel.edu

Distributed under the MIT license. See ``LICENSE`` for more information.

[https://github.com/cameronabrams](https://github.com/cameronabrams/)

## Contributing

1. Fork it (<https://github.com/cameronabrams/sandlercorrespondingstates/fork>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request
