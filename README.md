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
$ sandlercorrespondingstates state -T 400 -P 0.5 -n methane   
Tr: 2.100840336134454
Pr: 0.10869565217391305
v: 0.06643327431492282
Hdep: 4.020845568856118
Sdep: 0.022707829240921282
defaults: {'Z': 0.9987607943819014, 'Hdep': -0.00504729598256478, 'Sdep': -0.005427301443814838}
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
>>> for prop, value in result.items():
>>>     print(f"{prop}: {value}")
Tr: 2.100840336134454
Pr: 0.10869565217391305
v: 0.06643327431492282
Hdep: 4.020845568856118
Sdep: 0.022707829240921282
defaults: {'Z': 0.9987607943819014, 'Hdep': -0.00504729598256478, 'Sdep': -0.005427301443814838}
```

## Release History

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
