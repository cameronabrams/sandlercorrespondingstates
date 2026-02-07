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
T                      =   400 kelvin
P                      =   7.5 megapascal
v                      =  0.000434569 meter ** 3 / mole
u                      =  216.113 joule / mole
h                      =  3475.38 joule / mole
s                      = -25.7768 joule / kelvin / mole
Pv                     =  3259.27 joule / mole
Tc                     =  190.4 kelvin
Pc                     =    4.6 megapascal
Tr                     =  2.10084 dimensionless
Pr                     =  1.63043 dimensionless
Z                      =   0.98 dimensionless
h_departure            = -438.148 joule / mole
s_departure            = -1.12968 joule / kelvin / mole
neg_h_dep_over_Tc_read =   0.55 calorie / kelvin / mole
neg_s_dep_read         =   0.27 calorie / kelvin / mole
```

```sh
$ sandlercorrespondingstates delta -n methane -T1 400 -P1 1 -T2 500 -P2 2 --show-states
State-change calculations for methane using Corresponding states equation of state:

State 1:                                                     State 2:
T                      =   400 kelvin                        T                      =   500 kelvin
P                      =     1 megapascal                    P                      =     2 megapascal
v                      =  0.00332579 meter ** 3 / mole       v                      =  0.00207862 meter ** 3 / mole
u                      =  587.746 joule / mole               u                      =  4150.18 joule / mole
h                      =  3913.53 joule / mole               h                      =  8307.41 joule / mole
s                      = -7.97796 joule / kelvin / mole      s                      = -4.04963 joule / kelvin / mole
Pv                     =  3325.79 joule / mole               Pv                     =  4157.23 joule / mole
Tc                     =  190.4 kelvin                       Tc                     =  190.4 kelvin
Pc                     =    4.6 megapascal                   Pc                     =    4.6 megapascal
Tr                     =  2.10084 dimensionless              Tr                     =  2.62605 dimensionless
Pr                     =  0.217391 dimensionless             Pr                     =  0.434783 dimensionless
Z                      =      1 dimensionless                Z                      =      1 dimensionless
h_departure            =     -0 joule / mole                 h_departure            = -15.9327 joule / mole
s_departure            = -0.08368 joule / kelvin / mole      s_departure            = -0.2092 joule / kelvin / mole
neg_h_dep_over_Tc_read =      0 calorie / kelvin / mole      neg_h_dep_over_Tc_read =   0.02 calorie / kelvin / mole
neg_s_dep_read         =   0.02 calorie / kelvin / mole      neg_s_dep_read         =   0.05 calorie / kelvin / mole

Property changes:
ΔT  =    100 kelvin
ΔP  =      1 megapascal
Δh  =  4393.88 joule / mole
Δs  =  3.92832 joule / kelvin / mole
Δu  =  3562.43 joule / mole
Δv  = -0.00124717 meter ** 3 / mole
ΔPv =  831.446 joule / mole
ΔZ  =      0 dimensionless
```

### API

```python
>>> from sandlercorrespondingstates import CSSstate
>>> state = CSState(T=400, P=0.5).set_compound('methane')
>>> print(result.report(additional_vars=['Tc', 'Pc', 'Tr', 'Pr', 'Z', 'h_departure', 's_departure', 'neg_h_dep_over_Tc_read', 'neg_s_dep_read']))
T                      =   400 kelvin
P                      =   0.5 megapascal
v                      =  0.00665157 meter ** 3 / mole
u                      =  595.712 joule / mole
h                      =  3921.5 joule / mole
s                      = -2.08929 joule / kelvin / mole
Pv                     =  3325.79 joule / mole
Tc                     =  190.4 kelvin
Pc                     =    4.6 megapascal
Tr                     =  2.10084 dimensionless
Pr                     =  0.108696 dimensionless
Z                      =      1 dimensionless
h_departure            =  7.96634 joule / mole
s_departure            =  0.04184 joule / kelvin / mole
neg_h_dep_over_Tc_read =  -0.01 calorie / kelvin / mole
neg_s_dep_read         =  -0.01 calorie / kelvin / mole
```

## Release History

* 0.5.0
    * pint integration
* 0.4.1
    * help updated
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
