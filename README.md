# PHYPSO

A set of very simple hyperbolic solvers using techniques to use efficient computing kernels called from a python main program. 

### Table of content

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->


- [Installation](#installation)
    - [Download the source repository](#download-the-source-repository)
    - [Using docker](#using-docker)
    - [Using pip](#using-pip)
- [General Godunov solver](#general-godunov-solver)
    - [Usage](#usage)
- [Saint-Venant's equation solver](#saint-venants-equation-solver)
    - [Usage](#usage-1)
- [Burgers' equation solver](#burgers-equation-solver)
    - [Basic usage](#basic-usage)
        - [Examples:](#examples)
    - [Use `pythran` to accelerate Python](#use-pythran-to-accelerate-python)
        - [Howto](#howto)
        - [Acceleration](#acceleration)
    - [Use `f2py` to accelerate python](#use-f2py-to-accelerate-python)
        - [Howto](#howto-1)
        - [Acceleration](#acceleration-1)
- [Developers' corner](#developers-corner)
    - [Build the docker images](#build-the-docker-images)
    - [Test suite](#test-suite)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


## Installation 

### Download the source repository

Either use git:

```
git clone https://git.unistra.fr/m.boileau/phypso.git
```

or download [zip file](https://git.unistra.fr/m.boileau/phypso/repository/master/archive.zip).

### Using docker

To get a full environment ready for running, [Docker](https://www.docker.com/) is a solution.
From the host, clone the repository:

```
[host] git clone https://git.unistra.fr/m.boileau/phypso.git
```

then run a docker container using the `boileaum/phypso:env` image:

```
[host] docker run -ti -v $(pwd):/home/euler/phypso boileaum/phypso:env
```

> **Note:** to run docker from a Mac with support for matplotlib display, simply run:

> ```
> [host] ./run_docker_mac.sh
> ```

Now from the container:

```
[container] cd phypso
[container] make
```

### Using pip

The `docker/Dockerfile-deps` file provides the first required dependencies such as `python3` and `gfortran`.

Install the additional python libraries with `pip`:

```
pip install -r requirements.txt
```

From the project root directory, compile the C-executable, C-library, Cython and pythran versions of all programs using:

```
make
```

## General Godunov solver


### Usage

1. Go to `./godunov` subdirectory.
- Run using the python main program:

```
./godunov.py -h
usage: godunov.py [-h] [--problem {burgers,stvenant}] [--tmax final_time]
                  [--nmax number_of_pts] [--profile] [--plot]
                  [--kernel {python,numpy,pythran,bench}]

Solve hyperbolic problem

optional arguments:
  -h, --help            show this help message and exit
  --problem {burgers,stvenant}
                        select Problem type
  --tmax final_time     simulation final time
  --nmax number_of_pts  number of grid points
  --profile             activate profiling
  --plot                activate plotting
  --kernel {python,numpy,pythran,bench}
                        select kernel type
```

For example, run a benchmark to compare different kernel versions:

```
./godunov.py --kernel bench
```

## Saint-Venant's equation solver

### Usage

1. Go to `./stvenant` subdirectory.
- Run using the python main program:

```
./stvenant.py [-h] [--noplot]
```



## Burgers' equation solver


### Basic usage

1. Go to `./burgers` subdirectory.
- Run using the python main program:

```
./burgers.py -h
usage: burgers.py [-h] [--tmax final_time] [--nmax number_of_pts] [--profile]
                  [--plot] [--kernel {pythran,numba,python,numpy,fortran}]

Solve Burgers problem

optional arguments:
  -h, --help            show this help message and exit
  --tmax final_time     simulation final time
  --nmax number_of_pts  number of grid points
  --profile             activate profiling
  --plot                activate plotting
  --kernel {pythran,numba,python,numpy,fortran}
                        select kernel type
```

#### Examples:

- Using a basic python kernel (very unefficient):

```
$ ./burgers.py --nmax 1000 --profile
```

- Using the native numpy kernel (much faster but not optimal):

```
$ ./burgers.py --nmax 1000 --profile --kernel numpy
```

### Use `pythran` to accelerate Python


From [Pythran website](http://pythran.readthedocs.io/en/latest/):

> Pythran is a Python to c++ compiler for a subset of the Python language, with a focus on scientific computing. It takes a Python module annotated with a few interface description and turns it into a native python module with the same interface, but (hopefully) faster.

#### Howto

Pythran simply requires to explicit the signatures of python functions as comments in python file header:

```
$ head  godunov.py 
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#pythran export sol_exact(float, float)
#pythran export riemann(float, float, float)
#pythran export timeloop(float, int)
#pythran export xmin, xmax
"""
Godunov solver for the Burgers equation
"""
```

- Use pythran to compile the python submodule `godunov.py` as a `pythran_godunov.so` object file:

```
make pythran
```

- Execute as if it were standard python

``` 
./burgers.py --kernel pythran
```

#### Acceleration

- Example of execution output using the pythran module:

```
$ ./burgers.py --nmax 1000 --profile --kernel pythran
tmax = 1.0
nmax = 1000
kernel: pythran
Mean time [s] over 10 executions = 0.002362
L2 error = 0.018877
```

- Example of execution output using the native python module:

```
$ ./burgers.py --nmax 1000 --profile --kernel python
tmax = 1.0
nmax = 1000
kernel: python
Mean time [s] over 10 executions = 1.071154
L2 error = 0.018877
```

### Use `f2py` to accelerate python

#### Howto

- Compile the fortran file `godunov.f90` with `f2py`

```
make fortran
```

- Run with the `--kernel fortran` option to use the compiled module

```
./burgers.py --profile --nmax 1000 --kernel fortran
```

#### Acceleration

```
$ ./burgers.py --profile --nmax 1000 --kernel fortran
tmax = 1.0
nmax = 1000
kernel: fortran
Mean time [s] over 10 executions = 0.004342
L2 error = 0.018877
```

## Developers' corner


### Build the docker images

```
cd docker
docker-compose build
```

### Test suite

Tests are performed with [pytest](https://docs.pytest.org). To run the tests with verbose level and standard output:

```
pytest -vs
```


