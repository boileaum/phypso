# Burgers' equation solver

Everything in `./burgers/` subdirectory.

## Basic usage

### Get help with:

```
usage: burgers.py [-h] [--tmax final_time] [--nmax number_of_pts] [--profile]
                  [--plot] [--kernel {python,numpy,fortran}]

Solve Burgers problem

optional arguments:
  -h, --help            show this help message and exit
  --tmax final_time     simulation final time
  --nmax number_of_pts  number of grid points
  --profile             activate profiling
  --plot                activate plotting
  --kernel {python,numpy,fortran}
                        select kernel type
```

### Examples:

- Using a basic python kernel (very unefficient):

```
$ ./burgers.py --nmax 1000 --profile
```

- Using the native numpy kernel (much faster):

```
$ ./burgers.py --nmax 1000 --profile --kernel numpy
```

## Use pythran to accelerate Python


### Installation on Mac

- Install the [anaconda suite](https://www.anaconda.com/download/#macos)
- Install pythran with pip

```
pip install pythran
```


### Howto

- Use pythran to compile the python submodule `godunov.py` to produce a `godunov.so` object file:

```
pythran godunov.py
```

- Execute as if it where standard python

``` 
python3 burgers.py 1.0 --nmax 1000 --profile
```

### Acceleration

- Example of execution output using the pythran module:

```
$ ./burgers.py --profile --nmax 1000
tmax = 1.0
nmax = 1000
Mean time [s] over 10 executions = 0.0023488752000048406
L2 error = 0.036801482187378914
```

- Example of execution output using the native python module:

```
$ rm godunov.so 
$ ./burgers.py --profile --nmax 1000
tmax = 1.0
nmax = 1000
Mean time [s] over 10 executions = 0.8791365289000168
L2 error = 0.036801482187378914
```

## Use f2py to accelerate python

### Howto

- Compile the fortran file `godunov.f90` with `f2py`

```
$ make
```

- Run with the `--kernel fortran` option to use the compiled module

```
$ ./burgers.py --profile --nmax 1000 --kernel fortran
```

### Acceleration

```
$ ./burgers.py --profile --nmax 1000 --kernel fortran
tmax = 1.0
nmax = 1000
Mean time [s] over 10 executions = 0.007240815297700464
L2 error = 0.036801482187378914
```
