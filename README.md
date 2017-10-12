## Installation on Mac

- Install the [anaconda suite](https://www.anaconda.com/download/#macos)
- Install pythran with pip

```
pip install pythran
```

## Use pythran to accelerate Python

### Howto

- Compile the python module with pythran

```
pythran godunov.py
```

- Execute as standard python

``` 
python3 burgers.py 1.0 --nmax 1000 --profile
```

### Acceleration

- Example of execution output using the pythran module:

```
root@e01cd82d6730:/TP/burgers# ./burgers.py --profile --nmax 1000
tmax = 1.0
nmax = 1000
Mean time [s] over 10 executions = 0.0023488752000048406
L2 error = 0.036801482187378914
```

- Example of execution output using the python module:

```
root@e01cd82d6730:/TP/burgers# rm godunov.so 
root@e01cd82d6730:/TP/burgers# ./burgers.py --profile --nmax 1000
tmax = 1.0
nmax = 1000
Mean time [s] over 10 executions = 0.8791365289000168
L2 error = 0.036801482187378914
```