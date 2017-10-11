### Installation

```
brew reinstall boost-python --with-python3 --c++11
```

### Use pythran to accelerate Python:

```
pythran godunov.py
python3 burgers.py 1.0 --nmax 1000 --profilei
```
