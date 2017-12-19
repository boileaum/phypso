F2PY=f2py
FC=gfortran
FCFLAGS="-fcheck=bounds"
OPTFLAG=-O3
F2PYFLAGS=-DF2PY_REPORT_ON_ARRAY_COPY=1 --quiet --f90exec=$(FC) --opt=$(OPTFLAG) --f90flags=$(FCFLAGS)

.PHONY: clean

all: fortran pythran

fortran: godunov_burgers_fortran*.so
pythran: godunov_burgers_pythran.so

godunov_burgers_fortran*.so:
	$(F2PY) $(F2PYFLAGS) -m godunov_burgers_fortran -c godunov_burgers.f90

godunov_burgers_pythran.so:
	pythran godunov_burgers.py -o godunov_burgers_pythran.so

clean:
	rm -rf *.o *.mod *.so *.tmp *.pyf *.pyc build/* *.so.dSYM __pycache__
