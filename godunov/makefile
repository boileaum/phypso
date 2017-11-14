F2PY=f2py
FC=gfortran
FCFLAGS="-fcheck=bounds"
OPTFLAG=-O3
F2PYFLAGS=-DF2PY_REPORT_ON_ARRAY_COPY=1 --quiet --f90exec=$(FC) --opt=$(OPTFLAG) --f90flags=$(FCFLAGS)
#F2PYFLAGS=--quiet -lgomp --f90exec=$(FC) --opt=$(OPTFLAG) --f90flags=$(FCFLAGS)

.PHONY: clean

all: pythran

fortran:
	$(F2PY) $(F2PYFLAGS) -m godunov_fortran -c godunov.f90

pythran:
	pythran riemann_burgers.py -o riemann_burgers_pythran.so
	pythran riemann_stvenant.py -o riemann_stvenant_pythran.so

clean:
	rm -rf *.o *.mod *.so *.tmp *.pyf *.pyc build/* *.so.dSYM __pycache__
