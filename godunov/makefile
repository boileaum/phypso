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

pythran_burgers:
	cd burgers && pythran riemann.py -o riemann_pythran.so

pythran_stvenant:
	cd stvenant && pythran riemann.py -o riemann_pythran.so

pythran: pythran_burgers pythran_stvenant

clean_burgers:
	cd burgers && rm -rf *.o *.mod *.so *.tmp *.pyf *.pyc build/* *.so.dSYM __pycache__

clean_stvenant:
	cd stvenant && rm -rf *.o *.mod *.so *.tmp *.pyf *.pyc build/* *.so.dSYM __pycache__

clean: clean_burgers clean_stvenant
	rm -rf *.o *.mod *.so *.tmp *.pyf *.pyc build/* *.so.dSYM __pycache__
