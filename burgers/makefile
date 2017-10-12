F2PY=f2py
FC=gfortran
FCFLAGS="-fcheck=bounds"
OPTFLAG=-O3
F2PYFLAGS=-DF2PY_REPORT_ON_ARRAY_COPY=1 --quiet --f90exec=$(FC) --opt=$(OPTFLAG) --f90flags=$(FCFLAGS)
#F2PYFLAGS=--quiet -lgomp --f90exec=$(FC) --opt=$(OPTFLAG) --f90flags=$(FCFLAGS)

.PHONY: clean

all: godunov

godunov:
	$(F2PY) $(F2PYFLAGS) -m libgodunov -c godunov.f90

clean:
	rm -rf *.o *.out *.mod *.so *.tmp *.pyf *.pyc build/* *.so.dSYM
