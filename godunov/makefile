.PHONY: clean all pythran clean_burgers clean_stvenant

all: pythran

pythran: burgers/riemann_pythran.so stvenant/riemann_pythran.so

burgers/riemann_pythran.so: burgers/riemann.py
	cd burgers && pythran riemann.py -o riemann_pythran.so

stvenant/riemann_pythran.so: stvenant/riemann.py
	cd stvenant && pythran riemann.py -o riemann_pythran.so

clean_burgers:
	cd burgers && rm -rf *.o *.mod *.so *.tmp *.pyf *.pyc build/* *.so.dSYM __pycache__

clean_stvenant:
	cd stvenant && rm -rf *.o *.mod *.so *.tmp *.pyf *.pyc build/* *.so.dSYM __pycache__

clean: clean_burgers clean_stvenant
	rm -rf *.o *.mod *.so *.tmp *.pyf *.pyc build/* *.so.dSYM __pycache__
