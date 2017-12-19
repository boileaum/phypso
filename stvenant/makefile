COPTS= -fPIC -std=c99

.PHONY: clean all C libC pythran Cython

all: C libC Cython pythran
C: stvenant.exe
libC: libstvenant_c.so
pythran: riemann_pythran.so
Cython: cstvenant*.so

stvenant.o: stvenant.c
	gcc $(COPTS) -c stvenant.c 

stvenant.exe: stvenant.o
	gcc $(COPTS) stvenant.c -lm -o stvenant.exe

libstvenant_c.so: stvenant.o
	gcc $(COPTS) -shared -o libstvenant_c.so stvenant.o

cstvenant*.so: cstvenant.pyx cstvenant.pxd
	python3 setup.py build_ext --inplace

riemann_pythran.so: riemann.py
	pythran riemann.py -o riemann_pythran.so

clean:
	rm -rf *.o *.mod *.so *.tmp *.pyf *.pyc build *.so.dSYM __pycache__ cstvenant.c plotriem *.exe

	
