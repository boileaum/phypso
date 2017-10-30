COPTS= -fPIC
#COPTS= -Wno-unused-result -Wsign-compare -Wunreachable-code -DNDEBUG -fwrapv -O2 -Wall -Wstrict-prototypes -march=core2 -mtune=haswell -mssse3 -ftree-vectorize -fPIC -fPIE -fstack-protector-strong -O2 -pipe -march=core2 -mtune=haswell -mssse3 -ftree-vectorize -fPIC -fstack-protector-strong -O2 -pipe

.PHONY: clean

all: C libC Cython pythran

C:
	gcc stvenant.c -lm -o stvenant.exe

libC:
	gcc $(COPTS) -c stvenant.c 
	gcc $(COPTS) -shared -o libstvenant_c.so stvenant.o

Cython:
	python3 setup.py build_ext --inplace

pythran:
	pythran riemann.py -o riemann_pythran.so

clean:
	rm -rf *.o *.out *.mod *.so *.tmp *.pyf *.pyc build/* *.so.dSYM
