
.PHONY: clean

all: C libC

C:
	gcc stvenant.c -lm -o stvenant.exe

libC:
	gcc -fPIC -c stvenant.c 
	gcc -shared -o libstvenant_c.so stvenant.o

clean:
	rm -rf *.o *.out *.mod *.so *.tmp *.pyf *.pyc build/* *.so.dSYM
