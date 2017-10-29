cdef extern from "stvenant.c":
    void riemann(double *wL, double *wR, double xi, double *w)
