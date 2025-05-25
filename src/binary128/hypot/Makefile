FUNCTION_UNDER_TEST := hypotq

CFLAGS += -I../support -DMPFR_WANT_FLOAT128
LDFLAGS += -lquadmath

include ../support/Makefile.bivariate

check_worst.o: check_worst.c
	$(CC) $(OPENMP) $(CFLAGS) $(CORE_MATH_DEFINES) -I . -c -o $@ $<

$(FUNCTION_UNDER_TEST)_mpfr.o: ../support/fenv_mpfr.h

clean::
	rm -f check_worst
