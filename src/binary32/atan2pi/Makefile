FUNCTION_UNDER_TEST := atan2pif

include ../support/Makefile.bivariate

all:: check_special

check_special.o: check_special.c
	$(CC) $(CFLAGS) $(CORE_MATH_DEFINES) -c $(OPENMP) -o $@ $<

check_special: check_special.o $(CORE_MATH_OBJS) $(FUNCTION_UNDER_TEST)_mpfr.o
	$(CC) $(LDFLAGS) $(OPENMP) -o $@ $^ -lmpfr -lm

clean::
	rm -f check_special
