FUNCTION_UNDER_TEST := atanh

include ../support/Makefile.univariate

all:: check_special

check_special: check_special.o $(FUNCTION_UNDER_TEST).o $(FUNCTION_UNDER_TEST)_mpfr.o
	$(CC) $(LDFLAGS) $(OPENMP) -o $@ $^ -lmpfr -lm

check_special.o: check_special.c
	$(CC) $(CFLAGS) $(OPENMP) -o $@ $^ -c

clean::
	rm -f check_special
