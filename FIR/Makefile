CC			= gcc
CFLAGS	+= -Wall
LIBRARY	= -lm

all: runfilters testfilters
runfilters: run-filters.c
	$(CC) $(LIBRARY) -lm run-filters.c filters.c kiss_fft.c kiss_fftr.c -o run-filters
testfilters: test-filters.c
	$(CC) $(LIBRARY) test-filters.c kiss_fft.c kiss_fftr.c filters.c -o test-filters;