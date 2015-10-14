CC			= gcc
CFLAGS	+= -Wall
LIBRARY	= -L/usr/local/lib -lm -lfftw3

all: filters

filters: filters.c
	$(CC) $(CFLAGS) $(LIBRARY) -o filters filters.c

filterskiss: -g -lm filters.c kiss_fft.c kiss_fftr.c -o filters

#gcc -g -lm filters.c kiss_fft.c kiss_fftr.c -o filters
#gnuplot
#set terminal canvas
#set output 'plot.html'
#load 'filters.gnuplot'