CC = icpc

#CFLAGS= -g -O3 -xHost -mkl -fno-alias -qopenmp
CFLAGS= -g -O0 -xHost -mkl -fno-alias -qopenmp

OPTFLAGS = -qopt-report -qopt-report-file=$@.optrpt 

TARGETS = nbody0
TARGETOBJECTS = nbody0.o
TARGETS1 = nbody1
TARGETOBJECTS1 = nbody1.o

NSIZE=16384

.SUFFIXES: .o .c

all: $(TARGETS) $(TARGETS1)


$(TARGETS): $(TARGETOBJECTS)
	$(CC) -o $@ $(CFLAGS) $^

$(TARGETS1): $(TARGETOBJECTS1)
	$(CC) -o $@ $(CFLAGS) $^

.c.o:
	$(CC) $(CFLAGS) -c $(OPTFLAGS) -o $@ $<

clean: 
	rm -f $(TARGETOBJECTS) $(TARGETS) *.optrpt
