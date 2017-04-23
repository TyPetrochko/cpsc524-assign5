CC = icpc

CFLAGS= -g -O3 -xHost -mkl -fno-alias -qopenmp
#CFLAGS= -g -O0 -xHost -mkl -fno-alias -qopenmp

OPTFLAGS = -qopt-report -qopt-report-file=$@.optrpt 

TARGETS = nbody0
TARGETOBJECTS = nbody0.o

TARGETS1 = nbody1
TARGETOBJECTS1 = nbody1.o

TARGETS2 = nbody2
TARGETOBJECTS2 = nbody2.o

TARGETS2 = nbody3
TARGETOBJECTS2 = nbody3.o

NSIZE=16384

.SUFFIXES: .o .c

all: $(TARGETS) $(TARGETS1) $(TARGETS2) $(TARGETS3)


$(TARGETS): $(TARGETOBJECTS)
	$(CC) -o $@ $(CFLAGS) $^

$(TARGETS1): $(TARGETOBJECTS1)
	$(CC) -o $@ $(CFLAGS) $^

$(TARGETS2): $(TARGETOBJECTS2)
	$(CC) -o $@ $(CFLAGS) $^

$(TARGETS3): $(TARGETOBJECTS3)
	$(CC) -o $@ $(CFLAGS) $^

.c.o:
	$(CC) $(CFLAGS) -c $(OPTFLAGS) -o $@ $<

clean: 
	rm -f $(TARGETOBJECTS) $(TARGETS) *.optrpt
	rm -f $(TARGETOBJECTS1) $(TARGETS1) *.optrpt
	rm -f $(TARGETOBJECTS2) $(TARGETS2) *.optrpt
	rm -f $(TARGETOBJECTS3) $(TARGETS3) *.optrpt

