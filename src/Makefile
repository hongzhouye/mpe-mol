LIBXC_ROOT=/home/hzye2011/opt_local/libxc-4.1.1

CC = icc
CFLAGS = -O2
CLIBS = $(LIBXC_ROOT)/lib/libxc.a
CINCS = -I$(LIBXC_ROOT)/include
FC = gfortran
FCFLAGS = -O3 -fopenmp 
FCLIBS = -L/usr/lib -llapack -lblas

mpe.out: cienergy.o  fci.o  GetLSQ.o  GetPdP.o  GKSenergy.o \
	GKS.o  GKSPdP.o  main.o  RunMPE.o \
	fciutils.o  oep.o  partition.o  util90.o \
	libxc_wrapper.o libxc_utils.o
	$(FC) $(FCFLAGS) cienergy.o  fci.o  GetLSQ.o  GetPdP.o  GKSenergy.o \
	GKS.o  GKSPdP.o  main.o  RunMPE.o \
	fciutils.o  oep.o  partition.o  util90.o \
	libxc_wrapper.o libxc_utils.o \
	$(CLIBS) $(FCLIBS) -o mpe.out

cienergy.o: cienergy.F fciutils.o
	$(FC) $(FCFLAGS) -c cienergy.F $(FCLIBS)

fci.o: fci.F fciutils.o
	$(FC) $(FCFLAGS) -c fci.F $(FCLIBS)

GetLSQ.o: GetLSQ.F
	$(FC) $(FCFLAGS) -c GetLSQ.F $(FCLIBS)

GetPdP.o: GetPdP.F fciutils.o
	$(FC) $(FCFLAGS) -c GetPdP.F $(FCLIBS)

GKSenergy.o: GKSenergy.F
	$(FC) $(FCFLAGS) -c GKSenergy.F $(FCLIBS)

GKS.o: GKS.F
	$(FC) $(FCFLAGS) -c GKS.F $(FCLIBS)

GKSPdP.o: GKSPdP.F
	$(FC) $(FCFLAGS) -c GKSPdP.F $(FCLIBS)

main.o: main.F
	$(FC) $(FCFLAGS) -c main.F $(FCLIBS)

RunMPE.o: RunMPE.F
	$(FC) $(FCFLAGS) -c RunMPE.F $(FCLIBS)

fciutils.o: fciutils.f
	$(FC) $(FCFLAGS) -c fciutils.f $(FCLIBS)

oep.o: oep.f libxc_wrapper.c
	$(FC) $(FCFLAGS) -c oep.f $(FCLIBS)

partition.o: partition.f
	$(FC) $(FCFLAGS) -c partition.f $(FCLIBS)

util90.o: util90.f
	$(FC) $(FCFLAGS) -c util90.f $(FCLIBS)

libxc_utils.o: libxc_utils.c libxc_utils.h
	$(CC) $(CFLAGS) -c libxc_utils.c $(CINCS)

libxc_wrapper.o: libxc_wrapper.c libxc_utils.h
	$(CC) $(CFLAGS) -c libxc_wrapper.c $(CINCS)

clean:
	rm -f mpe.out *.o *.mod
