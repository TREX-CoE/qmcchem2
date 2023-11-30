default:
	cd .. ; OBJ= SRC= INCLUDES= make src/MAIN/qmc

INCLUDES=-I AO -I SAMPLING -I TOOLS -I JASTROW -I MAIN -I PROPERTIES -I ZMQ
SRC=ZMQ/f77_zmq_module.f90 qmckl_module.F90 trexio_module.F90
OBJ=IRPF90_temp/ZMQ/f77_zmq_module.o IRPF90_temp/qmckl_module.o IRPF90_temp/trexio_module.o

IRPF90_FLAGS ?= --align=32
IRPF90+=$(IRPF90_FLAGS) $(INCLUDES)


LIB+=$(LDFLAGS)
include irpf90.make
export

irpf90.make: $(filter-out IRPF90_temp/%, $(wildcard */*.irp.f)) $(wildcard *.irp.f) $(wildcard *.inc.f) Makefile 
	$(IRPF90)

IRPF90_temp/irp_touches.irp.o: IRPF90_temp/irp_touches.irp.F90
	$(FC) -O0 -c -g  IRPF90_temp/irp_touches.irp.F90 -o IRPF90_temp/irp_touches.irp.o

