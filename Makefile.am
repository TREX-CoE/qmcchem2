ACLOCAL_AMFLAGS = -I m4

#bin_PROGRAMS = \
#   ocaml/qmcchem

#   src/MAIN/qmc \
#   src/MAIN/qmc_create_walkers \
#

.PHONY: clean_ezfio clean_ocaml

EZFIO_DIR=$(abs_srcdir)/external/ezfio

IRPF90=$(abs_srcdir)/external/irpf90/bin/irpf90
QMCCHEM_PATH=$(abs_srcdir)
ARCHIVE=$(AR) crs

AM_LIBS=$(EZFIO_DIR)/lib/libezfio_irp.a
LIB=$(LIBS) $(AM_LIBS)


export FC FCFLAGS QMCCHEM_PATH IRPF90_FLAGS ARCHIVE LDFLAGS LIB

bin_EXTRA = ocaml/qmcchem

clean_irpf90:
	$(MAKE) -C src veryclean

src/irpf90.make:
	$(MAKE) -C src  irpf90.make

src/MAIN/qmc src/MAIN/qmc_create_walkers src/MAIN/qmcchem_info: \
    $(EZFIO_DIR)/lib/libezfio_irp.a src/irpf90.make
	$(MAKE) -C src

clean_ezfio:
	$(MAKE) -C $(EZFIO_DIR) clean

clean_ocaml:
	$(MAKE) -C ocaml clean

ocaml/ezfio.ml: $(EZFIO_DIR)/Ocaml/ezfio.ml
	cp $(EZFIO_DIR)/Ocaml/ezfio.ml ocaml/ezfio.ml

ocaml/Property.ml: src/PROPERTIES/properties_general.irp.f \
                   src/PROPERTIES/properties.irp.f  \
                   src/PROPERTIES/properties_energy.irp.f \
                   scripts/create_properties_python.py
	bash -c "source $(srcdir)/qmcchemrc ;\
         $(srcdir)/scripts/create_properties_python.py"

ocaml/qmcchem: ocaml/Property.ml ocaml/ezfio.ml
	$(MAKE) -C ocaml

ezfio_config/properties.config: src/irpf90.make
	python3 scripts/create_properties_ezfio.py

$(EZFIO_DIR)/Ocaml/ezfio.ml: $(EZFIO_DIR)/lib/libezfio_irp.a

$(EZFIO_DIR)/config/qmc.config: ezfio_config/qmc.config
	cp $(srcdir)/ezfio_config/qmc.config $(EZFIO_DIR)/config/qmc.config

$(EZFIO_DIR)/config/properties.config: ezfio_config/properties.config
	cp $(srcdir)/ezfio_config/properties.config $(EZFIO_DIR)/config/properties.config

$(EZFIO_DIR)/lib/libezfio_irp.a: $(EZFIO_DIR)/config/qmc.config $(EZFIO_DIR)/config/properties.config src/irpf90.make
	$(MAKE) -C $(EZFIO_DIR)

