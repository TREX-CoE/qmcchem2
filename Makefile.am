ACLOCAL_AMFLAGS = -I m4 

#bin_PROGRAMS = \
#   ocaml/qmcchem

#   src/MAIN/qmc \
#   src/MAIN/qmc_create_walkers \
#

CLEANFILES = \
 src/MAIN/qmc_create_walkers \
 src/MAIN/qmcchem_info \
 src/MAIN/vmc_test \
 src/MAIN/qmc \
 ocaml/Property.ml \
 ocaml/ezfio.ml \
 ocaml/qmcchem \
 ezfio_config/properties.config 
 

.PHONY: clean_ezfio clean_ocaml

EZFIO_DIR=$(abs_srcdir)/external/ezfio

IRPF90=$(abs_srcdir)/external/irpf90/bin/irpf90
QMCCHEM_PATH=$(abs_srcdir)
ARCHIVE=$(AR) crs

AM_LIBS=$(EZFIO_DIR)/lib/libezfio_irp.a
LIB=$(AM_LIBS) $(LIBS)


export FC FCFLAGS QMCCHEM_PATH IRPF90_FLAGS ARCHIVE LDFLAGS LIB IRPF90

bin_SCRIPTS = bin/qmcchem \
              bin/qmc \
              bin/qmc_create_walkers \
              bin/qmcchem_info


clean_irpf90:
	$(MAKE) -C src veryclean

src/irpf90.make: $(IRPSOURCES)
	$(MAKE) -C src  irpf90.make


src/MAIN/qmc_create_walkers:src/MAIN/qmc
src/MAIN/qmcchem_info:src/MAIN/qmc
src/MAIN/vmc_test:src/MAIN/qmc

src/MAIN/qmc: $(EZFIO_DIR)/lib/libezfio_irp.a src/irpf90.make
	$(MAKE) -C src

src/MAIN/vmc_test: $(EZFIO_DIR)/lib/libezfio_irp.a src/irpf90.make
	$(MAKE) -C src MAIN/vmc_test

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

ocaml/qmcchem: ocaml/Property.ml ocaml/ezfio.ml $(OCAMLSOURCES)
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


bin/qmcchem: ocaml/qmcchem
	cp ocaml/qmcchem bin/qmcchem

bin/qmc: src/MAIN/qmc
	cp src/MAIN/qmc bin/qmc

bin/qmcchem_info: src/MAIN/qmcchem_info
	cp src/MAIN/qmcchem_info bin/qmcchem_info

bin/qmc_create_walkers: src/MAIN/qmc_create_walkers
	cp src/MAIN/qmc_create_walkers bin/qmc_create_walkers

IRPSOURCES= \
  src/constants.F \
  src/det.irp.f \
  src/electrons.irp.f \
  src/ezfio_interface.irp.f \
  src/finish.irp.f \
  src/Makefile \
  src/mo.irp.f \
  src/mo_point.irp.f \
  src/nuclei.irp.f \
  src/point.irp.f \
  src/prepare_walkers.irp.f \
  src/properties.py \
  src/pseudo.irp.f \
  src/psi.irp.f \
  src/qmckl_module.F90 \
  src/qmckl.irp.f \
  src/simulation.irp.f \
  src/svd.irp.f \
  src/trexio_f.f90 \
  src/trexio.irp.f \
  src/types.F \
  src/wf.irp.f \
  src/AO/ao.irp.f \
  src/AO/ao_axis.irp.f \
  src/AO/ao_full.irp.f \
  src/AO/ao_oneD.irp.f \
  src/AO/ao_point.irp.f \
  src/JASTROW/jastrow_1b.irp.f \
  src/JASTROW/jastrow_core.irp.f \
  src/JASTROW/jastrow_full.irp.f \
  src/JASTROW/jastrow_mu.irp.f \
  src/JASTROW/jastrow_param.irp.f \
  src/JASTROW/jastrow_simple.irp.f \
  src/MAIN/.gitignore \
  src/MAIN/admc.irp.f \
  src/MAIN/admc.org \
  src/MAIN/admc.py \
  src/MAIN/admc_client.py \
  src/MAIN/admc_server.py \
  src/MAIN/bench_aos.irp.f \
  src/MAIN/bench_mos.irp.f \
  src/MAIN/qmc.irp.f \
  src/MAIN/qmc_create_walkers.irp.f \
  src/MAIN/qmcchem_info.irp.f \
  src/MAIN/vmc_test.irp.f \
  src/PROPERTIES/properties.irp.f \
  src/PROPERTIES/properties_ci.irp.f \
  src/PROPERTIES/properties_energy.irp.f \
  src/PROPERTIES/properties_general.irp.f \
  src/SAMPLING/block.irp.f \
  src/SAMPLING/brownian_step.irp.f \
  src/SAMPLING/dmc_step.irp.f \
  src/SAMPLING/fkmc_step.irp.f \
  src/SAMPLING/langevin_step.irp.f \
  src/SAMPLING/pdmc_step.irp.f \
  src/SAMPLING/reconfigure.irp.f \
  src/SAMPLING/srmc_step.irp.f \
  src/SAMPLING/vmc_step.irp.f \
  src/TOOLS/Util.irp.f \
  src/TOOLS/determinant.irp.f \
  src/TOOLS/info.irp.f \
  src/TOOLS/invert.irp.f \
  src/TOOLS/random.irp.f \
  src/TOOLS/sort.irp.f \
  src/ZMQ/f77_zmq.h \
  src/ZMQ/f77_zmq_module.f90 \
  src/ZMQ/qmc.irp.f \
  src/ZMQ/sockets.irp.f \
  src/ZMQ/utils.irp.f \
  src/ZMQ/worker.irp.f \
  src/ZMQ/zmq_ezfio.irp.f


OCAMLSOURCES= \
  ocaml/Block.ml \
  ocaml/Command_line.ml \
  ocaml/Command_line.mli \
  ocaml/Default.ml \
  ocaml/Input.ml \
  ocaml/Launcher.ml \
  ocaml/Message.ml \
  ocaml/QmcMd5.ml \
  ocaml/Qmcchem_config.ml \
  ocaml/Qmcchem_dataserver.ml \
  ocaml/Qmcchem_debug.ml \
  ocaml/Qmcchem_edit.ml \
  ocaml/Qmcchem_forwarder.ml \
  ocaml/Qmcchem_info.ml \
  ocaml/Qmcchem_md5.ml \
  ocaml/Qmcchem_result.ml \
  ocaml/Qmcchem_run.ml \
  ocaml/Qmcchem_stop.ml \
  ocaml/Qputils.ml \
  ocaml/Random_variable.ml \
  ocaml/Sample.ml \
  ocaml/Sample.mli \
  ocaml/Scheduler.ml \
  ocaml/Status.ml \
  ocaml/String_ext.ml \
  ocaml/Time.ml \
  ocaml/Watchdog.ml \
  ocaml/_tags \
  ocaml/c_bindings.c \
  ocaml/myocamlbuild.ml \
  ocaml/qmcchem.ml \
  ocaml/qptypes_generator.ml

EXTRA_DIST= README.md scripts/create_properties_python.py scripts/create_properties_ezfio.py $(OCAMLSOURCES) $(IRPSOURCES)
