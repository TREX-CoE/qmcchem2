.NOPARALLEL:

LIBS=
PKGS=
OCAMLCFLAGS="-g"
OCAMLOPTFLAGS="opt -O3 -remove-unused-arguments -rounds 16 -inline 100 -inline-max-unroll 100"
OCAMLBUILD=ocamlbuild -j 0 -cflags $(OCAMLCFLAGS) -lflags $(OCAMLCFLAGS) -ocamlopt $(OCAMLOPTFLAGS)
MLLFILES=$(wildcard *.mll) 
MLFILES=$(wildcard *.ml) ezfio.ml Qptypes.ml 
MLIFILES=$(wildcard *.mli) 
ALL_TESTS=$(patsubst %.ml,%.byte,$(wildcard test_*.ml))
ALL_EXE=$(patsubst %.ml,%.native,$(wildcard qp_*.ml)) qmcchem.native

.PHONY: default 


default: $(ALL_EXE) 

tests: $(ALL_TESTS) 


%.inferred.mli: $(MLFILES)
	$(OCAMLBUILD) $*.inferred.mli -use-ocamlfind  $(PKGS) 
	mv _build/$*.inferred.mli .

%.byte: $(MLFILES) $(MLIFILES)
	rm -f -- $*
	$(OCAMLBUILD) $*.byte  -use-ocamlfind  $(PKGS)
	ln -s $*.byte $*

%.native: $(MLFILES) $(MLIFILES)
	rm -f -- $*
	$(OCAMLBUILD) $*.native -use-ocamlfind $(PKGS)
	ln -s $*.native $*

qptypes_generator.byte: qptypes_generator.ml ezfio.ml
	$(OCAMLBUILD) qptypes_generator.byte  -use-ocamlfind 

Qptypes.ml: qptypes_generator.byte
	./qptypes_generator.byte > Qptypes.ml

clean: 
	rm -rf _build Qptypes.ml $(ALL_EXE) $(ALL_TESTS)

