BEGIN_PROVIDER [ character*(128), trexio_filename ]
 implicit none
 BEGIN_DOC
 ! Name of the TREXIO file
 END_DOC
 call ezfio_get_simulation_trexio_filename(trexio_filename)
END_PROVIDER

