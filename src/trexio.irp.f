BEGIN_PROVIDER [ logical, use_trexio ]
 implicit none
 BEGIN_DOC
 ! Is true, use TREXIO file
 END_DOC
 use_trexio = .False.
 call get_simulation_use_trexio(use_trexio)
END_PROVIDER

BEGIN_PROVIDER [ character*(128), trexio_filename ]
 implicit none
 BEGIN_DOC
 ! Name of the TREXIO file
 END_DOC
 call get_trexio_trexio_file(trexio_filename)

END_PROVIDER

BEGIN_PROVIDER [ integer*8, trexio_file ]
 use trexio
 implicit none
 BEGIN_DOC
 ! File handle for TREXIO
 END_DOC
 integer :: rc
 trexio_file = trexio_open(trim(trexio_filename), 'r', TREXIO_AUTO, rc)
 if (rc /= TREXIO_SUCCESS) then
   print *, irp_here
   print *, 'Unable to open TREXIO file ', trexio_filename
 endif

END_PROVIDER

