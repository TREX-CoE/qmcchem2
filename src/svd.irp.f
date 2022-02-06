 BEGIN_PROVIDER [ logical, use_svd ]
&BEGIN_PROVIDER [ integer, n_svd_coefs_full ]
 implicit none
 BEGIN_DOC
 ! If true, use SVD wave function
 END_DOC
 n_svd_coefs_full = -1
 call get_spindeterminants_n_svd_coefs(n_svd_coefs_full)
 use_svd = n_svd_coefs_full > 0
 if (.not.use_SVD) then
    n_svd_coefs_full = 1
 endif
END_PROVIDER

BEGIN_PROVIDER [ integer, n_svd_coefs ]
 implicit none
 BEGIN_DOC
 ! If true, use SVD wave function
 END_DOC
 integer :: i
 do i=1,n_svd_coefs_full
   if (psi_svd_coefs(i) < ci_threshold) then
     exit
   endif
!  print *,  i, psi_svd_coefs(i)
   n_svd_coefs = n_svd_coefs+1
 enddo
END_PROVIDER

 BEGIN_PROVIDER [ double precision, psi_svd_coefs, (               n_svd_coefs_full) ]
&BEGIN_PROVIDER [ double precision, psi_svd_alpha, (det_alpha_num, n_svd_coefs_full) ]   
&BEGIN_PROVIDER [ double precision, psi_svd_beta , (det_beta_num , n_svd_coefs_full) ]
   implicit none
   BEGIN_DOC
   ! !!!
   ! truncated SVD
   END_DOC
   call get_spindeterminants_psi_svd_coefs(psi_svd_coefs)
   call get_spindeterminants_psi_svd_alpha(psi_svd_alpha)
   call get_spindeterminants_psi_svd_beta(psi_svd_beta)
 END_PROVIDER

