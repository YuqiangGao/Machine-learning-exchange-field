! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   GAP (Gaussian Approximation Potental)
! HND X   
! HND X
! HND X   Portions of GAP were written by Albert Bartok-Partay, Gabor Csanyi, 
! HND X   Copyright 2006-2021.
! HND X
! HND X   Portions of GAP were written by Noam Bernstein as part of
! HND X   his employment for the U.S. Government, and are not subject
! HND X   to copyright in the USA.
! HND X
! HND X   GAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   GAP is distributed in the hope that it will be useful for non-commercial
! HND X   academic research, but WITHOUT ANY WARRANTY; without even the implied 
! HND X   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! HND X   ASL for more details.
! HND X
! HND X   You should have received a copy of the ASL along with this program
! HND X   (e.g. in a LICENSE.md file); if not, you can write to the original licensors,
! HND X   Gabor Csanyi or Albert Bartok-Partay. The ASL is also published at
! HND X   http://github.com/gabor1/ASL
! HND X
! HND X   When using this software, please cite the following reference:
! HND X
! HND X   A. P. Bartok et al Physical Review Letters vol 104 p136403 (2010)
! HND X
! HND X   When using the SOAP kernel or its variants, please additionally cite:
! HND X
! HND X   A. P. Bartok et al Physical Review B vol 87 p184115 (2013)  ! use error_module 
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!  descriptors for spin  written by Yuqiang Gao

#include "error.inc"

module descriptors_module
   use iso_fortran_env
   use error_module
   use system_module, only : dp, print, optional_default, system_timer, operator(//), split_string, string_to_int, split_string_simple, inoutput, OUTPUT, PRINT_VERBOSE, PRINT_NERD 
   use linkedlist_module
   use units_module 
   use periodictable_module
   use linearalgebra_module 
   use dictionary_module 
   use paramreader_module 
   use atoms_module 
   use atoms_types_module  
   use topology_module 
   use mpi_context_module
   use table_module 
#ifdef DESCRIPTORS_NONCOMMERCIAL
   use permutation_maker_module  
#endif
   use CInOutput_module 
   use clusters_module 
   use connection_module 
   use angular_functions_module  

   implicit none

   private
#ifdef GAP_VERSION
   integer, parameter :: gap_version = GAP_VERSION
#else
   integer, parameter :: gap_version = 0
#endif

   

  subroutine sosd_3b_ex_calc(this,at,descriptor_out,do_descriptor,do_grad_descriptor,args_str,error)
      type real_2d_array

         type(real_2d), dimension(:,:,:), allocatable :: x
      endtype real_2d_array

      type(sosd_3b), intent(in) :: this
      type(atoms), intent(in) :: at
      type(descriptor_data), intent(out) :: descriptor_out
      logical, intent(in), optional :: do_descriptor, do_grad_descriptor
      character(len=*), intent(in), optional :: args_str 
      integer, optional, intent(out) :: error

      type(Dictionary) :: params
      character(STRING_LENGTH) :: atom_mask_name
      logical :: has_atom_mask_name
      logical, dimension(:), pointer :: atom_mask_pointer

      type(cplx_1d), dimension(:), allocatable, save :: SphericalY_ij, SphericalE_ij,SphericalE_i,SphericalY_ik, SphericalE_ik
      type(cplx_2d), dimension(:), allocatable, save :: grad_SphericalY_ij
       complex(dp) :: spherical0

      !SPEED type(cplx_1d), dimension(:,:,:), allocatable :: fourier_so3
      !SPEED type(cplx_2d), dimension(:,:,:), allocatable :: grad_fourier_so3
      type(real_1d), dimension(:,:,:), allocatable, save :: fourier_so3_r, fourier_so3_i, global_fourier_so3_r, global_fourier_so3_i
      type(real_2d), dimension(:,:,:), allocatable, save :: grad_fourier_so3_r, grad_fourier_so3_i
      real(dp), allocatable :: t_g_r(:,:), t_g_i(:,:), t_f_r(:,:), t_f_i(:,:), t_g_f_rr(:,:), t_g_f_ii(:,:)
      integer :: alpha

      logical :: my_do_descriptor, my_do_grad_descriptor, do_two_l_plus_one
      integer :: d,x, i, j, n, a, b, k, l, m, i_pow, i_coeff, l_n_neighbours, n_i, &
         n_descriptors, n_cross, i_species, j_species, ia, jb, i_desc_i, &
         xml_version, sum_l_n_neighbours, i_pair, i_pair_i, n_index,l1,l2,m1,m2,nn , lr
      integer, dimension(3) :: shift_ij, shift_ik
      integer, dimension(:), allocatable :: i_desc
      integer, dimension(:,:), allocatable :: rs_index
      real(dp) :: r_ij, r_ik, arg_bess, arg_bess_e, mo_spher_bess_fi_ki_0, mo_spher_bess_fi_ki_l, mo_spher_bess_fi_ki_lm, mo_spher_bess_fi_ki_lmm, mo_spher_bess_fi_ki_l_e, mo_spher_bess_fi_ki_l_em, mo_spher_bess_fi_ki_l_emm, mo_spher_bess_fi_ki_lp, &
         exp_p, exp_m, f_cut, df_cut, norm_descriptor_i, radial_decay, dradial_decay, norm_radial_decay,tmpgrad_i, tmpgrad_r,  realnumber
      real(dp), dimension(3) :: u_ij, d_ij, e_i, e_j, m_j, u_ik, d_ik, e_k, m_k
      real(dp) :: e_ij, norm_mj,norm_mk,e_ik
      real(dp), dimension(:,:), allocatable, save :: radial_fun, radial_coefficient_j,radial_coefficient_k, grad_radial_fun, grad_radial_coefficient, grad_descriptor_i
      real(dp), dimension(:), allocatable, save ::  radial_e_fun 
      real(dp), dimension(:), allocatable, save :: descriptor_i
      real(dp), dimension(:), allocatable :: global_fourier_so3_r_array, global_fourier_so3_i_array
      type(real_2d_array), dimension(:), allocatable :: global_grad_fourier_so3_r_array, global_grad_fourier_so3_i_array
      integer, dimension(total_elements) :: species_map
      complex(dp), allocatable, save :: sphericalycartesian_all_t(:,:),sphericalycartesian_all_tk(:,:), gradsphericalycartesian_all_t(:,:,:),sphericalycartesian_all_e(:,:),sphericalycartesian_all_ek(:,:),sphericalycartesian_all_ei(:,:)
      complex(dp) :: g_tmp(3),c_tmp
      integer :: max_n_neigh
!$omp threadprivate(radial_fun, radial_coefficient, grad_radial_fun, grad_radial_coefficient)
!$omp threadprivate(sphericalycartesian_all_t, gradsphericalycartesian_all_t)
!$omp threadprivate(fourier_so3_r, fourier_so3_i)
!$omp threadprivate(SphericalY_ij,grad_SphericalY_ij)
!$omp threadprivate(descriptor_i, grad_descriptor_i)
!$omp threadprivate(grad_fourier_so3_r, grad_fourier_so3_i)

      INIT_ERROR(error)

      call system_timer('sosd_ex_calc')

      if(.not. this%initialised) then
         RAISE_ERROR("soap_calc: descriptor object not initialised", error)
      endif

      species_map = 0
      do i_species = 1, this%n_species
         if(this%species_Z(i_species) == 0) then
            species_map = 1
         else
            species_map(this%species_Z(i_species)) = i_species
         endif
      enddo

      my_do_descriptor = optional_default(.false., do_descriptor)
      my_do_grad_descriptor = optional_default(.false., do_grad_descriptor)

      if( .not. my_do_descriptor .and. .not. my_do_grad_descriptor ) return

      has_atom_mask_name = .false. ! allow atom mask column in the atom table
      atom_mask_pointer => null()  ! allow atom mask column in the atom table
      xml_version = 1423143769     ! This is the version number where the 2l+1 normalisation of soap vectors was introduced
      if(present(args_str)) then
         call initialise(params)
         
         call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
            help_string="Name of a logical property in the atoms object. For atoms where this property is " // &
            "true, descriptors are calculated.")

         call param_register(params, 'xml_version', '1423143769', xml_version, &
            help_string="Version of GAP the XML potential file was created")

         if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='soap_calc args_str')) then
            RAISE_ERROR("soap_calc failed to parse args_str='"//trim(args_str)//"'", error)
         endif
         
         call finalise(params)

         if( has_atom_mask_name ) then
            if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
               RAISE_ERROR("soap_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
            endif
         else
            atom_mask_pointer => null()
         endif

      endif

      if( this%cutoff_dexp > 0 ) then
         if( this%cutoff_rate == 0.0_dp ) then
            norm_radial_decay = 1.0_dp
         else
            norm_radial_decay = this%cutoff_rate / ( 1.0_dp + this%cutoff_rate )
         endif
      else
         norm_radial_decay = 1.0_dp
      endif

      do_two_l_plus_one = (xml_version >= 1423143769)

      allocate(rs_index(2,this%n_max*this%n_species))
      i = 0
      do i_species = 1, this%n_species
         do a = 1, this%n_max
            i = i + 1
            rs_index(:,i) = (/a,i_species/)
         enddo
      enddo

      call finalise(descriptor_out)

      d =sosd_3b_ex_dimensions(this,error)
     call print('the dimension ======'//d)

      if(associated(atom_mask_pointer)) then
         call descriptor_sizes(this,at,n_descriptors,n_cross, &
            mask=atom_mask_pointer,n_index=n_index,error=error)
      else
         call descriptor_sizes(this,at,n_descriptors,n_cross,n_index=n_index,error=error)
      endif

      allocate(descriptor_out%x(n_descriptors))
      allocate(i_desc(at%N))



      max_n_neigh = 0
      do n_i = 1, at%N
         max_n_neigh = max(max_n_neigh, n_neighbours(at, n_i))
      end do

!$omp parallel default(none) shared(this,my_do_grad_descriptor,d,max_n_neigh) private(i_species, a, l, n_i)
      allocate(descriptor_i(d))
            descriptor_i=0.0_dp
      if(my_do_grad_descriptor) allocate(grad_descriptor_i(d,3))

      allocate(radial_fun(0:this%l_max, this%n_max), radial_e_fun(0:this%l_max) , radial_coefficient_j(0:this%l_max, this%n_max),radial_coefficient_k(0:this%l_max, this%n_max)) 
      !SPEED allocate(fourier_so3(0:this%l_max,this%n_max,this%n_species), SphericalY_ij(0:this%l_max))
    !  allocate(fourier_so3_r(0:this%l_max,this%n_max,this%n_species), fourier_so3_i(0:this%l_max,this%n_max,this%n_species), SphericalE_ij(0:this%l_max), SphericalY_ij(0:this%l_max),SphericalE_ik(0:this%l_max), SphericalY_ik(0:this%l_max),SphericalE_i(0:this%l_max))
  allocate( SphericalE_ij(0:this%l_max), SphericalY_ij(0:this%l_max),SphericalE_ik(0:this%l_max), SphericalY_ik(0:this%l_max),SphericalE_i(0:this%l_max))


      if(my_do_grad_descriptor) then
         allocate(grad_radial_fun(0:this%l_max, this%n_max), grad_radial_coefficient(0:this%l_max, this%n_max))
         allocate(grad_SphericalY_ij(0:this%l_max))
      endif
    

      allocate(sphericalycartesian_all_t(0:this%l_max, -this%l_max:this%l_max))
         allocate(sphericalycartesian_all_e(0:this%l_max, -this%l_max:this%l_max))
                allocate(sphericalycartesian_all_tk(0:this%l_max, -this%l_max:this%l_max))
	  allocate(sphericalycartesian_all_ek(0:this%l_max, -this%l_max:this%l_max))
          allocate(sphericalycartesian_all_ei(0:this%l_max, this%l_max:this%l_max))
      if(my_do_grad_descriptor) then
          allocate(gradsphericalycartesian_all_t(0:this%l_max, -this%l_max:this%l_max, 3))
      end if

    !  do i_species = 1, this%n_species
    !     do a = 1, this%n_max
    !        do l = 0, this%l_max
    !           allocate(fourier_so3_r(l,a,i_species)%m(-l:l))
    !           allocate(fourier_so3_i(l,a,i_species)%m(-l:l))
    !           fourier_so3_r(l,a,i_species)%m(:) = 0.0_dp
    !           fourier_so3_i(l,a,i_species)%m(:) = 0.0_dp
    !        enddo
    !     enddo
    !  enddo
    

      do l = 0, this%l_max
         allocate(SphericalY_ij(l)%m(-l:l))
         allocate(SphericalE_ij(l)%m(-l:l))
                allocate(SphericalY_ik(l)%m(-l:l))
                allocate(SphericalE_ik(l)%m(-l:l))
        
         if(my_do_grad_descriptor) allocate(grad_SphericalY_ij(l)%mm(3,-l:l))
      enddo
         
      do l = 0, this%l_max
      allocate(SphericalE_i(l)%m(-l:l))
      enddo
      if (my_do_grad_descriptor) then
          !SPEED allocate( grad_fourier_so3(0:this%l_max,this%n_max,n_neighbours(at,i,max_dist=this%cutoff)) )
          allocate( grad_fourier_so3_r(0:this%l_max,this%n_max,max_n_neigh) )
          allocate( grad_fourier_so3_i(0:this%l_max,this%n_max,max_n_neigh) )
          do n_i=1, max_n_neigh
              do a = 1, this%n_max
                  do l = 0, this%l_max
                     !SPEED allocate(grad_fourier_so3(l,a,n_i)%mm(3,-l:l))
                     !SPEED grad_fourier_so3(l,a,n_i)%mm(:,:) = CPLX_ZERO
                     allocate(grad_fourier_so3_r(l,a,n_i)%mm(3,-l:l))
                     allocate(grad_fourier_so3_i(l,a,n_i)%mm(3,-l:l))
                  end do
              end do
          end do
      endif
!$omp end parallel

 
      i_desc = 0
      i_desc_i = 0
      do i = 1, at%N

         if( .not. any( at%Z(i) == this%Z ) .and. .not. any(this%Z == 0) ) cycle

         if(associated(atom_mask_pointer)) then
            if(.not. atom_mask_pointer(i)) cycle
         endif

         i_desc_i = i_desc_i + 1
         i_desc(i) = i_desc_i

         if(.not. this%global) then ! atomic SOAP
            if(my_do_descriptor) then
               allocate(descriptor_out%x(i_desc_i)%data(d))
               !slow, no need
               !descriptor_out%x(i_desc_i)%data = 0.0_dp
               allocate(descriptor_out%x(i_desc_i)%ci(n_index))
               descriptor_out%x(i_desc_i)%has_data = .false.
               descriptor_out%x(i_desc_i)%covariance_cutoff = 1.0_dp
            endif
            if(my_do_grad_descriptor) then
               l_n_neighbours = n_neighbours(at,i,max_dist=this%cutoff)

               allocate(descriptor_out%x(i_desc_i)%grad_data(d,3,0:l_n_neighbours))
               allocate(descriptor_out%x(i_desc_i)%ii(0:l_n_neighbours))
               allocate(descriptor_out%x(i_desc_i)%pos(3,0:l_n_neighbours))
               allocate(descriptor_out%x(i_desc_i)%has_grad_data(0:l_n_neighbours))
               ! slow, no need
               ! descriptor_out%x(i_desc_i)%grad_data = 0.0_dp
               descriptor_out%x(i_desc_i)%grad_data(:,:,0) = 0.0_dp
               descriptor_out%x(i_desc_i)%ii = 0
               descriptor_out%x(i_desc_i)%pos = 0.0_dp
               descriptor_out%x(i_desc_i)%has_grad_data = .false.

               allocate(descriptor_out%x(i_desc_i)%grad_covariance_cutoff(3,0:l_n_neighbours))
               descriptor_out%x(i_desc_i)%grad_covariance_cutoff = 0.0_dp
            endif
         endif
      enddo

      allocate( &
         global_fourier_so3_r_array((this%l_max+1)**2 * this%n_max * this%n_species), &
         global_fourier_so3_i_array((this%l_max+1)**2 * this%n_max * this%n_species), &
         global_grad_fourier_so3_r_array( count(i_desc/=0) ), &
         global_grad_fourier_so3_i_array( count(i_desc/=0) ) )




  
!$omp parallel do schedule(dynamic) default(none) shared(this, at, descriptor_out, my_do_descriptor, my_do_grad_descriptor, d, i_desc, species_map, rs_index, do_two_l_plus_one) &
!$omp shared(global_grad_fourier_so3_r_array, global_grad_fourier_so3_i_array, norm_radial_decay) &
!$omp private(i, j, i_species, j_species, a, b, l, m, n, n_i, r_ij, u_ij, d_ij, shift_ij, i_pow, i_coeff, ia, jb, alpha, i_desc_i) &
!$omp private(c_tmp) &
!$omp private(t_g_r, t_g_i, t_f_r, t_f_i, t_g_f_rr, t_g_f_ii) &
!$omp private(f_cut, df_cut, arg_bess, exp_p, exp_m, mo_spher_bess_fi_ki_l, mo_spher_bess_fi_ki_lp, mo_spher_bess_fi_ki_lm, mo_spher_bess_fi_ki_lmm, norm_descriptor_i) &
!$omp private(radial_decay, dradial_decay) &
!$omp reduction(+:global_fourier_so3_r_array,global_fourier_so3_i_array)
      do i = 1, at%N

         if(i_desc(i) == 0) then
            cycle
         else
            i_desc_i = i_desc(i)
         endif

         if(.not.this%global) then
            if(my_do_descriptor) then
               descriptor_out%x(i_desc_i)%ci(1) = i
               descriptor_out%x(i_desc_i)%has_data = .true.
            endif
            if(my_do_grad_descriptor) then
               descriptor_out%x(i_desc_i)%ii(0) = i
               descriptor_out%x(i_desc_i)%pos(:,0) = at%pos(:,i)
               descriptor_out%x(i_desc_i)%has_grad_data(0) = .true.
            endif
         endif

         !do a = 1, this%n_max
         !   radial_fun(0,a) = exp( -this%alpha * this%r_basis(a)**2 ) !* this%r_basis(a)
         !enddo
         !radial_coefficient(0,:) = matmul( radial_fun(0,:), this%transform_basis )
         radial_fun(0,:) = 0.0_dp
         radial_fun(0,1) = 1.0_dp
         radial_coefficient_j(0,:) = matmul( radial_fun(0,:), this%cholesky_overlap_basis)
         radial_coefficient_k(0,:) = matmul( radial_fun(0,:), this%cholesky_overlap_basis)
   
      !   do i_species = 1, this%n_species
      !      do a = 1, this%n_max
      !         !SPEED fourier_so3(0,a,i_species)%m(0) = radial_coefficient(0,a) * SphericalYCartesian(0,0,(/0.0_dp, 0.0_dp, 0.0_dp/))
      !         if( this%central_reference_all_species .or. this%species_Z(i_species) == at%Z(i) .or. this%species_Z(i_species) == 0 ) then
      !            fourier_so3_r(0,a,i_species)%m(:) = this%central_weight * real(radial_coefficient_j(0,a) * SphericalYCartesian(0,0,(/0.0_dp, 0.0_dp, 0.0_dp/)), dp)
      !            fourier_so3_i(0,a,i_species)%m(:) = this%central_weight * aimag(radial_coefficient_j(0,a) * SphericalYCartesian(0,0,(/0.0_dp, 0.0_dp, 0.0_dp/)))
      !         else
      !            fourier_so3_i(0,a,i_species)%m(:) = 0.0_dp
      !            fourier_so3_r(0,a,i_species)%m(:) = 0.0_dp
      !         endif

       !        do l = 1, this%l_max
       !           !SPEED fourier_so3(l,a,i_species) = CPLX_ZERO
       !           fourier_so3_r(l,a,i_species)%m(:) = 0.0_dp
       !           fourier_so3_i(l,a,i_species)%m(:) = 0.0_dp
       !        enddo
       !     enddo
       !  enddo
           i_pow=0
   call print('start calculate the descriptor'//i_pow)
 
! soap_calc 20 takes 0.0052 s
         n_i = 0
         do n = 1, n_neighbours(at,i)-1    !  留下一个是为了下一个k
            
            j = neighbour(at, i, n, distance = r_ij, cosines=u_ij, diff=d_ij, e_ij=e_ij, e_i=e_i, e_j=e_j,m_j=m_j, shift=shift_ij)  
               
              if( r_ij >= this%cutoff ) cycle
                           n_i = n_i + 1
                do nn=n+1,n_neighbours(at,i)
                   
                   k = neighbour(at, i, nn, distance = r_ik, cosines=u_ik, diff=d_ik, e_ij=e_ik, e_i=e_i, e_j=e_k,m_j=m_k, shift=shift_ik)  
                   
                      if( r_ik >= this%cutoff ) cycle

          

            i_species = species_map(at%Z(j))
            if( i_species == 0 ) cycle

            if(.not. this%global .and. my_do_grad_descriptor) then
               descriptor_out%x(i_desc_i)%ii(n_i) = j
               descriptor_out%x(i_desc_i)%pos(:,n_i) = at%pos(:,j) + matmul(at%lattice,shift_ij)
               descriptor_out%x(i_desc_i)%has_grad_data(n_i) = .true.
            endif

            f_cut = coordination_function(r_ij,this%cutoff, this%cutoff_transition_width)
            radial_decay = ( 1.0_dp + this%cutoff_rate ) / ( this%cutoff_rate + ( r_ij / this%cutoff_scale )**this%cutoff_dexp )
            radial_decay = norm_radial_decay * radial_decay

            f_cut = f_cut * radial_decay

            do a = 1, this%n_max
               arg_bess = 2.0_dp * this%alpha * r_ij * this%r_basis_3b(a) ! argument of modified spherical Bessel function of the first kind
		 arg_bess_e=2.0_dp * this%alpha_e
                mo_spher_bess_fi_ki_0= (exp(arg_bess_e)-exp(-arg_bess_e))/(2*arg_bess_e)
               exp_p = exp( -this%alpha*( r_ij + this%r_basis_3b(a) )**2 )
               exp_m = exp( -this%alpha*( r_ij - this%r_basis_3b(a) )**2 )

               do l = 0, this%l_max
                  if( l == 0 ) then
                      if(arg_bess_e == 0.0_dp) then
                          mo_spher_bess_fi_ki_l_e=1.0_dp
                       else
                          mo_spher_bess_fi_ki_l_e  = 0.5_dp * (exp(arg_bess_e) - exp(-arg_bess_e))/arg_bess_e 
			mo_spher_bess_fi_ki_l_em  = 0.5_dp * (exp(arg_bess_e) + exp(-arg_bess_e))/arg_bess_e 
                       endif
                     if(arg_bess == 0.0_dp) then
                        !mo_spher_bess_fi_ki_l = 1.0_dp
                        mo_spher_bess_fi_ki_l = exp( -this%alpha * (this%r_basis_3b(a)**2 + r_ij**2) ) !this is the coefficient  !modified spherical bessel function l=0
                       
			if(my_do_grad_descriptor) mo_spher_bess_fi_ki_lp = 0.0_dp                 ! first kind  l=0,x=0; 1
                     else                                                                           !  l/=0,x=0; 0
                        !mo_spher_bess_fi_ki_lm = cosh(arg_bess)/arg_bess
                        !mo_spher_bess_fi_ki_l = sinh(arg_bess)/arg_bess
                        mo_spher_bess_fi_ki_lm = 0.5_dp * (exp_m + exp_p) / arg_bess
                        mo_spher_bess_fi_ki_l  = 0.5_dp * (exp_m - exp_p) / arg_bess  ! function for l=0 mo sp bessel. i0(X)=sinhX/X
               
                      !  if(my_do_grad_descriptor) mo_spher_bess_fi_ki_lp = mo_spher_bess_fi_ki_lm - (2*l+1)*mo_spher_bess_fi_ki_l / arg_bess
                     endif
                  else
                      if(arg_bess_e == 0.0_dp) then
                              mo_spher_bess_fi_ki_l_e=0.0_dp
                       else
                         mo_spher_bess_fi_ki_l_emm = mo_spher_bess_fi_ki_l_em
                        mo_spher_bess_fi_ki_l_em = mo_spher_bess_fi_ki_l_e 
                       mo_spher_bess_fi_ki_l_e = mo_spher_bess_fi_ki_l_emm - (2*l-1)*mo_spher_bess_fi_ki_l_em / arg_bess_e 
                      endif
   
                     if(arg_bess == 0.0_dp) then
                        mo_spher_bess_fi_ki_l = 0.0_dp   !  l/=0,x=0; 0
                    
                      !  if(my_do_grad_descriptor) mo_spher_bess_fi_ki_lp = 0.0_dp
                     else
                        mo_spher_bess_fi_ki_lmm = mo_spher_bess_fi_ki_lm
                        mo_spher_bess_fi_ki_lm = mo_spher_bess_fi_ki_l
			
                      !  if(my_do_grad_descriptor) then
                      !     mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lp
                      !     mo_spher_bess_fi_ki_lp = mo_spher_bess_fi_ki_lm - (2*l+1)*mo_spher_bess_fi_ki_l / arg_bess
                      !  else
                           mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lmm - (2*l-1)*mo_spher_bess_fi_ki_lm / arg_bess   ! get il(r_basis)
                          
			! endif
                     endif
                  endif

                  !radial_fun(l,a) = exp( -this%alpha * (this%r_basis(a)**2 + r_ij**2) ) * mo_spher_bess_fi_ki_l !* this%r_basis(a)
                       
		 radial_e_fun(l) = mo_spher_bess_fi_ki_l_e 
                  
                radial_fun(l,a) = mo_spher_bess_fi_ki_l 
                 ! if(my_do_grad_descriptor) grad_radial_fun(l,a) = -2.0_dp * this%alpha * r_ij * mo_spher_bess_fi_ki_l + &
                 !    l*mo_spher_bess_fi_ki_l / r_ij + mo_spher_bess_fi_ki_lp * 2.0_dp * this%alpha * this%r_basis(a)
            enddo
         enddo
       



            radial_coefficient_j = matmul( radial_fun, this%transform_basis )
            radial_coefficient_j = radial_coefficient_j * f_cut

            sphericalycartesian_all_t = SphericalYCartesian_all(this%l_max, d_ij)   ! this is a function defined in angular_function.f95
	    SphericalYCartesian_all_e = SphericalYCartesian_all(this%l_max, e_j)
             SphericalYCartesian_all_ei = SphericalYCartesian_all(this%l_max, e_i)
              if(my_do_grad_descriptor) GradSphericalYCartesian_all_t = GradSphericalYCartesian_all(this%l_max, e_j)
            do l = 0, this%l_max
               do m = -l, l
                  SphericalY_ij(l)%m(m) = SphericalYCartesian_all_t(l,m)
		  SphericalE_ij(l)%m(m) = SphericalYCartesian_all_e(l,m)
                 
                !  spherical0=roundE_ij(l)%m(m)   
                !  call print(real(spherical0))
                !  call print(aimag(spherical0))
                  if(my_do_grad_descriptor) grad_SphericalY_ij(l)%mm(:,m) = GradSphericalYCartesian_all_t(l,m,:)
               enddo                                  
            enddo
             norm_mj=sqrt(m_j(1)**2+m_j(2)**2+m_j(3)**2)
             do l = 0, this%l_max
               do m = -l, l
                 SphericalE_i(l)%m(m) = SphericalYCartesian_all_ei(l,m)
                    if(my_do_grad_descriptor) grad_SphericalY_ij(l)%mm(:,m) = GradSphericalYCartesian_all_t(l,m,:)
               enddo                                  
            enddo

 !=======================================================  ================ = = 
              f_cut = coordination_function(r_ik,this%cutoff, this%cutoff_transition_width)
            radial_decay = ( 1.0_dp + this%cutoff_rate ) / ( this%cutoff_rate + ( r_ik / this%cutoff_scale )**this%cutoff_dexp )
            radial_decay = norm_radial_decay * radial_decay

            f_cut = f_cut * radial_decay

            do a = 1, this%n_max
               arg_bess = 2.0_dp * this%alpha * r_ik * this%r_basis_3b(a) ! argument of modified spherical Bessel function of the first kind
		 arg_bess_e=2.0_dp * this%alpha_e
                mo_spher_bess_fi_ki_0= (exp(arg_bess_e)-exp(-arg_bess_e))/(2*arg_bess_e)
               exp_p = exp( -this%alpha*( r_ik+ this%r_basis_3b(a) )**2 )
               exp_m = exp( -this%alpha*( r_ik- this%r_basis_3b(a) )**2 )

               do l = 0, this%l_max
                  if( l == 0 ) then
                      if(arg_bess_e == 0.0_dp) then
                          mo_spher_bess_fi_ki_l_e=1.0_dp
                       else
                          mo_spher_bess_fi_ki_l_e  = 0.5_dp * (exp(arg_bess_e) - exp(-arg_bess_e))/arg_bess_e 
			mo_spher_bess_fi_ki_l_em  = 0.5_dp * (exp(arg_bess_e) + exp(-arg_bess_e))/arg_bess_e 
                       endif
                     if(arg_bess == 0.0_dp) then
                        !mo_spher_bess_fi_ki_l = 1.0_dp
                        mo_spher_bess_fi_ki_l = exp( -this%alpha * (this%r_basis_3b(a)**2 + r_ik**2) ) !this is the coefficient  !modified spherical bessel function l=0
                       
			if(my_do_grad_descriptor) mo_spher_bess_fi_ki_lp = 0.0_dp                 ! first kind  l=0,x=0; 1
                     else                                                                           !  l/=0,x=0; 0
                        !mo_spher_bess_fi_ki_lm = cosh(arg_bess)/arg_bess
                        !mo_spher_bess_fi_ki_l = sinh(arg_bess)/arg_bess
                        mo_spher_bess_fi_ki_lm = 0.5_dp * (exp_m + exp_p) / arg_bess
                        mo_spher_bess_fi_ki_l  = 0.5_dp * (exp_m - exp_p) / arg_bess  ! function for l=0 mo sp bessel. i0(X)=sinhX/X
               
                      !  if(my_do_grad_descriptor) mo_spher_bess_fi_ki_lp = mo_spher_bess_fi_ki_lm - (2*l+1)*mo_spher_bess_fi_ki_l / arg_bess
                     endif
                  else
                      if(arg_bess_e == 0.0_dp) then
                              mo_spher_bess_fi_ki_l_e=0.0_dp
                       else
                         mo_spher_bess_fi_ki_l_emm = mo_spher_bess_fi_ki_l_em
                        mo_spher_bess_fi_ki_l_em = mo_spher_bess_fi_ki_l_e 
                       mo_spher_bess_fi_ki_l_e = mo_spher_bess_fi_ki_l_emm - (2*l-1)*mo_spher_bess_fi_ki_l_em / arg_bess_e 
                      endif
   
                     if(arg_bess == 0.0_dp) then
                        mo_spher_bess_fi_ki_l = 0.0_dp   !  l/=0,x=0; 0
                    
                      !  if(my_do_grad_descriptor) mo_spher_bess_fi_ki_lp = 0.0_dp
                     else
                        mo_spher_bess_fi_ki_lmm = mo_spher_bess_fi_ki_lm
                        mo_spher_bess_fi_ki_lm = mo_spher_bess_fi_ki_l
			
                      !  if(my_do_grad_descriptor) then
                      !     mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lp
                      !     mo_spher_bess_fi_ki_lp = mo_spher_bess_fi_ki_lm - (2*l+1)*mo_spher_bess_fi_ki_l / arg_bess
                      !  else
                           mo_spher_bess_fi_ki_l = mo_spher_bess_fi_ki_lmm - (2*l-1)*mo_spher_bess_fi_ki_lm / arg_bess   ! get il(r_basis)
                          
			! endif
                     endif
                  endif

                  !radial_fun(l,a) = exp( -this%alpha * (this%r_basis(a)**2 + r_ij**2) ) * mo_spher_bess_fi_ki_l !* this%r_basis(a)
                       
		 radial_e_fun(l) = mo_spher_bess_fi_ki_l_e 
                  
                radial_fun(l,a) = mo_spher_bess_fi_ki_l 
                 ! if(my_do_grad_descriptor) grad_radial_fun(l,a) = -2.0_dp * this%alpha * r_ij * mo_spher_bess_fi_ki_l + &
                 !    l*mo_spher_bess_fi_ki_l / r_ij + mo_spher_bess_fi_ki_lp * 2.0_dp * this%alpha * this%r_basis(a)
            enddo
         enddo
       



            radial_coefficient_k= matmul( radial_fun, this%transform_basis )
            radial_coefficient_k= radial_coefficient_k * f_cut

            sphericalycartesian_all_tk = SphericalYCartesian_all(this%l_max, d_ik)
	    SphericalYCartesian_all_ek = SphericalYCartesian_all(this%l_max, e_k)
                  if(my_do_grad_descriptor) GradSphericalYCartesian_all_t = GradSphericalYCartesian_all(this%l_max, e_k)
            do l = 0, this%l_max
               do m = -l, l
                  SphericalY_ik(l)%m(m) = SphericalYCartesian_all_tk(l,m)
		  SphericalE_ik(l)%m(m) = SphericalYCartesian_all_ek(l,m)

                !  spherical0=roundE_ij(l)%m(m)   
                !  call print(real(spherical0))
                !  call print(aimag(spherical0))
                  if(my_do_grad_descriptor) grad_SphericalY_ij(l)%mm(:,m) = GradSphericalYCartesian_all_t(l,m,:)
               enddo                                  
            enddo
             norm_mk=sqrt(m_k(1)**2+m_k(2)**2+m_k(3)**2)
!===         = = = ========================== =================== real part===============================================            
                     
               i_pow = 0
               do a=1, this%n_max
                 do b=1, a
                     do l1 = 0, this%l_max
                       do l2=0, this%l_max
                         do l=0, this%l_max
                           do lr=0, this%l_max

                          if((abs(l1-l2) .LE. l).AND.(l .LE. (l1+l2)))  then
                            i_pow = i_pow + 1
                           endif
                         do m=-lr,lr
                          do m1=-l1,l1
                           do m2= -l2,l2 
                           if(cg_calculate(l1,m1,l2,m2,l,m1+m2) .NE. 0)  then
                          
                           c_tmp  = radial_coefficient_j(lr,a) * SphericalY_ij(lr)%m(m)*CONJG(radial_coefficient_k(lr,b) * SphericalY_ik(lr)%m(m)) *(4/mo_spher_bess_fi_ki_0)**3*sqrt(real(1/(4*real(l)+2)))*radial_e_fun(l)*radial_e_fun(l1)*radial_e_fun(l2)* cg_calculate(l1,m1,l2,m2,l,m1+m2)*&
                                                                                                     SphericalE_i(l)%m(m1+m2)*SphericalE_ij(l1)%m(-m1)*(-11)**m1*SphericalE_ik(l2)%m(-m2)*(-12)**m2
                       ! notice::::   ======  Cnlm : the n refers to different orbitals and it is a series of equally descrete distance "a", as in radial_coefficient(0,a), here a=1,this%n_max
                               descriptor_i(i_pow) = descriptor_i(i_pow)  + real(c_tmp)
                                         endif
                   !  fourier_so3_i(l,a,b,i_species)%m(m) = fourier_so3_i(l,a,i_species)%m(m) + aimag(c_tmp)   !    complex      a  b  different  n  n'  sum  over  j

                   !  if(my_do_grad_descriptor) then
                   !     do m=-l,l    
                   !      g_tmp(:) = radial_coefficient(0,a) * SphericalY_ij(0)%m(0) *(4/mo_spher_bess_fi_ki_0)**2*sqrt((2*PI)/(2*l+1))*radial_e_fun(l)*SphericalE_ij(l)%m(-m) * (-1)**m*grad_SphericalY_ij(l)%mm(:,m)/norm_mj
      
                        ! grad_fourier_so3_r(l,a,n_i)%mm(:,m) = grad_fourier_so3_r(l,a,n_i)%mm(:,m) + real(c_tmp)
                        ! grad_fourier_so3_i(l,a,n_i)%mm(:,m) = grad_fourier_so3_i(l,a,n_i)%mm(:,m) + imag(c_tmp)
                    !    grad_fourier_so3_r(l,a,n_i)%mm(:,m) = real(g_tmp)
                    !    grad_fourier_so3_i(l,a,n_i)%mm(:,m) = aimag(g_tmp)   ! n_i list all of the neighbours; i_species sum over all neighbours with same species, split into different species.
                    !    enddo
                    ! endif ! my_do_grad_descriptor
                                              
                    enddo  !  m2
                  enddo  !m1  
               enddo  ! m
              enddo !lr  
               enddo   ! l1
            enddo   !  l2
                        enddo     !    l
           enddo     ! b
         enddo  ! a

!================================  imaginary part=============================
               do a=1, this%n_max
                 do b=1, a   
                     do l1 = 0, this%l_max
                       do l2=0, this%l_max
                         do l=0, this%l_max  
                           do lr=0, this%l_max
                          
                          if((abs(l1-l2) .LE. l).AND.(l .LE. (l1+l2)))  then
                            i_pow = i_pow + 1
                           endif
                           do m=-lr,lr
                            do m1=-l1,l1
                             do m2= -l2,l2 
                            if(cg_calculate(l1,m1,l2,m2,l,m1+m2) .NE. 0)  then
                       
                           c_tmp  = radial_coefficient_j(lr,a) * SphericalY_ij(lr)%m(m)*CONJG(radial_coefficient_k(lr,b) * SphericalY_ik(lr)%m(m)) *(4/mo_spher_bess_fi_ki_0)**3*sqrt(real(1/(4*real(l)+2)))*radial_e_fun(l)*radial_e_fun(l1)*radial_e_fun(l2)* cg_calculate(l1,m1,l2,m2,l,m1+m2)*&
                                                                                                     SphericalE_i(l)%m(m1+m2)*SphericalE_ij(l1)%m(-m1)*(-11)**m1*SphericalE_ik(l2)%m(-m2)*(-12)**m2
                          descriptor_i(i_pow) = descriptor_i(i_pow)  + imag(c_tmp)
                       ! notice::::   ======  Cnlm : the n refers to different orbitals and it is a series of equally descrete distance "a", as in radial_coefficient(0,a), here a=1,this%n_max
 
                   !  fourier_so3_i(l,a,b,i_species)%m(m) = fourier_so3_i(l,a,i_species)%m(m) + imag(c_tmp)   !    complex      a  b  different  n  n'  sum  over  j

                   !  if(my_do_grad_descriptor) then
                   !     do m=-l,l    
                   !      g_tmp(:) = radial_coefficient(0,a) * SphericalY_ij(0)%m(0) *(4/mo_spher_bess_fi_ki_0)**2*sqrt((2*PI)/(2*l+1))*radial_e_fun(l)*SphericalE_ij(l)%m(-m) * (-1)**m*grad_SphericalY_ij(l)%mm(:,m)/norm_mj
      
                        ! grad_fourier_so3_r(l,a,n_i)%mm(:,m) = grad_fourier_so3_r(l,a,n_i)%mm(:,m) + real(c_tmp)
                        ! grad_fourier_so3_i(l,a,n_i)%mm(:,m) = grad_fourier_so3_i(l,a,n_i)%mm(:,m) + imag(c_tmp)
                    !    grad_fourier_so3_r(l,a,n_i)%mm(:,m) = real(g_tmp)
                    !    grad_fourier_so3_i(l,a,n_i)%mm(:,m) = aimag(g_tmp)   ! n_i list all of the neighbours; i_species sum over all neighbours with same species, split into different species.
                    !    enddo
                    ! endif ! my_do_grad_descriptor
                   endif                           
                    enddo  !  m2
                  enddo  !m1    
                enddo  ! m
              enddo !lr
               enddo   ! l1
            enddo   !  l2
               enddo     !    l
           enddo     ! b
         enddo  ! a



        enddo ! j   this is loop over all of the neighbour of atom i  n_i
       enddo  ! k     

     !    if(this%global .and. my_do_grad_descriptor) then
     !       global_grad_fourier_so3_r_array(i_desc_i)%x = grad_fourier_so3_r
     !       global_grad_fourier_so3_i_array(i_desc_i)%x = grad_fourier_so3_i
     !       do n_i = lbound(grad_fourier_so3_r,3), ubound(grad_fourier_so3_r,3)
     !          do a = lbound(grad_fourier_so3_r,2), ubound(grad_fourier_so3_r,2)
     !             do l = lbound(grad_fourier_so3_r,1), ubound(grad_fourier_so3_r,1)
     !                global_grad_fourier_so3_r_array(i_desc_i)%x(l,a,n_i)%mm = grad_fourier_so3_r(l,a,n_i)%mm
     !                global_grad_fourier_so3_i_array(i_desc_i)%x(l,a,n_i)%mm = grad_fourier_so3_i(l,a,n_i)%mm
     !             enddo ! l
     !          enddo ! a
     !       enddo ! n_i
      !   endif

! I comment the lines below, I do not need and to avoid errors.
       !  if(this%global) then
       !     i_coeff = 0
       !     do ia = 1, this%n_species*this%n_max
       !        a = rs_index(1,ia)
       !        i_species = rs_index(2,ia)
       !        do l = 0, this%l_max
       !           global_fourier_so3_r_array(i_coeff+1:i_coeff+2*l+1) = global_fourier_so3_r_array(i_coeff+1:i_coeff+2*l+1) + fourier_so3_r(l,a,i_species) 
       !           global_fourier_so3_i_array(i_coeff+1:i_coeff+2*l+1) = global_fourier_so3_i_array(i_coeff+1:i_coeff+2*l+1) + fourier_so3_i(l,a,i_species) 
       !           i_coeff = i_coeff + 2*l+1
       !        enddo
       !     enddo
       !  endif
       !  descriptor_i=0.0_dp
       !  i_pow = 0
       !  do ia = 1, this%n_species*this%n_max
 !           a = rs_index(1,ia)
 !           i_species = rs_index(2,ia)
 !           do jb = 1, ia
 !              b = rs_index(1,jb)
 !              j_species = rs_index(2,jb)
 !    
 !              if(this%diagonal_radial .and. a /= b) cycle

!              do l1 = 0, this%l_max
 !                 do l2=0, this%l_max
 !                    do m1=-l1,l1
 !                      do m2= -l2,l2 
                    
 !                        do l=0,this%l_max
 !                         ! if(abs(l1-l2)<= l .and. l<= l1+l2 .and. -l<= (m1+m2) .and. (m1+m2) <= l ) then     
 !                              if(cg_calculate(l1,m1,l2,m2,l,m1+m2)==0)  cycle
 !                          i_pow = i_pow + 1
 !                  descriptor_i(i_pow) = cg_calculate(l1,m1,l2,m2,l,m1+m2)*real( 1/sqrt(4.0_dp*l+2.0_dp)*(4.0_dp/mo_spher_bess_fi_ki_0)*radial_e_fun(l)*SphericalE_i(l)%m(m1+m2)*&
 !                                   COMPLEX((fourier_so3_r(l1,a,i_species)%m(m1)*fourier_so3_r(l2,b,j_species)%m(m2) - fourier_so3_i(l1,a,i_species)%m(m1)*fourier_so3_i(l2,b,j_species)%m(m2)),(fourier_so3_r(l1,a,i_species)%m(m1)*fourier_so3_i(l2,b,j_species)%m(m2) + fourier_so3_i(l1,a,i_species)%m(m1)* fourier_so3_r(l2,b,j_species)%m(m2)) ))                
 !                      if(cg_calculate(l1,m1,l2,m2,l,m1+m2)==0)  call print('==cg = 0==')
 !                if( ia /= jb ) descriptor_i(i_pow) = descriptor_i(i_pow) * SQRT_TWO   
 !                          
 !                       ! endif
 !                    enddo ! l
 !                   enddo  !m2
 !                  enddo  !m1
 !                enddo !l2
 !              enddo !l1
 !           enddo !jb
 !        enddo !ia
     call print('this is i_pow======'//i_pow//'===d==='//d)
         descriptor_i(d) = 0.0_dp  ! sosd dimension
         norm_descriptor_i = sqrt(dot_product(descriptor_i,descriptor_i))

         if(.not. this%global .and. my_do_descriptor) then
            if(this%normalise .AND.  (norm_descriptor_i .GT. 0)) then
               descriptor_out%x(i_desc_i)%data = descriptor_i / norm_descriptor_i
            else
               descriptor_out%x(i_desc_i)%data = descriptor_i
            endif

            descriptor_out%x(i_desc_i)%data(d) = this%covariance_sigma0
         endif
          
        
!!!!!! grad
         if(my_do_grad_descriptor) then
! soap_calc 33 takes 0.047 s
	    allocate(t_g_r(this%n_max*3,0:this%l_max), t_g_i(this%n_max*3,0:this%l_max))
	    allocate(t_f_r(this%n_max*this%n_species, 1), t_f_i(this%n_max*this%n_species, 1))
	    allocate(t_g_f_rr(this%n_max*3, this%n_max*this%n_species), t_g_f_ii(this%n_max*3, this%n_max*this%n_species))
            !do n_i = 1, n_neighbours(at,i,max_dist=this%cutoff) 
           
         
            n_i = 0
            do n = 1, n_neighbours(at,i) 
               j = neighbour(at, i, n, distance = r_ij) 
               if( r_ij >= this%cutoff ) cycle

               n_i = n_i + 1

               if( species_map(at%Z(j)) == 0 ) cycle

               i_pow = 0
               grad_descriptor_i = 0.0_dp

               !SPEED do ia = 1, this%n_species*this%n_max
               !SPEED    a = rs_index(1,ia) 
               !SPEED    i_species = rs_index(2,ia)
               !SPEED    do jb = 1, ia
               !SPEED       b = rs_index(1,jb)
               !SPEED       j_species = rs_index(2,jb)
               !SPEED       do l = 0, this%l_max
               !SPEED          i_pow = i_pow + 1
               !SPEED          if(at%Z(j) == this%species_Z(i_species) .or. this%species_Z(i_species)==0) grad_descriptor_i(i_pow,:) = grad_descriptor_i(i_pow,:) + real( matmul(conjg(grad_fourier_so3(l,a,n_i)%mm),fourier_so3(l,b,j_species)%m) )
               !SPEED          if(at%Z(j) == this%species_Z(j_species) .or. this%species_Z(j_species)==0) grad_descriptor_i(i_pow,:) = grad_descriptor_i(i_pow,:) + real( matmul(grad_fourier_so3(l,b,n_i)%mm,conjg(fourier_so3(l,a,i_species)%m)) )
               !SPEED          !grad_descriptor_i(i_pow,:) = real( matmul(conjg(grad_fourier_so3(l,a,n_i)%mm),fourier_so3(l,b)%m) + matmul(grad_fourier_so3(l,b,n_i)%mm,conjg(fourier_so3(l,a)%m)) )
               !SPEED          if( ia /= jb ) grad_descriptor_i(i_pow, 1:3) = grad_descriptor_i(i_pow, 1:3) * SQRT_TWO
               !SPEED       enddo !l
               !SPEED    enddo !jb
               !SPEED enddo !ia

               !SPEED do ia = 1, this%n_species*this%n_max
               !SPEED    a = rs_index(1,ia) 
               !SPEED    i_species = rs_index(2,ia)
               !SPEED    do jb = 1, ia
               !SPEED       b = rs_index(1,jb)
               !SPEED       j_species = rs_index(2,jb)
               !SPEED       do l = 0, this%l_max
               !SPEED          i_pow = i_pow + 1
               !SPEED          if(at%Z(j) == this%species_Z(i_species) .or. this%species_Z(i_species)==0) grad_descriptor_i(i_pow,:) = grad_descriptor_i(i_pow,:) + &
               !SPEED             matmul(grad_fourier_so3_r(l,a,n_i)%mm,fourier_so3_r(l,b,j_species)%m) + matmul(grad_fourier_so3_i(l,a,n_i)%mm,fourier_so3_i(l,b,j_species)%m)
               !SPEED          if(at%Z(j) == this%species_Z(j_species) .or. this%species_Z(j_species)==0) grad_descriptor_i(i_pow,:) = grad_descriptor_i(i_pow,:) + &
               !SPEED             matmul(grad_fourier_so3_r(l,b,n_i)%mm,fourier_so3_r(l,a,i_species)%m) + matmul(grad_fourier_so3_i(l,b,n_i)%mm,fourier_so3_i(l,a,i_species)%m)
               !SPEED          if( ia /= jb ) grad_descriptor_i(i_pow, 1:3) = grad_descriptor_i(i_pow, 1:3) * SQRT_TWO
               !SPEED       enddo !l
               !SPEED    enddo !jb
               !SPEED enddo !ia
 
         do l=0, this%l_max
                  do a = 1, this%n_max
                     do alpha=1, 3
                	t_g_r(3*(a-1)+alpha, 1:2*l+1) = grad_fourier_so3_r(l,a,n_i)%mm(alpha,-l:l)
                	t_g_i(3*(a-1)+alpha, 1:2*l+1) = grad_fourier_so3_i(l,a,n_i)%mm(alpha,-l:l)
                     enddo
                  enddo
                  do ia = 1, this%n_species*this%n_max
                     a = rs_index(1,ia)
                     i_species = rs_index(2,ia)
                     
                     t_f_r(ia, 1:2*l+1) = fourier_so3_r(l,a,i_species)%m(-l:l)
                     t_f_i(ia, 1:2*l+1) = fourier_so3_i(l,a,i_species)%m(-l:l)
                  enddo
                  call dgemm('N','T',this%n_max*3, this%n_max*this%n_species, 2*l+1, 1.0_dp, &
                     t_g_r(1,1), size(t_g_r,1), t_f_r(1,1), size(t_f_r,1), 0.0_dp, t_g_f_rr(1,1), size(t_g_f_rr, 1))
                  call dgemm('N','T',this%n_max*3, this%n_max*this%n_species, 2*l+1, 1.0_dp, &
                     t_g_i(1,1), size(t_g_i,1), t_f_i(1,1), size(t_f_i,1), 0.0_dp, t_g_f_ii(1,1), size(t_g_f_ii, 1))
                  !t_g_f_rr = matmul(t_g_r,transpose(t_f_r))
                  !t_g_f_ii = matmul(t_g_i,transpose(t_f_i))
               
                  i_pow = l+1
                  do ia = 1, this%n_species*this%n_max
                     a = rs_index(1,ia)
                     i_species = rs_index(2,ia)
                     do jb = 1, ia !this%n_species*this%n_max !ia
                        b = rs_index(1,jb)
                        j_species = rs_index(2,jb)

                        if(this%diagonal_radial .and. a /= b) cycle
               
                        if(at%Z(j) == this%species_Z(i_species) .or. this%species_Z(i_species)==0) grad_descriptor_i(i_pow, 1:3) = grad_descriptor_i(i_pow, 1:3) + t_g_f_rr(3*(a-1)+1:3*a,jb) + t_g_f_ii(3*(a-1)+1:3*a,jb)
                        if(at%Z(j) == this%species_Z(j_species) .or. this%species_Z(j_species)==0) grad_descriptor_i(i_pow, 1:3) = grad_descriptor_i(i_pow, 1:3) + t_g_f_rr(3*(b-1)+1:3*b,ia) + t_g_f_ii(3*(b-1)+1:3*b,ia)
               

                        if(do_two_l_plus_one) grad_descriptor_i(i_pow, 1:3) = grad_descriptor_i(i_pow, 1:3) / sqrt(2.0_dp * l + 1.0_dp)
                        if( ia /= jb ) grad_descriptor_i(i_pow, 1:3) = grad_descriptor_i(i_pow, 1:3) * SQRT_TWO
                        i_pow = i_pow + this%l_max+1
                     enddo
                  enddo

                  !do a = 1, this%n_max
                  !   do b = 1, a
                  !      grad_descriptor_i(i_pow, 1:3) = t_g_f_rr(3*(a-1)+1:3*a,b) + t_g_f_ii(3*(a-1)+1:3*a,b) + &
                  !                                      t_g_f_rr(3*(b-1)+1:3*b,a) + t_g_f_ii(3*(b-1)+1:3*b,a)
                  !      if( a /= b ) grad_descriptor_i(i_pow, 1:3) = grad_descriptor_i(i_pow, 1:3) * SQRT_TWO
                  !      i_pow = i_pow + this%l_max+1
                  !   end do
                  !end do

               end do !l

                   
                        

 
           

               grad_descriptor_i(d, 1:3) = 0.0_dp
               

               if(.not. this%global) then
                  if( this%normalise ) then
                     descriptor_out%x(i_desc_i)%grad_data(:,:,n_i) = grad_descriptor_i/norm_descriptor_i
                     do k = 1, 3
                        descriptor_out%x(i_desc_i)%grad_data(:,k,n_i) = descriptor_out%x(i_desc_i)%grad_data(:,k,n_i) - descriptor_i * dot_product(descriptor_i,grad_descriptor_i(:,k)) / norm_descriptor_i**3
                      enddo  !  atom  i_desc_i    neighbour n_i
                  else
                     descriptor_out%x(i_desc_i)%grad_data(:,:,n_i) = grad_descriptor_i
                  endif

                  descriptor_out%x(i_desc_i)%grad_data(:,:,0) = descriptor_out%x(i_desc_i)%grad_data(:,:,0) - descriptor_out%x(i_desc_i)%grad_data(:,:,n_i)
               endif
            enddo !n i 
	    deallocate(t_f_r, t_f_i)
	    deallocate(t_g_r, t_g_i)
	    deallocate(t_g_f_rr, t_g_f_ii)




         endif 
!!!!! end grad
      enddo ! i    loop over all of the atom
!$omp end parallel do





      !SPEED if(allocated(fourier_so3)) then
      !SPEED    do i_species = 1, this%n_species
      !SPEED       do a = lbound(fourier_so3,2), ubound(fourier_so3,2)
      !SPEED          do l = lbound(fourier_so3,1), ubound(fourier_so3,1)
      !SPEED             deallocate(fourier_so3(l,a,i_species)%m)
      !SPEED          enddo
      !SPEED       enddo
      !SPEED    enddo
      !SPEED    deallocate(fourier_so3)
      !SPEED endif

!$omp parallel default(none) shared(this, max_n_neigh) private(i_species, a, l, n_i)
  !    if(allocated(fourier_so3_r)) then
  !       do i_species = lbound(fourier_so3_r,3), ubound(fourier_so3_r,3)
  !          do a = lbound(fourier_so3_r,2), ubound(fourier_so3_r,2)
  !             do l = lbound(fourier_so3_r,1), ubound(fourier_so3_r,1)
  !                deallocate(fourier_so3_r(l,a,i_species)%m)
  !             enddo!
  !          enddo
  !       enddo
  !       deallocate(fourier_so3_r)
  !    endif
  !    if(allocated(fourier_so3_i)) then
  !       do i_species = lbound(fourier_so3_i,3), ubound(fourier_so3_i,3)
  !          do a = lbound(fourier_so3_i,2), ubound(fourier_so3_i,2)
  !             do l = lbound(fourier_so3_i,1), ubound(fourier_so3_i,1)
  !                deallocate(fourier_so3_i(l,a,i_species)%m)
  !             enddo
  !          enddo
  !       enddo
  !       deallocate(fourier_so3_i)
  !    endif
      
  
   call print('this is line======10119')
         do l = lbound(SphericalY_ij,1), ubound(SphericalY_ij,1)
           do m=-l,l
          call print('this is line======10133'//SphericalY_ij(0)%m(m))
          enddo
    enddo

   
      if(allocated(SphericalY_ij)) then
          do l = lbound(SphericalY_ij,1), ubound(SphericalY_ij,1)
                    deallocate(SphericalY_ij(l)%m)
             enddo
            deallocate(SphericalY_ij)
        endif
     
  call print('this is line======10128')
                if(allocated(SphericalY_ik)) then
                        do l = lbound(SphericalY_ik,1), ubound(SphericalY_ik,1)
                           deallocate(SphericalY_ik(l)%m)
                       enddo
                      deallocate(SphericalY_ik)
               endif

  call print('this is line======10136')
      if(allocated(SphericalE_ij)) then
         do l = lbound(SphericalE_ij,1), ubound(SphericalE_ij,1)
            deallocate(SphericalE_ij(l)%m)
         enddo
         deallocate(SphericalE_ij)
      endif
  call print('this is line======10143')
         if(allocated(SphericalE_ik)) then
              do l = lbound(SphericalE_ik,1), ubound(SphericalE_ik,1)
                  deallocate(SphericalE_ik(l)%m)
               enddo
             deallocate(SphericalE_ik)
         endif
  call print('this is line======10150')
   
    if(allocated(SphericalE_i)) then
         do l = lbound(SphericalE_i,1), ubound(SphericalE_i,1)
            deallocate(SphericalE_i(l)%m)
         enddo
         deallocate(SphericalE_i)
      endif
  call print('this is line======10158')
      if(allocated(grad_SphericalY_ij)) then
         do l = lbound(grad_SphericalY_ij,1), ubound(grad_SphericalY_ij,1)
            deallocate(grad_SphericalY_ij(l)%mm)
         enddo
         deallocate(grad_SphericalY_ij)
      endif


   call print('this is line======10167')


      if (allocated(sphericalycartesian_all_t)) deallocate(sphericalycartesian_all_t)
       if (allocated(sphericalycartesian_all_e)) deallocate(sphericalycartesian_all_e)
      if (allocated(sphericalycartesian_all_tk)) deallocate(sphericalycartesian_all_tk)
       if (allocated(sphericalycartesian_all_ek)) deallocate(sphericalycartesian_all_ek)

       if (allocated(sphericalycartesian_all_ei)) deallocate(sphericalycartesian_all_ei)
      if (allocated(gradsphericalycartesian_all_t)) deallocate(gradsphericalycartesian_all_t)

      if(allocated(radial_fun)) deallocate(radial_fun)
      if(allocated(radial_e_fun)) deallocate(radial_e_fun)
      if(allocated(radial_coefficient_j)) deallocate(radial_coefficient_j)
            if(allocated(radial_coefficient_k)) deallocate(radial_coefficient_k)
      if(allocated(grad_radial_fun)) deallocate(grad_radial_fun)
      if(allocated(grad_radial_coefficient)) deallocate(grad_radial_coefficient)
      if(allocated(descriptor_i)) deallocate(descriptor_i)
      if(allocated(grad_descriptor_i)) deallocate(grad_descriptor_i)

        if (allocated(grad_fourier_so3_r)) then ! should really check for grad_fourier_so3_i also
            do n_i = 1, max_n_neigh
               do a = 1, this%n_max
                  do l = 0, this%l_max
                     !SPEED deallocate(grad_fourier_so3(l,a,n_i)%mm)
                     if(allocated(grad_fourier_so3_r(l,a,n_i)%mm)) deallocate(grad_fourier_so3_r(l,a,n_i)%mm)
                     if(allocated(grad_fourier_so3_i(l,a,n_i)%mm)) deallocate(grad_fourier_so3_i(l,a,n_i)%mm)
                  enddo
               enddo
            enddo
        endif
        !SPEED deallocate(grad_fourier_so3)
        if (allocated(grad_fourier_so3_r)) deallocate(grad_fourier_so3_r)
        if (allocated(grad_fourier_so3_i)) deallocate(grad_fourier_so3_i)
!$omp end parallel

   !   if(this%global)  I delete all of the code about the global descriptor for one atom specie object. 


   call print('this is line======100201')

  
      if(allocated(global_fourier_so3_r_array)) deallocate(global_fourier_so3_r_array)
      if(allocated(global_fourier_so3_i_array)) deallocate(global_fourier_so3_i_array)
   
      if(allocated(global_grad_fourier_so3_r_array)) then
         do i_desc_i = lbound(global_grad_fourier_so3_r_array,1), ubound(global_grad_fourier_so3_r_array,1)
            if(allocated(global_grad_fourier_so3_r_array(i_desc_i)%x)) then
               do n_i = lbound(global_grad_fourier_so3_r_array(i_desc_i)%x,3), ubound(global_grad_fourier_so3_r_array(i_desc_i)%x,3)
                  do a = lbound(global_grad_fourier_so3_r_array(i_desc_i)%x,2), ubound(global_grad_fourier_so3_r_array(i_desc_i)%x,2)
                     do l = lbound(global_grad_fourier_so3_r_array(i_desc_i)%x,1), ubound(global_grad_fourier_so3_r_array(i_desc_i)%x,1)
                        if(allocated(global_grad_fourier_so3_r_array(i_desc_i)%x(l,a,n_i)%mm)) deallocate(global_grad_fourier_so3_r_array(i_desc_i)%x(l,a,n_i)%mm)
                     enddo ! l
                  enddo ! a
               enddo ! n_i
               deallocate(global_grad_fourier_so3_r_array(i_desc_i)%x)
            endif
         enddo ! i_desc_i
         deallocate(global_grad_fourier_so3_r_array)
      endif
      if(allocated(global_grad_fourier_so3_i_array)) then
         do i_desc_i = lbound(global_grad_fourier_so3_i_array,1), ubound(global_grad_fourier_so3_i_array,1)
            if(allocated(global_grad_fourier_so3_i_array(i_desc_i)%x)) then
               do n_i = lbound(global_grad_fourier_so3_i_array(i_desc_i)%x,3), ubound(global_grad_fourier_so3_i_array(i_desc_i)%x,3)
                  do a = lbound(global_grad_fourier_so3_i_array(i_desc_i)%x,2), ubound(global_grad_fourier_so3_i_array(i_desc_i)%x,2)
                     do l = lbound(global_grad_fourier_so3_i_array(i_desc_i)%x,1), ubound(global_grad_fourier_so3_i_array(i_desc_i)%x,1)
                        if(allocated(global_grad_fourier_so3_i_array(i_desc_i)%x(l,a,n_i)%mm)) deallocate(global_grad_fourier_so3_i_array(i_desc_i)%x(l,a,n_i)%mm)
                     enddo ! l
                  enddo ! a
               enddo ! n_i
               deallocate(global_grad_fourier_so3_i_array(i_desc_i)%x)
            endif
         enddo ! i_desc_i
         deallocate(global_grad_fourier_so3_i_array)
      endif
  
      if(allocated(rs_index)) deallocate(rs_index)
      if(allocated(i_desc)) deallocate(i_desc)
  
      call system_timer('sosd_ex_calc')

   endsubroutine sosd_3b_ex_calc
   
   

endmodule descriptors_module
