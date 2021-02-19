!  This file is part of prob_module.  
!  

module prob_module
    use xtb_metadynamic
    
    
    implicit none

    integer :: prob_c
    logical :: prob_flag = .false.
    integer :: prob_metadyn_count = 0

    real(8) :: prob_a, prob_l, prob_inv_l
   
    

    contains

    subroutine prob_metadynamic(metavar,nat,at,xyz,ebias,g) !prob code
       use xtb_mctc_accuracy, only : wp
       use xtb_type_setvar
       use xtb_lsrmsd

       implicit none
       type(metadyn_setvar),intent(in) :: metavar
       integer, intent(in)    :: nat
       integer, intent(in)    :: at(nat)
       real(wp),intent(in)    :: xyz(3,nat)
       real(wp),intent(inout) :: ebias
       real(wp),intent(inout) :: g(3,nat)
       real(wp),allocatable   :: xyzref(:,:),grad(:,:),xyzdup(:,:)
       real(wp) :: U(3,3), x_center(3), y_center(3)
       real(wp) :: etmp,rmsdval,e
       integer  :: i,j,k,iref,iat
    
       real(wp) :: prob_kpush, prob_alp, prob_freq

       prob_inv_l =  1.0d0/prob_l

       prob_kpush = 0.084_wp
       prob_alp = 0.04_wp
       prob_freq = 4.0_wp

       if(metavar%nstruc < 1 ) return
      
      ! prob_flag = .false.
      if (.false.) then 
         if (metavar%nat == 0) then
            allocate( xyzref(3,nat), grad(3,nat),source = 0.0_wp )
            !$omp parallel default(none) &
            !$omp shared(metavar,nat,xyz, prob_a, prob_l, prob_inv_l) & ! omp shared(metavar,nat,xyz) 
            !$omp private(grad,xyzref,U,x_center,y_center,rmsdval,e,etmp) &
            !$omp reduction(+:ebias,g)
            !$omp do schedule(dynamic)
            do iref = 1, metavar%nstruc
               grad = 0.0_wp
               xyzref = metavar%xyz(:,:,iref)
               call rmsd(nat,xyz,xyzref,1,U,x_center,y_center,rmsdval, &
                        .true.,grad)
               ! e = metavar%factor(iref) * exp(-metavar%width(iref) * rmsdval**2)    !<-----------------
               e =  prob_a*prob_inv_l*exp(-1.73205_wp * rmsdval * prob_inv_l) * (prob_l + 1.73205_wp*rmsdval)
               
               ebias = ebias + e
               ! etmp = -2.0_wp * metavar%width(iref) * e * rmsdval
               etmp = 1.73205_wp * prob_a * prob_inv_l * exp(-1.73205_wp * rmsdval * prob_inv_l) - 1.73205_wp * prob_inv_l * e
               g = g + etmp*grad
            enddo
            !$omp enddo
            !$omp end parallel
         else
            allocate( xyzref(3,metavar%nat), xyzdup(3,metavar%nat), grad(3,metavar%nat) )
            !$omp parallel default(none) &
            !$omp shared(metavar,nat,xyz, prob_a, prob_l, prob_inv_l) & ! omp shared(metavar,nat,xyz)
            !$omp private(grad,xyzref,xyzdup,U,x_center,y_center,rmsdval,e,etmp,i,iat) &
            !$omp reduction(+:ebias,g)
            !$omp do schedule(dynamic)
            do iref = 1, metavar%nstruc
               grad = 0.0_wp
               do i = 1, metavar%nat
                  iat = metavar%atoms(i)
                  xyzref(:,i) = metavar%xyz(:,iat,iref)
                  xyzdup(:,i) = xyz(:,iat)
               enddo
               call rmsd(metavar%nat,xyzdup,xyzref,1,U,x_center,y_center,rmsdval, &
                        .true.,grad)
               ! e = metavar%factor(iref) * exp(-metavar%width(iref) * rmsdval**2)  !<----------------
               e =  prob_a*prob_inv_l*exp(-1.73205_wp * rmsdval * prob_inv_l) * (prob_l + 1.73205_wp*rmsdval)
                        
               ebias = ebias + e
               ! etmp = -2.0_wp * metavar%width(iref) * e * rmsdval
               etmp = 1.73205_wp * prob_a * prob_inv_l * exp(-1.73205_wp * rmsdval * prob_inv_l) - 1.73205_wp * prob_inv_l * e

               do i = 1, metavar%nat
                  iat = metavar%atoms(i)
                  g(:,iat) = g(:,iat) + etmp*grad(:,i)
               enddo
            enddo
            !$omp enddo
            !$omp end parallel
            deallocate( xyzdup )
         endif
         deallocate( xyzref, grad )
      else

         if (metavar%nat == 0) then
            allocate( xyzref(3,nat), grad(3,nat),source = 0.0_wp )
            !$omp parallel default(none) &
            !$omp shared(metavar,nat,xyz,prob_kpush,prob_alp,prob_freq) &
            !$omp private(grad,xyzref,U,x_center,y_center,rmsdval,e,etmp) &
            !$omp reduction(+:ebias,g)
            !$omp do schedule(dynamic)
            do iref = 1, metavar%nstruc
               grad = 0.0_wp
               xyzref = metavar%xyz(:,:,iref)
               call rmsd(nat,xyz,xyzref,1,U,x_center,y_center,rmsdval, &
                        .true.,grad)
               ! e = metavar%factor(iref) * exp(-metavar%width(iref) * rmsdval**2) 
               e = (prob_kpush * cos(prob_freq*rmsdval)) * exp( -prob_alp*rmsdval**2 )
               ebias = ebias + e
               ! etmp = -2.0_wp * metavar%width(iref) * e * rmsdval 
               etmp =  -2.0_wp*prob_alp*prob_kpush*rmsdval*Cos(prob_freq*rmsdval) * exp(-prob_alp*rmsdval**2) - & 
                        prob_freq*prob_kpush*Sin(prob_freq*rmsdval) * exp(-prob_alp*rmsdval**2)
               ! etmp = (-2*prob_alp*prob_kpush*rmsdval*Cos(prob_freq*rmsdval))/ &
               !    E**(prob_alp*rmsdval**2) - &
               !   (prob_freq*prob_kpush*Sin(prob_freq*rmsdval))/ &
               !    E**(prob_alp*rmsdval**2) 
               g = g + etmp*grad
               ! print *, "nat zero called"
            enddo
            !$omp enddo
            !$omp end parallel
         else
            allocate( xyzref(3,metavar%nat), xyzdup(3,metavar%nat), grad(3,metavar%nat) )
            !$omp parallel default(none) &
            !$omp shared(metavar,nat,xyz) &
            !$omp private(grad,xyzref,xyzdup,U,x_center,y_center,rmsdval,e,etmp,i,iat) &
            !$omp reduction(+:ebias,g)
            !$omp do schedule(dynamic)
            do iref = 1, metavar%nstruc
               grad = 0.0_wp
               do i = 1, metavar%nat
                  iat = metavar%atoms(i)
                  xyzref(:,i) = metavar%xyz(:,iat,iref)
                  xyzdup(:,i) = xyz(:,iat)
               enddo
               call rmsd(metavar%nat,xyzdup,xyzref,1,U,x_center,y_center,rmsdval, &
                        .true.,grad)
               e = metavar%factor(iref) * exp(-metavar%width(iref) * rmsdval**2)
               ebias = ebias + e
               etmp = -2.0_wp * metavar%width(iref) * e * rmsdval
               do i = 1, metavar%nat
                  iat = metavar%atoms(i)
                  g(:,iat) = g(:,iat) + etmp*grad(:,i)
               enddo
            enddo
            !$omp enddo
            !$omp end parallel
            deallocate( xyzdup )
            print *, "nat nonzero called!!!!!"
         endif
         deallocate( xyzref, grad )         

      end if
       
      ! prob_metadyn_count = prob_metadyn_count + 1
      ! !  write(stdout,*) "Calling prob_metadynamic()..." 
      !  print *, 'Call prob_metadynamic()...', prob_metadyn_count
    end subroutine prob_metadynamic
    

    subroutine prob_func()  
      use xtb_mctc_accuracy, only : wp
      use xtb_type_setvar
      use xtb_lsrmsd

      use forpy_mod
      implicit none

        integer :: ierror
        type(module_py) :: opt
        type(ndarray) :: nd_c, nd_b_ub, nd_A_ub, nd_x
        type(object) :: retval, attr
        type(tuple) :: args
        
        real(dp) :: c(2) = [-300._dp, -500._dp]
        real(dp) :: b_ub(3) = [170._dp, 150._dp, 180._dp]
        
        real(dp) :: A_ub(3,2)
        real(dp), dimension(:), pointer :: x
        real(dp) :: objective_fun_value
    
        type(module_py) :: pymod
        type(list) :: paths
        type(tuple) :: argsval
        type(object) :: returnval
        
        write(*,*) "prob_func call! success!"

        A_ub(1,1) = 1.0_dp 
        A_ub(2,1) = 1.0_dp 
        A_ub(3,1) = 0.0_dp 
        
        A_ub(1,2) = 2.0_dp
        A_ub(2,2) = 1.0_dp
        A_ub(3,2) = 3.0_dp
        
        ierror = forpy_initialize()
        ierror = import_py(opt, "scipy.optimize")
        
        ierror = ndarray_create(nd_c, c)
        ierror = ndarray_create(nd_b_ub, b_ub)
        ierror = ndarray_create(nd_A_ub, A_ub)
        
        ierror = tuple_create(args, 3)
        ierror = args%setitem(0, nd_c)
        ierror = args%setitem(1, nd_A_ub)
        ierror = args%setitem(2, nd_b_ub)
        
        ierror = call_py(retval, opt, "linprog", args)
        
        ierror = retval%getattribute(attr, "x")
        ierror = cast(nd_x, attr)
        ierror = nd_x%get_data(x)
        call attr%destroy
        
        ierror = retval%getattribute(attr, "fun")
        ierror = cast(objective_fun_value, attr)
        call attr%destroy
        
        ierror = get_sys_path(paths)
        ierror = paths%append(".")
    
        ierror = import_py(pymod, "pymodule")
        ierror = tuple_create(argsval, 1)
        ierror = argsval%setitem(0, nd_c)
        ierror = call_py(returnval, pymod, "initialize", argsval)
    
        call pymod%destroy
        call paths%destroy
        call argsval%destroy
        call returnval%destroy
 
        print *, "Solution: x = ", x
        print *, "Valueee of objective function: fun = ", objective_fun_value
        
    
        call retval%destroy
        call args%destroy
        call nd_c%destroy
        call nd_b_ub%destroy
        call nd_A_ub%destroy
        call nd_x%destroy
        call opt%destroy
        
        call forpy_finalize()
        
     end subroutine prob_func 
    
    ! subroutine load_metadynamic(metavar,nat,at,xyz)
    !    use xtb_mctc_io, only : stdout
    !    use xtb_mctc_accuracy, only : wp
    !    use xtb_fixparam
    !    use xtb_readin
    !    implicit none
    !    type(metadyn_setvar),intent(inout) :: metavar
    !    integer, intent(in)    :: nat
    !    integer, intent(in)    :: at(nat)
    !    real(wp),intent(in)    :: xyz(3,nat)
    
    !    integer  :: nstruc
    
    !    if (.not.allocated(metavar%fname)) return
    
    !    nstruc = metavar%maxsave
    !    call readlog(metavar%fname,nat,at,metavar%xyz,nstruc)
    !    metavar%nstruc = nstruc
    !    write(stdout,'(a,1x,i0,1x,a)') &
    !       "metadynamics with", nstruc, "initial structures loaded"
    ! end subroutine load_metadynamic
    
    ! subroutine load_rmsdbias(metavar,nat,at,xyz)
    !    use xtb_mctc_io, only : stdout
    !    use xtb_mctc_accuracy, only : wp
    !    use xtb_mctc_convert, only : aatoau
    !    use xtb_mctc_systools, only : getline
    !    use xtb_fixparam
    !    use xtb_readin
    !    implicit none
    !    type(metadyn_setvar),intent(inout) :: metavar
    !    integer, intent(in)    :: nat
    !    integer, intent(in)    :: at(nat)
    !    real(wp),intent(in)    :: xyz(3,nat)
    
    !    integer :: ii, nn, mm, cc, stat, unit
    !    integer, parameter :: initial_size = 16
    !    real(wp) :: kk, aa, xx, yy, zz
    !    character(len=:), allocatable :: line
    !    character(len=4) :: sym
    !    real(wp), allocatable :: tmp_xyz(:, :, :), tmp_realloc3(:, :, :)
    !    real(wp), allocatable :: tmp_par(:, :), tmp_realloc2(:, :)
    
    !    if (.not.allocated(metavar%fname)) return
    
    !    write(stdout, '("#", *(1x, g0))') &
    !       "Reading bias information from", metavar%fname
    
    !    allocate(tmp_xyz(3, nat, initial_size), tmp_par(2, initial_size))
    
    !    open(unit, file=metavar%fname, iostat=stat)
    !    cc = 0
    !    rdxyz: do while(stat == 0)
    !       call getline(unit, line, stat)
    !       if (stat /= 0) exit rdxyz
    !       read(line, *, iostat=stat) nn
    !       if (stat /= 0) exit rdxyz
    !       if (nn /= nat) then
    !          stat = 1
    !          exit rdxyz
    !       end if
    
    !       call getline(unit, line, stat)
    !       read(line, *, iostat=stat) kk, aa
    !       if (stat /= 0) exit rdxyz
    
    !       nn = size(tmp_xyz, 3)
    !       if (cc >= nn) then
    !          mm = nn + nn/2 + 1
    !          nn = min(nn, mm)
    !          call move_alloc(tmp_xyz, tmp_realloc3)
    !          call move_alloc(tmp_par, tmp_realloc2)
    !          allocate(tmp_xyz(3, nat, mm), tmp_par(2, mm))
    !          tmp_xyz(:, :, :nn) = tmp_realloc3(:, :, :nn)
    !          tmp_par(:, :nn) = tmp_realloc2(:, :nn)
    !          deallocate(tmp_realloc3, tmp_realloc2)
    !       end if
    !       cc = cc + 1
    
    !       tmp_par(:, cc) = [kk, aa]
    
    !       do ii = 1, nat
    !          call getline(unit, line, stat)
    !          if (stat /= 0) exit rdxyz
    !          read(line, *, iostat=stat) sym, xx, yy, zz
    !          if (stat /= 0) exit rdxyz
    !          tmp_xyz(:, ii, cc) = [xx, yy, zz] * aatoau
    !       end do
    !    end do rdxyz
    !    if (is_iostat_end(stat)) stat = 0
    !    if (stat /= 0) return
    !    close(unit, iostat=stat)
    
    !    write(stdout, '("#", *(1x, g0))') &
    !       "Read bias potential for", cc, "structures"
    
    !    call metavar%allocate(nat, cc)
    !    metavar%nstruc = cc
    !    metavar%xyz(:, :, :) = tmp_xyz(:, :, :cc)
    !    metavar%factor(:) = tmp_par(1, :cc)
    !    metavar%width(:) = tmp_par(2, :cc)
    
    ! end subroutine load_rmsdbias
    
    ! subroutine set_metadynamic(metavar,nat,at,xyz)
    !    use xtb_mctc_io, only : stdout
    !    use xtb_mctc_accuracy, only : wp
    !    use xtb_fixparam
    !    use xtb_readin
    !    implicit none
    !    type(metadyn_setvar),intent(inout) :: metavar
    !    integer, intent(in)    :: nat
    !    integer, intent(in)    :: at(nat)
    !    real(wp),intent(in)    :: xyz(3,nat)
    
    !    integer  :: nstruc
    
    !    nstruc = metavar%maxsave
    !    metavar%nstruc = nstruc
    !    metavar%xyz(:,:,nstruc) = xyz
    !    metavar%width=1.0
    !    write(stdout,'(a,1x,i0,1x,a)') &
    !       "metadynamics with", nstruc, "initial structures loaded"
    ! end subroutine set_metadynamic
    
    end module prob_module
    
