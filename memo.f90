!=================================================================================
! MPI_ALLREDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, COMM, IERROR)
!         <type> SENDBUF(*), RECVBUF(*)
!         INTEGER COUNT, DATATYPE, OP, COMM, IERROR
! comm内の全プロセスのsendbufに、opで指定した演算を施して、全プロセスのrecvbufへ送る。
! http://www.cv.titech.ac.jp/~hiro-lab/study/mpi_reference/chapter3.html#mpi_allreduce
!=================================================================================
   subroutine hecmw_allreduce_i(hecMESH, val, n, ntag)
     use hecmw_util
     implicit none
     integer(kind=kint):: n, ntag
     integer(kind=kint), dimension(n) :: val
     type (hecmwST_local_mesh) :: hecMESH
 #ifndef HECMW_SERIAL
     integer(kind=kint):: ierr
     integer(kind=kint), dimension(:), allocatable :: VALM
  
     allocate (valm(n))
     valm= 0
     if (ntag .eq. hecmw_sum) then
       call mpi_allreduce                                              &
         &       (val, valm, n, mpi_integer, mpi_sum,                       &
         &        hecmesh%MPI_COMM, ierr)
     endif

     !----------------------------------------------------------------------------
     ! val              : 各プロセスからの入力
     ! valm             : 演算処理'mpi_sum'を施し出力した結果
     ! mpi_integer      : データ型
     ! mpi_sum          : 演算処理(和)
     ! mpi_max          : 演算処理(最大値)
     ! mpi_min          : 演算処理(最小値)
     ! hecmesh%MPI_COMM : MPI用コミュニケータ
     ! ierr             : エラーコード
     !----------------------------------------------------------------------------
  
     if (ntag .eq. hecmw_max) then
       call mpi_allreduce                                              &
         &       (val, valm, n, mpi_integer, mpi_max,                       &
         &        hecmesh%MPI_COMM, ierr)
     endif
  
     if (ntag .eq. hecmw_min) then
       call mpi_allreduce                                              &
         &       (val, valm, n, mpi_integer, mpi_min,                       &
         &        hecmesh%MPI_COMM, ierr)
     endif
  
  
     val= valm
     deallocate (valm)
 #endif
   end subroutine hecmw_allreduce_i
!=================================================================================
   subroutine hecmw_allreduce_i1 (hecMESH, s, ntag)
     use hecmw_util
     implicit none
     integer(kind=kint)::  ntag, s
     type (hecmwST_local_mesh) :: hecMESH
 #ifndef HECMW_SERIAL
     integer(kind=kint), dimension(1) :: val
  
     val(1) = s
     call hecmw_allreduce_i(hecmesh, val, 1, ntag )
     s = val(1)
 #endif
   end subroutine hecmw_allreduce_i1
!=================================================================================
   subroutine hecmw_mpc_scale(hecMESH)
    implicit none
    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    integer(kind=kint) :: i, j, k
    real(kind=kreal) :: WVAL

    !$omp parallel default(none),private(i,j,k,WVAL),shared(hecMESH)
    !$omp do
    do i = 1, hecMESH%mpc%n_mpc
      k = hecMESH%mpc%mpc_index(i-1)+1
      ! mpc_index 拘束グループ用一次元インデックス
      WVAL = 1.d0 / hecMESH%mpc%mpc_val(k)
      ! mpc_val 拘束自由度係数
      hecMESH%mpc%mpc_val(k) = 1.d0
      do j = hecMESH%mpc%mpc_index(i-1)+2, hecMESH%mpc%mpc_index(i)
        hecMESH%mpc%mpc_val(j) = hecMESH%mpc%mpc_val(j) * WVAL
      enddo
      hecMESH%mpc%mpc_const(i) = hecMESH%mpc%mpc_const(i) * WVAL
      ! mpc_const 実数(正体不明)
    enddo
    !$omp end do
    !$omp end parallel

  end subroutine hecmw_mpc_scale
!=================================================================================
  subroutine forward_increment_lagrange(cstep,ndof,mmat,hecMESH,fstrSOLID,infoCTChange,wkarray,uc)
    integer, intent(in)                    :: cstep
    integer, intent(in)                    :: ndof
    real(kind=kreal), intent(in)           :: mmat(:)
    type( hecmwst_local_mesh ), intent(in) :: hecmesh
    type(fstr_solid), intent(inout)        :: fstrSOLID
    type(fstr_info_contactchange)          :: infoCTChange
    real(kind=kreal), intent(out)          :: wkarray(:)
    real(kind=kreal), intent(out)          :: uc(:)
    integer :: i, j, k, m, grpid, slave, nn, iSS, sid, etype, iter
    real(kind=kreal) :: fdum, conv, dlambda, shapefunc(l_max_surface_node), lambda(3)
 
    call fstr_scan_contact_state_exp( cstep, hecmesh, fstrsolid, infoctchange )
    if( .not. infoctchange%active ) return
 
    uc = 0.0d0
 
    iter = 0
    do
      wkarray = 0.0d0
      do i=1,size(fstrsolid%contacts)
        do j= 1, size(fstrsolid%contacts(i)%slave)
          if( fstrsolid%contacts(i)%states(j)%state == contactfree ) cycle
          ! contactfreeなら以降のループ内処理はしない
          if( fstrsolid%contacts(i)%states(j)%distance>epsilon(1.d0) ) then
            fstrsolid%contacts(i)%states(j)%state = contactfree
            cycle
          endif
          if( iter==0 ) then
            fstrsolid%contacts(i)%states(j)%multiplier(:) =0.d0
            fstrsolid%contacts(i)%states(j)%wkdist =0.d0
            cycle
          endif
          slave = fstrsolid%contacts(i)%slave(j)
 
          sid = fstrsolid%contacts(i)%states(j)%surface
          nn = size( fstrsolid%contacts(i)%master(sid)%nodes )
          etype = fstrsolid%contacts(i)%master(sid)%etype
          call getshapefunc( etype, fstrsolid%contacts(i)%states(j)%lpos(:), shapefunc )
          wkarray( slave ) = -fstrsolid%contacts(i)%states(j)%multiplier(1)
          do k=1,nn
            iss = fstrsolid%contacts(i)%master(sid)%nodes(k)
            wkarray( iss ) = wkarray( iss ) + shapefunc(k) * fstrsolid%contacts(i)%states(j)%multiplier(1)
          enddo
        enddo
      enddo
 
      if(iter > 0)then
        do i=1,size(fstrsolid%contacts)
          do j= 1, size(fstrsolid%contacts(i)%slave)
            if( fstrsolid%contacts(i)%states(j)%state == contactfree ) cycle
            slave = fstrsolid%contacts(i)%slave(j)
            sid = fstrsolid%contacts(i)%states(j)%surface
            nn = size( fstrsolid%contacts(i)%master(sid)%nodes )
            etype = fstrsolid%contacts(i)%master(sid)%etype
            call getshapefunc( etype, fstrsolid%contacts(i)%states(j)%lpos(:), shapefunc )
            fstrsolid%contacts(i)%states(j)%wkdist = -wkarray( slave )/mmat( (slave-1)*ndof+1 )
            do k=1,nn
              iss = fstrsolid%contacts(i)%master(sid)%nodes(k)
              fstrsolid%contacts(i)%states(j)%wkdist = fstrsolid%contacts(i)%states(j)%wkdist  &
                   + shapefunc(k) * wkarray(iss) / mmat( (iss-1)*ndof+1 )
            enddo
          enddo
        enddo
      endif
 
      conv = 0.d0
      wkarray = 0.d0
      do i=1,size(fstrsolid%contacts)
        do j= 1, size(fstrsolid%contacts(i)%slave)
          if( fstrsolid%contacts(i)%states(j)%state == contactfree ) cycle
          slave = fstrsolid%contacts(i)%slave(j)
          sid = fstrsolid%contacts(i)%states(j)%surface
          nn = size( fstrsolid%contacts(i)%master(sid)%nodes )
          etype = fstrsolid%contacts(i)%master(sid)%etype
          call getshapefunc( etype, fstrsolid%contacts(i)%states(j)%lpos(:), shapefunc )
          fdum = 1.d0/mmat( (slave-1)*ndof+1 )
          do k=1,nn
            iss = fstrsolid%contacts(i)%master(sid)%nodes(k)
            fdum = fdum + shapefunc(k)*shapefunc(k)/mmat( (iss-1)*ndof+1 )
          enddo
          dlambda= (fstrsolid%contacts(i)%states(j)%distance-fstrsolid%contacts(i)%states(j)%wkdist) /fdum
          conv = conv + dlambda*dlambda;
          fstrsolid%contacts(i)%states(j)%multiplier(1) = fstrsolid%contacts(i)%states(j)%multiplier(1) + dlambda
          if( fstrsolid%contacts(i)%fcoeff>0.d0 ) then
            if( fstrsolid%contacts(i)%states(j)%state == contactslip ) then
              fstrsolid%contacts(i)%states(j)%multiplier(2) =             &
              fstrsolid%contacts(i)%fcoeff * fstrsolid%contacts(i)%states(j)%multiplier(1)
            else    ! stick
              !      fstrSOLID%contacts(i)%states(j)%multiplier(2) =
            endif
          endif
          lambda = fstrsolid%contacts(i)%states(j)%multiplier(1)* fstrsolid%contacts(i)%states(j)%direction
          wkarray((slave-1)*ndof+1:(slave-1)*ndof+3) = lambda(:)
          do k=1,nn
            iss = fstrsolid%contacts(i)%master(sid)%nodes(k)
            wkarray((iss-1)*ndof+1:(iss-1)*ndof+3) = wkarray((iss-1)*ndof+1:(iss-1)*ndof+3) -lambda(:)*shapefunc(k)
          enddo
        enddo
      enddo
      if( dsqrt(conv)<1.d-8 ) exit
      iter = iter+1
    enddo
 
    do i=1,hecmesh%n_node*ndof
      uc(i) = wkarray(i)/mmat(i)
    enddo
  end subroutine forward_increment_lagrange
!=================================================================================
  subroutine fstr_updatestate( hecMESH, fstrSOLID, tincr)
    use m_fstr
    use m_static_lib
    use m_elastoplastic
    use mcreep
    use mviscoelastic
    type(hecmwst_local_mesh) :: hecmesh
    type(fstr_solid) :: fstrSOLID
    real(kind=kreal) :: tincr
    integer(kind=kint) :: itype, is, iE, ic_type, icel, ngauss, i
 
    ! 温度の割当があるときだけ
    if( associated( fstrsolid%temperature ) ) then
      do i = 1, hecmesh%n_node
        fstrsolid%last_temp(i) = fstrsolid%temperature(i)
      end do
    endif
 
    do itype = 1, hecmesh%n_elem_type
      is = hecmesh%elem_type_index(itype-1) + 1
      ie = hecmesh%elem_type_index(itype  )
      ic_type= hecmesh%elem_type_item(itype)
      if( ic_type == 301 ) ic_type = 111
      if( hecmw_is_etype_link(ic_type) ) cycle
      if( hecmw_is_etype_patch(ic_type) ) cycle
 
      ngauss = numofquadpoints( ic_type )
      do icel = is, ie
        if( iselastoplastic( fstrsolid%elements(icel)%gausses(1)%pMaterial%mtype ) ) then
          do i = 1, ngauss
            call updateepstate( fstrsolid%elements(icel)%gausses(i) ) ! 弾塑性
          enddo
        elseif( fstrsolid%elements(icel)%gausses(1)%pMaterial%mtype == norton ) then
          if( tincr>0.d0 ) then
            do i = 1, ngauss
              call updateviscostate( fstrsolid%elements(icel)%gausses(i) ) ! 粘性
            enddo
          endif
        elseif( isviscoelastic( fstrsolid%elements(icel)%gausses(1)%pMaterial%mtype ) ) then
          if( tincr > 0.d0 ) then
            do i = 1, ngauss
              call updateviscoelasticstate( fstrsolid%elements(icel)%gausses(i) ) ! 粘弾性
            enddo
          endif
        endif
 
        do i = 1, ngauss
          fstrsolid%elements(icel)%gausses(i)%strain_bak = fstrsolid%elements(icel)%gausses(i)%strain
          fstrsolid%elements(icel)%gausses(i)%stress_bak = fstrsolid%elements(icel)%gausses(i)%stress
        enddo
      enddo
    enddo
  end subroutine fstr_updatestate
!=================================================================================