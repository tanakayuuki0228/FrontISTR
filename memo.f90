!=================================================================================
MPI_ALLREDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, COMM, IERROR)
        <type> SENDBUF(*), RECVBUF(*)
        INTEGER COUNT, DATATYPE, OP, COMM, IERROR
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