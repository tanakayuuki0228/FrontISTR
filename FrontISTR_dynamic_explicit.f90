!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module contains subroutines for nonlinear explicit dynamic analysis

module fstr_dynamic_nlexplicit
    use m_fstr
    use m_static_lib
    use m_dynamic_output
    use m_fstr_EIG_setMASS
    use m_dynamic_mat_ass_bc_ac
    use m_dynamic_mat_ass_bc
    use m_dynamic_mat_ass_bc_vl
    use m_dynamic_mat_ass_load
    use m_fstr_Update
    use m_fstr_Restart
    use m_dynamic_mat_ass_couple
    use m_fstr_rcap_io
    use mContact
  
  contains
  
    !C================================================================C
    !C-- subroutine  fstr_solve_LINEAR_DYNAMIC
    !C================================================================C
    subroutine fstr_solve_dynamic_nlexplicit(hecMESH,hecMAT,fstrSOLID,fstrEIG   &
        ,fstrDYN,fstrRESULT,fstrPARAM,infoCTChange &
        ,fstrCPL, restrt_step_num )
      implicit none
      type(hecmwST_local_mesh)             :: hecMESH !読み込んだメッシュデータを格納する構造体
      !---------------------------------------------------------------------------
      !構造体型 hecmwST_local_mesh  /FrontISTR/hecmw1/doc/0803_001d_hecmw_PC_cluster_201_io.pdf p.53
      !---------------------------------------------------------------------------
      type(hecmwST_matrix)                 :: hecMAT
      type(fstr_eigen)                     :: fstrEIG
      type(fstr_solid)                     :: fstrSOLID
      type(hecmwST_result_data)            :: fstrRESULT
      type(fstr_param)                     :: fstrPARAM
      type(fstr_dynamic)                   :: fstrDYN
      type(hecmwST_matrix_lagrange)        :: hecLagMAT !< type hecmwST_matrix_lagrange
      type(fstr_info_contactChange)        :: infoCTChange !< fstr_info_contactChange
      type(fstr_couple)                    :: fstrCPL !for COUPLE
      type(hecmwST_matrix), pointer :: hecMATmpc
      integer(kind=kint), allocatable :: mark(:)
      integer(kind=kint) :: nnod, ndof, nn, numnp
      integer(kind=kint) :: i, j, ids, ide, kk
      integer(kind=kint) :: kkk0, kkk1
      integer(kind=kint) :: ierror
      integer(kind=kint) :: iiii5, iexit
      integer(kind=kint) :: revocap_flag
      real(kind=kreal), allocatable :: prevB(:)
      real(kind=kreal) :: a1, a2, a3, b1, b2, b3, c1, c2
      real(kind=kreal) :: bsize, res
      real(kind=kreal) :: time_1, time_2
      integer(kind=kint) :: restrt_step_num
      real(kind=kreal), parameter :: PI = 3.14159265358979323846D0
  
      a1 = 0.0d0; a2 = 0.0d0; a3 = 0.0d0; b1 = 0.0d0; b2 = 0.0d0; b3 = 0.0d0
      c1 = 0.0d0; c2 = 0.0d0
  
      !メッシュデータの読み込み？
      call hecmw_mpc_mat_init_explicit(hecMESH, hecMAT, hecMATmpc)

      ! \home\frontISTR2\FrontISTR\hecmw1\src\solver\mpc\hecmw_mpc_prepost.f90

      !---------------------------------------------------------------------------
      !サブルーチン hecmw_mpc_mat_init_explicit !計算前の初期化
      !---------------------------------------------------------------------------
      ! subroutine hecmw_mpc_mat_init_explicit(hecMESH, hecMAT, hecMATmpc)
      !   implicit none
      !   type (hecmwST_local_mesh), intent(inout), target :: hecMESH
      !   type (hecmwST_matrix), intent(in), target :: hecMAT
      !   type (hecmwST_matrix), pointer :: hecMATmpc
      !   integer(kind=kint) :: totalmpc, MPC_METHOD

      !---------------------------------------------------------------------------
      !target属性，pointer属性
      !---------------------------------------------------------------------------
      !  形状無指定配列でも上下限を維持，
      !  手続きの定義において，仮引数の配列にpointer属性を付け，実引数となる配列にtarget属性を付ける
      !---------------------------------------------------------------------------

      !   totalmpc = hecMESH%mpc%n_mpc !拘束グループ数

      !---------------------------------------------------------------------------
      !type(派生型，構造体)
      !---------------------------------------------------------------------------
      !  複数の型の変数を同時に保持，保持する型は派生型でもOK
      !  成分は派生型変数名の後ろに，%，成分名を書く
      !---------------------------------------------------------------------------

      !   call hecmw_allreduce_I1 (hecMESH, totalmpc, hecmw_sum)

      !-----------------------------------------------------------------------------
      ! サブルーチン hecmw_allreduce_I1
      !-----------------------------------------------------------------------------
      ! 拘束グループ数'totalmpc'をコミュニケータ内で合計してtotalmpcに返している
      !-----------------------------------------------------------------------------
      !-----------------------------------------------------------------------------
      ! 引数 hecmw_sum
      !-----------------------------------------------------------------------------
      ! hecmw_allreduce内のmpi_allreduceで演算を'和'に指定するやつ
      ! どこかのモジュールか何かで定義されているのでしょう，たぶん...
      !-----------------------------------------------------------------------------

      !   if (totalmpc == 0) then !拘束グループ無しの場合
      !     hecMATmpc => hecMAT !hecMATmpcの指示先がhecMATと同じになる
      !     return
      !   endif
    
      !   call hecmw_mpc_scale(hecMESH) !hecmw_mpc_prepost.f90 line455 何が行われているかは不明
    
      !   ! Force MPC_METHOD=3
      !   MPC_METHOD = 3
      !   call hecmw_mat_set_mpc_method(hecMAT, MPC_METHOD) 
      
      !---------------------------------------------------------------------------
      ! サブルーチン hecmat_mat_sett_mpc_method
      !---------------------------------------------------------------------------
      ! hecMATのIarrayの13番目(idx_i_mpc_method)をMPC_METHODにする
      ! mpcの何かを設定？
      !---------------------------------------------------------------------------
    
      !   allocate(hecMATmpc) ! hecMATmpcの格納領域を割当
      !   call hecmw_mat_init(hecMATmpc) ! hecMATmpcを初期化
    
      ! hecMATの情報をhecMATmpcへ
      !   hecMATmpc%N = hecMAT%N !内点数 
      !   hecMATmpc%NP = hecMAT%NP !総節点数
      !   hecMATmpc%NDOF = hecMAT%NDOF !節点自由度数
      !   allocate(hecMATmpc%B(size(hecMAT%B))) !右辺ベクトルの割当
      !   allocate(hecMATmpc%X(size(hecMAT%X))) !解ベクトルの割当
      ! end subroutine hecmw_mpc_mat_init_explicit
      !---------------------------------------------------------------------------
  
      ! メッシュ情報引渡？
      hecMAT%NDOF=hecMESH%n_dof
      nnod=hecMESH%n_node
      ndof=hecMAT%NDOF
      nn=ndof*ndof
  
      if( fstrPARAM%fg_couple == 1) then ! fg_cople 連成解析のパラメータ？
        if( fstrPARAM%fg_couple_type==5 .or. &
            fstrPARAM%fg_couple_type==6 ) then
          allocate( prevB(hecMAT%NP*ndof)      ,stat=ierror ) ! prevBの割当，失敗するとierrorにエラー番号が入る
          prevB = 0.0d0
          if( ierror /= 0 ) then !割当が失敗
            write(idbg,*) 'stop due to allocation error <fstr_solve_NONLINEAR_DYNAMIC, prevB>'
            write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
            call flush(idbg) !同じデータに複数の処理がアクセスできるようにするやつ
            call hecmw_abort( hecmw_comm_get_comm()) !処理中断
          endif
        endif
      endif
  
      fstrSOLID%dunode(:) =0.d0
  
      a1 = 1.d0/fstrDYN%t_delta**2
      a2 = 1.d0/(2.d0*fstrDYN%t_delta)
  
      call setMASS(fstrSOLID,hecMESH,hecMAT,fstrEIG) !Set up lumped mass matrix
      call hecmw_mpc_trans_mass(hecMESH, hecMAT, fstrEIG%mass) ! mpc用に何かしている
  
      allocate(mark(hecMAT%NP * hecMAT%NDOF))
      call hecmw_mpc_mark_slave(hecMESH, hecMAT, mark) ! mpc用に何かしている
  
      do j = 1 ,ndof*nnod !節点自由度×節点数
        fstrDYN%VEC1(j) = (a1 + a2 *fstrDYN%ray_m) * fstrEIG%mass(j)
        ! VEC1 temporary quantity, ray_m rayleigh damping parameter Rm
        ! lumped質量行列の値
        if(mark(j) == 1) fstrDYN%VEC1(j) = 1.d0
        if(dabs(fstrDYN%VEC1(j)) < 1.0e-20) then ! dabs() 倍精度の絶対値
          !レーリー減衰 計算の安定性を確認？
          if( hecMESH%my_rank == 0 ) then
            write(*,*) 'stop due to fstrDYN%VEC(j) = 0 ,  j = ', j
            write(imsg,*) 'stop due to fstrDYN%VEC(j) = 0 ,  j = ', j
          end if
          call hecmw_abort( hecmw_comm_get_comm())
        endif
      end do
  
      deallocate(mark)
  
      !C-- output of initial state
      if( restrt_step_num == 1 ) then
        do j = 1 ,ndof*nnod
          fstrDYN%DISP(j,3) = fstrDYN%DISP(j,1) - fstrDYN%VEL (j,1)/(2.d0*a2)  + fstrDYN%ACC (j,1)/ (2.d0*a1)
          fstrDYN%DISP(j,2) = fstrDYN%DISP(j,1) - fstrDYN%VEL (j,1)/ a2 + fstrDYN%ACC (j,1)/ (2.d0*a1) * 4.d0
        end do
        ! DISP(j,3)がt-⊿t，DISP(j,2)がt-2⊿tの変位，2次までテーラー展開
  
        ! 初期値出力
        call fstr_dynamic_Output(hecMESH, fstrSOLID, fstrDYN, fstrPARAM)
        call dynamic_output_monit(hecMESH, fstrPARAM, fstrDYN, fstrEIG, fstrSOLID)
      end if
  
      ! 接触解析用
      if( associated( fstrSOLID%contacts ) )  then
        call initialize_contact_output_vectors(fstrSOLID,hecMAT)
        call forward_increment_Lagrange(1,ndof,fstrDYN%VEC1,hecMESH,fstrSOLID,infoCTChange,&
          & fstrDYN%DISP(:,2),fstrSOLID%ddunode)
      endif
  
      do i= restrt_step_num, fstrDYN%n_step
  
        fstrDYN%i_step = i
        fstrDYN%t_curr = fstrDYN%t_delta * i
  
        !C-- mechanical boundary condition ! 力学的境界条件←応力規定境界
        ! 右辺ベクトルに応力規定境界分を入れる
        call dynamic_mat_ass_load (hecMESH, hecMAT, fstrSOLID, fstrDYN, fstrPARAM)
        do j=1, hecMESH%n_node*  hecMESH%n_dof
          hecMAT%B(j)=hecMAT%B(j)-fstrSOLID%QFORCE(j)
        end do
  
        !C ********************************************************************************
        !C for couple analysis 連成解析
        ! TYPE =    1: 片方向連成          (FrontISTRはデータ受信から開始)
        !           2: 片方向連成          (FrontISTRはデータ送信から開始)
        !           3: Staggered双方向連成 (FrontISTRはデータ受信から開始)
        !           4: Staggered双方向連成 (FrontISTRはデータ送信から開始)
        !           5: 分離反復双方向連成  (FrontISTRはデータ受信から開始)
        !           6: 分離反復双方向連成  (FrontISTRはデータ送信から開始)
        if( fstrPARAM%fg_couple == 1 ) then
          if( fstrPARAM%fg_couple_type==5 .or. &
              fstrPARAM%fg_couple_type==6 ) then
            do j = 1, hecMAT%NP * ndof
              prevB(j) = hecMAT%B(j)
            enddo
          endif
        endif
        do
          if( fstrPARAM%fg_couple == 1 ) then
            ! REVOCAP_Coupler(連成解析ツール)を使う場合
            if( fstrPARAM%fg_couple_type==1 .or. &
              fstrPARAM%fg_couple_type==3 .or. &
              fstrPARAM%fg_couple_type==5 ) call fstr_rcap_get( fstrCPL )
            if( fstrPARAM%fg_couple_first /= 0 ) then
              bsize = dfloat( i ) / dfloat( fstrPARAM%fg_couple_first )
              if( bsize > 1.0 ) bsize = 1.0
              do kkk0 = 1, fstrCPL%coupled_node_n
                kkk1 = 3 * kkk0
                fstrCPL%trac(kkk1-2) = bsize * fstrCPL%trac(kkk1-2)
                fstrCPL%trac(kkk1-1) = bsize * fstrCPL%trac(kkk1-1)
                fstrCPL%trac(kkk1  ) = bsize * fstrCPL%trac(kkk1  )
              enddo
            endif
            if( fstrPARAM%fg_couple_window > 0 ) then
              j = i - restrt_step_num + 1
              kk = fstrDYN%n_step - restrt_step_num + 1
              bsize = 0.5*(1.0-cos(2.0*PI*dfloat(j)/dfloat(kk)))
              do kkk0 = 1, fstrCPL%coupled_node_n
                kkk1 = 3 * kkk0
                fstrCPL%trac(kkk1-2) = bsize * fstrCPL%trac(kkk1-2)
                fstrCPL%trac(kkk1-1) = bsize * fstrCPL%trac(kkk1-1)
                fstrCPL%trac(kkk1  ) = bsize * fstrCPL%trac(kkk1  )
              enddo
            endif
            ! 節点情報をhecMATに足し込む
            call dynamic_mat_ass_couple( hecMESH, hecMAT, fstrSOLID, fstrCPL )
          endif
          !C ********************************************************************************
  
          call hecmw_mpc_trans_rhs(hecMESH, hecMAT, hecMATmpc)
  
          ! a1 = 1.d0/fstrdyn%t_delta**2
          ! a2 = 1.d0/(2.d0*fstrdyn%t_delta)
          do j = 1 ,ndof*nnod
            hecMATmpc%B(j) = hecMATmpc%B(j) + 2.d0*a1* fstrEIG%mass(j) * fstrDYN%DISP(j,1)  &
              + (- a1 + a2 * fstrDYN%ray_m) * fstrEIG%mass(j) * fstrDYN%DISP(j,3)
          end do
  
          !C
          !C-- geometrical boundary condition ! 幾何学的境界条件←変位規定境界
          ! 右辺ベクトルに変位規定境界分を入れる
          ! ニューマークβ法を使ってる．γ=0だと計算停止?
          call dynamic_explicit_ass_bc(hecMESH, hecMATmpc, fstrSOLID, fstrDYN) 
          call dynamic_explicit_ass_vl(hecMESH, hecMATmpc, fstrSOLID, fstrDYN)
          call dynamic_explicit_ass_ac(hecMESH, hecMATmpc, fstrSOLID, fstrDYN)
          !call dynamic_mat_ass_bc   (hecMESH, hecMATmpc, fstrSOLID, fstrDYN, fstrPARAM, hecLagMAT)
          !call dynamic_mat_ass_bc_vl(hecMESH, hecMATmpc, fstrSOLID, fstrDYN, fstrPARAM, hecLagMAT)
          !call dynamic_mat_ass_bc_ac(hecMESH, hecMATmpc, fstrSOLID, fstrDYN, fstrPARAM, hecLagMAT)
  
          ! Finish the calculation
          do j = 1 ,ndof*nnod
            ! lumpedなので逆数をかけるだけ
            hecMATmpc%X(j) = hecMATmpc%B(j) / fstrDYN%VEC1(j)
            if(dabs(hecMATmpc%X(j)) > 1.0d+5) then
              if( hecMESH%my_rank == 0 ) then
                print *, 'Displacement increment too large, please adjust your step size!',i,hecMATmpc%X(j)
                write(imsg,*) 'Displacement increment too large, please adjust your step size!',i,hecMATmpc%B(j),fstrDYN%VEC1(j)
              end if
              call hecmw_abort( hecmw_comm_get_comm())
            end if
          end do
          call hecmw_mpc_tback_sol(hecMESH, hecMAT, hecMATmpc)
  
          !C *****************************************************
          !C for couple analysis
          if( fstrPARAM%fg_couple == 1 ) then
            if( fstrPARAM%fg_couple_type>1 ) then
              do j=1, fstrCPL%coupled_node_n
                if( fstrCPL%dof == 3 ) then
                  kkk0 = j*3
                  kkk1 = fstrCPL%coupled_node(j)*3
  
                  fstrCPL%disp (kkk0-2) = hecMAT%X(kkk1-2)
                  fstrCPL%disp (kkk0-1) = hecMAT%X(kkk1-1)
                  fstrCPL%disp (kkk0  ) = hecMAT%X(kkk1  )
  
                  fstrCPL%velo (kkk0-2) = -b1*fstrDYN%ACC(kkk1-2,1) - b2*fstrDYN%VEL(kkk1-2,1) + &
                    b3*( hecMAT%X(kkk1-2) - fstrDYN%DISP(kkk1-2,1) )
                  fstrCPL%velo (kkk0-1) = -b1*fstrDYN%ACC(kkk1-1,1) - b2*fstrDYN%VEL(kkk1-1,1) + &
                    b3*( hecMAT%X(kkk1-1) - fstrDYN%DISP(kkk1-1,1) )
                  fstrCPL%velo (kkk0  ) = -b1*fstrDYN%ACC(kkk1,1) - b2*fstrDYN%VEL(kkk1,1) + &
                    b3*( hecMAT%X(kkk1) - fstrDYN%DISP(kkk1,1) )
                  fstrCPL%accel(kkk0-2) = -a1*fstrDYN%ACC(kkk1-2,1) - a2*fstrDYN%VEL(kkk1-2,1) + &
                    a3*( hecMAT%X(kkk1-2) - fstrDYN%DISP(kkk1-2,1) )
                  fstrCPL%accel(kkk0-1) = -a1*fstrDYN%ACC(kkk1-1,1) - a2*fstrDYN%VEL(kkk1-1,1) + &
                    a3*( hecMAT%X(kkk1-1) - fstrDYN%DISP(kkk1-1,1) )
                  fstrCPL%accel(kkk0  ) = -a1*fstrDYN%ACC(kkk1,1) - a2*fstrDYN%VEL(kkk1,1) + &
                    a3*( hecMAT%X(kkk1) - fstrDYN%DISP(kkk1,1) )
                else
                  kkk0 = j*2
                  kkk1 = fstrCPL%coupled_node(j)*2
  
                  fstrCPL%disp (kkk0-1) = hecMAT%X(kkk1-1)
                  fstrCPL%disp (kkk0  ) = hecMAT%X(kkk1  )
  
                  fstrCPL%velo (kkk0-1) = -b1*fstrDYN%ACC(kkk1-1,1) - b2*fstrDYN%VEL(kkk1-1,1) + &
                    b3*( hecMAT%X(kkk1-1) - fstrDYN%DISP(kkk1-1,1) )
                  fstrCPL%velo (kkk0  ) = -b1*fstrDYN%ACC(kkk1,1) - b2*fstrDYN%VEL(kkk1,1) + &
                    b3*( hecMAT%X(kkk1) - fstrDYN%DISP(kkk1,1) )
                  fstrCPL%accel(kkk0-1) = -a1*fstrDYN%ACC(kkk1-1,1) - a2*fstrDYN%VEL(kkk1-1,1) + &
                    a3*( hecMAT%X(kkk1-1) - fstrDYN%DISP(kkk1-1,1) )
                  fstrCPL%accel(kkk0  ) = -a1*fstrDYN%ACC(kkk1,1) - a2*fstrDYN%VEL(kkk1,1) + &
                    a3*( hecMAT%X(kkk1) - fstrDYN%DISP(kkk1,1) )
                endif
              end do
              call fstr_rcap_send( fstrCPL )
            endif
  
            select case ( fstrPARAM%fg_couple_type )
              case (4)
                call fstr_rcap_get( fstrCPL )
              case (5)
                call fstr_get_convergence( revocap_flag )
                if( revocap_flag==0 ) then
                  do j = 1, hecMAT%NP * ndof
                    hecMAT%B(j) = prevB(j)
                  enddo
                  cycle
                endif
              case (6)
                call fstr_get_convergence( revocap_flag )
                if( revocap_flag==0 ) then
                  do j = 1, hecMAT%NP * ndof
                    hecMAT%B(j) = prevB(j)
                  enddo
                  call fstr_rcap_get( fstrCPL )
                  cycle
                else
                  if( i /= fstrDYN%n_step ) call fstr_rcap_get( fstrCPL )
                endif
            end select
          endif
          exit
        enddo
  
        !C *****************************************************
        !C-- contact corrector
        !C
        do j = 1 ,ndof*nnod
          fstrSOLID%unode(j)  = fstrDYN%DISP(j,1)
          fstrSOLID%dunode(j)  = hecMAT%X(j)-fstrDYN%DISP(j,1)
        enddo
        if( associated( fstrSOLID%contacts ) )  then
          !call fstr_scan_contact_state( 1, fstrDYN%t_delta, kcaSLAGRANGE, hecMESH, fstrSOLID, infoCTChange )
          call forward_increment_Lagrange(1,ndof,fstrDYN%VEC1,hecMESH,fstrSOLID,infoCTChange,&
            & fstrDYN%DISP(:,2),fstrSOLID%ddunode)
          do j = 1 ,ndof*nnod
            hecMAT%X(j)  = hecMAT%X(j) + fstrSOLID%ddunode(j)
          enddo
        endif
  
        !C-- new displacement, velocity and acceleration
        do j = 1 ,ndof*nnod
          fstrDYN%ACC (j,1) = a1*(hecMAT%X(j) - 2.d0*fstrDYN%DISP(j,1) + fstrDYN%DISP(j,3))
          fstrDYN%VEL (j,1) = a2*(hecMAT%X(j) - fstrDYN%DISP(j,3))
          fstrSOLID%unode(j)  = fstrDYN%DISP(j,1)
          fstrSOLID%dunode(j)  = hecMAT%X(j)-fstrDYN%DISP(j,1)
          fstrDYN%DISP(j,3) = fstrDYN%DISP(j,1)
          fstrDYN%DISP(j,1) = hecMAT%X(j)
          hecMAT%X(j)  = fstrSOLID%dunode(j)
        end do
  
        ! ----- update strain, stress, and internal force
        call fstr_UpdateNewton( hecMESH, hecMAT, fstrSOLID, fstrDYN%t_curr, fstrDYN%t_delta, 1 )
  
        do j = 1 ,ndof*nnod
          fstrSOLID%unode(j) = fstrSOLID%unode(j) + fstrSOLID%dunode(j)
        end do
        call fstr_UpdateState( hecMESH, fstrSOLID, fstrDYN%t_delta )
  
        if( fstrDYN%restart_nout > 0 ) then
          if ( mod(i,fstrDYN%restart_nout).eq.0 .or. i.eq.fstrDYN%n_step ) then
            call fstr_write_restart_dyna_nl(i,hecMESH,fstrSOLID,fstrDYN,fstrPARAM)
          end if
        end if
        !
        !C-- output new displacement, velocity and acceleration
        call fstr_dynamic_Output(hecMESH, fstrSOLID, fstrDYN, fstrPARAM)
        call dynamic_output_monit(hecMESH, fstrPARAM, fstrDYN, fstrEIG, fstrSOLID)
  
      enddo
  
      if( fstrPARAM%fg_couple == 1) then
        if( fstrPARAM%fg_couple_type==5 .or. &
            fstrPARAM%fg_couple_type==6 ) then
          deallocate( prevB      ,stat=ierror )
          if( ierror /= 0 ) then
            write(idbg,*) 'stop due to deallocation error <fstr_solve_NONLINEAR_DYNAMIC, prevB>'
            write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
            call flush(idbg)
            call hecmw_abort( hecmw_comm_get_comm())
          endif
        endif
      endif
  
      call hecmw_mpc_mat_finalize_explicit(hecMESH, hecMAT, hecMATmpc)
  
    end subroutine fstr_solve_dynamic_nlexplicit
  
    !< This subroutine implements Forward increment Lagrange multiplier method( NJ Carpenter et al. Int.J.Num.Meth.Eng.,32(1991),103-128 )
    subroutine forward_increment_Lagrange(cstep,ndof,mmat,hecMESH,fstrSOLID,infoCTChange,wkarray,uc)
      integer, intent(in)                    :: cstep
      integer, intent(in)                    :: ndof
      real(kind=kreal), intent(in)           :: mmat(:)
      type( hecmwST_local_mesh ), intent(in) :: hecMESH       !< type mesh
      type(fstr_solid), intent(inout)        :: fstrSOLID
      type(fstr_info_contactChange)          :: infoCTChange
      real(kind=kreal), intent(out)          :: wkarray(:)
      real(kind=kreal), intent(out)          :: uc(:)
      integer :: i, j, k, m, grpid, slave, nn, iSS, sid, etype, iter
      real(kind=kreal) :: fdum, conv, dlambda, shapefunc(l_max_surface_node), lambda(3)
  
      call fstr_scan_contact_state_exp( cstep, hecMESH, fstrSOLID, infoCTChange )
      if( .not. infoCTChange%active ) return
  
      uc = 0.0d0
  
      iter = 0
      do
        wkarray = 0.0d0
        do i=1,size(fstrSOLID%contacts)
          do j= 1, size(fstrSOLID%contacts(i)%slave)
            if( fstrSOLID%contacts(i)%states(j)%state == CONTACTFREE ) cycle
            if( fstrSOLID%contacts(i)%states(j)%distance>epsilon(1.d0) ) then
              fstrSOLID%contacts(i)%states(j)%state = CONTACTFREE
              cycle
            endif
            if( iter==0 ) then
              fstrSOLID%contacts(i)%states(j)%multiplier(:) =0.d0
              fstrSOLID%contacts(i)%states(j)%wkdist =0.d0
              cycle
            endif
            slave = fstrSOLID%contacts(i)%slave(j)
  
            sid = fstrSOLID%contacts(i)%states(j)%surface
            nn = size( fstrSOLID%contacts(i)%master(sid)%nodes )
            etype = fstrSOLID%contacts(i)%master(sid)%etype
            call getShapeFunc( etype, fstrSOLID%contacts(i)%states(j)%lpos(:), shapefunc )
            wkarray( slave ) = -fstrSOLID%contacts(i)%states(j)%multiplier(1)
            do k=1,nn
              iSS = fstrSOLID%contacts(i)%master(sid)%nodes(k)
              wkarray( iSS ) = wkarray( iSS ) + shapefunc(k) * fstrSOLID%contacts(i)%states(j)%multiplier(1)
            enddo
          enddo
        enddo
  
        if(iter > 0)then
          do i=1,size(fstrSOLID%contacts)
            do j= 1, size(fstrSOLID%contacts(i)%slave)
              if( fstrSOLID%contacts(i)%states(j)%state == CONTACTFREE ) cycle
              slave = fstrSOLID%contacts(i)%slave(j)
              sid = fstrSOLID%contacts(i)%states(j)%surface
              nn = size( fstrSOLID%contacts(i)%master(sid)%nodes )
              etype = fstrSOLID%contacts(i)%master(sid)%etype
              call getShapeFunc( etype, fstrSOLID%contacts(i)%states(j)%lpos(:), shapefunc )
              fstrSOLID%contacts(i)%states(j)%wkdist = -wkarray( slave )/mmat( (slave-1)*ndof+1 )
              do k=1,nn
                iSS = fstrSOLID%contacts(i)%master(sid)%nodes(k)
                fstrSOLID%contacts(i)%states(j)%wkdist = fstrSOLID%contacts(i)%states(j)%wkdist  &
                     + shapefunc(k) * wkarray(iSS) / mmat( (iSS-1)*ndof+1 )
              enddo
            enddo
          enddo
        endif
  
        conv = 0.d0
        wkarray = 0.d0
        do i=1,size(fstrSOLID%contacts)
          do j= 1, size(fstrSOLID%contacts(i)%slave)
            if( fstrSOLID%contacts(i)%states(j)%state == CONTACTFREE ) cycle
            slave = fstrSOLID%contacts(i)%slave(j)
            sid = fstrSOLID%contacts(i)%states(j)%surface
            nn = size( fstrSOLID%contacts(i)%master(sid)%nodes )
            etype = fstrSOLID%contacts(i)%master(sid)%etype
            call getShapeFunc( etype, fstrSOLID%contacts(i)%states(j)%lpos(:), shapefunc )
            fdum = 1.d0/mmat( (slave-1)*ndof+1 )
            do k=1,nn
              iSS = fstrSOLID%contacts(i)%master(sid)%nodes(k)
              fdum = fdum + shapefunc(k)*shapefunc(k)/mmat( (iSS-1)*ndof+1 )
            enddo
            dlambda= (fstrSOLID%contacts(i)%states(j)%distance-fstrSOLID%contacts(i)%states(j)%wkdist) /fdum
            conv = conv + dlambda*dlambda;
            fstrSOLID%contacts(i)%states(j)%multiplier(1) = fstrSOLID%contacts(i)%states(j)%multiplier(1) + dlambda
            if( fstrSOLID%contacts(i)%fcoeff>0.d0 ) then
              if( fstrSOLID%contacts(i)%states(j)%state == CONTACTSLIP ) then
                fstrSOLID%contacts(i)%states(j)%multiplier(2) =             &
                fstrSOLID%contacts(i)%fcoeff * fstrSOLID%contacts(i)%states(j)%multiplier(1)
              else    ! stick
                !      fstrSOLID%contacts(i)%states(j)%multiplier(2) =
              endif
            endif
            lambda = fstrSOLID%contacts(i)%states(j)%multiplier(1)* fstrSOLID%contacts(i)%states(j)%direction
            wkarray((slave-1)*ndof+1:(slave-1)*ndof+3) = lambda(:)
            do k=1,nn
              iSS = fstrSOLID%contacts(i)%master(sid)%nodes(k)
              wkarray((iSS-1)*ndof+1:(iSS-1)*ndof+3) = wkarray((iSS-1)*ndof+1:(iSS-1)*ndof+3) -lambda(:)*shapefunc(k)
            enddo
          enddo
        enddo
        if( dsqrt(conv)<1.d-8 ) exit
        iter = iter+1
      enddo
  
      do i=1,hecMESH%n_node*ndof
        uc(i) = wkarray(i)/mmat(i)
      enddo
    end subroutine forward_increment_Lagrange
  
  end module fstr_dynamic_nlexplicit
  
  !===============================================================================
  ! メモ
  !===============================================================================
  ! コンパイルする時↓
  !-------------------------------------------------------------------------------
  ! mpiifort frontistr_test.f90 -I/home/frontISTR2/FrontISTR/build/hecmw1 -L/home/frontISTR2/FrontISTR/build/hecmw1 -lhecmw -qmkl