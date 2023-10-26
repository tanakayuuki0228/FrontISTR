!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module contains subroutines for nonlinear explicit dynamic analysis

!////////////////////////////////////////////////////////////////////////////////////////////
! frountISTRの動解析用プログラム
! 本ファイルの中身は、以下のファイルをほぼコピペしたもの：
!   FronISTRインストールディレクトリ/FrontISTR/fistr1/src/analysis/dynamic/transit/fstr_dynamic_nlexplicit.f90
!
! 変更点(tanaka)：
! ・module名を fstr_test_dynamic_nlexplicit に変更
! ・subroutine phase_test_init を追加
! ・subroutine phase_test_output を追加
! ・real(kind=kreal), parameter Gc,L0 を追加
! ・subroutine phase_test_solve を追加
! ・use m_dynamic_test_output を追加
! ・use m_fstr_test_Update を追加
! ・連成解析,接触解析部分をコメントアウト
! ・hecMAT を hecMAT_DISP に変更, hecMAT_PHASE を追加 (方程式の係数行列,右辺ベクトルなどを持つ構造体)
! ・use m_test_phase, type fstrPHASE を追加 (phase,Gc,L0 を持つ構造体)
! ・subroutine phase_test_update_Kmatrix2 を追加
! ・subroutine phase_test_update_RHS2 を追加
!////////////////////////////////////////////////////////////////////////////////////////////

module fstr_test_dynamic_nlexplicit ! 2023/08/22(tanaka)
  use m_fstr
  use m_static_lib
  ! use m_dynamic_test_output ! 2023/08/24(tanaka)xx
  use m_fstr_EIG_setMASS
  use m_dynamic_mat_ass_bc_ac
  use m_dynamic_mat_ass_bc
  use m_dynamic_mat_ass_bc_vl
  use m_dynamic_mat_ass_load
  use m_fstr_test_Update ! 2023/08/24(tanaka)
  use m_fstr_Restart
  use m_dynamic_mat_ass_couple
  use m_fstr_rcap_io
  use mContact
  use m_test_phase ! 2023/09/13(tanaka)
  use m_fstr_test_Q_Update ! 2023/10/20(tanaka)
  use m_test_fstr_NodalStress ! 2023/10/25

contains

  !C================================================================C
  !C-- subroutine  fstr_solve_LINEAR_DYNAMIC
  !C================================================================C
  subroutine fstr_solve_dynamic_nlexplicit(hecMESH,hecMAT_DISP,hecMAT_PHASE,fstrSOLID,fstrHEAT,fstrEIG   &
      ,fstrDYN,fstrRESULT,fstrPARAM,infoCTChange &
      ,fstrCPL, restrt_step_num ) ! 2023/09/12(tanaka)

    implicit none

    type(hecmwST_local_mesh)             :: hecMESH
    type(hecmwST_matrix)                 :: hecMAT_DISP,hecMAT_PHASE ! 2023/09/12(tanaka)
    type(fstr_eigen)                     :: fstrEIG
    type(fstr_solid)                     :: fstrSOLID
    type(fstr_heat)                      :: fstrHEAT ! 2023/09/19(tanaka)
    type(hecmwST_result_data)            :: fstrRESULT
    type(fstr_param)                     :: fstrPARAM
    type(fstr_dynamic)                   :: fstrDYN
    type(hecmwST_matrix_lagrange)        :: hecLagMAT !< type hecmwST_matrix_lagrange
    type(fstr_info_contactChange)        :: infoCTChange !< fstr_info_contactChange
    type(fstr_couple)                    :: fstrCPL !for COUPLE
    type(hecmwST_matrix), pointer :: hecMATmpc
    type(fstr_phase) :: fstrPHASE ! 2023/09/12(tanaka)
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
    integer(kind=kint) :: res_output_interval ! エネルギー以外の結果出力間隔 2023/10/12(tanaka)
    integer(kind=kint) :: ene_output_interval ! エネルギーの出力間隔
    logical :: ene_output_flag ! 出力エネルギー計算のフラグ ene_output_flag = mod(istep,ene_output_interval)==0 とする


    a1 = 0.0d0; a2 = 0.0d0; a3 = 0.0d0; b1 = 0.0d0; b2 = 0.0d0; b3 = 0.0d0
    c1 = 0.0d0; c2 = 0.0d0

    call hecmw_mpc_mat_init_explicit(hecMESH, hecMAT_DISP, hecMATmpc) ! 2023/09/12(tanaka)

    hecMAT_DISP%NDOF=hecMESH%n_dof ! 2023/09/12(tanaka)
    nnod=hecMESH%n_node
    ndof=hecMAT_DISP%NDOF ! 2023/09/12(tanaka)
    nn=ndof*ndof

    ! write(*,*) hecMESH%n_node ! デバッグ用 2023/08/22(tanaka)
    ! write(*,*) hecMAT%NP ! デバッグ用 2023/08/22(tanaka)
    ! write(*,*) hecMAT%N ! デバッグ用 2023/08/22(tanaka)

    ! コメントアウト ! 2023/09/12(tanaka)
    ! if( fstrPARAM%fg_couple == 1) then
    !   if( fstrPARAM%fg_couple_type==5 .or. &
    !       fstrPARAM%fg_couple_type==6 ) then
    !     allocate( prevB(hecMAT_DISP%NP*ndof)      ,stat=ierror ) ! 2023/09/12(tanaka)
    !     prevB = 0.0d0
    !     if( ierror /= 0 ) then
    !       write(idbg,*) 'stop due to allocation error <fstr_solve_NONLINEAR_DYNAMIC, prevB>'
    !       write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
    !       call flush(idbg)
    !       call hecmw_abort( hecmw_comm_get_comm())
    !     endif
    !   endif
    ! endif

    fstrSOLID%dunode(:) = 0.d0

    a1 = 1.d0/fstrDYN%t_delta**2        ! 1/(Δt^2)
    a2 = 1.d0/( 2.d0*fstrDYN%t_delta )  ! 1/(2Δt)

    call setMASS(fstrSOLID,hecMESH,hecMAT_DISP,fstrEIG) ! 2023/09/12(tanaka)
    call hecmw_mpc_trans_mass(hecMESH, hecMAT_DISP, fstrEIG%mass) ! 2023/09/12(tanaka)

    ! phase field方程式関連の初期設定 (配列のallocate・初期値・規定境界条件)
    call phase_test_init(hecMESH,fstrPHASE,fstrHEAT,fstrDYN%n_step) ! 2023/08/22(tanaka) ! 2023/09/12(tanaka)

    ! phase-fieldの方程式を解く ! 2023/09/19(tanaka)
    ene_output_flag = .true. ! 初期エネルギーは出力する
    call phase_test_solve(hecMESH, fstrSOLID, fstrDYN, fstrHEAT, hecMAT_PHASE, fstrPHASE , ene_output_flag)

    ! phase-fieldを解で更新する ! 2023/09/19(tanaka)
    do j = 1 ,nnod
      fstrPHASE%phase(j) = hecMAT_PHASE%X(j)
    end do

    allocate(mark(hecMAT_DISP%NP * hecMAT_DISP%NDOF)) ! 2023/09/12(tanaka)
    call hecmw_mpc_mark_slave(hecMESH, hecMAT_DISP, mark) ! 2023/09/12(tanaka)

    do j = 1 ,ndof*nnod
      fstrDYN%VEC1(j) = (a1 + a2 *fstrDYN%ray_m) * fstrEIG%mass(j)
      if(mark(j) == 1) fstrDYN%VEC1(j) = 1.d0
      if(dabs(fstrDYN%VEC1(j)) < 1.0e-20) then
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
        !/////////////////////////////////////////////////////////
        ! Taylor展開: (a2:=1/(2Δt), a1:=1/(Δt^2))
        !   U(t- Δt) = U(t) - V(t)Δt    + (1/2)A(t)( Δt)^2 + ...
        !   U(t-2Δt) = U(t) - V(t)(2Δt) + (1/2)A(t)(2Δt)^2 + ... 
        fstrDYN%DISP(j,3) = fstrDYN%DISP(j,1) - fstrDYN%VEL(j,1) / (2.d0*a2)  + fstrDYN%ACC(j,1) / (2.d0*a1)
        fstrDYN%DISP(j,2) = fstrDYN%DISP(j,1) - fstrDYN%VEL(j,1) / a2         + fstrDYN%ACC(j,1) / (2.d0*a1) * 4.d0
      end do

      ! call fstr_dynamic_Output(hecMESH, fstrSOLID, fstrDYN, fstrPARAM)
      ! call dynamic_output_monit(hecMESH, fstrPARAM, fstrDYN, fstrEIG, fstrSOLID)
      call phase_test_output(hecMESH, fstrSOLID, fstrDYN, fstrPHASE, 0, res_output_interval, ene_output_interval) ! 2023/08/22(tanaka)
    end if

    ! コメントアウト ! 2023/09/12(tanaka)
    ! if( associated( fstrSOLID%contacts ) )  then
    !   call initialize_contact_output_vectors(fstrSOLID,hecMAT_DISP) ! 2023/09/12(tanaka)
    !   call forward_increment_Lagrange(1,ndof,fstrDYN%VEC1,hecMESH,fstrSOLID,infoCTChange,&
    !     & fstrDYN%DISP(:,2),fstrSOLID%ddunode)
    ! endif

    do i= restrt_step_num, fstrDYN%n_step
      fstrDYN%i_step = i
      fstrDYN%t_curr = fstrDYN%t_delta * i

      !C-- mechanical boundary condition
      call dynamic_mat_ass_load (hecMESH, hecMAT_DISP, fstrSOLID, fstrDYN, fstrPARAM) ! 2023/09/12(tanaka)
      do j=1, hecMESH%n_node*  hecMESH%n_dof
        hecMAT_DISP%B(j)=hecMAT_DISP%B(j)-fstrSOLID%QFORCE(j) ! 2023/09/12(tanaka)
      end do
      call hecmw_update_R(hecMESH,hecMAT_DISP%B,hecMESH%n_node, hecMESH%n_dof) ! 2023/10/03(tanaka)
      !C ********************************************************************************
      ! コメントアウト ! 2023/09/12(tanaka)
      ! !C for couple analysis
      ! if( fstrPARAM%fg_couple == 1 ) then
      !   if( fstrPARAM%fg_couple_type==5 .or. &
      !       fstrPARAM%fg_couple_type==6 ) then
      !     do j = 1, hecMAT_DISP%NP * ndof ! 2023/09/12(tanaka)
      !       prevB(j) = hecMAT_DISP%B(j) ! 2023/09/12(tanaka)
      !     enddo
      !   endif
      ! endif
      do
        ! コメントアウト ! 2023/09/12(tanaka)
        ! if( fstrPARAM%fg_couple == 1 ) then
        !   if( fstrPARAM%fg_couple_type==1 .or. &
        !     fstrPARAM%fg_couple_type==3 .or. &
        !     fstrPARAM%fg_couple_type==5 ) call fstr_rcap_get( fstrCPL )
        !   if( fstrPARAM%fg_couple_first /= 0 ) then
        !     bsize = dfloat( i ) / dfloat( fstrPARAM%fg_couple_first )
        !     if( bsize > 1.0 ) bsize = 1.0
        !     do kkk0 = 1, fstrCPL%coupled_node_n
        !       kkk1 = 3 * kkk0
        !       fstrCPL%trac(kkk1-2) = bsize * fstrCPL%trac(kkk1-2)
        !       fstrCPL%trac(kkk1-1) = bsize * fstrCPL%trac(kkk1-1)
        !       fstrCPL%trac(kkk1  ) = bsize * fstrCPL%trac(kkk1  )
        !     enddo
        !   endif
        !   if( fstrPARAM%fg_couple_window > 0 ) then
        !     j = i - restrt_step_num + 1
        !     kk = fstrDYN%n_step - restrt_step_num + 1
        !     bsize = 0.5*(1.0-cos(2.0*PI*dfloat(j)/dfloat(kk)))
        !     do kkk0 = 1, fstrCPL%coupled_node_n
        !       kkk1 = 3 * kkk0
        !       fstrCPL%trac(kkk1-2) = bsize * fstrCPL%trac(kkk1-2)
        !       fstrCPL%trac(kkk1-1) = bsize * fstrCPL%trac(kkk1-1)
        !       fstrCPL%trac(kkk1  ) = bsize * fstrCPL%trac(kkk1  )
        !     enddo
        !   endif
        !   call dynamic_mat_ass_couple( hecMESH, hecMAT_DISP, fstrSOLID, fstrCPL ) ! 2023/09/12(tanaka)
        ! endif
        !C ********************************************************************************

        call hecmw_mpc_trans_rhs(hecMESH, hecMAT_DISP, hecMATmpc) ! 2023/09/12(tanaka)

        do j = 1 ,ndof*nnod
          hecMATmpc%B(j) = hecMATmpc%B(j) + 2.d0*a1* fstrEIG%mass(j) * fstrDYN%DISP(j,1)  &
            + (- a1 + a2 * fstrDYN%ray_m) * fstrEIG%mass(j) * fstrDYN%DISP(j,3)
        end do

        !C
        !C-- geometrical boundary condition

        call dynamic_explicit_ass_bc(hecMESH, hecMATmpc, fstrSOLID, fstrDYN)
        call dynamic_explicit_ass_vl(hecMESH, hecMATmpc, fstrSOLID, fstrDYN)
        call dynamic_explicit_ass_ac(hecMESH, hecMATmpc, fstrSOLID, fstrDYN)
        !call dynamic_mat_ass_bc   (hecMESH, hecMATmpc, fstrSOLID, fstrDYN, fstrPARAM, hecLagMAT)
        !call dynamic_mat_ass_bc_vl(hecMESH, hecMATmpc, fstrSOLID, fstrDYN, fstrPARAM, hecLagMAT)
        !call dynamic_mat_ass_bc_ac(hecMESH, hecMATmpc, fstrSOLID, fstrDYN, fstrPARAM, hecLagMAT)

        ! Finish the calculation
        do j = 1 ,ndof*nnod
          hecMATmpc%X(j) = hecMATmpc%B(j) / fstrDYN%VEC1(j)
          if(dabs(hecMATmpc%X(j)) > 1.0d+5) then
            if( hecMESH%my_rank == 0 ) then
              print *, 'Displacement increment too large, please adjust your step size!',i,hecMATmpc%X(j)
              write(imsg,*) 'Displacement increment too large, please adjust your step size!',i,hecMATmpc%B(j),fstrDYN%VEC1(j)
            end if
            call hecmw_abort( hecmw_comm_get_comm())
          end if
        end do
        call hecmw_mpc_tback_sol(hecMESH, hecMAT_DISP, hecMATmpc) ! 2023/09/12(tanaka)

        !C *****************************************************
        ! コメントアウト ! 2023/09/12(tanaka)
        ! !C for couple analysis
        ! if( fstrPARAM%fg_couple == 1 ) then
        !   if( fstrPARAM%fg_couple_type>1 ) then
        !     do j=1, fstrCPL%coupled_node_n
        !       if( fstrCPL%dof == 3 ) then
        !         kkk0 = j*3
        !         kkk1 = fstrCPL%coupled_node(j)*3

        !         fstrCPL%disp (kkk0-2) = hecMAT_DISP%X(kkk1-2) ! 2023/09/12(tanaka)
        !         fstrCPL%disp (kkk0-1) = hecMAT_DISP%X(kkk1-1) ! 2023/09/12(tanaka)
        !         fstrCPL%disp (kkk0  ) = hecMAT_DISP%X(kkk1  ) ! 2023/09/12(tanaka)

        !         fstrCPL%velo (kkk0-2) = -b1*fstrDYN%ACC(kkk1-2,1) - b2*fstrDYN%VEL(kkk1-2,1) + &
        !           b3*( hecMAT_DISP%X(kkk1-2) - fstrDYN%DISP(kkk1-2,1) ) ! 2023/09/12(tanaka)
        !         fstrCPL%velo (kkk0-1) = -b1*fstrDYN%ACC(kkk1-1,1) - b2*fstrDYN%VEL(kkk1-1,1) + &
        !           b3*( hecMAT_DISP%X(kkk1-1) - fstrDYN%DISP(kkk1-1,1) ) ! 2023/09/12(tanaka)
        !         fstrCPL%velo (kkk0  ) = -b1*fstrDYN%ACC(kkk1,1) - b2*fstrDYN%VEL(kkk1,1) + &
        !           b3*( hecMAT_DISP%X(kkk1) - fstrDYN%DISP(kkk1,1) ) ! 2023/09/12(tanaka)
        !         fstrCPL%accel(kkk0-2) = -a1*fstrDYN%ACC(kkk1-2,1) - a2*fstrDYN%VEL(kkk1-2,1) + &
        !           a3*( hecMAT_DISP%X(kkk1-2) - fstrDYN%DISP(kkk1-2,1) ) ! 2023/09/12(tanaka)
        !         fstrCPL%accel(kkk0-1) = -a1*fstrDYN%ACC(kkk1-1,1) - a2*fstrDYN%VEL(kkk1-1,1) + &
        !           a3*( hecMAT_DISP%X(kkk1-1) - fstrDYN%DISP(kkk1-1,1) ) ! 2023/09/12(tanaka)
        !         fstrCPL%accel(kkk0  ) = -a1*fstrDYN%ACC(kkk1,1) - a2*fstrDYN%VEL(kkk1,1) + &
        !           a3*( hecMAT_DISP%X(kkk1) - fstrDYN%DISP(kkk1,1) ) ! 2023/09/12(tanaka)
        !       else
        !         kkk0 = j*2
        !         kkk1 = fstrCPL%coupled_node(j)*2

        !         fstrCPL%disp (kkk0-1) = hecMAT_DISP%X(kkk1-1) ! 2023/09/12(tanaka)
        !         fstrCPL%disp (kkk0  ) = hecMAT_DISP%X(kkk1  ) ! 2023/09/12(tanaka)

        !         fstrCPL%velo (kkk0-1) = -b1*fstrDYN%ACC(kkk1-1,1) - b2*fstrDYN%VEL(kkk1-1,1) + &
        !           b3*( hecMAT_DISP%X(kkk1-1) - fstrDYN%DISP(kkk1-1,1) ) ! 2023/09/12(tanaka)
        !         fstrCPL%velo (kkk0  ) = -b1*fstrDYN%ACC(kkk1,1) - b2*fstrDYN%VEL(kkk1,1) + &
        !           b3*( hecMAT_DISP%X(kkk1) - fstrDYN%DISP(kkk1,1) ) ! 2023/09/12(tanaka)
        !         fstrCPL%accel(kkk0-1) = -a1*fstrDYN%ACC(kkk1-1,1) - a2*fstrDYN%VEL(kkk1-1,1) + &
        !           a3*( hecMAT_DISP%X(kkk1-1) - fstrDYN%DISP(kkk1-1,1) ) ! 2023/09/12(tanaka)
        !         fstrCPL%accel(kkk0  ) = -a1*fstrDYN%ACC(kkk1,1) - a2*fstrDYN%VEL(kkk1,1) + &
        !           a3*( hecMAT_DISP%X(kkk1) - fstrDYN%DISP(kkk1,1) ) ! 2023/09/12(tanaka)
        !       endif
        !     end do
        !     call fstr_rcap_send( fstrCPL )
        !   endif

        !   select case ( fstrPARAM%fg_couple_type )
        !     case (4)
        !       call fstr_rcap_get( fstrCPL )
        !     case (5)
        !       call fstr_get_convergence( revocap_flag )
        !       if( revocap_flag==0 ) then
        !         do j = 1, hecMAT_DISP%NP * ndof ! 2023/09/12(tanaka)
        !           hecMAT_DISP%B(j) = prevB(j) ! 2023/09/12(tanaka)
        !         enddo
        !         cycle
        !       endif
        !     case (6)
        !       call fstr_get_convergence( revocap_flag )
        !       if( revocap_flag==0 ) then
        !         do j = 1, hecMAT_DISP%NP * ndof ! 2023/09/12(tanaka)
        !           hecMAT_DISP%B(j) = prevB(j) ! 2023/09/12(tanaka)
        !         enddo
        !         call fstr_rcap_get( fstrCPL )
        !         cycle
        !       else
        !         if( i /= fstrDYN%n_step ) call fstr_rcap_get( fstrCPL )
        !       endif
        !   end select
        ! endif
        exit
      enddo

      !C *****************************************************
      !C-- contact corrector
      !C
      do j = 1 ,ndof*nnod
        fstrSOLID%unode(j)  = fstrDYN%DISP(j,1)
        fstrSOLID%dunode(j)  = hecMAT_DISP%X(j)-fstrDYN%DISP(j,1) ! 2023/09/12(tanaka)
      enddo
      ! コメントアウト ! 2023/09/12(tanaka)
      ! if( associated( fstrSOLID%contacts ) )  then
      !   !call fstr_scan_contact_state( 1, fstrDYN%t_delta, kcaSLAGRANGE, hecMESH, fstrSOLID, infoCTChange )
      !   call forward_increment_Lagrange(1,ndof,fstrDYN%VEC1,hecMESH,fstrSOLID,infoCTChange,&
      !     & fstrDYN%DISP(:,2),fstrSOLID%ddunode)
      !   do j = 1 ,ndof*nnod
      !     hecMAT_DISP%X(j)  = hecMAT_DISP%X(j) + fstrSOLID%ddunode(j) ! 2023/09/12(tanaka)
      !   enddo
      ! endif

      !C-- new displacement, velocity and acceleration
      do j = 1 ,ndof*nnod
        fstrDYN%ACC (j,1) = a1*(hecMAT_DISP%X(j) - 2.d0*fstrDYN%DISP(j,1) + fstrDYN%DISP(j,3)) ! 2023/09/12(tanaka)
        fstrDYN%VEL (j,1) = a2*(hecMAT_DISP%X(j) - fstrDYN%DISP(j,3)) ! 2023/09/12(tanaka)
        fstrSOLID%unode(j)  = fstrDYN%DISP(j,1)
        fstrSOLID%dunode(j)  = hecMAT_DISP%X(j)-fstrDYN%DISP(j,1) ! 2023/09/12(tanaka)
        fstrDYN%DISP(j,3) = fstrDYN%DISP(j,1)
        fstrDYN%DISP(j,1) = hecMAT_DISP%X(j) ! 2023/09/12(tanaka)
        hecMAT_DISP%X(j)  = fstrSOLID%dunode(j) ! 2023/09/12(tanaka)
      end do

      ! ----- update strain, stress, and internal force
      ! 積分点の応力,ひずみを更新
      call fstr_UpdateStressStrain( hecMESH, hecMAT_DISP, fstrSOLID, fstrDYN%t_curr, fstrDYN%t_delta, 1, fstrPHASE )

      do j = 1 ,ndof*nnod
        fstrSOLID%unode(j) = fstrSOLID%unode(j) + fstrSOLID%dunode(j)
      end do
      call fstr_UpdateState( hecMESH, fstrSOLID, fstrDYN%t_delta )

      ! ! 節点応力,節点ひずみを更新
      ! call fstr_dynamic_Output(hecMESH, fstrSOLID, fstrDYN, fstrPARAM)

      ! phase-fieldの方程式を解く ! 2023/09/19(tanaka)
      ene_output_flag = mod(i,ene_output_interval) == 0
      call phase_test_solve(hecMESH, fstrSOLID, fstrDYN, fstrHEAT, hecMAT_PHASE, fstrPHASE, ene_output_flag)

      ! phase-fieldを解で更新する ! 2023/09/19(tanaka)
      do j = 1 ,nnod
        fstrPHASE%phase(j) = hecMAT_PHASE%X(j)
      end do

      ! 内力ベクトルQFORCEを更新
      call fstr_UpdateQFORCE( hecMESH, hecMAT_DISP, fstrSOLID, fstrDYN%t_curr, fstrDYN%t_delta, 1 , fstrPHASE )

      if( fstrDYN%restart_nout > 0 ) then
        if ( mod(i,fstrDYN%restart_nout).eq.0 .or. i.eq.fstrDYN%n_step ) then
          call fstr_write_restart_dyna_nl(i,hecMESH,fstrSOLID,fstrDYN,fstrPARAM)
        end if
      end if
      !
      !C-- output new displacement, velocity and acceleration
      ! call fstr_dynamic_Output(hecMESH, fstrSOLID, fstrDYN, fstrPARAM)
      ! call dynamic_output_monit(hecMESH, fstrPARAM, fstrDYN, fstrEIG, fstrSOLID)
      if(mod(i,res_output_interval)==0) then
        call phase_test_output(hecMESH, fstrSOLID, fstrDYN, fstrPHASE, i, res_output_interval, ene_output_interval) ! 2023/08/22(tanaka)
      endif
    enddo

    ! コメントアウト ! 2023/09/12(tanaka)
    ! if( fstrPARAM%fg_couple == 1) then
    !   if( fstrPARAM%fg_couple_type==5 .or. &
    !       fstrPARAM%fg_couple_type==6 ) then
    !     deallocate( prevB      ,stat=ierror )
    !     if( ierror /= 0 ) then
    !       write(idbg,*) 'stop due to deallocation error <fstr_solve_NONLINEAR_DYNAMIC, prevB>'
    !       write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
    !       call flush(idbg)
    !       call hecmw_abort( hecmw_comm_get_comm())
    !     endif
    !   endif
    ! endif

    call hecmw_mpc_mat_finalize_explicit(hecMESH, hecMAT_DISP, hecMATmpc) ! 2023/09/12(tanaka)

  end subroutine fstr_solve_dynamic_nlexplicit

  ! !< This subroutine implements Forward increment Lagrange multiplier method( NJ Carpenter et al. Int.J.Num.Meth.Eng.,32(1991),103-128 )
  ! subroutine forward_increment_Lagrange(cstep,ndof,mmat,hecMESH,fstrSOLID,infoCTChange,wkarray,uc)
  !   integer, intent(in)                    :: cstep
  !   integer, intent(in)                    :: ndof
  !   real(kind=kreal), intent(in)           :: mmat(:)
  !   type( hecmwST_local_mesh ), intent(in) :: hecMESH       !< type mesh
  !   type(fstr_solid), intent(inout)        :: fstrSOLID
  !   type(fstr_info_contactChange)          :: infoCTChange
  !   real(kind=kreal), intent(out)          :: wkarray(:)
  !   real(kind=kreal), intent(out)          :: uc(:)
  !   integer :: i, j, k, m, grpid, slave, nn, iSS, sid, etype, iter
  !   real(kind=kreal) :: fdum, conv, dlambda, shapefunc(l_max_surface_node), lambda(3)

  !   call fstr_scan_contact_state_exp( cstep, hecMESH, fstrSOLID, infoCTChange )
  !   if( .not. infoCTChange%active ) return

  !   uc = 0.0d0

  !   iter = 0
  !   do
  !     wkarray = 0.0d0
  !     do i=1,size(fstrSOLID%contacts)
  !       do j= 1, size(fstrSOLID%contacts(i)%slave)
  !         if( fstrSOLID%contacts(i)%states(j)%state == CONTACTFREE ) cycle
  !         if( fstrSOLID%contacts(i)%states(j)%distance>epsilon(1.d0) ) then
  !           fstrSOLID%contacts(i)%states(j)%state = CONTACTFREE
  !           cycle
  !         endif
  !         if( iter==0 ) then
  !           fstrSOLID%contacts(i)%states(j)%multiplier(:) =0.d0
  !           fstrSOLID%contacts(i)%states(j)%wkdist =0.d0
  !           cycle
  !         endif
  !         slave = fstrSOLID%contacts(i)%slave(j)

  !         sid = fstrSOLID%contacts(i)%states(j)%surface
  !         nn = size( fstrSOLID%contacts(i)%master(sid)%nodes )
  !         etype = fstrSOLID%contacts(i)%master(sid)%etype
  !         call getShapeFunc( etype, fstrSOLID%contacts(i)%states(j)%lpos(:), shapefunc )
  !         wkarray( slave ) = -fstrSOLID%contacts(i)%states(j)%multiplier(1)
  !         do k=1,nn
  !           iSS = fstrSOLID%contacts(i)%master(sid)%nodes(k)
  !           wkarray( iSS ) = wkarray( iSS ) + shapefunc(k) * fstrSOLID%contacts(i)%states(j)%multiplier(1)
  !         enddo
  !       enddo
  !     enddo

  !     if(iter > 0)then
  !       do i=1,size(fstrSOLID%contacts)
  !         do j= 1, size(fstrSOLID%contacts(i)%slave)
  !           if( fstrSOLID%contacts(i)%states(j)%state == CONTACTFREE ) cycle
  !           slave = fstrSOLID%contacts(i)%slave(j)
  !           sid = fstrSOLID%contacts(i)%states(j)%surface
  !           nn = size( fstrSOLID%contacts(i)%master(sid)%nodes )
  !           etype = fstrSOLID%contacts(i)%master(sid)%etype
  !           call getShapeFunc( etype, fstrSOLID%contacts(i)%states(j)%lpos(:), shapefunc )
  !           fstrSOLID%contacts(i)%states(j)%wkdist = -wkarray( slave )/mmat( (slave-1)*ndof+1 )
  !           do k=1,nn
  !             iSS = fstrSOLID%contacts(i)%master(sid)%nodes(k)
  !             fstrSOLID%contacts(i)%states(j)%wkdist = fstrSOLID%contacts(i)%states(j)%wkdist  &
  !                  + shapefunc(k) * wkarray(iSS) / mmat( (iSS-1)*ndof+1 )
  !           enddo
  !         enddo
  !       enddo
  !     endif

  !     conv = 0.d0
  !     wkarray = 0.d0
  !     do i=1,size(fstrSOLID%contacts)
  !       do j= 1, size(fstrSOLID%contacts(i)%slave)
  !         if( fstrSOLID%contacts(i)%states(j)%state == CONTACTFREE ) cycle
  !         slave = fstrSOLID%contacts(i)%slave(j)
  !         sid = fstrSOLID%contacts(i)%states(j)%surface
  !         nn = size( fstrSOLID%contacts(i)%master(sid)%nodes )
  !         etype = fstrSOLID%contacts(i)%master(sid)%etype
  !         call getShapeFunc( etype, fstrSOLID%contacts(i)%states(j)%lpos(:), shapefunc )
  !         fdum = 1.d0/mmat( (slave-1)*ndof+1 )
  !         do k=1,nn
  !           iSS = fstrSOLID%contacts(i)%master(sid)%nodes(k)
  !           fdum = fdum + shapefunc(k)*shapefunc(k)/mmat( (iSS-1)*ndof+1 )
  !         enddo
  !         dlambda= (fstrSOLID%contacts(i)%states(j)%distance-fstrSOLID%contacts(i)%states(j)%wkdist) /fdum
  !         conv = conv + dlambda*dlambda;
  !         fstrSOLID%contacts(i)%states(j)%multiplier(1) = fstrSOLID%contacts(i)%states(j)%multiplier(1) + dlambda
  !         if( fstrSOLID%contacts(i)%fcoeff>0.d0 ) then
  !           if( fstrSOLID%contacts(i)%states(j)%state == CONTACTSLIP ) then
  !             fstrSOLID%contacts(i)%states(j)%multiplier(2) =             &
  !             fstrSOLID%contacts(i)%fcoeff * fstrSOLID%contacts(i)%states(j)%multiplier(1)
  !           else    ! stick
  !             !      fstrSOLID%contacts(i)%states(j)%multiplier(2) =
  !           endif
  !         endif
  !         lambda = fstrSOLID%contacts(i)%states(j)%multiplier(1)* fstrSOLID%contacts(i)%states(j)%direction
  !         wkarray((slave-1)*ndof+1:(slave-1)*ndof+3) = lambda(:)
  !         do k=1,nn
  !           iSS = fstrSOLID%contacts(i)%master(sid)%nodes(k)
  !           wkarray((iSS-1)*ndof+1:(iSS-1)*ndof+3) = wkarray((iSS-1)*ndof+1:(iSS-1)*ndof+3) -lambda(:)*shapefunc(k)
  !         enddo
  !       enddo
  !     enddo
  !     if( dsqrt(conv)<1.d-8 ) exit
  !     iter = iter+1
  !   enddo

  !   do i=1,hecMESH%n_node*ndof
  !     uc(i) = wkarray(i)/mmat(i)
  !   enddo
  ! end subroutine forward_increment_Lagrange

  subroutine phase_test_init(hecMESH,fstrPHASE,fstrHEAT,n_step) ! 2023/08/22(tanaka)
    !////////////////////////////////////////////////////////////////////////////////////////////
    ! Phase field(PF)の初期化と初期値格納,パラメーターの取得
    ! 
    ! PF : fstrPHASE%phaseをallocate & 1(未損傷状態)で初期化
    !
    ! その他の設定:
    ! ・ひずみ履歴場を使用する場合 : fstrPHASE%strain_energy_historyをallocate & 0で初期化
    ! ・引張-圧縮分解を使用する場合: fstrPHASE%stress_plus          をallocate & 0で初期化
    !
    ! 本subroutineは、以下のファイルのheat_initを参考に作成：
    !   FronISTRインストールディレクトリ/FrontISTR/fistr1/src/analysis/heat/heat_init.f90
    !////////////////////////////////////////////////////////////////////////////////////////////
    implicit none

    type(hecmwST_local_mesh), intent(in)  :: hecMESH
    type(fstr_phase),         intent(out) :: fstrPHASE
    type(fstr_heat),          intent(in)  :: fstrHEAT
    integer                 , intent(in)  :: n_step ! 総時間ステップ数 fstrDYN%n_step
    integer(kind=kint) :: i, j, ierror, ib, ii, ios, fn

    ! 外部ファイルからパラメーターを取得 (装置番号10が使用可能かどうか確認する方がよい)
    fn = 10
    open(fn,file='parameter.dat',action='read')
    read(fn,'(E14.7)',iostat=ios) fstrPHASE%Gc        ! 臨界エネルギー解放率
    read(fn,'(E14.7)',iostat=ios) fstrPHASE%L0        ! 微小な正則化長さパラメータ
    read(fn,'(1i)',   iostat=ios) fstrPHASE%history   ! ひずみ履歴場(1:考慮)
    read(fn,'(1i)',   iostat=ios) fstrPHASE%asymmetry ! 引張-圧縮分解(1:スペクトル分解)
                                          ! ↑変数名は%TensComp_Splitなどの方が良い？
    close(fn)
    ! fstrPHASE%Gc = 3.d0
    ! fstrPHASE%L0 = 0.25E-3
    ! 
    fstrPHASE%Omega = fstrPHASE%Gc / ( 2.d0 * fstrPHASE%L0 )


    !///// 割付と初期化 /////
    ! phase-field(PF)
    allocate( fstrPHASE%phase(hecMESH%n_node), stat=ierror )
    if( ierror /= 0 ) then
      write(idbg,*) 'stop due to allocation error of phase field'
      call hecmw_abort( hecmw_comm_get_comm() )
    end if
    !___ 1(未損傷状態)で初期化 ___
    fstrPHASE%phase(:) = 1.d0

    ! ひずみ履歴場
    if( fstrPHASE%history==1 ) then
      !↓ 積分点数4は、どの要素タイプでも対応できるように変数で入力した方が良い
      allocate( fstrPHASE%strain_energy_history(hecMESH%n_elem,NumOfQuadPoints(hecmesh%elem_type_item(1))), stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to allocation error of strain_energy_history'
        call hecmw_abort( hecmw_comm_get_comm() )
      end if

      !___ 0で初期化 ___
      fstrPHASE%strain_energy_history(:,:) = 0.d0
    endif

    ! 引張-圧縮分解(ひずみのスペクトル分解(Miehe et al. (2010)の方法))
    ! σ = (σ^+) + (σ^-)と分解したときの(σ^+)格納用配列
    if( fstrPHASE%asymmetry==1 ) then
      !↓ 4,4は、どの要素タイプでも対応できるように変数で入力した方が良い
      allocate( fstrPHASE%stress_plus(hecMESH%n_elem,NumOfQuadPoints(hecmesh%elem_type_item(1)),hecMESH%n_dof*2), stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to allocation error of stress_plus'
        call hecmw_abort( hecmw_comm_get_comm() )
      end if

      !___ 0で初期化 ___
      fstrPHASE%stress_plus(:,:,:) = 0.d0
    endif

    ! 出力するエネルギーの配列を割付
    ! energy_global( :, 1 ) 時刻
    ! energy_global( :, 2 ) ひずみエネルギー
    ! energy_global( :, 3 ) 運動エネルギー
    ! energy_global( :, 1 ) 破壊エネルギー
    allocate( fstrPHASE%energy_global( fstrDYN%n_step+1, 4 ))
    !///// 割付と初期化(終) /////


    !///// PFの初期条件・規定境界条件 /////

    ! 初期条件: 解析制御ファイルでPFの初期値が指定されている場合は格納 (解析制御ファイルでは初期温度として指定)
    if( associated(g_initialcnd) ) then
      do j=1,size(g_initialcnd)
        if( g_initialcnd(j)%cond_name=="temperature" ) then
          ! ↓配列演算できる?
          do i= 1, hecmesh%n_node
            fstrPHASE%phase(i)= g_initialcnd(j)%realval(i)
            ! write(*,*) hecmesh%global_node_id(i),fstr_PHASE%phase(i) ! デバッグ用
          enddo
          exit
        endif
      enddo
      ! write(ilog,*) ' Initial condition of temperatures: OK'
    endif

    ! PF規定境界条件: を反映 (解析制御ファイル !FIXTEMP で指定)
    if( fstrHEAT%T_FIX_tot>0 ) then
      do ib = 1, fstrHEAT%T_FIX_tot
        ii = fstrHEAT%T_FIX_node(ib)
        fstrPHASE%phase(ii) = fstrHEAT%T_FIX_VAL(ib)
      enddo
    endif
   
  end subroutine phase_test_init

  subroutine phase_test_output(hecMESH, fstrSOLID, fstrDYN, fstrPHASE, istep, res_output_interval, ene_output_interval) ! 2023/08/22(tanaka)
    !////////////////////////////////////////////////////////////////////////////////////////////
    ! phase-fieldの結果出力
    ! 出力ファイル名 test.res.info.'プロセス番号(分割メッシュ番号)'.'時間ステップ番号'
    ! 出力内容 :
    ! ・内点数、内点の節点番号、要素数、要素番号(、外点数、外点の節点番号)
    ! 出力ファイル名 test.res.result.'プロセス番号(分割メッシュ番号)'.'時間ステップ番号'
    ! 出力内容 :
    ! ・phase、変位、速度、加速度、節点応力、節点ひずみ、積分点応力、積分点ひずみ
    ! ※output_ctrl.dat(入力ファイル)で出力する物理量を指定可能
    !
    ! 節点応力orひずみが出力指定された場合は,その都度[fstr_test_NodalStress2D]で計算してから出力
    !
    !////////////////////////////////////////////////////////////////////////////////////////////
    implicit none

    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(fstr_solid),         intent(in) :: fstrSOLID
    type(fstr_dynamic),       intent(in) :: fstrDYN
    type(fstr_phase),         intent(in) :: fstrPHASE
    integer,                  intent(in) :: istep

    integer :: res_output_interval ! エネルギー以外の物理量の出力間隔
    integer :: ene_output_interval ! エネルギーの出力間隔

    integer(kind=kint) :: i, j, ios, fn

    integer(kind=kint), save :: num_str
    integer(kind=kint), save :: outputinfo(8)
    character(256),     save :: filename

    character(10) :: filenum
    character(2) :: ranknum

    fn = 10

    if( istep==0 )then
      ! 全体制御ファイルからファイル名(絶対パス込み)の取得
      open(fn, file='hecmw_ctrl.dat', action='read' )
      do i= 1, 10
        read(fn, '(a)', iostat=ios ) filename
        if( index(filename,'!RESULT')==1 )then
          read(fn, '(a)', iostat=ios ) filename
          exit
        endif
        if( ios<0 )then
          exit
        endif
      enddo
      close(fn)

      ! 物理量の出力指定を取得
      open(fn,file='output_ctrl.dat',action='read')

      do i = 1,8
        read(fn, '(1i)', iostat=ios ) outputinfo(i)
      enddo
      
      ! 結果ファイルの出力間隔を取得
      ! エネルギー以外
      read(10, '(6i)', iostat=ios ) res_output_interval
      write(*, '("res_output_interval has been set to ", i0)' ) res_output_interval
      ! エネルギー
      read(10, '(6i)', iostat=ios ) ene_output_interval
      write(*, '("ene_output_interval has been set to ", i0)' ) ene_output_interval
      
      close(10)

      ! 出力する応力(ひずみ)の成分数
      if( hecMESH%n_dof==2 )then
        num_str = 3
      elseif( hecMESH%n_dof==3 )then
        num_str = 6
      endif
    endif
    
    ! 結果を書き込むファイル名作成
    write(filenum,'(i10)') istep
    write(ranknum,'(i2)') hecMESH%my_rank

    ! 出力データの情報を出力
    if( istep==0 )then
      ! バイナリ・ストリーム出力でファイルopen
      open(fn, file=trim(adjustl(filename))//'.info.'//trim(adjustl(ranknum)), &
               action='write', access='stream', status='replace', form='unformatted' )
               
      ! 一括書込
      write(fn) hecMESH%n_dof, &                                  ! 1節点の変位自由度 (2 or 3)
                num_str, &                                        ! 応力・ひずみの成分数 (3 or 6)
                NumOfQuadPoints(hecmesh%elem_type_item(1)), &     ! 1要素あたりの積分点数 (4 or ?? )
                fstrDYN%t_delta, &                                ! dt
                hecMESH%nn_internal, &                            ! 内点数
                hecMESH%global_node_id(1:hecMESH%nn_internal), &  ! 内点数・内点の節点番号
                hecmesh%import_index(hecmesh%n_neighbor_pe), &    ! 外点数
                hecMESH%global_node_id(hecmesh%import_item), &    ! 外点の節点番号
                hecMESH%n_elem, &                                 ! 要素数
                hecMESH%global_elem_id(1:hecMESH%n_elem), &       ! 要素番号
                outputinfo, &                                     ! outputinfo(1:8)
                hecMESH%node(1:3*hecMESH%nn_internal)             ! 節点座標(2D計算でも3Dの情報を出力)
      close(fn)
    endif

    ! 結果の書き込み
    open(fn, file=trim(adjustl(filename))//'.result.'//&
                  trim(adjustl(ranknum))//'.'//trim(adjustl(filenum)),&
             action='write', access='stream', status='replace', form='unformatted' )

    !___ Nodal Phase field ___
    if( outputinfo(1)==1 ) write(fn) fstrPHASE%phase(1:hecMESH%nn_internal)

    !___ Nodal Displacement ___
    if( outputinfo(2)==1 ) write(fn) fstrDYN%DISP( 1:hecMESH%n_dof*hecMESH%nn_internal, 1 )
    
    !___ Nodal Particle Velocity ___
    if( outputinfo(3)==1 ) write(fn) fstrDYN%VEL( 1:hecMESH%n_dof*hecMESH%nn_internal, 1 )
    
    !___ Nodal Particle Acceleration ___
    if( outputinfo(4)==1 ) write(fn) fstrDYN%ACC( 1:hecMESH%n_dof*hecMESH%nn_internal, 1 )
    

    !___ Nodal Stresses or Nodal Strains ___
    if( outputinfo(5)==1 .or. outputinfo(6) )then ! 2023/10/25(tanaka)

      ! 節点応力と節点ひずみを更新
      call fstr_test_NodalStress2D( hecMESH, fstrSOLID )
    
      !___ Nodal Stresses ___
      if( outputinfo(5)==1 ) write(fn) fstrSOLID%STRESS( 1:num_str*hecMESH%nn_internal )

      !___ Nodal Strains ___
      if( outputinfo(6)==1 ) write(fn) fstrSOLID%STRAIN( 1:num_str*hecMESH%nn_internal )
    
    endif

    ! デバッグ用
    ! if(hecMESH%global_node_id(i)>1100) then
    !   write(*,*) hecMESH%global_node_id(i)
    !   write(*,*) fstrPHASE%phase(i)
    ! endif
    ! write(*,*) hecmesh%n_neighbor_pe
    ! write(*,*) hecmesh%export_index(hecmesh%n_neighbor_pe)
    ! write(*,*) hecmesh%import_index(hecmesh%n_neighbor_pe)
    ! write(*,*) hecmesh%export_item
    ! write(*,*) hecmesh%import_item
    ! write(*,*) hecMESH%global_node_id(hecmesh%export_item)
    ! enddo

    !___ Stresses @ integral points ___
    if( outputinfo(7)==1 )then
      do i = 1,hecMESH%n_elem
        do j = 1,NumOfQuadPoints(hecmesh%elem_type_item(1))
            write(fn) fstrSOLID%elements(i)%gausses(j)%stress(:)
        enddo
      enddo
    endif

    !___ Strains @ integral points ___
    if( outputinfo(8)==1 )then
      do i = 1,hecMESH%n_elem
        do j = 1,NumOfQuadPoints(hecmesh%elem_type_item(1))
            write(fn) fstrSOLID%elements(i)%gausses(j)%strain(:)
        enddo
      enddo
    endif

    close(fn)

    !___ 全ひずみエネルギー・全運動エネルギー・全破壊エネルギー ___
    if(istep==fstrDYN%n_step) then
      open(fn, file=trim(adjustl(filename))//'.result.'//&
                    trim(adjustl(ranknum))//'.'//'energy',&
              action='write', status='replace' )
      ! open(fn, file=trim(adjustl(filename))//'.result.'//&
      !               trim(adjustl(ranknum))//'.'//'energy',&
      !          action='write', access='stream', status='replace', form='unformatted' )
      write(fn,'E15.8') fstrPHASE%energy_global
      ! write(fn) fstrPHASE%energy_global
      close(fn)
    endif
    
  end subroutine phase_test_output

  subroutine phase_test_solve(hecMESH, fstrSOLID, fstrDYN, fstrHEAT, hecMAT_PHASE, fstrPHASE, ene_output_flag)
    !////////////////////////////////////////////////////////////////////////////////////////////
    ! phase-fieldの方程式を解くサブルーチン(+領域全体のひずみエネルギーの計算)
    !////////////////////////////////////////////////////////////////////////////////////////////
    use hecmw_matrix_misc
    use hecmw_solver_direct_clustermkl
    use m_heat_get_amplitude

    implicit none

    type(hecmwST_local_mesh),intent(in)    :: hecMESH
    type(fstr_solid),        intent(in)    :: fstrSOLID
    type(fstr_dynamic),      intent(in)    :: fstrDYN
    type(fstr_heat),         intent(in)    :: fstrHEAT
    type(hecmwST_matrix),    intent(inout) :: hecMAT_PHASE
    type(fstr_phase),        intent(inout) :: fstrPHASE
    logical,                 intent(in)    :: ene_output_flag ! 出力エネルギー計算のフラグ

    integer(kind=kint) :: ndof, itype ,is, iE, ic_type, icel, iiS, nn, j, i, bup_n_dof, ib, id, ii, isect, cid

    ! 1要素を構成する節点の情報格納用配列 (max_ncon:element connectivityの最大数)
    integer(kind=kint) :: nodLOCAL(fstrSOLID%max_ncon)          ! 節点番号
    real(kind=kreal)   :: phase(fstrSOLID%max_ncon)             ! PF
    real(kind=kreal)   :: vel(hecMESH%n_dof*fstrSOLID%max_ncon) ! 節点速度(2次元に限定しないように要修正)　←　修正2023/10/26(tanaka)
    real(kind=kreal)   :: ecoord(3,fstrSOLID%max_ncon)          ! 節点の三次元座標

    real(kind=kreal)   :: K(fstrSOLID%max_ncon,fstrSOLID%max_ncon),RHS(fstrSOLID%max_ncon)
    real(kind=kreal)   :: strain_energy_elem(hecMESH%n_elem)   ! 要素ひずみエネルギー
    real(kind=kreal)   :: kinetic_energy_elem(hecMESH%n_elem)  ! 要素運動エネルギー
    real(kind=kreal)   :: fracture_energy_elem(hecMESH%n_elem) ! 要素破壊エネルギー
    real(kind=kreal)   :: rho ! 質量密度
    ! real(kind=kreal)   :: dummy_history(4)
    ! real(kind=kreal)   :: dummy_stress(4,4)


    ! hecMAT_PHASE の初期化
    call hecmw_mat_clear( hecMAT_PHASE )
    call hecmw_mat_clear_b( hecMAT_PHASE )
    hecMAT_PHASE%X(:) = 0.d0

    ! 領域全体のエネルギーの初期化
    fstrPHASE%strain_energy_global   = 0.d0
    fstrPHASE%kinetic_energy_global  = 0.d0
    fstrPHASE%fracture_energy_global = 0.d0

    ! 節点自由度
    ndof = hecMAT_PHASE%NDOF

    ! itypeが何を表すのかメモ書きが欲しい
    ! itypeは要素タイプのインデックス( 1 ~ 解析で使用している要素タイプの数 )
    ! 複数の要素タイプを解析に用いる場合, do itype = 1, hecMESH%n_elem_type のようなループが要る
    ! 以下では要素番号の最初と最後,要素タイプを取得 → 次の要素ループで使用
    itype = 1 ! 現段階では要素タイプは1つしか使わないので1だけ
    is = hecMESH%elem_type_index(itype-1)+1
    iE = hecMESH%elem_type_index(itype  )
    ic_type= hecMESH%elem_type_item(itype)

    !element loop (OpenMP並列)
    ! OMPの共有属性設定
    !$omp parallel default(none), &
      !$omp&  private(icel,iiS,j,nn,nodLOCAL,i,ecoord,phase,vel, &
      !$omp&          strain_energy_elem,kinetic_energy_elem,fracture_energy_elem,K,RHS,isect,cid,rho), &
      !$omp&  shared(ene_output_flag,iS,iE,hecMESH,fstrSOLID,fstrDYN,ndof,hecMAT_PHASE,ic_type,fstrPHASE,&
      !$omp&         dummy_history,dummy_stress)
    !$omp do

    ! 要素ループ
    do icel = is, iE
      iiS = hecMESH%elem_node_index(icel-1)
      nn = hecMESH%elem_node_index(icel)-iiS ! １要素の節点数, 単一要素タイプならfstrSOLID%max_nconでもよい
      !if( nn>150 ) stop "elemental nodes > 150!"

      ! 要素icelを構成する節点に関するループ
      do j = 1, nn
        nodLOCAL(j) = hecMESH%elem_node_item(iiS+j)
        
        ! 節点nodLOCAL(j)の三次元座標
        do i = 1, 3
          ecoord(i,j) = hecMESH%node(3*nodLOCAL(j)+i-3) ! 2次元の場合でも座標は3次元分ある
        enddo

        ! 節点のphase-field 2023/08/24(tanaka)
        phase(j) = fstrPHASE%phase( nodLOCAL(j) )

        do i = 1, hecMESH%n_dof
          vel(i) = fstrDYN%VEL( hecMESH%n_dof*(nodLOCAL(j)-1) + i, 1 ) ! 節点速度は次元の分だけ割り当てられている
        enddo
      enddo

      isect = hecmesh%section_ID(icel)
      cid   = hecmesh%section%sect_mat_ID_item(isect)
      rho   = fstrsolid%materials(cid)%variables(m_density)

      ! 2023/09/13(tanaka)
      ! 要素の係数行列 K を計算　(+出力エネルギーの計算)
      call phase_test_update_Kmatrix( ic_type,nn,ecoord,fstrSOLID%elements(icel)%gausses(:),fstrSOLID%elements(icel)%iset, &
                                    phase,vel,fstrPHASE%L0,fstrPHASE%Omega,K,hecMESH%n_dof,rho,strain_energy_elem(icel),&
                                    kinetic_energy_elem(icel),fracture_energy_elem(icel),ene_output_flag,&
                                    fstrPHASE%history,fstrPHASE%asymmetry,fstrPHASE%strain_energy_history(icel,:),&
                                    fstrPHASE%stress_plus(icel,:,:) )

      ! 要素の右辺ベクトル RHS を計算
      call phase_test_update_RHS2( ic_type,nn,ecoord,fstrPHASE%Omega,RHS,hecMESH%n_dof )

      ! 要素の係数行列 K を全体の係数行列 hecMAT_PHASE%AU,hecMAT_PHASE%D にアセンブル
      call hecmw_mat_ass_elem(hecMAT_PHASE, nn, nodLOCAL, K)

      ! 要素の右辺ベクトル RHS を全体の右辺ベクトル hecMAT_PHASE%B にアセンブル
      do j = 1, nn
        !$omp atomic
        hecMAT_PHASE%B(nodLOCAL(j)) = hecMAT_PHASE%B(nodLOCAL(j)) + RHS(j)
      enddo
    enddo ! icel
    !$omp end do
    !$omp end parallel

    if(ene_output_flag) then
      ! 出力エネルギーの計算
      ! 時刻格納
      fstrPHASE%energy_global(fstrDYN%i_step+1,1) = fstrDYN%t_curr
      ! 要素のエネルギーを領域のエネルギーに足し合わせ
      fstrPHASE%energy_global(fstrDYN%i_step+1,2) = sum(strain_energy_elem(hecMESH%elem_internal_list(:)))
      fstrPHASE%energy_global(fstrDYN%i_step+1,3) = sum(kinetic_energy_elem(hecMESH%elem_internal_list(:)))
      fstrPHASE%energy_global(fstrDYN%i_step+1,4) = sum(fracture_energy_elem(hecMESH%elem_internal_list(:)))
      ! ! すべての分散領域のエネルギーを足し合わせ
      ! call hecmw_allreduce_r1 (hecMESH, fstrPHASE%strain_energy_global, hecmw_sum)
      ! call hecmw_allreduce_r1 (hecMESH, fstrPHASE%kinetic_energy_global, hecmw_sum)
      ! call hecmw_allreduce_r1 (hecMESH, fstrPHASE%fracture_energy_global, hecmw_sum)
    endif

    ! phase-field規定境界条件を係数行列,右辺ベクトルに反映 (解析制御ファイル !FIXTEMP で指定)
    if(fstrHEAT%T_FIX_tot>0) then
      do ib = 1, fstrHEAT%T_FIX_tot
        ii = fstrHEAT%T_FIX_node(ib)
        id = fstrHEAT%T_FIX_ampl(ib)
        call hecmw_mat_ass_bc(hecMAT_PHASE, ii, 1, fstrHEAT%T_FIX_VAL(ib))
      enddo
    endif

    ! cluster_sparse_solver で求解
    hecMAT_PHASE%Iarray(97) = 1 !Need numerical factorization
    ! bup_n_dof = hecMESH%n_dof
    ! hecMESH%n_dof = 1
    call hecmw_solve_direct_clustermkl(hecMESH,hecMAT_PHASE)
    ! hecMESH%n_dof = bup_n_dof

    ! 解ベクトルをMPI通信で更新
    call hecmw_update_R(hecMESH,hecMAT_PHASE%X,hecMESH%n_node, ndof) ! 2023/09/22(tanaka)

    ! デバッグ用
    ! do j=1,hecMESH%n_node
    !   write(*,*) hecMAT_PHASE%X(j)
    ! enddo
  end subroutine phase_test_solve

  subroutine phase_test_update_Kmatrix( ic_type,nn,ecoord,gausses,iset,phase,vel,L0,Omega,K,ndof,&
    rho,strain_energy_elem,kinetic_energy_elem,fracture_energy_elem,ene_output_flag,&
    history,asymmetry,strain_energy_history,stress_plus )
    !////////////////////////////////////////////////////////////////////////////////////////////
    ! phase-field方程式の要素係数行列 K (と要素のエネルギー) を計算するサブルーチン
    !////////////////////////////////////////////////////////////////////////////////////////////
    use m_MatMatrix

    implicit none

    type(tGaussStatus), intent(in)  :: gausses(:)
    integer(kind=kint), intent(in)  :: ic_type, nn ! ic_typeは要素タイプ,nnは1要素あたりの節点数
    integer(kind=kint), intent(in)  :: iset,ndof
    integer(kind=kint), intent(in)  :: history ! 1でひずみ履歴場
    integer(kind=kint), intent(in)  :: asymmetry ! 1でひずみのスペクトル分解
    real(kind=kreal),   intent(in)  :: ecoord(3,nn)
    real(kind=kreal),   intent(in)  :: L0,Omega,rho
    real(kind=kreal),   intent(out) :: K(nn,nn)
    real(kind=kreal),   intent(out) :: strain_energy_elem ! 要素のひずみエネルギー
    real(kind=kreal),   intent(out) :: kinetic_energy_elem ! 要素の運動エネルギー
    real(kind=kreal),   intent(out) :: fracture_energy_elem ! 要素の破壊エネルギー
    real(kind=kreal),   intent(in)  :: phase(nn)
    real(kind=kreal),   intent(in)  :: vel(2*nn)
    logical,            intent(in)  :: ene_output_flag ! 出力エネルギー計算のフラグ
    real(kind=kreal)                :: strain_energy_history(NumOfQuadPoints(ic_type))
    real(kind=kreal)                :: stress_plus(NumOfQuadPoints(ic_type),ndof*2)

    real(kind=kreal)   :: H(nn) ! 形状関数
    real(kind=kreal)   :: D(ndof*2,ndof*2) ! Dマトリックス
    real(kind=kreal)   :: HH(nn,nn) ! H'H  (「'」は転置の意味)
    real(kind=kreal)   :: BB(nn,nn) ! B'B  (「'」は転置の意味)
    real(kind=kreal)   :: det, wg
    real(kind=kreal)   :: gderiv(nn,ndof)
    real(kind=kreal)   :: localCoord(ndof)
    real(kind=kreal)   :: cdsys(3,3)
    integer(kind=kint) :: j, i, LX
    real(kind=kreal)   :: phase_gauss ! 積分点のphase-field
    real(kind=kreal)   :: phase_deriv_gauss(ndof) ! 積分点のphase-field勾配
    real(kind=kreal)   :: vel_gauss(ndof) ! 積分点の速度
    real(kind=kreal)   :: strain_energy_gauss ! 積分点のひずみエネルギー密度
    real(kind=kreal)   :: strain_energy_gauss_out ! 積分点のひずみエネルギー密度(出力エネルギー計算用)

    ! 引張圧縮非対称性
    real(kind=kreal)   :: eig_val(ndof) ! ひずみテンソルの固有値(dsyevの引数w)
    real(kind=kreal)   :: eig_vec(ndof,ndof) ! ひずみの固有ベクトル(dsyevの引数a),固有ベクトルは列ベクトル(:,n)で帰ってくる
    real(kind=kreal)   :: eig_diag(ndof,ndof) ! 固有値を対角に並べた行列
    character          :: jobz(1) = 'V' ! Nで固有値だけ,Vで固有ベクトルも返す
    character          :: uplo(1) = 'U' ! 入力する行列のformat,Uで上三角,Lで下三角をa(eig_vec)に入れる
    double precision   :: work(ndof*3) ! dsyevで使う配列,ndof*3にあまり意味はない(推奨サイズになるようにしている)
    integer            :: lwork ! workのサイズ,2次元なら6,3次元なら9にする
    integer            :: info ! call dsyev(...) のあとinfo=0なら正常
    real(kind=kreal)   :: strain_plus(ndof*2) ! 2次元のときはstrain_plus(4)は常にゼロ

    ! 要素の係数行列 K, 要素のひずみエネルギー, 運動エネルギー, 破壊エネルギーを初期化
    K(:,:) = 0.d0
    strain_energy_elem = 0.d0
    kinetic_energy_elem = 0.d0
    fracture_energy_elem = 0.d0

    ! 積分点ループ
    do LX=1, NumOfQuadPoints(ic_type)
      ! 積分点の局所座標を取得
      call getQuadPoint( ic_type, LX, localCoord(:) )
      ! 形状関数のグローバル座標微分とヤコビアンを計算
      call getGlobalDeriv( ic_type, nn, localcoord, ecoord, det, gderiv )
      ! 形状関数を計算
      call getShapeFunc( ic_type, localcoord, H(:) )
      ! Dマトリックスを計算
      call MatlMatrix( gausses(LX), iset, D, 1.d0, 1.d0, cdsys, 0.d0 )

      ! H'H を計算
      do j=1,nn
        do i=1,nn
          HH(i,j) = H(i) * H(j)
        enddo
      enddo

      ! B'B を計算
      do j=1,nn
        do i=1,nn
          BB(i,j) = dot_product( gderiv(i,:), gderiv(j,:) )
        enddo
      enddo

      ! 積分点のphase-field
      phase_gauss = dot_product( H(:), phase(:) ) ! ガウス点のphase-field ! 2023/09/14(tanaka)

      ! 重みxヤコビアン
      wg=getWeight( ic_type, LX )*det

      ! 引張圧縮非対称性による分岐
      if(asymmetry==1) then
      !ひずみテンソルのスペクトル分解で非対称性を入れる場合
        ! 2次元のstrain(+)
        if(ndof == 2) then
          ! ひずみテンソルの上三角部分をeig_vec(dsyevの入力行列)に格納
          eig_vec(1,1) = gausses(LX)%strain(1)
          eig_vec(2,2) = gausses(LX)%strain(2)
          eig_vec(1,2) = gausses(LX)%strain(3)

          ! dsyevで固有値,固有ベクトル計算
          lwork = 6 ! workの配列サイズ(推奨値にしている)
          call dsyev(jobz,uplo,ndof,eig_vec,ndof,eig_val,work,lwork,info)

          ! 固有値をeig_diagに格納
          eig_diag(:,:) = 0.d0
          eig_diag(1,1) = ( eig_val(1) + dabs(eig_val(1)) ) / 2
          eig_diag(2,2) = ( eig_val(2) + dabs(eig_val(2)) ) / 2

          ! strain(+) = eig_vec eig_diag (eig_vec)' 固有空間から元の空間に戻す
          eig_diag(:,:) = matmul(matmul(eig_vec(:,:),eig_diag(:,:)),transpose(eig_vec(:,:)))

          ! strain(+)をstrain_plusに格納
          strain_plus(:) = 0.d0
          strain_plus(1) = eig_diag(1,1)
          strain_plus(2) = eig_diag(2,2)
          strain_plus(3) = eig_diag(1,2)

        ! 3次元のstrain(+)
        elseif(ndof == 3) then
          ! ひずみテンソルの上三角部分をeig_vec(dsyevの入力行列)に格納
          eig_vec(1,1) = gausses(LX)%strain(1)
          eig_vec(2,2) = gausses(LX)%strain(2)
          eig_vec(3,3) = gausses(LX)%strain(3)
          eig_vec(1,2) = gausses(LX)%strain(4)
          eig_vec(2,3) = gausses(LX)%strain(5)
          eig_vec(1,3) = gausses(LX)%strain(6)

          ! dsyevで固有値,固有ベクトル計算
          lwork = 9 ! workの配列サイズ(推奨値にしている)
          call dsyev(jobz,uplo,ndof,eig_vec,ndof,eig_val,work,lwork,info)

          ! 固有値が正の成分を足し合わせてstrain(+)を計算
          strain_plus(:) = 0.d0
          do j=1,ndof
            strain_plus(1) = strain_plus(1) + ( eig_val(j) + dabs(eig_val(j)) ) / 2 * eig_vec(j,1) * eig_vec(j,1)
            strain_plus(2) = strain_plus(2) + ( eig_val(j) + dabs(eig_val(j)) ) / 2 * eig_vec(j,2) * eig_vec(j,2)
            strain_plus(3) = strain_plus(3) + ( eig_val(j) + dabs(eig_val(j)) ) / 2 * eig_vec(j,3) * eig_vec(j,3)
            strain_plus(4) = strain_plus(4) + ( eig_val(j) + dabs(eig_val(j)) ) / 2 * eig_vec(j,1) * eig_vec(j,2)
            strain_plus(5) = strain_plus(5) + ( eig_val(j) + dabs(eig_val(j)) ) / 2 * eig_vec(j,2) * eig_vec(j,3)
            strain_plus(6) = strain_plus(6) + ( eig_val(j) + dabs(eig_val(j)) ) / 2 * eig_vec(j,3) * eig_vec(j,1)
          enddo
        endif

        ! stress(+)を計算
        stress_plus(LX,:) = matmul( D(:,:), strain_plus(:) )

        ! 積分点のひずみエネルギー密度(+)
        ! strain_energy = 0.5 * stress(+)(ij)strain(+)(ij)
        strain_energy_gauss = 0.5d0 * dot_product( strain_plus(1:ndof), stress_plus(LX,1:ndof) ) + & ! xx,yy,(zz)成分
        dot_product( strain_plus(ndof+1:ndof*(ndof+1)/2), stress_plus(LX,ndof+1:ndof*(ndof+1)/2) ) ! xy,(yz,zx)成分

        if(ene_output_flag) then
          ! 劣化関数をかけたひずみエネルギーの計算
          ! t = dt*(istep+1) の変位と t = dt*istep のphaseから計算する
          ! strain_energy_out = phase**2.d0 * strain_energy(+) + strain_energy(-)
          strain_energy_gauss_out = phase_gauss**2.d0 * strain_energy_gauss + & ! phase**2.d0 * (+)
          0.5d0 * dot_product( gausses(LX)%strain(1:ndof)-strain_plus(1:ndof), & ! (-)のxx,yy,(zz)成分
          gausses(LX)%stress(1:ndof)-stress_plus(LX,1:ndof) ) + & ! (-)のxx,yy,(zz)成分
          dot_product( gausses(LX)%strain(ndof+1:ndof*(ndof+1)/2)-strain_plus(ndof+1:ndof*(ndof+1)/2), & ! xy,(yz,zx)成分
          gausses(LX)%stress(ndof+1:ndof*(ndof+1)/2)-stress_plus(LX,ndof+1:ndof*(ndof+1)/2) ) ! xy,(yz,zx)成分
        endif
        
      else
        ! 非対称性を入れない場合
        ! 積分点のひずみエネルギー密度(+)
        ! strain_energy_gauss = 0.5 * stress(ij)strain(ij)
        strain_energy_gauss = 0.5d0 * dot_product( gausses(LX)%strain(1:ndof), gausses(LX)%stress(1:ndof) ) + &! (-)のxx,yy,(zz)成分
        dot_product( gausses(LX)%strain(ndof+1:ndof*(ndof+1)/2), gausses(LX)%stress(ndof+1:ndof*(ndof+1)/2) )! (-)のxy,(yz,zx)成分

        if(ene_output_flag) then
          ! 劣化関数をかけたひずみエネルギーの計算
          ! t = dt*(istep+1) の変位と t = dt*istep のphaseから計算する
          ! strain_energy_out = phase**2.d0 * strain_energy
          strain_energy_gauss_out = phase_gauss**2.d0 * strain_energy_gauss ! phase**2.d0 * (+)
        endif
    endif

      ! ひずみ履歴場
      if(history==1) then
        strain_energy_gauss = max(strain_energy_gauss,strain_energy_history(LX))
        strain_energy_history(LX) = strain_energy_gauss
      endif

      ! 係数行列 K を計算
      K(:,:) = K(:,:) + 2.d0 * ( ( strain_energy_gauss + Omega ) * HH(:,:) + Omega * L0**2.d0 * BB(:,:) ) * wg


      ! 出力エネルギー(のため)の計算
      if(ene_output_flag) then
        ! 積分点のphase-field勾配
        do j=1,ndof
          phase_deriv_gauss(j) = dot_product( gderiv(j,:), phase(:) )
        enddo

        ! 積分点の速度
        vel_gauss(:) = 0.d0
        do j=1,nn
          do i=1,ndof
            vel_gauss(i) = vel_gauss(i) + H(j) * vel(ndof*(j-1)+i)
          enddo
        enddo

        ! 要素のひずみエネルギー(t = dt*(istep+1) の変位と t = dt*istep のphaseから計算)
         strain_energy_elem = strain_energy_elem + strain_energy_gauss_out * wg
        
        ! 要素の運動エネルギー(t = dt*(istep+1) の速度から計算)
        kinetic_energy_elem = kinetic_energy_elem + 0.5d0 * rho * dot_product( vel_gauss(:),vel_gauss(:) ) * wg

        ! 要素の破壊エネルギー(t = dt*istep のphaseから計算)
        fracture_energy_elem = fracture_energy_elem + Omega * &
        ( phase_gauss**2.d0 + L0**2.d0 * dot_product( phase_deriv_gauss(:),phase_deriv_gauss(:) ) ) * wg
      endif
    enddo
  end subroutine phase_test_update_Kmatrix

  subroutine phase_test_update_RHS2( ic_type,nn,ecoord,Omega,RHS,ndof )
    !////////////////////////////////////////////////////////////////////////////////////////////
    ! phase-field方程式の要素右辺ベクトルを計算するサブルーチン
    !////////////////////////////////////////////////////////////////////////////////////////////
    implicit none

    integer(kind=kint), intent(in)  :: ic_type, nn, ndof
    real(kind=kreal),   intent(in)  :: ecoord(3,nn)
    real(kind=kreal),   intent(in)  :: Omega
    real(kind=kreal),   intent(out) :: RHS(nn)

    integer(kind=kint)            :: LX
    real(kind=kreal)              :: H(nn)
    real(kind=kreal)              :: localCoord(ndof)
    real(kind=kreal)              :: gderiv(nn,ndof)
    real(kind=kreal)              :: I(nn)
    real(kind=kreal)              :: det, wg

    ! 要素の右辺ベクトル RHS を初期化
    RHS(:) = 0.d0

    ! 積分点ループ
    do LX=1, NumOfQuadPoints(ic_type)
      call getQuadPoint( ic_type, LX, localCoord(:) )
      call getShapeFunc( ic_type, localcoord, H(:) )
      call getGlobalDeriv( ic_type, nn, localcoord, ecoord, det, gderiv )
      !
      wg=getWeight( ic_type, LX )*det

      ! 係数行列 RHS を計算
      RHS(1:nn) = RHS(1:nn) + 2.d0 * Omega * H(1:nn) * wg
    enddo
  end subroutine phase_test_update_RHS2


end module fstr_test_dynamic_nlexplicit
