!-------------------------------------------------------------------------------
 ! Copyright (c) 2019 FrontISTR Commons
 ! This software is released under the MIT License, see LICENSE.txt
 !-------------------------------------------------------------------------------
  
 module fstr_solver_dynamic
  
    use m_fstr
    use fstr_dynamic_nlexplicit
    use fstr_dynamic_nlimplicit
    use fstr_frequency_analysis  !Frequency analysis module
   
  contains
   
    !C================================================================C
    !C================================================================C
    subroutine fstr_solve_dynamic(hecMESH,hecMAT,fstrSOLID,fstrEIG   &
        ,fstrDYNAMIC,fstrRESULT,fstrPARAM &
        ,fstrCPL,fstrFREQ, hecLagMAT  &
        ,conMAT )
      use m_fstr_setup
      implicit none
      type(hecmwst_local_mesh)             :: hecmesh
      type(hecmwst_matrix)                 :: hecMAT
      type(fstr_eigen)                     :: fstrEIG
      type(fstr_solid)                     :: fstrSOLID
      type(hecmwst_result_data)            :: fstrRESULT
      type(fstr_param)                     :: fstrPARAM
      type(fstr_dynamic)                   :: fstrDYNAMIC
      type(fstr_couple)                    :: fstrCPL !for COUPLE
      type(fstr_freqanalysis)              :: fstrFREQ
      type(hecmwst_matrix_lagrange)        :: hecLagMAT
      type(fstr_info_contactchange)        :: infoCTChange
      type(hecmwst_matrix)                 :: conMAT
      integer(kind=kint) :: i, j, num_monit, ig, is, iE, ik, in, ing, iunitS, iunit, ierror, flag, limit
      character(len=HECMW_FILENAME_LEN) :: fname, header
      integer(kind=kint) :: restrt_step_num, ndof
      integer(kind=kint) :: restrt_step(1)
   
      num_monit = 0
   
      ! dtの確認
      if(dabs(fstrdynamic%t_delta) < 1.0e-20) then
        if( hecmesh%my_rank == 0 ) then
          write(imsg,*) 'stop due to fstrDYNAMIC%t_delta = 0'
        end if
        call hecmw_abort( hecmw_comm_get_comm())
      end if
      call fstr_dynamic_alloc( hecmesh, fstrdynamic )
   
      !C-- file open for local use
      !C
      if(fstrdynamic%idx_resp == 1) then   ! time history analysis
        call hecmw_ctrl_is_subdir( flag, limit )
        if( flag == 0)then
          header = ""
        else
          header = adjustl("MONITOR/")
          call hecmw_ctrl_make_subdir( adjustl("MONITOR/test.txt"), limit )
        endif
        ig = fstrdynamic%ngrp_monit
        is = hecmesh%node_group%grp_index(ig-1)+1
        ie = hecmesh%node_group%grp_index(ig)
        do ik=is,ie
          in = hecmesh%node_group%grp_item(ik)
          if (in > hecmesh%nn_internal) cycle
          num_monit = num_monit+1
          ing = hecmesh%global_node_id(in)
          iunits = 10*(num_monit-1)
   
          iunit = iunits + fstrdynamic%dynamic_IW4
          write(fname,'(a,i0,a)') trim(header)//'dyna_disp_',ing,'.txt'
          if(fstrdynamic%restart_nout < 0 ) then
            open(iunit,file=fname, position = 'append', iostat=ierror)
          else
            open(iunit,file=fname, status = 'replace', iostat=ierror)
          endif
          if( ierror /= 0 ) then
            write(*,*) 'stop due to file opening error:',trim(fname)
            call hecmw_abort( hecmw_comm_get_comm())
          end if
   
          iunit = iunits + fstrdynamic%dynamic_IW5
          write(fname,'(a,i0,a)') trim(header)//'dyna_velo_',ing,'.txt'
          if(fstrdynamic%restart_nout < 0 ) then
            open(iunit,file=fname, position = 'append', iostat=ierror)
          else
            open(iunit,file=fname, status = 'replace', iostat=ierror)
          endif
          if( ierror /= 0 ) then
            write(*,*) 'stop due to file opening error',trim(fname)
            call hecmw_abort( hecmw_comm_get_comm())
          end if
   
          iunit = iunits + fstrdynamic%dynamic_IW6
          write(fname,'(a,i0,a)') trim(header)//'dyna_acce_',ing,'.txt'
          if(fstrdynamic%restart_nout < 0 ) then
            open(iunit,file=fname, position = 'append', iostat=ierror)
          else
            open(iunit,file=fname, status = 'replace', iostat=ierror)
          endif
          if( ierror /= 0 ) then
            write(*,*) 'stop due to file opening error',trim(fname)
            call hecmw_abort( hecmw_comm_get_comm())
          end if
   
          iunit = iunits + fstrdynamic%dynamic_IW7
          write(fname,'(a,i0,a)') trim(header)//'dyna_force_',ing,'.txt'
          if(fstrdynamic%restart_nout < 0 ) then
            open(iunit,file=fname, position = 'append', iostat=ierror)
          else
            open(iunit,file=fname, status = 'replace', iostat=ierror)
          endif
          if( ierror /= 0 ) then
            write(*,*) 'stop due to file opening error',trim(fname)
            call hecmw_abort( hecmw_comm_get_comm())
          end if
          iunit = iunits + fstrdynamic%dynamic_IW8
          write(fname,'(a,i0,a)') trim(header)//'dyna_strain_',ing,'.txt'
          if(fstrdynamic%restart_nout < 0 ) then
            open(iunit,file=fname, position = 'append', iostat=ierror)
          else
            open(iunit,file=fname, status = 'replace', iostat=ierror)
          endif
          if( ierror /= 0 ) then
            write(*,*) 'stop due to file opening error',trim(fname)
            call hecmw_abort( hecmw_comm_get_comm())
          end if
   
          iunit = iunits + fstrdynamic%dynamic_IW9
          write(fname,'(a,i0,a)') trim(header)//'dyna_stress_',ing,'.txt'
          if(fstrdynamic%restart_nout < 0 ) then
            open(iunit,file=fname, position = 'append', iostat=ierror)
          else
            open(iunit,file=fname, status = 'replace', iostat=ierror)
          endif
          if( ierror /= 0 ) then
            write(*,*) 'stop due to file opening error',trim(fname)
            call hecmw_abort( hecmw_comm_get_comm())
          end if
        enddo
      endif
   
      !C
      !C-- initial condition
      !C
      fstrdynamic%DISP = 0.d0
      fstrdynamic%VEL  = 0.d0
      fstrdynamic%ACC  = 0.d0
      fstrsolid%unode(:)  =0.d0
      fstrsolid%QFORCE(:) =0.d0
      fstrdynamic%kineticEnergy=0.d0
      fstrdynamic%strainEnergy=0.d0
      fstrdynamic%totalEnergy=0.d0
      call fstr_updatestate( hecmesh, fstrsolid, 0.d0 )
   
      ! ---- restart
   
      restrt_step_num = 1
      fstrdynamic%i_step = 0
      infoctchange%contactNode_previous = 0
   
      if(associated(g_initialcnd))then
        ndof = hecmat%NDOF
        do j = 1, size(g_initialcnd)
          if(g_initialcnd(j)%cond_name == "velocity")then
            do i= 1, hecmesh%n_node
              ing = g_initialcnd(j)%intval(i)
              if(ing <= 0) cycle
              fstrdynamic%VEL(ndof*i-(ndof-ing),1) = g_initialcnd(j)%realval(i)
            end do
          elseif(g_initialcnd(j)%cond_name == "acceleration")then
            do i = 1, hecmesh%n_node
              ing = g_initialcnd(j)%intval(i)
              if(ing <= 0) cycle
              fstrdynamic%ACC(ndof*i-(ndof-ing),1) = g_initialcnd(j)%realval(i)
            enddo
          endif
        enddo
      endif
   
      if(fstrdynamic%restart_nout >= 0 ) then
        call dynamic_bc_init   (hecmesh, hecmat, fstrsolid, fstrdynamic)
        call dynamic_bc_init_vl(hecmesh, hecmat, fstrsolid, fstrdynamic)
        call dynamic_bc_init_ac(hecmesh, hecmat, fstrsolid, fstrdynamic)
      endif
   
      !restart
      if(fstrdynamic%restart_nout < 0 ) then
        if( .not. associated( fstrsolid%contacts ) ) then
          call fstr_read_restart_dyna_nl(restrt_step_num,hecmesh,fstrsolid,fstrdynamic,fstrparam)
        else
          call fstr_read_restart_dyna_nl(restrt_step_num,hecmesh,fstrsolid,fstrdynamic,fstrparam,&
            infoctchange%contactNode_previous)
        endif
        restrt_step_num = restrt_step_num + 1
        fstrdynamic%restart_nout = - fstrdynamic%restart_nout
        hecmat%Iarray(98) = 1
      end if
   
      if(fstrdynamic%idx_resp == 1) then   ! time history analysis
   
        if(fstrdynamic%idx_eqa == 1) then     ! implicit dynamic analysis
          if(.not. associated( fstrsolid%contacts ) ) then
            call fstr_solve_dynamic_nlimplicit(1, hecmesh,hecmat,fstrsolid,fstreig   &
              ,fstrdynamic,fstrresult,fstrparam &
              ,fstrcpl, restrt_step_num )
          elseif( fstrparam%contact_algo == kcaslagrange ) then
            call fstr_solve_dynamic_nlimplicit_contactslag(1, hecmesh,hecmat,fstrsolid,fstreig   &
              ,fstrdynamic,fstrresult,fstrparam &
              ,fstrcpl,heclagmat,restrt_step_num,infoctchange   &
              ,conmat )
          endif
   
        else if(fstrdynamic%idx_eqa == 11) then  ! explicit dynamic analysis
          call fstr_solve_dynamic_nlexplicit(hecmesh,hecmat,fstrsolid,fstreig   &
            ,fstrdynamic,fstrresult,fstrparam,infoctchange &
            ,fstrcpl, restrt_step_num )
        endif
   
      else if(fstrdynamic%idx_resp == 2) then  ! frequency response analysis
   
        if( fstrparam%nlgeom ) then
          if( hecmesh%my_rank .eq. 0 ) then
            write(imsg,*) 'stop: steady-state harmonic response analysis is not yet available !'
          end if
          call hecmw_abort( hecmw_comm_get_comm())
        end if
   
        if( hecmesh%my_rank .eq. 0 ) then
          call fstr_solve_frequency_analysis(hecmesh, hecmat, fstrsolid, fstreig, fstrdynamic, &
            fstrresult, fstrparam, fstrcpl, fstrfreq, heclagmat, &
            restrt_step_num)
        end if
      end if
   
      !C-- file close for local use
      if(fstrdynamic%idx_resp == 1) then   ! time history analysis
        do i=1,num_monit
          iunits = 10*(i-1)
          close(iunits + fstrdynamic%dynamic_IW4)
          close(iunits + fstrdynamic%dynamic_IW5)
          close(iunits + fstrdynamic%dynamic_IW6)
          close(iunits + fstrdynamic%dynamic_IW7)
          close(iunits + fstrdynamic%dynamic_IW8)
          close(iunits + fstrdynamic%dynamic_IW9)
        enddo
      endif
      !C-- end of finalization
   
      call fstr_dynamic_finalize( fstrdynamic )
      call hecmat_finalize( hecmat )
   
      deallocate( fstreig%mass )
   
    end subroutine fstr_solve_dynamic
   
  end module fstr_solver_dynamic