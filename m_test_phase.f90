module m_test_phase
  use hecmw

  implicit none

  public

  ! phase-fieldの構造体 2023/09/12(tanaka)
  type fstr_phase
    real(kind=kreal), pointer     :: phase(:)
    real(kind=kreal)              :: Gc
    real(kind=kreal)              :: L0
    real(kind=kreal)              :: Omega
    real(kind=kreal), pointer     :: strain_energy_history(:,:) ! ひずみ履歴場
    real(kind=kreal), pointer     :: stress_plus(:,:,:)
    real(kind=kreal), allocatable :: energy_global(:,:) ! 時刻,領域全体のひずみエネルギー,運動エネルギー,破壊エネルギー
    integer(kind=kint)            :: history ! 1でひずみ履歴場あり
    integer(kind=kint)            :: asymmetry ! 1でひずみテンソルをスペクトル分解
  end type fstr_phase

end module m_test_phase