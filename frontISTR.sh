#!/bin/sh
###################################################################
# Simple automatic build script for FrontISTR
###################################################################
# Copyright (c) 2022 Yosuke ISHII
# The major design pattern of this script was abstracted
# from Michio Ogawa's script, which is subject to the same license
###################################################################
# Copyright (c) 2017-2021 Michio Ogawa
# This software is released under the MIT License, see LICENSE.txt
###################################################################
# Requirements
#  - curl
#  - git
#  - cmake
#__________________________________________________________________________
# Usage
#  sh fistr_build.sh
#__________________________________________________________________________
#  FrontISTR (fistr1, hec2rcap, hecmw_part1, hecmw_vis1, neu2fstr, rconv, rmerge)
#  の実行ファイルの格納される場所： ${BUILD_ROOT}/FrontISTR/bin
#
#  FrontISTRとHEC-MWのライブラリ一式の格納される場所：${BUILD_ROOT}/FrontISTR/lib
#  FrontISTRとHEC-MWのインクルードファイルは一式の格納される場所：${BUILD_ROOT}/FrontISTR/include
#
#  その他の各ソフトの実行ファイル・インクルードファイル・ライブラリの格納場所：
#     ${BUILD_ROOT}/local/bin
#     ${BUILD_ROOT}/local/include
#     ${BUILD_ROOT}/local/lib
#__________________________________________________________________________

#__________________________________________________________________________
# 自分の環境に合わせて編集
#__________________________________________________________________________

# ビルドするディレクトリ
BUILD_ROOT=`pwd`

# make時の並列数(自分のPCのコア数と一緒にするのがいいらしい)
MAKE_PAR=2

#__________________________________________________________________________
# 編集ここまで(一番下のmain関数内も適宜コメントアウトしたりすること)
#__________________________________________________________________________


# FrontISTR以外の各ソフトのbin,lib,includeを格納するディレクトリ
LIB_ROOT=${BUILD_ROOT}/local

# Misc. settings
CURL_FLAGS="-# -S --connect-timeout 10 --max-time 60 --retry 2"



########################################
# install Intel OneAPI compiler
########################################
install_compiler() {
    sudo mkdir ${BUILD_ROOT}/tmp
    cd ${BUILD_ROOT}/tmp

    sudo mkdir /usr/local/share/keyrings

    curl -sS https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | sudo gpg --dearmor -o /usr/local/share/keyrings/intel-oneapi-keyring.gpg
    sudo echo "deb [signed-by=/usr/local/share/keyrings/intel-oneapi-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/intel-oneapi.list
    
    sudo apt update
    sudo apt install intel-hpckit intel-basekit -y

    source /opt/intel/oneapi/setvars.sh > /dev/null
    source ~/.bashrc
}

########################################
# Set compiler (One API)
########################################
set_compiler() {
    CC=/opt/intel/oneapi/compiler/2023.1.0/linux/bin/intel64/icc
    CXX=/opt/intel/oneapi/compiler/2023.1.0/linux/bin/intel64/icpc
    FC=/opt/intel/oneapi/compiler/2023.1.0/linux/bin/intel64/ifort
    MPICC=/opt/intel/oneapi/mpi/2021.9.0/bin/mpiicc
    MPICXX=/opt/intel/oneapi/mpi/2021.9.0/bin/mpiicpc
    MPIFC=/opt/intel/oneapi/mpi/2021.9.0/bin/mpiifort
    CFLAGS="-O3 -xHost -mkl -warn all"; CXXFLAGS="-O3 -xHost -mkl -warn all"; FCFLAGS="-O3 -mkl -xHost -warn all"
    OMP="-qopenmp"
}

########################################
# cmake-3.20.2
########################################
CMAKE="cmake-3.20.2-linux-x86_64"
get_cmake() {
	if [ ! -d ${CMAKE} ]; then
		echo ">>>>> Getting " ${CMAKE} " <<<<<"
    curl ${CURL_FLAGS} -L -O https://cmake.org/files/LatestRelease/${CMAKE}.tar.gz
	else
		echo "Already download ${CMAKE}"
	fi
}
extract_cmake() {
	echo "extract latest binary cmake"
	tar xvf ${CMAKE}.tar.gz
	PATH=`pwd`/${CMAKE}/bin:$PATH
}

########################################
# metis-5.1.0
########################################
METIS="metis-5.1.0"
get_metis() {
  if [ ! -f ${METIS}.tar.gz ]; then
		echo ">>>>> Getting " ${METIS} " <<<<<"
    curl ${CURL_FLAGS} -L -O \
			http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/${METIS}.tar.gz
  else
    echo "Already downloaded ${METIS}.tar.gz"
  fi
}
build_metis() {
  if [ -f ${LIB_ROOT}/lib/libmetis.a ]; then
    echo "skip building ${METIS}"
    return
  fi
  if [ -f ${METIS}.tar.gz ]; then
    tar xvf ${METIS}.tar.gz
    cd ${METIS}
    make config prefix=${LIB_ROOT} cc=${CC} openmp=1
    make -j${MAKE_PAR}
    make install
    cd ${BUILD_ROOT}
  else
    echo "No ${METIS} archive"
  fi
}

########################################
# parmetis-4.0.3
# ATTN : License and edit build_mumps()
########################################
PARMETIS="parmetis-4.0.3"
get_parmetis() {
  if [ ! -f ${PARMETIS}.tar.gz ]; then
		echo ">>>>> Getting " ${PARMETIS} " <<<<<"
    curl ${CURL_FLAGS} -L -O http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/${PARMETIS}.tar.gz
  else
    echo "Already downloaded ${PARMETIS}.tar.gz"
  fi
}
build_parmetis() {
  if [ -f ${LIB_ROOT}/lib/libparmetis.a ]; then
    echo "skip to build ${PARMETIS}"
    return
  fi
  if [ -f ${PARMETIS}.tar.gz ]; then
    tar xvf ${PARMETIS}.tar.gz
    cd ${PARMETIS}
    make config prefix=${LIB_ROOT} cc=${MPICC} cxx=${MPICXX} openmp=1
    make -j${MAKE_PAR}
    make install
    cd ${BUILD_ROOT}
  else
    echo "No ${PARMETIS} archive"
  fi
}

########################################
# MUMPS-5.4.0.0
########################################
MUMPS="mumps"
get_mumps() {
  if [ ! -d ${MUMPS} ]; then
		echo ">>>>> Getting " ${MUMPS} " <<<<<"
    git clone https://github.com/scivision/${MUMPS}.git
    cd ${MUMPS}
    git checkout -b 5.4.0.0
  else
    echo "Already downloaded ${MUMPS}"
  fi
}
build_mumps() {
  if [ -f ${LIB_ROOT}/lib/libpord.a \
         -a -f ${LIB_ROOT}/lib/libdmumps.a \
         -a -f ${LIB_ROOT}/lib/libmumps_common.a ]; then
    echo "skip to build ${MUMPS}"
    return
  fi
  if [ -d ${MUMPS} ]; then
    cd ${MUMPS}
    mkdir build
    cd build
    cmake \
      -DCMAKE_INSTALL_PREFIX=${LIB_ROOT} \
      -DCMAKE_C_COMPILER=${MPICC} \
      -DCMAKE_CXX_COMPILER=${MPICXX} \
      -DCMAKE_Fortran_COMPILER=${MPIFC} \
      -Dopenmp=ON \
      -Dmetis=ON \
	    ..
    make -j${MAKE_PAR}
    make install
    cd ${BUILD_ROOT}
  else
    echo "No ${MUMPS} archive"
  fi
}

########################################
# Trilinos 13.4.0
########################################
TRILINOS="Trilinos"
get_trilinos() {
  if [ ! -d ${TRILINOS} ]; then
		echo "##### Getting " ${TRILINOS} "#####"
    git clone https://github.com/trilinos/${TRILINOS}.git
    cd ${TRILINOS}
    git checkout -b trilinos-release-13-4-0
  else
    echo "Already downloaded ${TRILINOS}"
  fi
}
build_trilinos() {
  if [ -f ${LIB_ROOT}/TrilinosRepoVersion.txt ]; then
    echo "skip to build ${TRILINOS}"
    return
  fi
  if [ -d ${TRILINOS} ]; then
    cd ${TRILINOS}
    mkdir build
    cd build
    cmake \
      -DCMAKE_INSTALL_PREFIX=${LIB_ROOT} \
      -DCMAKE_C_COMPILER=${MPICC} \
      -DCMAKE_CXX_COMPILER=${MPICXX} \
      -DCMAKE_Fortran_COMPILER=${MPIFC} \
      -DTPL_ENABLE_MPI=ON \
      -DTPL_ENABLE_LAPACK=ON \
      -DTPL_ENABLE_SCALAPACK=ON \
      -DTPL_ENABLE_METIS=ON \
      -DTPL_ENABLE_MUMPS=ON \
      -DTrilinos_ENABLE_ML=ON \
      -DTrilinos_ENABLE_Zoltan=ON \
      -DTrilinos_ENABLE_OpenMP=ON \
      -DTrilinos_ENABLE_Amesos=ON \
      -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF \
      -DTPL_ENABLE_MKL=ON \
      -DTPL_ENABLE_PARDISO_MKL=ON \
      -DMKL_INCLUDE_DIRS="${MKLROOT}/include" \
      -DMKL_LIBRARY_DIRS="${MKLROOT}/lib/intel64" \
      -DPARDISO_MKL_INCLUDE_DIRS="${MKLROOT}/include" \
      -DPARDISO_MKL_LIBRARY_DIRS="${MKLROOT}/lib/intel64" \
      -DAmesos_ENABLE_PARDISO_MKL=ON \
      -DBLAS_LIBRARY_DIRS="${MKLROOT}/lib/intel64" \
      -DLAPACK_LIBRARY_DIRS="${MKLROOT}/lib/intel64" \
      -DSCALAPACK_LIBRARY_DIRS="${MKLROOT}/lib/intel64" \
      -DBLAS_LIBRARY_NAMES="mkl_intel_lp64;mkl_intel_thread;mkl_core" \
      -DLAPACK_LIBRARY_NAMES="mkl_intel_lp64;mkl_intel_thread;mkl_core" \
      -DSCALAPACK_LIBRARY_NAMES="mkl_scalapack_lp64;mkl_blacs_intelmpi_lp64" \
      ..
    make -j${MAKE_PAR}
    make install
    cd ${BUILD_ROOT}
  else
    echo "No ${TRILINOS} archive"
  fi
}

########################################
# REVOCAP_Refiner-1.1.04
########################################
REFINER="REVOCAP_Refiner-1.1.04"
get_refiner() {
  if [ ! -f ${REFINER}.tar.gz -o -d REVOCAP_Refiner ]; then
		echo ">>>>> Getting " ${REFINER} " <<<<<"
    curl -L https://www.frontistr.com/download/link.php?${REFINER}.tar.gz -o ${REFINER}.tar.gz
    tar xvf ${REFINER}.tar.gz
    #git clone -b v1.1.04 https://github.com/FrontISTR/REVOCAP_Refiner
    echo "refiner"
  else
    echo "Already downloaded ${REFINER}.tar.gz"
  fi
}
build_refiner() {
  if [ -f ${LIB_ROOT}/lib/libRcapRefiner.a ]; then
    echo "skip to build ${REFINER}"
    return
  fi
  if [ -f ${REFINER}.tar.gz ]; then
    tar xvf ${REFINER}.tar.gz
    cd ${REFINER}
    cp MakefileConfig.LinuxIntelCompiler MakefileConfig.in
    make
    cp lib/x86_64-linux-intel/libRcapRefiner.a ${LIB_ROOT}/lib
    cp Refiner/rcapRefiner.h ${LIB_ROOT}/include
    cd ${BUILD_ROOT}
  else
    echo "No ${REFINER} archvie"
  fi
}

########################################
# FrontISTR
########################################
FRONTISTR="FrontISTR"
get_fistr() {
  if [ ! -d ${FRONTISTR} ]; then
		echo ">>>>> Getting " ${FRONTISTR} " <<<<<"
    git clone https://gitlab.com/FrontISTR-Commons/${FRONTISTR}.git
  else
    echo "Already downloaded ${FRONTISTR}"
  fi
}

build_fistr_from_cmake() {
  #_________________________________
  # cmake => make => make install
  #_________________________________
  if [ -d ${FRONTISTR} ]; then
    cd ${FRONTISTR}
    mkdir build; cd build

    #tmp_CC_dir=`which icx`
    #CC_dir=${tmp_CC_dir%/bin/*}

    cmake \
      -DCMAKE_INSTALL_PREFIX=${HOME}/local \
      -DCMAKE_PREFIX_PATH=${LIB_ROOT} \
      -DCMAKE_C_COMPILER=${CC} \
      -DCMAKE_CXX_COMPILER=${CXX} \
      -DCMAKE_Fortran_COMPILER=${FC} \
      -DBLAS_LIBRARIES="${MKLROOT}/lib/intel64/libmkl_intel_lp64.so;${MKLROOT}/lib/intel64/libmkl_intel_thread.so;${MKLROOT}/lib/intel64/libmkl_core.so" \
      -DLAPACK_LIBRARIES="${MKLROOT}/lib/intel64/libmkl_intel_lp64.so;${MKLROOT}/lib/intel64/libmkl_intel_thread.so;${MKLROOT}/lib/intel64/libmkl_core.so" \
      -DSCALAPACK_LIBRARIES="${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.so;${MKLROOT}/lib/intel64/libmkl_intel_lp64.so;${MKLROOT}/lib/intel64/libmkl_intel_thread.so;${MKLROOT}/lib/intel64/libmkl_core.so;${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.so;iomp5;pthread;m;dl" \
      -DWITH_MKL=1 \
      ..    
    make -j${MAKE_PAR}
    make install
  else
    echo "No ${FRONTISTR} archive"
  fi
}

build_fistr_from_setup_sh() {
  #___________________________________________________________________________
  # cmakeを使用せずに、Makefile.confの編集 => setup.sh => make => make install
  #___________________________________________________________________________
  if [ -d ${FRONTISTR} ]; then
    cd ${FRONTISTR}

    # ${FRONTISTR}ディレクトリ直下にMakefile.confを作成する関数を呼ぶ(自力で全部書く!!!)
    set_Makefile_conf

    # ${FRONTISTR}ディレクトリ直下にあるsetup.shを実行してMakefile作成(REVOCAP Couplerは使用しないので --with-revocapは付けない)
    ./setup.sh -p --with-metis --with-parmetis --with-tools --with-refiner --with-mkl --with-mumps --with-ml --with-lapack

    wait    
    #make -j${MAKE_PAR} ← 並列させるとエラーが出る？
    make
    make install

    # 生成されたライブラリとインクロードファイルを一か所にまとめる。
    cd ${BUILD_ROOT}/${FRONTISTR}
    mkdir -p lib include
    cp ./hecmw1/lib/* ./lib
    cp ./fistr1/lib/* ./lib
    cp ./hecmw1/include/* ./include
    cp ./fistr1/include/* ./include

    cd ${BUILD_ROOT}
  else
    echo "No ${FRONTISTR} archive"
  fi
}

set_Makefile_conf() {
  #____________________________________________________________________________________
  # ${FRONTISTR}ディレクトリ直下にMakefile.confの作成 
  #____________________________________________________________________________________

  # ヒアドキュメントはタブがそのまま反映されてしまうので左寄せで書く。
    #
    # おそらくMPIのディレクトリは指定しなくてもOK
    #
    #//////////////////////////////////////////////////////
    # mpiifortの絶対パスを取得(シングルクォーテーション'ではなくバッククォート`なので注意)
    # MPIIFORT_dir=`which mpiifort`
    # mpiのインストールディレクトリ：MPIIFORT_dirにおいて、binの上
    # MPI_dir=${MPIIFORT_dir%/bin/*}
    #
    #MPIDIR         = ${MPI_dir}
    #MPIBINDIR      = \$(MPIDIR)/bin
    #MPILIBDIR      = \$(MPIDIR)/lib
    #MPIINCDIR      = \$(MPIDIR)/include
    #MPILIBS        = -lmpicxx -lmpifort
    #//////////////////////////////////////////////////////

  cat -<<EOF > ./Makefile.conf
##################################################
#                                                #
#     Setup Configuration File for FrontISTR     #
#                                                #
##################################################

# MPI
MPIDIR         = \opt\intel\oneapi\mpi\2021.9.0
MPIBINDIR      = \$(MPIDIR)/bin
MPILIBDIR      = \$(MPIDIR)/lib
MPIINCDIR      = \$(MPIDIR)/include
MPILIBS        = -lmpicxx -lmpifort

# for install option only
PREFIX         = ${BUILD_ROOT}/FrontISTR
BINDIR         = \$(PREFIX)/bin
LIBDIR         = \$(PREFIX)/lib
INCLUDEDIR     = \$(PREFIX)/include

# Metis
METISDIR       = ${LIB_ROOT}/bin
METISLIBDIR    = ${LIB_ROOT}/lib
METISINCDIR    = ${LIB_ROOT}/include
HECMW_METIS_VER= 5

# ParMetis
PARMETISDIR    = ${LIB_ROOT}/bin
PARMETISLIBDIR = ${LIB_ROOT}/lib
PARMETISINCDIR = ${LIB_ROOT}/include

# Refiner
REFINERDIR     = ${LIB_ROOT}/bin
REFINERINCDIR  = ${LIB_ROOT}/include
REFINERLIBDIR  = ${LIB_ROOT}/lib

# Coupler
REVOCAPDIR     = ${LIB_ROOT}/bin
REVOCAPINCDIR  = ${LIB_ROOT}/include
REVOCAPLIBDIR  = ${LIB_ROOT}/lib

# MUMPS
MUMPSDIR       = ${LIB_ROOT}/bin
MUMPSINCDIR    = ${LIB_ROOT}/include
MUMPSLIBDIR    = ${LIB_ROOT}/lib
MUMPSLIBS      = -ldmumps -lmumps_common -lpord -L\${HOME}/local/lib

# MKL PARDISO
MKLDIR     = \$(HOME)/
MKLINCDIR  = \$(MKLDIR)/include
MKLLIBDIR  = \$(MKLDIR)/lib

# ML
MLDIR          = ${LIB_ROOT}/bin
MLINCDIR       = ${LIB_ROOT}/include
MLLIBDIR       = ${LIB_ROOT}/lib
MLLIBS         = -lml -lamesos -ltrilinosss -lzoltan -lepetra -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchoscore -ldmumps -lmumps_common -lpord -lmetis

# C compiler settings
CC             = mpiicc
CFLAGS         = -qopenmp -qmkl -diag-disable=10441
LDFLAGS        = -lm -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64
OPTFLAGS       = -O3

# C++ compiler settings
CPP            = mpiicpc
CPPFLAGS       = -qopenmp -qmkl -diag-disable=10441
CPPLDFLAGS     = -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64
CPPOPTFLAGS    = -O3

# Fortran compiler settings
F90            = mpiifort
F90FLAGS       = -nofor-main -qopenmp -qmkl
F90LDFLAGS     = -L\${HOME}/local/lib -lstdc++ -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64
F90OPTFLAGS    = -O2
F90FPP         = -fpp
F90LINKER      = mpiifort 

MAKE           = make
AR             = ar rUuv
MV             = mv -f
CP             = cp -f
RM             = rm -f
MKDIR          = mkdir -p
EOF
}


################################################################
# Main関数
# 既にインストールされているものはコメントアウトしてスキップすること.
################################################################
mkdir -p ${LIB_ROOT}/bin ${LIB_ROOT}/lib ${LIB_ROOT}/include
export PATH="${LIB_ROOT}/bin:$PATH"


# コンパイラ（Intel oneAPI）のインストール
#install_compiler
#wait
set_compiler

# cmakeのインストール
#get_cmake
#extract_cmake
#wait

# 各ソフトウェアのイインストールに必要なファイル等のダウンロード
get_metis &
get_parmetis &
get_refiner &
get_mumps &
get_trilinos &
get_fistr &
wait

# 各ソフトウェアのインストール
build_metis &
build_parmetis &
build_refiner &
wait

build_mumps
build_trilinos

#___ FrontISTRのインストール ___
build_fistr_from_cmake

# [2022/10/27] 
# build_fistr_from_cmakeは、cmakeを使ってMakefileを作成してインストールするが、この方法ではなぜかHEC-MWのライブラリ「libfhecmw.a」が生成されない。
# 精密には、本来生成されるべき「libhecmw.a(C言語用ライブラリ)」と「libfhecmw.a(Fortran用ライブラリ)」が、「libhecmw.a」にまとめられて生成されているようである。
# この場合、自作プログラムでライブラリをリンクしても、うまく動作しないことを確認した（コンパイルは通るが実行時にmpi関連のエラーが出る）。
# なので、cmakeを使用せずに、FrontISTRディレクトリ直下にある「Makefile.conf」を力技で設定 => 「setup.sh」でMakefileを作成 => make => make installする関数を作成した。
#
# [2022/10/31] ↑のcmakeの方で問題無いことを確認（mpi関連のエラーはHecmwのソースコード自体に原因があることがわかった）。なので、build_fistr_from_cmakeでインストールすればよい。
build_fistr_from_setup_sh