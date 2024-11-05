#!/bin/bash

echo "Let's start windows-setup"
# export PATH="/c/msys64/mingw64/bin:/c/Program Files/Git/bin:$PATH"
export PATH="/mingw64/bin:$PATH"

# # install GLPK for Windows from binary
# wget https://downloads.sourceforge.net/project/winglpk/winglpk/GLPK-4.65/winglpk-4.65.zip
# upzip winglpk-4.65.zip
# mkdir /mingw64/local
# cp -r winglpk-4.65/* /mingw64/

# # install ipopt from binary
# wget https://github.com/coin-or/Ipopt/releases/download/releases%2F3.14.16/Ipopt-3.14.16-win64-msvs2019-md.zip
# unzip Ipopt-3.14.16-win64-msvs2019-md.zip
# mkdir /mingw64/local
# cp -r Ipopt-3.14.16-win64-msvs2019-md/* /mingw64/

# # # install ipopt via https://coin-or.github.io/Ipopt/INSTALL.html
# # # install Mumps
# # git clone https://github.com/coin-or-tools/ThirdParty-Mumps.git
# # cd ThirdParty-Mumps
# # ./get.Mumps
# # ./configure --prefix=/mingw64/
# # make
# # make install
# # cd ..

# # # install ipopt from source
# # git clone https://github.com/coin-or/Ipopt.git
# # cd Ipopt
# # mkdir build
# # cd build
# # ../configure --prefix=/mingw64/
# # make
# # make install
# # cd ..
# # cd ..

# # get FindIPOPT_DIR from casadi, which is better written
# git clone --depth 1 --branch 3.6.5 https://github.com/casadi/casadi.git

# # install ifopt from source
# git clone https://github.com/ethz-adrl/ifopt.git
# cd ifopt
# # move FindIPOPT.cmake around
# mv ifopt_ipopt/cmake/FindIPOPT.cmake ifopt_ipopt/cmake/FindIPOPT.cmakeold
# cp ../casadi/cmake/FindIPOPT.cmake ifopt_ipopt/cmake/
# cp ../casadi/cmake/canonicalize_paths.cmake ifopt_ipopt/cmake/
# cmake -A x64 -B build \
#   -DCMAKE_VERBOSE_MAKEFILE=ON \
#   -DCMAKE_INSTALL_PREFIX="/mingw64/local" \
#   -DCMAKE_PREFIX_PATH="/mingw64" \
#   -G "Visual Studio 17 2022"

# cmake --build build --config Release
# cmake --install build --config Release
# cd ..

# # mkdir build
# # cd build
# # cmake .. -DCMAKE_VERBOSE_MAKEFILE=ON \
# #   -DCMAKE_INSTALL_PREFIX="/mingw64/local" \
# #   -DCMAKE_PREFIX_PATH="/mingw64" \
# #   -DIPOPT_LIBRARIES="/mingw64/lib/libipopt.dll.a" \
# #   -DIPOPT_INCLUDE_DIRS="/mingw64/include/coin-or" \
# #   -G "Unix Makefiles"

# # make VERBOSE=1
# # make install
# # cd ..
# # cd ..

# eigen_dir=$(cygpath -w /mingw64/share/eigen3/cmake)
# echo $eigen_dir
# echo "Eigen3_DIR=$eigen_dir" >> $GITHUB_ENV
# ifopt_dir=$(cygpath -w /mingw64/local/share/ifopt/cmake)
# echo `ls /mingw64/local/share/ifopt/cmake`
# echo $ifopt_dir
# echo "ifopt_DIR=$ifopt_dir" >> $GITHUB_ENV
eigen_dir=$(cygpath -w /mingw64/share/eigen3/cmake)
echo $eigen_dir
echo "Eigen3_DIR=$eigen_dir" >> $GITHUB_ENV

# glpk_include_dir=$(cygpath -w /mingw64/include)
# glpk_library=$(cygpath -w /mingw64/lib)
# echo "GLPK_INCLUDE_DIR=$glpk_include_dir" >> $GITHUB_ENV
# echo "GLPK_LIBRARY=$glpk_library" >> $GITHUB_ENV