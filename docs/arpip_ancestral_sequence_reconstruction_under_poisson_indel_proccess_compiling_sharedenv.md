---
layout: page 
title: Compilation on a shared environment
---
### Preparing the environment


We define a unique location for installing the libraries (lib + include). This location should be accessible with read/write permissions by the user.


```
#!bash
export SharedLibraryPath=path/to/folder
export SharedIncludePath=${SharedLibraryPath}/include

```

For example:


```
SharedLibraryPath=$HOME/local
SharedIncludePath=$HOME/local/include
```




### Compiling and installing the dependencies

**cmake** in the case which the minimum version is not provided
```
Wget https://cmake.org/files/v3.16/cmake-3.16.3.tar.gz
tar zxvf cmake-3.16.3.tar.gz
cd cmake-3.16.3
./bootstrap --prefix=${SharedLibraryPath}
make -j$(nproc)
make install
export PATH=$PATH:$SharedLibraryPath/bin
export PATH=$SharedLibraryPath/bin:$PATH

```



**bpp-core** http://biopp.univ-montp2.fr/

```
#!bash
git clone https://github.com/BioPP/bpp-core
cd bpp-core
git checkout tags/v2.4.1 -b v241
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${SharedLibraryPath} ..
make install
```

**bpp-seq** http://biopp.univ-montp2.fr/

```
#!bash
git clone https://github.com/BioPP/bpp-seq
cd bpp-seq
git checkout tags/v2.4.1 -b v241
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${SharedLibraryPath} ..
make install
```

**bpp-phyl**  http://biopp.univ-montp2.fr/

```
#!bash
git clone https://github.com/BioPP/bpp-phyl
cd bpp-phyl
git checkout tags/v2.4.1 -b v241
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${SharedLibraryPath} ..
make install
```
**boost - C++ Libraries** http://www.boost.org/

```
#!bash
wget https://boostorg.jfrog.io/artifactory/main/release/1.79.0/source/boost_1_79_0.tar.gz
tar xvf boost_1_79_0.tar.gz
cd boost_1_79_0
./bootstrap.sh --libdir=${SharedLibraryPath}/lib --includedir=${SharedIncludePath}
./b2 --libdir=${SharedLibraryPath}/lib --includedir=${SharedIncludePath}
./b2 install --prefix=${SharedLibraryPath}
```

**glog - Google Logging Library** https://github.com/google/glog/

```
#!bash
git clone -b v0.5.0 https://github.com/google/glog
cd glog
cmake -H. -Bbuild -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX=${SharedLibraryPath}
cmake --build build --target install

```
**gtest - Google Test Library** https://github.com/google/googletest/

```
#!bash
git clone https://github.com/google/googletest.git -b release-1.11.0
cd googletest        
cmake -H. -Bbuild -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX=${SharedLibraryPath}
cmake --build build --target install

```


### Compiling ARPIP


*Dynamic linking*
```
#!bash
git clone https://github.com/acg-team/bpp-arpip/
cd bpp-arpip
cmake --target ARPIP -- -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=${SharedIncludePath} -DCMAKE_PREFIX_PATH=${SharedLibraryPath} CMakeLists.txt
make
```
*Static linking*
```
#!bash
git clone https://github.com/acg-team/bpp-arpip/
cd bpp-arpip
cmake --target ARPIP -- -DCMAKE_BUILD_TYPE=Release-static -DCMAKE_PREFIX_PATH=${SharedIncludePath} -DCMAKE_PREFIX_PATH=${SharedLibraryPath} CMakeLists.txt
make
```