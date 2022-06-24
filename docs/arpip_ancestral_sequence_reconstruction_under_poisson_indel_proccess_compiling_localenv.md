---
layout: page
title: Compilation on a local environment
---

The user shoud have reading/writing rights on the system folders (i.e. /usr/local, /usr/local/include).


## Compiling and installing the dependencies


**bpp-core** http://biopp.univ-montp2.fr/

```
#!bash
git clone https://github.com/BioPP/bpp-core
cd bpp-core
git checkout tags/v2.4.1 -b v241
mkdir build
cd build
cmake ..
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
cmake ..
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
cmake  ..
make install
```

**boost - C++ Libraries** http://www.boost.org/

```
#!bash
wget https://boostorg.jfrog.io/artifactory/main/release/1.79.0/source/boost_1_79_0.tar.gz
tar xvf boost_1_79_0.tar.gz
cd boost_1_79_0
./bootstrap.sh --prefix=/usr/
./b2 
./b2 install 
```


**glog - Google Logging Library** https://github.com/google/glog

```
#!bash
git clone -b v0.5.0 https://github.com/google/glog
cd glog
cmake -H. -Bbuild -G "Unix Makefiles"
cmake --build build --target install

```

**gtest - Google Test Library** https://github.com/google/googletest/

```
#!bash
git clone https://github.com/google/googletest.git -b release-1.11.0
cd googletest        # Main directory of the cloned repository.
mkdir build          # Create a directory to hold the build output.
cd build
cmake ..             # Generate native build scripts for GoogleTest.
```



---

## Compiling ARPIP


*Dynamic linking*
```
#!bash
git clone https://github.com/acg-team/bpp-arpip/
cd bpp-arpip
cmake --target ARPIP -- -DCMAKE_BUILD_TYPE=Release CMakeLists.txt
make
```

*Static linking* 
```
#!bash
git clone https://github.com/acg-team/bpp-arpip/
cd bpp-arpip
cmake --target ARPIP -- -DCMAKE_BUILD_TYPE=Release-static CMakeLists.txt
make
```

