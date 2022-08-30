---
layout: page
title: Preparing your system
---

## Development tools

Compiling ARPIP and its dependencies requires the basic development tools. Most of the times these tools are already present in the system.


### On Ubuntu 14.04/16.04/20.04

```
#!bash
sudo apt-get install build-essential git

```

### On Centos 6.9

```
#!bash
wget -O /etc/yum.repos.d/slc6-devtoolset.repo http://linuxsoft.cern.ch/cern/devtoolset/slc6-devtoolset.repo
yum install --nogpgcheck devtoolset-2 -y

```

Activate the environment either with:

```
#!bash
scl enable devtoolset-2 'bash'

```

or with appending the following line to your `.bahsrc` file:


```
#!bash
source scl_source enable devtoolset-2
```


### On Mac Os (10.13)

[Download and Install XCode](https://developer.apple.com/xcode/)


---

## Cmake (all the platforms)

ARPIP depends on CMake version 3.16.3


[**cmake 3.16.3**](http://cmake.org/)

```
#!bash
wget https://cmake.org/files/v3.16/cmake-3.16.3-Linux-x86_64.sh
chmod +x cmake-3.16.3-Linux-x86_64.sh
sudo ./cmake-3.16.3-Linux-x86_64.sh
export PATH=$PATH:path/to/cmake-3.16.3-Linux-x86_64/bin
```
