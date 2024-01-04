# Galois Field

## Software requirements
CMake and at least C++17 

## Overview
Implementation of the Galois field GF(2^293) in a normal basis

## Features
- **Adding operation**
- **Multiplication operation** ( using lamda matrix )
- **Trace calculation**
- **Inverse calculation** ( using algorithm Itoh-Tsujii )
- **Power calculation** ( using Horner`s method )
- **Square calculation** ( just cyclic shift >>> 1 )

## Build and Usage

Follow these steps to build and run:

1. Clone the GaloisField repository to your local machine:
```bash
https://github.com/avept/GaloisField_Normal_Basis.git
```

2. Navigate to the project directory:
```bash 
cd GaloisField
```

3. Create a build directory:
```bash
mkdir build
cd build
```
   
4. Generate the build system using CMake:
```bash
cmake ..
```

5. Build the library:
```bash
make -j4
```
