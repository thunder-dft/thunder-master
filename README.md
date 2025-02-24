# Fireball  - Set Up Guide

## 1. Required Environment 
To use Fireball, you need to install the following software:

- Intel MPI (Version 2021.3+)
- Intel Compiler  (Version 2021.3+)
- Intel MKL (Version 2021.3+)

The following optional software can greatly enhance Fireball in solving the eigenvalue problem:
- ElPA (Version 2024.03.001+)

The following optional software can support accelerating Ewald by GPU:
- Kokkos (Version 4.0+)

## 1.1 Installing Kokkos

**Note:** It is important to compile Kokkos for each specific GPU architecture you are using to avoid performance loss. 

**Software Requirements:**
- Ensure that the complete CUDA toolkit is installed on your system. The CUDA toolkit installed via Python is incomplete and cannot be used to compile Kokkos.
- It is recommended to download the CUDA toolkit runfile from the official NVIDIA website and install it manually.
- Be aware that improper combinations of CUDA and g++ versions may lead to errors. We have successfully tested CUDA 12.8 with g++ 11.5.

To install Kokkos, follow these steps:

1. **Download Kokkos**: Clone the Kokkos repository from GitHub using the following command:
   ```shell
   git clone https://github.com/kokkos/kokkos.git
   ```

2. **Build Kokkos**: Create a build directory and run CMake to configure the build. Then compile using make:
   ```shell
   cd kokkos
   mkdir build
   cd build
   cmake .. -DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_CUDA=ON -DCMAKE_INSTALL_PREFIX=/pure_kokkos_path
   make -j
   make install
   ```

Ensure that the `KOKKOS_PATH` in your machine file points to the correct installation path of Kokkos.

### Compiling Kokkos for Different GPU Architectures

When compiling Kokkos, you need to specify the architecture of the GPU you are targeting. Here are some examples for different architectures:

- **NVIDIA Volta (e.g., V100)**:
  ```shell
  cmake .. -DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_VOLTA70=ON
  ```

- **NVIDIA Ampere (e.g., A100)**:
  ```shell
  cmake .. -DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_AMPERE80=ON
  ```

- **AMD GFX906 (e.g., MI50)**:
  ```shell
  cmake .. -DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_HIP=ON -DKokkos_ARCH_AMD_GFX906=ON
  ```

Make sure to replace the architecture flag with the one that matches your GPU. This ensures that Kokkos is optimized for your specific hardware.

**Reminder:** Ensure that the version of Kokkos you are installing supports your specific GPU architecture. Refer to the Kokkos documentation for compatibility details.

## 1.2 Installing Fortran Language Compatibility Layer (FLCL)

**Note:** Ensure that the Fortran compiler used for compiling FLCL is the same as the one used for compiling Fireball to maintain compatibility.

To install FLCL, follow these steps:

1. **Download FLCL**: Clone the FLCL repository from GitHub using the following command:
   ```shell
    git clone https://github.com/kokkos/kokkos-fortran-interop.git
   ```

2. **Build FLCL**: Navigate into the FLCL directory, create a build directory, and run CMake to configure the build. Then compile using make:
   ```shell
   cd flcl
   mkdir build
   cd build
   cmake .. -DCMAKE_Fortran_COMPILER=ifx # Use the same Fortran compiler as used for compiling Fireball
   -DKokkos_DIR=/Kokkos_ROOT/lib/cmake/Kokkos -DFLCL_BUILD_EXAMPLES=OFF
   make -j
   ```

Ensure that the `FLCL_PATH` in your machine file points to the correct installation path of FLCL.

## 2. How to Compile

### 2.1 Modify the Options File
Edit the `OPTIONS` file to select the appropriate settings. Below is a default configuration suitable for parallel computation with the explanation of some important options.

```plaintext
# Choose the COMPILER
COMPILER = ifort

# Choose whether use kokkos
ENABLE_KOKKOS = FALSE

# Do you use DEBUG or Optimization Mode? (DEBUG/OPT)
MODE = OPT

# Turn on OpenMP or MPI (enter TRUE .or. FALSE)?
OPENMP = TRUE
MPI = TRUE
OPENMP_PROCEDURE = TRUE # OpenMP for vna3c
OPENMP_EWALD = TRUE #OpenMP for Ewald

# Choose the MACHINE file
MACHINE = CU10

# Choose HORSFIELD or McWEDA
THEORY = McWEDA

# Choose HARRIS, DOGS, or KS (Kohn-Sham)
SCF = DOGS

# Choose the DIAGONALIZATION Method, options are SCALAPACK\LAPACK\ELPA
DIAGONALIZATION = SCALAPACK

# Choose the XC functional
XC = LDA
```

### 2.2 Modify the Machine Files
Edit the files in the `machine` folder to suit your current environment and your `OPTIONS` file.  Please refer to CU10 to write your Machine File. Following are some important variables.

```plaintext

# For ScaLAPACK or LAPACK
 MKLROOT = /Your/path/mkl
# For Kokkos
 KOKKOS_DEVICES="HIP,OpenMP"
 KOKKOS_ARCH = "AMD_GFX906"
 KOKKOS_PATH=/Your/path/kokkos
 FLCL_PATH=/Your/path/flcl
 LINKFLAGS_kokkos = /Your/path/libkokkoscontainers.a  /Your/path/lib64/libkokkoscore.a /Your/path/lib64/libkokkossimd.a  ${FLCL_PATH}/lib64/libflcl.a 

```

### 2.3 Compile Using the PCCOMPILE File
For compiling serial Fireball
```shell
make fireball-ase.x
```
