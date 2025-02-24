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
To install Kokkos, follow these steps:

1. **Download Kokkos**: Clone the Kokkos repository from GitHub using the following command:
   ```shell
   git clone https://github.com/kokkos/kokkos.git
   ```

2. **Checkout the desired version**: Navigate into the Kokkos directory and checkout the desired version (e.g., 4.0):
   ```shell
   cd kokkos
   git checkout tags/4.0
   ```

3. **Build Kokkos**: Create a build directory and run CMake to configure the build. Then compile using make:
   ```shell
   mkdir build
   cd build
   cmake .. -DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_CUDA=ON
   make -j
   ```

4. **Install Kokkos**: Optionally, you can install Kokkos to a specific directory:
   ```shell
   make install DESTDIR=/your/installation/path
   ```

Ensure that the `KOKKOS_PATH` in your machine file points to the correct installation path of Kokkos.

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
# For ELPA
 ELPA_PATH=/Your/path/ELPA
 elpaversion=elpa-2024.03.001
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
For compiling parallel Fireball(OpenMP/ELPA/ScaLAPACK/Kokkos)
```shell
make fireball-ase-elpa.x
```

## 3. How to Run the Code
Write your structure to be calculated by following the example in the `fire-opt.py` Python file. Make sure to modify the MPI control parameters to manage the number of resources used. Below is an example:

```python
import os
os.system('./clean.sh')
from thunder_ase.fireball import Fireball
import numpy as np
import ase
from ase import units
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.optimize import FIRE

max_step = 1
cell_min = 4

atoms = ase.io.read('Ga2O3_monoclinic.cif')
super_matrix = np.eye(3) * [int(cell_min / i)+1 for i in atoms.cell.lengths()]
atoms = ase.build.make_supercell(atoms, super_matrix)

Fdata_path = '/mytest/Fdata-McWEDA-0.15-3SN.Os3.35p3.80-3SNP.Gas4.85p5.60d5.60.Ins5.45p6.20d6.20'
kwargs = {
    'ipi': 1,
    }

fireball = Fireball(command=' mpirun -np 4  /mytest/fireball-ase-elpa.x',
    Fdata_path=Fdata_path,
    **kwargs)
atoms.calc = fireball

dyn = FIRE(atoms, trajectory='opt.traj', logfile='opt.log')
fireball.dynamics(dyn,fmax=0.2, steps=max_step)
```

Run the code using a command like:

```shell
python fire-opt.py
```
