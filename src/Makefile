# @Copyright 2022
# FIREBALL Committee
# Hong Kong University, Hong Kong Quantum AI Laboratory, Ltd. - James P. Lewis, Chair
# Universidad Autonoma de Madrid - Jose Ortega
# Academy of Sciences of the Czech Republic - Pavel Jelinek
# Arizona State University - Otto F. Sankey

# Previous and/or current contributors:
# Auburn University - Jian Jun Dong
# Caltech - Brandon Keith
# Czech Institute of Physics - Prokop Hapala
# Czech Institute of Physics - Vladimír Zobač
# Dublin Institute of Technology - Barry Haycock
# Pacific Northwest National Laboratory - Kurt Glaesemann
# University of Texas at Austin - Alex Demkov
# Ohio University - Dave Drabold
# Washington University - Pete Fedders
# West Virginia University - Ning Ma and Hao Wang
# also Gary Adams, Juergen Frisch, John Tomfohr, Kevin Schmidt,
#      and Spencer Shellman
#
# RESTRICTED RIGHTS LEGEND
# Use, duplication, or disclosure of this software and its documentation
# by the Government is subject to restrictions as set forth in subdivision
# { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
# clause at 52.227-7013. 
#
# Choose the options in the OPTIONS file
include include/OPTIONS

# Note: These are the machines files used by the Lewis Research Group. You need 
# to create a MACHINES file specific to your architecture.  You can use one of
# the machines file below as a starting point.  
include MACHINES/GENERIC
include MACHINES/HKQAI

##############################################################################
###
### BEGIN LISTING OBJECTS
###
##############################################################################
# ****************************************************************************
# BEGIN                    GENERICALLY COMPILED ROUTINES
# ****************************************************************************
DIRS = a.GLOBAL b.FUNCTIONS c.SYSTEM e.FDATA f.MPI g.XC_FUNCTIONALS h.DASSEMBLERS i.GRID h.SOLVESH j.ASSEMBLERS k.DASSEMBLERS l.SCF m.MD n.NAMD o.OUTPUT p.THEORY u.UTIL v.BEGIN x.DRIVERS

# Objects from a.GLOBAL
GLOBAL = M_precision.o M_assemble_blocks.o M_welcome.o

# Objects from b.FUNCTIONS
FUNCTIONS = clebsch_gordan.o clm.o factorial.o simpson.o distance.o magnitude.o smoother.o Ylm.o a_cross_b.o Depsilon_2c.o Depsilon_3c.o invert3x3.o mult3x1.o mult3x3.o Dsmoother.o

# Objects from c.SYSTEM
SYSTEM_BEGIN =  M_species.o M_atom_functions.o M_atomPP_functions_KB.o M_atomPP_ion_functions_KB.o M_integrals_2c.o
SYSTEM_LOOPRC = M_species.o M_atom_functions.o M_atomPP_functions_KB.o M_integrals_2c.o
SYSTEM_CREATE = M_species.o M_atom_functions.o M_integrals_2c.o M_integrals_3c.o M_atomPP_functions_KB.o M_atomPP_ion_functions_KB.o
SYSTEM_LIGHTNING = M_species.o M_kpoints.o M_configuraciones.o M_atom_functions.o M_neighbors.o M_neighbors_PP.o M_rotations.o M_grid.o rotate_PP.o
SYSTEM = M_species.o M_kpoints.o M_configuraciones.o M_neighbors.o M_neighbors_PP.o M_atom_functions.o M_rotations.o M_integrals_2c.o M_grid.o M_Drotations.o M_Drotations_PP.o M_build_forces.o rotate_PP.o M_vdW.o

BEGIN = M_psi.o M_rcatms_Coulomb.o M_rcatms_Coulomb_DMOL_excited.o M_rcatms_DMOL_excited.o M_rcatms.o M_vnn.o
BEGIN_LOOPRC = M_psi.o M_rcatms_Coulomb.o M_rcatms.o M_looprc.o

FUNCTIONS_CREATE = initialize_noMPI.o finalize_noMPI.o
FUNCTIONS_LIGHTNING = make_cells.o cl.o epsilon.o

MD = M_gear_quench.o
NAMD = M_non_adiabatic_gamma.o M_non_adiabatic_coupling.o M_density_MDET.o
OUTPUT = writeout_energies.o writeout_xyz.o writeout_xsf.o
UTIL = absorption_gamma.o dos_gamma.o M_isosurfaces.o
SOLVESH = M_diagonalization_Lowdin.o M_kspace_Lowdin.o M_density_matrix.o

# LDA is default - we change later if we want something other than LDA
XC_FUNCTIONALS = M_xc_option3_1c.o M_xc_option3_2c.o

# McWEDA is default - we change later if we want to do Horsfield
CREATE = M_overlap.o M_Goverlap.o M_kinetic.o M_Coulomb.o M_dipole_z.o M_vna_Harris.o M_vna_DOGS.o M_vnl_KB.o M_rho_2c_McWEDA_Harris.o M_rhoS_2c_McWEDA_Harris.o M_rho_3c_McWEDA_Harris.o M_rhoS_3c_McWEDA_Harris.o M_vxc_McWEDA_Harris.o M_vxc_McWEDA_DOGS.o M_bcna_Harris.o M_bcna_DOGS.o
# ****************************************************************************
# END                     GENERICALLY COMPILED ROUTINES
# ****************************************************************************

ifneq (,$(findstring McWEDA, ${THEORY}))
# ****************************************************************************
# BEGIN                      MCWEDA COMPILED ROUTINES
# 
# Note that McWEDA is the default for compiling. 
# Set THEORY = HORSFIELD if you want to compile the Horsfield subroutines. 
# ****************************************************************************
# ****************************************************************************
# Here we choose how to do the SCF
# No SCF - choose HARRIS in include/OPTIONS
# with SCF - choose DOGS in include/OPTIONS
# SCF on density - choose KS in include/OPTIONS
# ****************************************************************************
# ********************************* HARRIS ***********************************
ifneq (,$(findstring HARRIS, ${SCF}))
ASSEMBLERS = M_assemble_2c_Harris.o M_assemble_3c_Harris.o M_assemble_ewald_Harris.o M_assemble_rho_McWEDA_Harris.o M_assemble_usr_McWEDA_Harris.o M_assemble_vxc_McWEDA_Harris.o M_assemble_PP_2c.o M_assemble_PP_3c.o
DASSEMBLERS = M_Dassemble_rho_McWEDA_Harris.o M_Dassemble_2c_Harris.o M_Dassemble_PP_2c.o M_Dassemble_PP_3c.o M_Dassemble_ewald_Harris.o M_Dassemble_usr_McWEDA_Harris.o M_Dassemble_3c_Harris.o M_Dassemble_vxc_3c_McWEDA_Harris.o M_Dassemble_vxc_McWEDA_Harris.o 
FDATA = M_Fdata_1c_McWEDA_Harris.o M_Fdata_2c.o M_Fdata_3c.o
SCF = M_charges_Lowdin.o Qmixer_Anderson_Harris.o
endif

# **************************i******** DOGS ***********************************
ifneq (,$(findstring DOGS, ${SCF}))
ASSEMBLERS = M_assemble_2c_DOGS.o M_assemble_3c_DOGS.o M_assemble_PP_2c.o M_assemble_PP_3c.o M_assemble_rho_McWEDA_DOGS.o M_assemble_vxc_McWEDA_DOGS.o M_assemble_usr_McWEDA_DOGS.o M_assemble_ewald_DOGS.o
DASSEMBLERS = M_Dassemble_2c_DOGS.o M_Dassemble_3c_DOGS.o M_Dassemble_PP_2c.o M_Dassemble_PP_3c.o M_Dassemble_rho_McWEDA_DOGS.o M_Dassemble_vxc_3c_McWEDA_DOGS.o M_Dassemble_vxc_McWEDA_DOGS.o M_Dassemble_usr_McWEDA_DOGS.o M_Dassemble_ewald_DOGS.o
FDATA = M_Fdata_1c_McWEDA_DOGS.o M_Fdata_2c.o M_Fdata_3c.o
SCF = M_charges_Lowdin.o Qmixer_Anderson_DOGS.o
endif

# ************************************ KS ************************************
ifneq (,$(findstring KS, ${SCF}))
ASSEMBLERS = M_assemble_2c_grid.o M_assemble_3c_Harris.o M_assemble_ewald_DOGS.o M_assemble_rho_McWEDA_DOGS.o M_assemble_vxc_McWEDA_grid.o M_assemble_PP_2c.o M_assemble_PP_3c.o M_assemble_usr_grid.o
FDATA = M_Fdata_1c_McWEDA_DOGS.o M_Fdata_2c.o M_Fdata_3c.o
SCF = M_charges_Mulliken.o Qmixer_Anderson_grid.o
endif
# ****************************************************************************
# END                         MCWEDA COMPILED ROUTINES
# ****************************************************************************
endif


ifneq (,$(findstring HORSFIELD, ${THEORY}))
# ****************************************************************************
# BEGIN                      HORSFIELD COMPILED ROUTINES
# 
# Note that McWEDA is the default for compiling. 
# Set THEORY = HORSFIELD if you want to compile the Horsfield subroutines. 
# ****************************************************************************
CREATE = M_overlap.o M_Goverlap.o M_kinetic.o M_Coulomb.o M_dipole_z.o M_vna_Harris.o M_vna_DOGS.o M_vnl_KB.o M_vxc_Horsfield_Harris.o M_vxc_Horsfield_DOGS.o M_bcna_Harris.o M_bcna_DOGS.o lda_ceperley_alder.o M_xc3c_Horsfield_Harris.o M_xc3c_Horsfield_DOGS.o
# ****************************************************************************
# Here we choose how to do the SCF
# No SCF - choose HARRIS in include/OPTIONS
# with SCF - choose DOGS in include/OPTIONS
# SCF on density - choose KS in include/OPTIONS
# ****************************************************************************
# ********************************* HARRIS ***********************************
ifneq (,$(findstring HARRIS, ${SCF}))
ASSEMBLERS = M_assemble_2c_Harris.o M_assemble_3c_Harris.o M_assemble_ewald_Harris.o M_assemble_usr_Horsfield_Harris.o M_assemble_vxc_Horsfield_Harris.o M_assemble_PP_2c.o M_assemble_PP_3c.o
DASSEMBLERS = M_Dassemble_2c_Harris.o M_Dassemble_3c_Harris.o M_Dassemble_ewald_Harris.o M_Dassemble_usr_Horsfield_Harris.o M_Dassemble_vxc_Horsfield_Harris.o M_Dassemble_PP_2c.o M_Dassemble_PP_3c.o 
FDATA = M_Fdata_1c_Horsfield_Harris.o M_Fdata_2c.o M_Fdata_3c.o
SCF = M_charges_Lowdin.o Qmixer_Anderson_Harris.o
endif

# **************************i******** DOGS ***********************************
ifneq (,$(findstring DOGS, ${SCF}))
ASSEMBLERS = M_assemble_2c_DOGS.o M_assemble_3c_DOGS.o M_assemble_ewald_DOGS.o M_assemble_usr_Horsfield_DOGS.o M_assemble_vxc_Horsfield_DOGS.o M_assemble_PP_2c.o M_assemble_PP_3c.o
DASSEMBLERS = M_Dassemble_2c_DOGS.o M_Dassemble_3c_DOGS.o M_Dassemble_ewald_DOGS.o M_Dassemble_usr_Horsfield_DOGS.o M_Dassemble_vxc_Horsfield_DOGS.o M_Dassemble_PP_2c.o M_Dassemble_PP_3c.o 
FDATA = M_Fdata_1c_Horsfield_DOGS.o M_Fdata_2c.o M_Fdata_3c.o
SCF = M_charges_Lowdin.o Qmixer_Anderson_DOGS.o
endif

# ************************************ KS ************************************
ifneq (,$(findstring KS, ${SCF}))
#ASSEMBLERS = M_assemble_2c_grid.o M_assemble_3c_Harris.o M_assemble_ewald_DOGS.o M_assemble_rho_McWEDA_DOGS.o M_assemble_vxc_McWEDA_grid.o M_assemble_PP_2c.o M_assemble_PP_3c.o M_assemble_usr_grid.o
#FDATA = M_Fdata_1c_McWEDA_DOGS.o M_Fdata_2c.o M_Fdata_3c.o
#SCF = M_charges_Mulliken.o Qmixer_Anderson_grid.o
endif
# ****************************************************************************
# END                         HORSFIELD COMPILED ROUTINES
# ****************************************************************************
endif

# ****************************************************************************
# BEGIN                   DIAGONALIZATION COMPILED ROUTINES
# ****************************************************************************
# **************************** GAMMA - DEFAULT PACKAGES **********************
SOLVESH = M_diagonalization_Lowdin_gamma.o M_kspace_Lowdin.o M_density_matrix.o

# ************************* GAMMA, KPOINTS, MDET, or KS **********************
ifneq (,$(findstring GAMMAD, ${DIAGONALIZATION}))
SOLVESH = M_diagonalization_Lowdin_gammad.o M_kspace_Lowdin.o M_density_matrix.o
endif
ifneq (,$(findstring KPOINTS, ${DIAGONALIZATION}))
SOLVESH = M_diagonalization_Lowdin_kpoints.o M_kspace_Lowdin.o M_density_matrix.o
UTIL = absorption_gamma.o dos_kpoints.o M_isosurfaces.o 
endif
ifneq (,$(findstring MDET, ${DIAGONALIZATION}))
SOLVESH = M_diagonalization_Lowdin_gamma.o M_kspace_Lowdin.o M_density_matrix_MDET.o
SYSTEM = M_species.o M_kpoints_MDET.o M_configuraciones.o M_neighbors.o M_neighbors_PP.o M_atom_functions.o
endif
ifneq (,$(findstring KS, ${DIAGONALIZATION}))
GRID = M_project_grid.o M_fft_grid.o writeout_xsf.o
SOLVESH = M_diagonalization_kpoints.o M_kspace_grid.o M_density_matrix_grid.o M_project_grid.o M_fft_grid.o
endif

# ********************************** XC **************************************
ifneq (,$(findstring BLYP, ${XC}))
XC_FUNCTIONALS = M_xc_option9_1c.o M_xc_option9_2c.o
endif

OBJECTS_LOOPRC = $(GLOBAL) $(FUNCTIONS) $(SYSTEM_BEGIN) $(FUNCTIONS_BEGIN) $(XC_FUNCTIONALS) M_vxc_McWEDA_Harris.o $(BEGIN_LOOPRC)
OBJECTS_BEGIN = $(GLOBAL) $(FUNCTIONS) $(SYSTEM_BEGIN) $(FUNCTIONS_BEGIN) $(XC_FUNCTIONALS) M_vxc_McWEDA_Harris.o $(BEGIN)
OBJECTS_CREATE = $(GLOBAL) $(FUNCTIONS) $(SYSTEM_CREATE) $(FUNCTIONS_CREATE) $(XC_FUNCTIONALS) $(CREATE)
OBJECTS_LIGHTNING = $(GLOBAL) $(FUNCTIONS) $(SYSTEM_LIGHTNING) $(FUNCTIONS_LIGHTNING) $(FDATA) lda_ceperley_alder.o $(SOLVESH) $(ASSEMBLERS) $(SCF) $(OUTPUT) $(UTIL)
OBJECTS_FIREBALL = $(GLOBAL) $(FUNCTIONS) $(SYSTEM) $(FUNCTIONS_LIGHTNING) $(FDATA) lda_ceperley_alder.o $(SOLVESH) $(ASSEMBLERS) $(DASSEMBLERS) $(SCF) $(MD) $(OUTPUT) $(UTIL)
OBJECTS_MDET = $(GLOBAL) $(FUNCTIONS) $(SYSTEM)  $(FUNCTIONS_LIGHTNING) $(FDATA) lda_ceperley_alder.o $(SOLVESH) $(ASSEMBLERS) $(DASSEMBLERS) $(SCF) $(MD) $(NAMD) $(OUTPUT)
OBJECTS_LIGHTNING_GRID = $(GLOBAL) $(FUNCTIONS) $(SYSTEM_LIGHTNING) $(FUNCTIONS_LIGHTNING) $(FDATA) lda_ceperley_alder.o $(SOLVESH) $(GRID) $(ASSEMBLERS) $(SCF) $(OUTPUT) $(UTIL)
##############################################################################
###
### END LISTING OBJECTS
###
##############################################################################

# ********************************* LINKING **********************************
begin-looprc.x : $(OBJECTS_LOOPRC) begin-looprc.o
	$(F90) -o begin-looprc.x $(FFLAGS) $(OBJECTS_LOOPRC) begin-looprc.o $(LFLAGS)

begin.x : $(OBJECTS_BEGIN) begin.o
	$(F90) -o begin.x $(FFLAGS) $(OBJECTS_BEGIN) begin.o $(LFLAGS)

ifneq (,$(findstring McWEDA, ${THEORY}))
create.x : $(OBJECTS_CREATE) create_McWEDA.o
	$(F90) -o create.x $(FFLAGS) $(OBJECTS_CREATE) create_McWEDA.o $(LFLAGS)
endif
ifneq (,$(findstring HORSFIELD, ${THEORY}))
create.x : $(OBJECTS_CREATE) create_Horsfield.o
	$(F90) -o create.x $(FFLAGS) $(OBJECTS_CREATE) create_Horsfield.o $(LFLAGS)
endif

lightning.x: $(OBJECTS_LIGHTNING) lightning.o
	$(F90) -o lightning.x $(FFLAGS) $(OMP) $(OBJECTS_LIGHTNING) lightning.o $(LFLAGS)

lightning-grid.x: $(OBJECTS_LIGHTNING) lightning-grid.o
	$(F90) -o lightning-grid.x $(FFLAGS) $(OMP) $(OBJECTS_LIGHTNING) lightning-grid.o $(LFLAGS) $(LFLAGS_FFT)
	
fireball.x: $(OBJECTS_FIREBALL) fireball.o
	$(F90) -o fireball.x $(FFLAGS) $(OMP) $(OBJECTS_FIREBALL) fireball.o $(LFLAGS) 

fireball-MDET.x: $(OBJECTS_MDET) fireball-MDET.o
	$(F90) -o fireball-MDET.x $(FFLAGS) $(OBJECTS_MDET) fireball-MDET.o $(LFLAGS)

.PHONY: clean veryclean

clean:
	rm -f -r *.o *.mod *.mod.*
	rm -f -r core .nfs* rii_files begin-looprc.x.ip* begin.x.ip* \
			create.x.ip* fireball.x.ip* fireball-MDET.x.ip* \
			lightning.x.ip* *genmod.f90 .#* ldtmp* *.vo *~ *.il
	rm -f -r include/*.o include/*.mod include/*.mod.*
	rm -f -r a.GLOBAL/*.o a.GLOBAL/*.mod a.GLOBAL/*.mod.*
	rm -f -r b.FUNCTIONS/*.o b.FUNCTIONS/*.mod b.FUNCTIONS/*.mod.*
	rm -f -r c.SYSTEM/*.o c.SYSTEM/*.mod c.SYSTEM/*.mod.*
	rm -f -r d.FUNCTIONS_EXTRA/*.o d.FUNCTIONS_EXTRA/*.mod d.FUNCTIONS_EXTRA/*.mod.*
	rm -f -r e.FDATA/*.o e.FDATA/*.mod e.FDATA/*.mod.*
	rm -f -r f.MPI/*.o f.MPI/*.mod f.MPI/*.mod.*
	rm -f -r g.XC_FUNCTIONALS/*.o g.XC_FUNCTIONALS/*.mod g.XC_FUNCTIONALS/*.mod.*
	rm -f -r h.SOLVESH/*.o h.SOLVESH/*.mod h.SOLVESH/*.mod.*
	rm -f -r j.ASSEMBLERS/*.o j.ASSEMBLERS/*.mod j.ASSEMBLERS/*.mod.*
	rm -f -r k.DASSEMBLERS/*.o k.DASSEMBLERS/*.mod k.DASSEMBLERS/*.mod.*
	rm -f -r l.SCF/*.o l.SCF/*.mod l.SCF/*.mod.*
	rm -f -r m.MD/*.o m.MD/*.mod m.MD/*.mod.*
	rm -f -r n.NAMD/*.o n.NAMD/*.mod n.NAMD/*.mod.*
	rm -f -r o.OUTPUT/*.o o.OUTPUT/*.mod o.OUTPUT/*.mod.*
	rm -f -r p.THEORY/*.o p.THEORY/*.mod p.THEORY/*.mod.*

veryclean: clean
	rm -f begin-looprc.x begin.x create.x lightning.x lightning-grid.x \
			fireball.x fireball-MDET.x 

all:
	make begin-looprc.x
	make begin.x
	make create.x
	make fireball.x
	make fireball-MDET.x
	make lightning.x
	make lightning-grid.x
	
######################################################################
###
### DEPENDENCIES
###
######################################################################
# *****************************************************************************
# a.GLOBAL
# *****************************************************************************
M_precision.o : a.GLOBAL/M_precision.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_assemble_blocks.o : a.GLOBAL/M_assemble_blocks.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_welcome.o : a.GLOBAL/M_welcome.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^		

# *****************************************************************************
# b.FUNCTIONS
# *****************************************************************************
a_cross_b.o : b.FUNCTIONS/a_cross_b.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
clebsch_gordan.o : b.FUNCTIONS/clebsch_gordan.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
cl.o : b.FUNCTIONS/cl.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
clm.o : b.FUNCTIONS/clm.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
distance.o : b.FUNCTIONS/distance.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
erf0.o : b.FUNCTIONS/erf0.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
epsilon.o : b.FUNCTIONS/epsilon.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
Depsilon_2c.o : b.FUNCTIONS/Depsilon_2c.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
Depsilon_3c.o : b.FUNCTIONS/Depsilon_3c.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
factorial.o : b.FUNCTIONS/factorial.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
invert3x3.o : b.FUNCTIONS/invert3x3.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
magnitude.o : b.FUNCTIONS/magnitude.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
mult3x1.o : b.FUNCTIONS/mult3x1.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
mult3x3.o : b.FUNCTIONS/mult3x3.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
simpson.o : b.FUNCTIONS/simpson.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
smoother.o : b.FUNCTIONS/smoother.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
Dsmoother.o : b.FUNCTIONS/Dsmoother.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
Ylm.o : b.FUNCTIONS/Ylm.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^

# *****************************************************************************
# c.SYSTEM
# *****************************************************************************
M_atom_functions.o : c.SYSTEM/M_atom_functions.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_atomPP_functions_KB.o : c.SYSTEM/M_atomPP_functions_KB.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_atomPP_ion_functions_KB.o : c.SYSTEM/M_atomPP_ion_functions_KB.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_build_forces.o : c.SYSTEM/M_build_forces.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_configuraciones.o : c.SYSTEM/M_configuraciones.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Drotations.o : c.SYSTEM/M_Drotations.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Drotations_PP.o : c.SYSTEM/M_Drotations_PP.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_integrals_2c.o : c.SYSTEM/M_integrals_2c.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_integrals_3c.o : c.SYSTEM/M_integrals_3c.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_kpoints.o : c.SYSTEM/M_kpoints.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_kpoints_MDET.o : c.SYSTEM/M_kpoints_MDET.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_neighbors.o : c.SYSTEM/M_neighbors.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_neighbors_PP.o : c.SYSTEM/M_neighbors_PP.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_species.o : c.SYSTEM/M_species.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_rotations.o : c.SYSTEM/M_rotations.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_vdW.o : c.SYSTEM/M_vdW.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
make_cells.o : c.SYSTEM/make_cells.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
rotate_PP.o : c.SYSTEM/rotate_PP.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^

# *****************************************************************************
# e.FDATA
# *****************************************************************************
M_Fdata_1c_McWEDA_Harris.o : e.FDATA/M_Fdata_1c_McWEDA_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Fdata_1c_McWEDA_DOGS.o : e.FDATA/M_Fdata_1c_McWEDA_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Fdata_1c_Horsfield_Harris.o : e.FDATA/M_Fdata_1c_Horsfield_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Fdata_1c_Horsfield_DOGS.o : e.FDATA/M_Fdata_1c_Horsfield_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Fdata_2c.o : e.FDATA/M_Fdata_2c.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Fdata_3c.o : e.FDATA/M_Fdata_3c.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^

# *****************************************************************************
# f.MPI
# *****************************************************************************
finalize_MPI.o : f.MPI/finalize_MPI.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
finalize_noMPI.o : f.MPI/finalize_noMPI.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
initialize_MPI.o : f.MPI/initialize_MPI.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
initialize_noMPI.o : f.MPI/initialize_noMPI.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^

# *****************************************************************************
# g.XC_FUNCTIONALS
# *****************************************************************************
lda_ceperley_alder.o : g.XC_FUNCTIONALS/lda_ceperley_alder.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_xc_option3_1c.o : g.XC_FUNCTIONALS/M_xc_option3_1c.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_xc_option3_2c.o : g.XC_FUNCTIONALS/M_xc_option3_2c.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_xc_option9_1c.o : g.XC_FUNCTIONALS/M_xc_option9_1c.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_xc_option9_2c.o : g.XC_FUNCTIONALS/M_xc_option9_2c.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^

# *****************************************************************************
# h.SOLVESH
# *****************************************************************************
M_density_matrix.o : h.SOLVESH/M_density_matrix.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_density_matrix_grid.o : h.SOLVESH/M_density_matrix_grid.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_density_matrix_MDET.o : h.SOLVESH/M_density_matrix_MDET.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_kspace_grid.o : h.SOLVESH/M_kspace_grid.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_kspace_Lowdin.o : h.SOLVESH/M_kspace_Lowdin.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_diagonalization_Lowdin_gamma.o : h.SOLVESH/M_diagonalization_Lowdin_gamma.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_diagonalization_Lowdin_gammad.o : h.SOLVESH/M_diagonalization_Lowdin_gammad.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^		
M_diagonalization_Lowdin_kpoints.o : h.SOLVESH/M_diagonalization_Lowdin_kpoints.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_diagonalization_kpoints.o : h.SOLVESH/M_diagonalization_kpoints.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^

# *****************************************************************************
# i.GRID
# *****************************************************************************
M_grid.o : i.GRID/M_grid.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_project_grid.o : i.GRID/M_project_grid.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_fft_grid.o : i.GRID/M_fft_grid.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^

# *****************************************************************************
# j.ASSEMBLERS
# *****************************************************************************
M_assemble_2c_grid.o : j.ASSEMBLERS/M_assemble_2c_grid.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_assemble_2c_Harris.o : j.ASSEMBLERS/M_assemble_2c_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_assemble_2c_DOGS.o : j.ASSEMBLERS/M_assemble_2c_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_assemble_3c_Harris.o : j.ASSEMBLERS/M_assemble_3c_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_assemble_3c_DOGS.o : j.ASSEMBLERS/M_assemble_3c_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_assemble_ewald_Harris.o : j.ASSEMBLERS/M_assemble_ewald_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_assemble_ewald_DOGS.o : j.ASSEMBLERS/M_assemble_ewald_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_assemble_usr_grid.o : j.ASSEMBLERS/M_assemble_usr_grid.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_assemble_usr_Horsfield_Harris.o : j.ASSEMBLERS/M_assemble_usr_Horsfield_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_assemble_usr_Horsfield_DOGS.o : j.ASSEMBLERS/M_assemble_usr_Horsfield_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_assemble_usr_McWEDA_Harris.o : j.ASSEMBLERS/M_assemble_usr_McWEDA_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_assemble_usr_McWEDA_DOGS.o : j.ASSEMBLERS/M_assemble_usr_McWEDA_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_assemble_vxc_Horsfield_Harris.o : j.ASSEMBLERS/M_assemble_vxc_Horsfield_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_assemble_vxc_Horsfield_DOGS.o : j.ASSEMBLERS/M_assemble_vxc_Horsfield_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_assemble_rho_McWEDA_Harris.o : j.ASSEMBLERS/M_assemble_rho_McWEDA_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_assemble_rho_McWEDA_DOGS.o : j.ASSEMBLERS/M_assemble_rho_McWEDA_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_assemble_vxc_McWEDA_Harris.o : j.ASSEMBLERS/M_assemble_vxc_McWEDA_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_assemble_vxc_McWEDA_DOGS.o : j.ASSEMBLERS/M_assemble_vxc_McWEDA_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_assemble_vxc_McWEDA_grid.o : j.ASSEMBLERS/M_assemble_vxc_McWEDA_grid.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_assemble_PP_2c.o : j.ASSEMBLERS/M_assemble_PP_2c.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_assemble_PP_3c.o : j.ASSEMBLERS/M_assemble_PP_3c.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^

# *****************************************************************************
# k.DASSEMBLERS
# *****************************************************************************
M_Dassemble_blocks.o : k.DASSEMBLERS/M_Dassemble_blocks.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Dassemble_2c_Harris.o : k.DASSEMBLERS/M_Dassemble_2c_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Dassemble_2c_DOGS.o : k.DASSEMBLERS/M_Dassemble_2c_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Dassemble_PP_2c.o : k.DASSEMBLERS/M_Dassemble_PP_2c.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Dassemble_PP_3c.o : k.DASSEMBLERS/M_Dassemble_PP_3c.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Dassemble_3c_Harris.o : k.DASSEMBLERS/M_Dassemble_3c_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Dassemble_3c_DOGS.o : k.DASSEMBLERS/M_Dassemble_3c_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Dassemble_vxc_Horsfield_Harris.o : k.DASSEMBLERS/M_Dassemble_vxc_Horsfield_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Dassemble_vxc_Horsfield_DOGS.o : k.DASSEMBLERS/M_Dassemble_vxc_Horsfield_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Dassemble_rho_McWEDA_Harris.o : k.DASSEMBLERS/M_Dassemble_rho_McWEDA_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Dassemble_rho_McWEDA_DOGS.o : k.DASSEMBLERS/M_Dassemble_rho_McWEDA_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Dassemble_vxc_McWEDA_Harris.o : k.DASSEMBLERS/M_Dassemble_vxc_McWEDA_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Dassemble_vxc_McWEDA_DOGS.o : k.DASSEMBLERS/M_Dassemble_vxc_McWEDA_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Dassemble_vxc_3c_McWEDA_Harris.o : k.DASSEMBLERS/M_Dassemble_vxc_3c_McWEDA_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Dassemble_vxc_3c_McWEDA_DOGS.o : k.DASSEMBLERS/M_Dassemble_vxc_3c_McWEDA_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Dassemble_usr_Horsfield_Harris.o : k.DASSEMBLERS/M_Dassemble_usr_Horsfield_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Dassemble_usr_Horsfield_DOGS.o : k.DASSEMBLERS/M_Dassemble_usr_Horsfield_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Dassemble_usr_McWEDA_Harris.o : k.DASSEMBLERS/M_Dassemble_usr_McWEDA_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Dassemble_usr_McWEDA_DOGS.o : k.DASSEMBLERS/M_Dassemble_usr_McWEDA_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Dassemble_ewald_Harris.o : k.DASSEMBLERS/M_Dassemble_ewald_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Dassemble_ewald_DOGS.o : k.DASSEMBLERS/M_Dassemble_ewald_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^

# *****************************************************************************
# l.SCF
# *****************************************************************************
M_charges_Lowdin.o : l.SCF/M_charges_Lowdin.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_charges_Mulliken.o : l.SCF/M_charges_Mulliken.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
Qmixer_Anderson_Harris.o : l.SCF/Qmixer_Anderson_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
Qmixer_Anderson_DOGS.o : l.SCF/Qmixer_Anderson_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
Qmixer_Anderson_grid.o : l.SCF/Qmixer_Anderson_grid.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^

# *****************************************************************************
# m.MD
# *****************************************************************************
M_gear_quench.o : m.MD/M_gear_quench.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^

# *****************************************************************************
# n.NAMD
# *****************************************************************************
M_non_adiabatic_gamma.o : n.NAMD/M_non_adiabatic_gamma.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_non_adiabatic_coupling.o : n.NAMD/M_non_adiabatic_coupling.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_density_MDET.o : n.NAMD/M_density_MDET.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^

# *****************************************************************************
# o.OUTPUT
# *****************************************************************************
writeout_energies.o : o.OUTPUT/writeout_energies.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
writeout_xsf.o : o.OUTPUT/writeout_xsf.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
writeout_xyz.o : o.OUTPUT/writeout_xyz.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^

# *****************************************************************************
# p.THEORY
# *****************************************************************************
M_bcna_Harris.o : p.THEORY/M_bcna_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_bcna_DOGS.o : p.THEORY/M_bcna_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Coulomb.o : p.THEORY/M_Coulomb.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_dipole_z.o : p.THEORY/M_dipole_z.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_kinetic.o : p.THEORY/M_kinetic.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_overlap.o : p.THEORY/M_overlap.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_Goverlap.o : p.THEORY/M_Goverlap.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_rho_2c_McWEDA_Harris.o : p.THEORY/M_rho_2c_McWEDA_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_rhoS_2c_McWEDA_Harris.o : p.THEORY/M_rhoS_2c_McWEDA_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_rho_3c_McWEDA_Harris.o : p.THEORY/M_rho_3c_McWEDA_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_rhoS_3c_McWEDA_Harris.o : p.THEORY/M_rhoS_3c_McWEDA_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_vna_Harris.o : p.THEORY/M_vna_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_vna_DOGS.o : p.THEORY/M_vna_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_vnl_KB.o : p.THEORY/M_vnl_KB.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_vxc_Horsfield_Harris.o :  p.THEORY/M_vxc_Horsfield_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^		
M_vxc_Horsfield_DOGS.o :  p.THEORY/M_vxc_Horsfield_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^		
M_vxc_McWEDA_Harris.o :  p.THEORY/M_vxc_McWEDA_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_vxc_McWEDA_DOGS.o : p.THEORY/M_vxc_McWEDA_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_xc3c_Horsfield_Harris.o : p.THEORY/M_xc3c_Horsfield_Harris.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_xc3c_Horsfield_DOGS.o : p.THEORY/M_xc3c_Horsfield_DOGS.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^

# *****************************************************************************
# u.UTIL
# *****************************************************************************
absorption_gamma.o : u.UTIL/absorption_gamma.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
dos_gamma.o : u.UTIL/dos_gamma.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
dos_kpoints.o : u.UTIL/dos_kpoints.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_isosurfaces.o : u.UTIL/M_isosurfaces.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^

# *****************************************************************************
# v.BEGIN
# *****************************************************************************
M_looprc.o : v.BEGIN/M_looprc.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_psi.o : v.BEGIN/M_psi.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_rcatms.o : v.BEGIN/M_rcatms.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_rcatms_Coulomb.o : v.BEGIN/M_rcatms_Coulomb.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_rcatms_Coulomb_DMOL_excited.o : v.BEGIN/M_rcatms_Coulomb_DMOL_excited.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_rcatms_DMOL_excited.o : v.BEGIN/M_rcatms_DMOL_excited.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
M_vnn.o : v.BEGIN/M_vnn.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^

# *****************************************************************************
# x.DRIVERS
# ***************************************************************************** 
begin-looprc.o : x.DRIVERS/begin-looprc.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
begin.o : x.DRIVERS/begin.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
create_Horsfield.o : x.DRIVERS/create_Horsfield.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^ 
create_McWEDA.o : x.DRIVERS/create_McWEDA.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^ 
fireball.o : x.DRIVERS/fireball.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
fireball-MDET.o : x.DRIVERS/fireball-MDET.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
lightning.o : x.DRIVERS/lightning.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
lightning-grid.o : x.DRIVERS/lightning-grid.f90
		$(F90) $(DEFS) -I.. -I../.. $(FFLAGS) -c $^
