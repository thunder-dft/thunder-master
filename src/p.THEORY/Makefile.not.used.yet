THEORY=M_Coulomb.o  M_dipole_z.o  M_kinetic.o  M_overlap.o
VXC_MCWEDA = M_vxc_McWEDA_Harris.o
VXC_HORSFIELD = M_vxc_Horsfield_Harris.o
VXC_DOGS = M_vxc_McWEDA_DOGS.o

all: libthd_VXC_MCWEDA.a libthd_VXC_HORSFIELD.a libthd_THEORY.a

.PHONY: clean

clean:
	rm -rf *.a *.o *.mod *genmod.o *~

libthd_THEORY.a: $(THEORY)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

libthd_VXC_MCWEDA.a: $(VXC_MCWEDA)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

libthd_VXC_DOGS.a: $(VXC_DOGS)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

libthd_VXC_HORSFIELD.a: $(VXC_HORSFIELD)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

%.o: %.f90
	${F90} ${FFLAGS} -c $< -I../include
