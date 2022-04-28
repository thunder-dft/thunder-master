NAMD=M_adiabatic.o M_non_adiabatic_gamma.o

all: libthd_NAMD.a 

.PHONY: clean

clean:
	rm -rf *.a *.o *.mod *genmod.o *~

libthd_NAMD.a: $(NAMD)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

%.o: %.f90
	${F90} ${FFLAGS} -c $< -I../include
