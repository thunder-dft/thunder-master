OUTPUT = writeout_energies.o writeout_xsf.o writeout_xyz.o

all: libthd_OUTPUT.a

.PHONY: clean

clean:
	rm -rf *.a *.o *.mod *__genmod.* *~

libthd_OUTPUT.a: $(OUTPUT)
	ar rcs $@ $^ && mv $@ ../libs 

%.o: %.f90
	${F90} ${FFLAGS} -c $< -I../include
