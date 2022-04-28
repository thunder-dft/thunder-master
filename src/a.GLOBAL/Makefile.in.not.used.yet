GLOBAL = M_assemble_blocks.o M_precision.o welcome.o

all: libthd_GLOBAL.a

.PHONY: clean

clean:
	rm -rf *.a *.o *.mod *genmod.f90 *~

libthd_GLOBAL.a: $(GLOBAL)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

%.o: %.f90
	${F90} ${FFLAGS} -c $< -I../include
