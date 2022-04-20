UTIL=M_isosurfaces.o absorption_gamma.o  dos_gamma.o

all: libthd_UTIL.a

.PHONY: clean

clean:
	rm -rf *.a *.o *.mod *genmod.f90 *~

libthd_UTIL.a: $(UTIL)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

%.o: %.f90
	${F90} ${FFLAGS} -c $< -I../include
