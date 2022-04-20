FUNCTIONS_EXTRA = cl.o epsilon.o

all: libthd_FUNCTIONS_EXTRA.a

.PHONY: clean

clean:
	rm -rf *.a *.o *.mod *genmod.f90 *~

libthd_FUNCTIONS_EXTRA.a: $(FUNCTIONS_EXTRA)
	ar rcs $@ $^ && mv $@ ../libs

%.o: %.f90
	${F90} ${FFLAGS} -c $< -I../include
