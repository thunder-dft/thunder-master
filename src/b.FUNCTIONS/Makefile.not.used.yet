FUNCTIONS=a_cross_b.o  clebsch_gordan.o  clm.o  Depsilon_2c.o  Depsilon_3c.o  distance.o  Dsmoother.o  erf0.o  factorial.o  invert3x3.o  magnitude.o  mult3x1.o  mult3x3.o  simpson.o  smoother.o  Ylm.o

all: libthd_FUNCTIONS.a

.PHONY: clean

clean:
	rm -rf *.a *.o *.mod *genmod.f90 *~

libthd_FUNCTIONS.a: $(FUNCTIONS)
	ar rcs $@ $^ && mv $@ ../libs 

%.o: %.f90
	${F90} ${FFLAGS} -c $< -I../include
