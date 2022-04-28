OBJECTS_MPI = finalize_MPI.o initialize_MPI.o
OBJECTS_noMPI = finalize_noMPI.o initialize_noMPI.o


all: libthd_MPI.a libthd_noMPI.a


.PHONY: clean

clean:
	rm -rf *.a *.o *.mod *genmod.f90 *~

libthd_MPI.a: $(OBJECTS_MPI)
	ar rcs $@ $^ && mv $@ ../libs 

libthd_noMPI.a: $(OBJECTS_noMPI)
	ar rcs $@ $^ && mv $@ ../libs 

%.o: %.f90
	${F90} ${FFLAGS} -c $<
