# LDA is default - we change later if we want something other than LDA
XC_FUNCTIONALS = M_xc_option3_1c.o M_xc_option3_2c.o lda_ceperley_alder.o

# ****************************** XC ********************************
ifneq (,$(findstring BLYP, ${XC}))
XC_FUNCTIONALS = M_xc_option9_1c.o M_xc_option9_2c.o lda_ceperley_alder.o
endif

all: libthd_XC_FUNCTIONALS.a

.PHONY: clean

clean:
	rm -rf *.a *.o *.mod *genmod.f90 *~

libthd_XC_FUNCTIONALS.a: $(XC_FUNCTIONALS)
	ar rcs $@ $^ && mv $@ ../libs && cp *.mod ../include

%.o: %.f90
	${F90} ${FFLAGS} -c $< -I../include


