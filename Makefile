# Makefile for rho_ops
#
#CPP = c++ -DHAVE_LIBXC /usr/lib/x86_64-linux-gnu/libxc.so
CPP = c++
CPPFLAGS = -Wall -Wl,--allow-multiple-definition -O3 -fopenmp -llapack -lblas -lstdc++ 
Cln = /bin/rm -rf
NAME=rho_ops
###########################################
###########################################
SCR= main.cpp Corr_indicators.cpp utils_IO.cpp gauss_quad.cpp sphere_lebedev_rule.cpp legendre_quadrature.cpp hcubature.c Integrals_quadrature.cpp Integrals_atomic.cpp DMN_ops_p_class.cpp Input_commands.cpp Integrals_DMN.cpp Mathematical_Functions.cpp MOp_class.cpp String_ops.cpp gitver.cpp DMN_ops_class.cpp D_read_calc_rho.cpp Integrals.cpp MO_class.cpp NO_class.cpp NO_DMN_class.cpp gnuplot.cpp
SCR0= mescal.cpp mescal_utils.cpp mescal_IO.cpp
OBJECTS= main.o Corr_indicators.o utils_IO.o gauss_quad.o sphere_lebedev_rule.o legendre_quadrature.o hcubature.o Integrals_quadrature.o Integrals_atomic.o DMN_ops_p_class.o Input_commands.o Integrals_DMN.o Mathematical_Functions.o MOp_class.o String_ops.o gitver.o DMN_ops_class.o D_read_calc_rho.o Integrals.o MO_class.o NO_class.o NO_DMN_class.o gnuplot.o
OBJECTS0= mescal.o mescal_utils.o mescal_IO.o
LIB= libcuba.a libmescal.a 

all:
	./check_mescal.sh
	./gitversion.sh
	make mescal 
	make RHO_OPS

mescal: $(OBJECTS0) $(SCR0) Mescal.h
	ar rvs libmescal.a $(OBJECTS0)

RHO_OPS: $(OBJECTS) $(SCR) $(LIB) cuba.h Makefile 
	$(CPP) $(CPPFLAGS) $(OBJECTS) $(LIB) -o RHO_OPS

%.o: %.cpp   
	$(CPP) $(CPPFLAGS) -c $*.cpp 
%.o: %.c
	$(CPP) -c $*.c 

clean:
	$(Cln) *.o
	$(Cln) *RHO_OPS
	$(Cln) *~
	$(Cln) libmescal.a
	$(Cln) *mescal*cpp
	$(Cln) *Mescal.h
	$(Cln) gitver.cpp gitver.o gitver.h
	$(Cln) $(NAME).tar.gz 
tar:
	mkdir $(NAME)
	cp *.cpp *.h *.c Makefile ./$(NAME)
	cp -r mescal ./$(NAME)
	tar -pczf $(NAME).tar.gz ./$(NAME)
	rm -r ./$(NAME)
	cp ./*.tar.gz ../
tex:
	latex Manual_RHO.OPS.tex
	latex Manual_RHO.OPS.tex
	dvipdf Manual_RHO.OPS.dvi
