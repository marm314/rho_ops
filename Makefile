# Makefile for rho_ops
#
CPP = icc
#CPPFLAGS = -Wall -fopenmp -llapack -lstdc++
CPPFLAGS = -Wall -fopenmp -mkl -lstdc++
Cln = /bin/rm -rf
NAME=rho_ops
###########################################
###########################################
SCR= main.cpp sphere_lebedev_rule.cpp legendre_quadrature.cpp hcubature.c Integrals_quadrature.cpp Critical_Points.cpp DMN_ops_p_class.cpp Input_commands.cpp Integrals_DMN.cpp Mathematical_Functions.cpp MOp_class.cpp String_ops.cpp DMN_ops_class.cpp D_read_calc_rho.cpp Integrals.cpp main.cpp MO_class.cpp NO_class.cpp NO_DMN_class.cpp gnuplot.cpp
OBJECTS= main.o sphere_lebedev_rule.o legendre_quadrature.o hcubature.o Integrals_quadrature.o Critical_Points.o DMN_ops_p_class.o Input_commands.o Integrals_DMN.o Mathematical_Functions.o MOp_class.o String_ops.o DMN_ops_class.o D_read_calc_rho.o Integrals.o MO_class.o NO_class.o NO_DMN_class.o gnuplot.o
LIB=libcuba.a

all: 
	make rhoops

rhoops: $(OBJECTS) $(SCR)  $(LIB) cuba.h Makefile 
	$(CPP) $(CPPFLAGS) $(OBJECTS) $(LIB) -o $(NAME)_EXE

%.o: %.cpp   
	$(CPP) $(CPPFLAGS) -c $*.cpp 
%.o: %.c
	$(CPP) -c $*.c 

clean:
	$(Cln) *.o
	$(Cln) *$(NAME)_EXE
	$(Cln) *~
	$(Cln) $(NAME).tar.gz 
tar:
	mkdir $(NAME)
	cp *.cpp *.h *.c Makefile ./$(NAME)
	tar -pczf $(NAME).tar.gz ./$(NAME)
	rm -r ./$(NAME)
	cp ./*.tar.gz ../
