# Makefile for rho_ops
#
#CPP = icc
#CPPFLAGS = -Wall -O3 -fopenmp -lstdc++
CPP = g++
CPPFLAGS = -Wall -Wl,--allow-multiple-definition -O3 -fopenmp -llapack -lstdc++
Cln = /bin/rm -rf
NAME=rho_ops
###########################################
###########################################
SCR= main.cpp gauss_quad.cpp sphere_lebedev_rule.cpp legendre_quadrature.cpp hcubature.c Integrals_quadrature.cpp Critical_Points.cpp DMN_ops_p_class.cpp Input_commands.cpp Integrals_DMN.cpp Mathematical_Functions.cpp MOp_class.cpp String_ops.cpp gitver.cpp DMN_ops_class.cpp D_read_calc_rho.cpp Integrals.cpp MO_class.cpp NO_class.cpp NO_DMN_class.cpp gnuplot.cpp
OBJECTS= main.o gauss_quad.o sphere_lebedev_rule.o legendre_quadrature.o hcubature.o Integrals_quadrature.o Critical_Points.o DMN_ops_p_class.o Input_commands.o Integrals_DMN.o Mathematical_Functions.o MOp_class.o String_ops.o gitver.o DMN_ops_class.o D_read_calc_rho.o Integrals.o MO_class.o NO_class.o NO_DMN_class.o gnuplot.o
LIB=libcuba.a

all:
	./gitversion.sh 
	make RHO_OPS

RHO_OPS: $(OBJECTS) $(SCR)  $(LIB) cuba.h Makefile 
	$(CPP) $(CPPFLAGS) $(OBJECTS) $(LIB) -o RHO_OPS

%.o: %.cpp   
	$(CPP) $(CPPFLAGS) -c $*.cpp 
%.o: %.c
	$(CPP) -c $*.c 

clean:
	$(Cln) *.o
	$(Cln) *RHO_OPS
	$(Cln) *~
	$(Cln) gitver.cpp gitver.o gitver.h
	$(Cln) $(NAME).tar.gz 
tar:
	mkdir $(NAME)
	cp *.cpp *.h *.c Makefile ./$(NAME)
	tar -pczf $(NAME).tar.gz ./$(NAME)
	rm -r ./$(NAME)
	cp ./*.tar.gz ../
tex:
	latex Manual_RHO.OPS.tex
	latex Manual_RHO.OPS.tex
	dvipdf Manual_RHO.OPS.dvi
