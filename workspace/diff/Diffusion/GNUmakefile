PROG	 	= Nucdif
PACKAGE  	= $(LIBDIR)/libqprop.a

ESQLBASICS    	= 
C++BASICS     	=   Nucdif.cpp  density.cpp grid.cpp fluid.cpp hamop.cpp 

PROGBASICS      = Nucdif

LOCALLIBS	=$(LIBDIR)/libqprop.a 

include         ../GNUmakefile.tmpl

#	User definable options

MYLIBS          = -lqprop  -lm 

COPTI           = -O8 #-pg
C++OPTI         = -O8 #-pg
C++SPECIAL	= -funroll-all-loops -finline-functions -fexpensive-optimizations -ffast-math -Wno-deprecated#-DTEST
LDFLAGS         = -Xlinker -defsym -Xlinker MAIN__=main\
                    
