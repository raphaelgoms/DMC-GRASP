FLAGS = -lpthread -g
LIBS=-Wall -I libGRASP -I libDTS -I libFuncao -I libUtil -I libBFGS -I /usr/local/include -I .
CC=g++

libGRASPObjs = \
xmeans.o \
cgrasp.o \
main.o

libFuncaoObjs = \
ap.o \
Funcao.o \
Rosenbrock2.o \
Zakharov.o \
SumSquares.o \
Branin.o \
Easom.o \
GoldsteinPrice.o \
Hartmann.o \
Shekel.o \
Shubert.o \
Beale.o \
Booth.o \
Bohachevsky.o \
Hump.o \
Matyas.o \
Schwefel.o \
Colville.o \
Perm.o \
Perm0.o \
PowerSum.o \
Griewank.o \
Rastrigin.o \
Trid.o \
Powell.o \
DixonPrice.o \
Ackley.o \
Levy.o \
Sphere.o 

libUtilObjs = \
mt19937ar.o \
simclist.o \
Util.o 

libGRASPObjsPre = $(addprefix libGRASP/,$(libGRASPObjs) )
libFuncaoObjsPre = $(addprefix libFuncao/,$(libFuncaoObjs) )
libUtilObjsPre  = $(addprefix libUtil/, $(libUtilObjs)  )

OBJECTS=$(libGRASPObjsPre) $(libFuncaoObjsPre) $(libUtilObjsPre) 

all:  CGrasp

CGrasp: $(OBJECTS)
	$(CC) -g -o CGrasp $(OBJECTS) -L /usr/X11/lib $(LIBS) $(FLAGS)

.c.o: $<
	$(CC) -g -c $< -o $@ $(LIBS) 

.cpp.o: $<
	$(CC) -g -c $< -o $@ $(LIBS)

clean:
	@find -iname "CGrasp" -exec rm {} \;
	@find -iname "*.o" -exec rm {} \;
	@find -iname "*.so" -exec rm {} \;
	@find -iname "*~" -exec rm {} \;
	@find -iname "*.swp" -exec rm {} \;

