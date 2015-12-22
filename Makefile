CXXFLAGS   += -Wall -std=c++11
LDFLAGS     = -larmadillo -L/usr/lib64/atlas -lsatlas -O3
BIN         = out
SOURCES     = rho-linearV.cc io.cc
OBJ         = ${SOURCES:%.cc=%.o}

all: ${BIN}

${BIN}: ${OBJ}
	${CXX} ${OBJ} ${LDFLAGS} ${COPTFLAGS} -o $@

main.o: plot.h hamil.h io.h

rho-linearV.o: io.h

gcd.o: gcd.h

hamil.o: hamil.h gcd.h

plot.o: plot.h

io.o: io.h

clean:
	${RM} ${OBJ} plot.o

