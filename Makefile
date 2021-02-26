PROGRAM = betkeWeil
CC = g++
CFLAGS = -O2

USERDIR = /mnt/c/Users/ferin

CAPDBINLIB = ${USERDIR}/lib/capd/bin
CAPDLIBS = `${CAPDBINLIB}/capd-config --cflags --libs`

INCLUDE = -I.

HEADERFILES = main.h betkeWeil.h jetTools.h tools.h

${PROGRAM}: ${PROGRAM}.cpp pugixml.o ${HEADERFILES}
	${CC} ${CFLAGS} ${INCLUDE} ${PROGRAM}.cpp -o ${PROGRAM} pugixml.o ${CAPDLIBS}

pugixml.o: pugixml.hpp
	${CC} ${CFLAGS} ${INCLUDE} -o pugixml.o -c pugixml.cpp

clean:
	rm *~ ${PROGRAM}

tar:
	tar cvfz ${PROGRAM}.tgz ${PROGRAM}.cpp ${HEADERFILES} Makefile
