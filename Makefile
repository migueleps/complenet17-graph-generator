EXEC_NAME=genk
CC=g++

CFLAGS= -w -std=c++11 -lm -pthread

SRC =                   \
        Fase.cpp        \
        Label.cpp       \
        IGtrie.cpp      \
        Isomorphism.cpp \
        Timer.cpp       \
        Timer2.cpp      \
        DynamicGraph.cpp\
        GraphMatrix.cpp \
        GraphUtils.cpp  \
        Random.cpp      \
        nauty/nauty.c   \
        nauty/nautil.c  \
        nauty/naugraph.c \
        GenSubgInv.cpp

OBJ =  ${SRC:.cpp=.o}

all: ${EXEC_NAME}

${EXEC_NAME}: ${OBJ}
	${CC} ${CFLAGS} ${CLIBS} -o ${EXEC_NAME} ${OBJ}

%.o: %.cpp
	${CC} ${CFLAGS} -c -o $@ $+

clean:
	rm ${EXEC_NAME} *.o *~ *# -rf
