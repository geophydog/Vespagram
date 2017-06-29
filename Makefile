OBJ = sacio.o vespa.o
vespa : $(OBJ)
	cc -g -o vespa $(OBJ) -lm

$(OBJ) : sacio.h

clean :
	rm -rf vespa $(OBJ)
