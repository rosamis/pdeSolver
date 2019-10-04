    CC     = gcc -std=c11 -g
    
    CFLAGS = 
    LFLAGS = -lm

      PROG = pdeSolver
      OBJS = utils.o \
             pdeSolver.o\
             SistemasLineares.o \
             $(PROG).o

.PHONY: faxina clean distclean purge all doc

%.o: %.c %.h utils.h pdeSolver.h
	$(CC) -c $(CFLAGS) $<

$(PROG):  $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)

faxina:
	@rm -f *~ *.o #pdeSolver

clean:   faxina
	@rm -f *.o core a.out
	@rm -f $(PROG)


doc: $(OBJ) 
	doxygen ./Doxyfile
#doc: *.c trabalho-1.doxy *.h
#	@doxygen trabalho-1.doxy
