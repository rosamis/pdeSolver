    CC     = gcc -std=c11 -g -O3 -mavx2 -march=native
    
      PROG = pdeSolver
      OBJS = utils.o \
             pdeSolver.o\
             SistemasLineares.o \
             $(PROG).o
             
      LIKWID = /home/soft/likwid
LIKWID_FLAGS = -I$(LIKWID)/include -D LIKWID_PERFMON
 LIKWID_LIBS = -L$(LIKWID)/lib -llikwid

      CFLAGS = $(LIKWID_FLAGS) 
      LFLAGS = $(LIKWID_LIBS) -lm
#    AVX_FLAGS = -march=native -mavx -O3 -falign-functions=32 -falign-loops=32 -fstrict-aliasing
#   AVX_LOG_FLAGS = -fopt-info-vec -fopt-info-vec-missed 


.PHONY: faxina clean distclean purge all doc likwid avx normaL2 gauss debug

all likwid:  $(PROG)

debug:  CFLAGS += -DDEBUG

#avx:           CFLAGS += $(AVX_FLAGS) $(AVX_LOG_FLAGS)
#avx likwid:    CFLAGS += -DLIKWID_PERFMON
#avx likwid:    LFLAGS += -llikwid

gauss: CFLAGS += -Dgauss
normaL2: CFLAGS += -DnormaL2
gaussTEMPO: CFLAGS += -DgaussTEMPO
normaL2TEMPO: CFLAGS += -DnormaL2TEMPO

gauss: likwid
normaL2: likwid

gaussTEMPO: $(PROG)
normaL2TEMPO: $(PROG)

#likwid avx debug: $(PROG)


%.o: %.c %.h utils.h pdeSolver.h
	$(CC) -c $(CFLAGS) $<

$(PROG):  $(PROG).o

$(PROG):  $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)

faxina:
	@rm -f *~ *.o #pdeSolver

clean:   faxina
	@rm -f *.o core a.out
	@rm -f $(PROG)


#doc: $(OBJ) 
#	doxygen ./Doxyfile
#doc: *.c trabalho-1.doxy *.h
#	@doxygen trabalho-1.doxy
