CC=			gcc
CFLAGS=		-g -Wall -Wc++-compat -O3
CPPFLAGS=
INCLUDES=
OBJS=		libsais.o libsais64.o kalloc.o kthread.o misc.o io.o rld0.o bre.o rle.o rope.o mrope.o \
			dawg.o fm-index.o ssa.o sais-ss.o build.o search.o bwa-sw.o
PROG=		ropebwt3
LIBS=		-lpthread -lz -lm

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address -ldl
endif

ifneq ($(omp),0)
	CPPFLAGS=-DLIBSAIS_OPENMP
	CFLAGS+=-fopenmp
	LIBS+=-fopenmp
endif

.SUFFIXES:.c .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

ropebwt3:$(OBJS) main.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

rld0.o:rld0.c rld0.h bre.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DRLD_HAVE_BRE $(INCLUDES) $< -o $@

clean:
		rm -fr *.o a.out $(PROG) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(CPPFLAGS) -- *.c)

# DO NOT DELETE

bre.o: bre.h
build.o: rb3priv.h fm-index.h rld0.h mrope.h rope.h io.h rle.h bre.h ketopt.h
build.o: kthread.h
bwa-sw.o: rb3priv.h fm-index.h rld0.h mrope.h rope.h io.h align.h kalloc.h
bwa-sw.o: dawg.h khashl-km.h ksort.h
dawg.o: dawg.h kalloc.h libsais.h io.h rb3priv.h khashl-km.h
fm-index.o: rb3priv.h fm-index.h rld0.h mrope.h rope.h io.h rle.h kthread.h
fm-index.o: kalloc.h khashl-km.h
io.o: rb3priv.h io.h kseq.h
kalloc.o: kalloc.h
kthread.o: kthread.h
libsais.o: libsais.h
libsais64.o: libsais.h libsais64.h
main.o: rb3priv.h fm-index.h rld0.h mrope.h rope.h io.h ketopt.h
misc.o: rb3priv.h
mrope.o: mrope.h rope.h rle.h
rld0.o: rld0.h
rle.o: rle.h
rope.o: rle.h rope.h
sais-ss.o: rb3priv.h libsais.h libsais64.h
search.o: fm-index.h rb3priv.h rld0.h mrope.h rope.h io.h align.h ketopt.h
search.o: kthread.h kalloc.h
ssa.o: rb3priv.h fm-index.h rld0.h mrope.h rope.h io.h kalloc.h kthread.h
ssa.o: ketopt.h ksort.h
