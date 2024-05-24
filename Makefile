CC=			nvcc
CFLAGS=		-O2 -Xcompiler -Wunused-result
OBJS=		kthread.o misc.o bseq.o sdust.o index.o launch_data.o map.o#sketch.o
PROG=		minimap
PROG_EXTRA=	sdust minimap-lite
LIBS=		-lm -lz -lpthread

.SUFFIXES:.cu .o

.cu.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

extra:all $(PROG_EXTRA)

minimap:main.o libminimap.a
		$(CC) $(CFLAGS) $< -o $@ -L. -lminimap $(LDFLAGS) $(LIBS)

minimap-lite:example.o libminimap.a
		$(CC) $(CFLAGS) $< -o $@ -L. -lminimap $(LDFLAGS) $(LIBS)

libminimap.a:$(OBJS)
		$(AR) -csru $@ $(OBJS)

sdust:sdust.cu kdq.h kvec.h kseq.h sdust.h
		$(CC) -D_SDUST_MAIN $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@ $(LDFLAGS) -lz

clean:
		rm -fr gmon.out *.o a.out $(PROG) $(PROG_EXTRA) *~ *.a *.dSYM session*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

bseq.o: bseq.h kseq.h
example.o: minimap.h bseq.h kseq.h
index.o: minimap.h bseq.h kvec.h khash.h
main.o: minimap.h bseq.h
#cuda_map.o: bseq.h kvec.h minimap.h sdust.h ksort.h
map.o: bseq.h kvec.h minimap.h sdust.h ksort.h
misc.o: minimap.h bseq.h ksort.h
sdust.o: kdq.h kvec.h sdust.h
#sketch.o: kvec.h minimap.h bseq.h
