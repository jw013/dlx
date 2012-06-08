# GNU make

CFLAGS = -std=c99 -pedantic -Wextra -Wall -Wshadow -Wstrict-prototypes -Wcast-align\
	 -Wstrict-aliasing
CTAGS = ctags
MAKEDEPFLAG = -M

OBJ = dlx.o dlx_read.o

all: ${OBJ} dlx_test

dlx_test: dlx_read_test.c ${OBJ}
	${CC} -o $@ ${CFLAGS} $^

%.o: %.c
	${CC} ${CFLAGS} -c $<

depend: *.c
	${CC} ${CFLAGS} ${MAKEDEPFLAG} $^ > $@

tags: *.c
	${CTAGS} $^

clean: 
	-rm -f ${OBJ} dlx_test

.PHONY: clean

include depend
