# LAST EDIT: Fri Sep 20 15:09:27 2002 by Andreas Loebel (opt0.zib.de) 
# $Id: Makefile,v 1.3 2002/09/20 13:09:40 bzfloebe Exp $

CFLAGS = -O6 -g -Wall -lgomp -fopenmp
CC     = gcc
VPATH  = ../src
 

LIB_MOD    = common arith inout log
valid_MOD  = valid porta four_mot portsort largecalc mp $(LIB_MOD)
xporta_MOD = xporta porta four_mot portsort largecalc mp $(LIB_MOD)
ALL_MOD    = $(valid_MOD) $(xporta_MOD)

libporta_MOD = libporta $(ALL_MOD)
test_MOD = test libporta porta four_mot portsort largecalc mp common arith inout log

all: xporta valid test

valid: $(valid_MOD:%=%.d) $(valid_MOD:%=%.o) 
	$(CC) -fopenmp -o bin/$@ $(valid_MOD:%=%.o)

xporta: $(xporta_MOD:%=%.d) $(xporta_MOD:%=%.o) 
#	$(CC) -ggdb -fopenmp -o bin/$@ $(xporta_MOD:%=%.o)
	$(CC) -fopenmp -o bin/$@ $(xporta_MOD:%=%.o)

#libporta: $(libporta_MOD:%=%.o) 
#	$(CC) -fopenmp -o bin/$@ $(libporta_MOD:%=%.o)

test: $(test_MOD:%=%.o)
	g++ -fopenmp -o bin/$@ $(test_MOD:%=%.o)
	#c++ -o bin/$@ $(test_MOD:%=%.o) -ggdb -pthread -lstdc++ -L ../efence -lefence


clean:
	rm -f *.o bin/xporta bin/valid

%.d: %.c
	$(SHELL) -ec 'gcc -MM $(CPPFLAGS) $< \
		| sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
		[ -s $@ ] || rm -f $@'

include $(ALL_MOD:%=%.d)

