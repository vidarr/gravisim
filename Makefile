CC=gcc
CFLAGS=-g -Wstrict-prototypes -ansi
LNFLAGS= -lm

USE_GSL=1

ifeq ($(USE_GSL),1)
	LNFLAGS+=-lgsl -lgslcblas
	CFLAGS+=-DUSE_GSL
endif

.PHONY: all

all: gravisim

gravisim: gravisim.c
	$(CC) $(LNFLAGS) $(CFLAGS) $(DEBUG_FLAGS) -o gravisim $^
