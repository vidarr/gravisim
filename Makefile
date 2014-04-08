CC=gcc
CFLAGS=-g -Wstrict-prototypes -ansi
LNFLAGS= -lm
.PHONY: all

all: gravisim

gravisim: gravisim.c
	$(CC) $(LNFLAGS) $(CFLAGS) $(DEBUG_FLAGS) -o gravisim $^
