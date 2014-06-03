CC=gcc
CFLAGS=-g -Wstrict-prototypes -ansi
LNFLAGS= -lm
LNFLAGS+=-lgsl -lgslcblas
SDL_LNFLAGS+=`sdl-config --cflags --libs`

SRC_DIR=src
BIN_DIR=bin
TEST_DIR=test

.PHONY: all clean

all: $(BIN_DIR)/gravisim

$(BIN_DIR): 
	[ ! -e $(BIN_DIR) ] && mkdir $(BIN_DIR)

$(BIN_DIR)/gravisim.o: $(SRC_DIR)/gravisim.c $(SRC_DIR)/solver_gauss_legendre_2.h
	$(CC) $(CFLAGS) $(DEBUG_FLAGS) -c $(SRC_DIR)/gravisim.c -o $@

$(BIN_DIR)/solver_gauss_legendre_2.o: $(SRC_DIR)/solver_gauss_legendre_2.c $(SRC_DIR)/solver_gauss_legendre_2.h
	$(CC) $(CFLAGS) $(DEBUG_FLAGS) -c $(SRC_DIR)/solver_gauss_legendre_2.c -o $@

$(BIN_DIR)/gravisim: $(BIN_DIR) $(BIN_DIR)/gravisim.o $(BIN_DIR)/solver_gauss_legendre_2.o $(SRC_DIR)/solver_gauss_legendre_2.h
	$(CC) $(LNFLAGS) $(CFLAGS) $(DEBUG_FLAGS) -o $@ $(BIN_DIR)/gravisim.o $(BIN_DIR)/solver_gauss_legendre_2.o

$(BIN_DIR)/viewer.o: $(SRC_DIR)/viewer.c
	$(CC) $(CLFAGS) $(DEBUG_FLAGS) -c $(SRC_DIR)/viewer.c -o $@

$(BIN_DIR)/viewer: $(BIN_DIR) $(BIN_DIR)/viewer.o
	$(CC) $(LNFLAGS) $(SDL_LNFLAGS) -o $@ $(BIN_DIR)/viewer.o

$(TEST_DIR)/exponential: $(TEST_DIR)/exponential.c
	$(CC) $(LNFLAGS) $(CFLAGS) $(DEBUG_FLAGS) -o $@ $(TEST_DIR)/exponential.c

clean:
	[ -e $(BIN_DIR) ] && rm -r $(BIN_DIR)
