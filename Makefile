ifdef CLANG
CC = clang++
else
CC = g++
endif

CPPFLAGS = -std=c++17 -Wall -Wextra -Werror

INCLUDE_PATH = -Iparlaylib/include/

ifdef CILKPLUS
CC = clang++
CPPFLAGS += -DPARLAY_CILKPLUS -DCILK -fcilkplus
else ifdef OPENCILK
CPPFLAGS += -DPARLAY_OPENCILK -DCILK -fopencilk
else ifdef SERIAL
CPPFLAGS += -DPARLAY_SEQUENTIAL
else
CPPFLAGS += -pthread
endif

ifdef DEBUG
CPPFLAGS += -DDEBUG -Og
else ifdef PERF
CC = g++
CPPFLAGS += -Og -mcx16 -march=native -g
else ifdef MEMCHECK
CPPFLAGS += -Og -mcx16 -DPARLAY_SEQUENTIAL
else
CPPFLAGS += -O3 -mcx16 -march=native
endif

ifdef STDALLOC
CPPFLAGS += -DPARLAY_USE_STD_ALLOC
endif

all: sssp

sssp:	sssp.cc sssp.h dijkstra.h graph.h utils.h
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) sssp.cc -o sssp

symmetrize:	symmetrize.cc graph.h
	$(CC) $(CPPFLAGS) $(INCLUDE_PATH) symmetrize.cc -o symmetrize

clean:
	rm sssp symmetrize

