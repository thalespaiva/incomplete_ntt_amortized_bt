CC = g++ -std=c++17
FLAGS = -Wall -Wno-unused-variable -Wno-unused-function -Wno-write-strings -Wno-unused-result -funroll-all-loops -march=native -lm
OPT_FLAGS = -O3 -fwhole-program -flto $(FLAGS)
DEBUG_FLAGS = -O0 -g $(FLAGS)
INCLUDE_DIRS = ./src/third-party/hexl/hexl/include/ ./include/
LIB_DIRS = ./src/third-party/hexl/build/hexl/lib 
LD_LIBS = hexl
INCLUDE_FLAGS = $(addprefix -I, $(INCLUDE_DIRS))
LIB_FLAGS = $(addprefix -L, $(LIB_DIRS)) $(addprefix -l, $(LD_LIBS))  -fopenmp
SRC = polynomial.cpp lwe.cpp rlwe.cpp gsw.cpp misc.cpp third-party/misc_tp.cpp \
	  lmk.cpp forward_ntt.cpp clwe_sampling.cpp operation_counters.cpp pre_comp_data.cpp

all: main_amortized

main_amortized: $(addprefix ./src/, $(SRC))  main_amortized_csv.cpp ./src/amortized_bootstrap_general.cpp
	$(CC) -g -o main_amortized $^ $(OPT_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS)

hexl: hexl/build
	cmake --build ./src/third-party/hexl/build

hexl/build:
	cmake -S ./src/third-party/hexl/ -B ./src/third-party/hexl/build

clean: 
	rm --f main_amortized main_rlwe main_non_amortized main_intt
