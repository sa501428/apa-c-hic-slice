#!/bin/bash

# Compile with C++11 support and all necessary warnings
g++ -std=c++11 -Wall -Wextra -o apa4 main.cpp apa.cpp apa_matrix.cpp vector_tools.cpp bedpe_builder.cpp -lz

# Example usage:
# ./apa4 intra 5000 50000 100 data.hicslice forward.bed reverse.bed output.txt
# ./apa4 inter 0 0 100 data.hicslice forward.bed reverse.bed output.txt

