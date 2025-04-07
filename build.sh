#!/bin/bash

# Compile the program
g++ -o apa4 main.cpp apa.cpp bedpe_builder.cpp -lz

# Example usage:
# ./apa4 forward.bed reverse.bed 5000 50000 data.hicslice
# ./apa4 -both-intra-inter forward.bed reverse.bed 5000 50000 data.hicslice
# ./apa4 -only-inter forward.bed reverse.bed 5000 50000 data.hicslice

