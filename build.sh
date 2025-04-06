g++ -o bed_slice_analyzer bed_slice_analyzer.cpp bedpe_builder.cpp hic_slice_reader.cpp -lz

# ./bed_slice_analyzer forward.bed reverse.bed 5000 50000 data.hicslice
# ./bed_slice_analyzer -both-intra-inter forward.bed reverse.bed 5000 50000 data.hicslice

