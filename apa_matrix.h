#ifndef APA_MATRIX_H
#define APA_MATRIX_H

#include <vector>
#include <string>
#include "vector_tools.h"

// Structure to hold APA matrix
class APAMatrix {
public:
    APAMatrix(int size);
    
    void add(int relX, int relY, float value);
    void normalize(const std::vector<float>& rowSums, const std::vector<float>& colSums);
    void save(const std::string& filename) const;

private:
    std::vector<std::vector<float>> matrix;
    int width;
};

#endif 