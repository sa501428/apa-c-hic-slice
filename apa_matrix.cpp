#include "apa_matrix.h"
#include <fstream>
#include <iomanip>
#include <stdexcept>

APAMatrix::APAMatrix(int size) : width(size) {
    if (size <= 0) {
        throw std::runtime_error("APAMatrix size must be positive");
    }
    matrix.resize(size, std::vector<float>(size, 0.0f));
}

void APAMatrix::add(int relX, int relY, float value) {
    if (relX >= 0 && relX < width && relY >= 0 && relY < width) {
        matrix[relX][relY] += value;
    }
}

void APAMatrix::normalize(const std::vector<float>& rowSums, const std::vector<float>& colSums) {
    std::vector<std::vector<float>> normalized(width, std::vector<float>(width, 0.0f));
    
    for (int r = 0; r < width; ++r) {
        for (int c = 0; c < width; ++c) {
            float normVal = rowSums[r] * colSums[c];
            if (normVal > 0) {
                normalized[r][c] = matrix[r][c] / normVal;
            }
            // If normVal is 0, leave normalized[r][c] as 0
        }
    }
    
    matrix = std::move(normalized);
}

void APAMatrix::save(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Cannot open output file: " + filename);
    }

    // Write as space-separated values with high precision
    out << std::fixed << std::setprecision(6);
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < width; j++) {
            if (j > 0) out << " ";
            out << matrix[i][j];
        }
        out << "\n";
    }
} 