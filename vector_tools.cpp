#include "vector_tools.h"

namespace vector_tools {
    float getAverage(const std::vector<float>& vec) {
        float sum = 0.0f;
        int count = 0;
        for (float val : vec) {
            if (val > 0) {
                sum += val;
                count++;
            }
        }
        return count > 0 ? sum / count : 1.0f;  // Return 1.0 instead of 0 to prevent division by zero
    }

    void scaleByAverage(std::vector<float>& vec) {
        float avg = getAverage(vec);
        if (avg > 0) {
            for (float& val : vec) {
                val /= avg;
            }
        }
    }
} 