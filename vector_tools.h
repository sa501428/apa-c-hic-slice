#ifndef VECTOR_TOOLS_H
#define VECTOR_TOOLS_H

#include <vector>

namespace vector_tools {
    // Get average of non-zero values
    float getAverage(const std::vector<float>& vec);

    // Scale vector by its average
    void scaleByAverage(std::vector<float>& vec);
}

#endif 