#include "ext_math_c.h"

void ext_mat_ones(ARM_MATRIX_INSTANCE * pDstA)
{
    uint16_t n = pDstA->numRows, m = pDstA->numCols;
    FLOAT_T * pA = pDstA->pData;
    ARM_FILL(1.0, pA, n*m);
}
