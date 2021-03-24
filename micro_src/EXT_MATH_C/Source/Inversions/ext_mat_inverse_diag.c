#include "ext_math_c.h"

void ext_mat_inverse_diag(
    const ARM_MATRIX_INSTANCE * pSrcA,
          ARM_MATRIX_INSTANCE * pDstB)
{
    uint16_t n = pSrcA->numRows;

    FLOAT_T *pA = pSrcA->pData;
    FLOAT_T *pB = pDstB->pData;

    uint16_t i;

    for(i = 0; i < n; ++i)
    {
        *pB = 1.0 / *pA;
        pA += n+1;
        pB += n+1;
    }
}
