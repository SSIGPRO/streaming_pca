#include "ext_math_c.h"

void ext_mat_inverse_triu(
    const ARM_MATRIX_INSTANCE * pSrcA,
          ARM_MATRIX_INSTANCE * pDstB)
{
    uint16_t n = pSrcA->numRows;

    FLOAT_T *pA = pSrcA->pData;
    FLOAT_T *pB = pDstB->pData;
    FLOAT_T * pIn;
    FLOAT_T * pOut;
    FLOAT_T sum;
    uint16_t i, j, k;

    // Evaluate inverse
    for(i = 0; i < n; ++i) // row loop
    {
        // Evaluate row at the right of the diagonal
        for(j = 0; j < n-i; ++j) // column loop
        {
            pIn = pA + j; // select input column
            pOut = pB;

            sum = (j == 0) ? -1.0 : 0.0f;

            for(k = 0; k < j; ++k)
            {
                sum += *pOut * *pIn;
                ++pOut;
                pIn += n;
            }

            *pOut = -sum / *pIn;
        }
        pA += n+1;
        pB += n+1;
    }
}
