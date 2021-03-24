#include "ext_math_c.h"

void ext_mat_xtx_mult(
  const ARM_MATRIX_INSTANCE * pSrcA,
        ARM_MATRIX_INSTANCE * pDstB)
{
    FLOAT_T *pInA = pSrcA->pData, *pOutB = pDstB->pData;
    uint16_t n = pSrcA->numRows, m = pSrcA->numCols;
    uint16_t i, j, k;

    for(i = 0; i < m; ++i)
    {
        FLOAT_T *pInA1 = pInA;
        FLOAT_T *pInA2 = pInA;
        FLOAT_T *pOutB1 = pOutB, *pOutB2 = pOutB;
        for(j = 0; j < m-i; ++j)
        {
            FLOAT_T sum = 0.0;
            FLOAT_T *pIn1 = pInA1, *pIn2 = pInA2;
            for(k = 0; k < n; ++k)
            {
                sum += *pIn1**pIn2;
                pIn1 += m;
                pIn2 += m;
            }
            *pOutB1++ = sum;
            if(j > 0)*pOutB2 = sum;
            pOutB2 += m;

            ++pInA2;
        }
        ++pInA;
        pOutB += m+1;
    }
}
