#include "ext_math_c.h"


void ext_mat_mult_simple(
    ARM_MATRIX_INSTANCE * pSrcA, /* n x m input matrix */
    ARM_MATRIX_INSTANCE * pSrcB, /* m x l input matrix */
    ARM_MATRIX_INSTANCE * pDstC) /* n x l output matrix */
{
    FLOAT_T *pInA = pSrcA->pData, *pB = pSrcB->pData, *pOut = pDstC->pData;
    uint16_t n = pSrcA->numRows, m = pSrcA->numCols, l = pSrcB->numCols;
    uint16_t i, j, k;

    for(i = 0; i < n; ++i)
    {
        FLOAT_T *pInB = pB;
        for(j = 0; j < l; ++j)
        {
            FLOAT_T sum = 0.0;
            FLOAT_T *pIn1 = pInA, *pIn2 = pInB;
            for(k = 0; k < m; ++k)
            {
                sum += *pIn1++**pIn2;
                pIn2 += l;
            }
            *pOut++ = sum;
            ++pInB;
        }
        pInA += m;
    }
}
