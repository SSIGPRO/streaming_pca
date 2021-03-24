#include "ext_math_c.h"


void ext_mat_mult(
    ARM_MATRIX_INSTANCE * pSrcA, /* n x m input matrix */
    ARM_MATRIX_INSTANCE * pSrcB, /* m x l input matrix */
    ARM_MATRIX_INSTANCE * pDstC) /* n x l output matrix */
{
    FLOAT_T *pInA = pSrcA->pData, *pB = pSrcB->pData, *pOutC = pDstC->pData;
    uint16_t n = pSrcA->numRows, m = pSrcA->numCols, l = pSrcB->numCols;

    uint16_t i, j, k;
    uint16_t numCols4 = l - (l % 2);
    uint16_t numCols4space = numCols4 >> 1;

    for(i = 0; i < n; ++i)
    {
        FLOAT_T *pInB_0 = pB;
        FLOAT_T *pInB_1 = pInB_0+numCols4space;
        FLOAT_T *pOutC_0 = pOutC;
        FLOAT_T *pOutC_1 = pOutC_0+numCols4space;

        for(j = 0; j < numCols4space; ++j)
        {
            FLOAT_T sum0 = 0.0, sum1 = 0.0;
            FLOAT_T *pIn1 = pInA;
            FLOAT_T *pIn2_0 = pInB_0;
            FLOAT_T *pIn2_1 = pInB_1;

            for(k = 0; k < m; ++k)
            {
                sum0 += *pIn1**pIn2_0;
                sum1 += *pIn1**pIn2_1;

                pIn2_0 += l;
                pIn2_1 += l;

                ++pIn1;
            }
            *pOutC_0++ = sum0;
            *pOutC_1++ = sum1;
            ++pInB_0;
            ++pInB_1;
        }
        for(j = 0; j < l - numCols4; ++j)
        {
            FLOAT_T sum = 0.0;
            FLOAT_T *pIn1 = pInA;
            FLOAT_T *pIn2 = pInB_1;

            for(k = 0; k < m; ++k)
            {
                sum += *pIn1++**pIn2;
                pIn2 += l;
            }
            *pOutC_1++ = sum;
            ++pInB_1;
        }
        pInA += m;
        pOutC += l;
    }
}
