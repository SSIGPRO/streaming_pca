#include "ext_math_c.h"


void ext_mat_mult2(
    ARM_MATRIX_INSTANCE * pSrcA, /* n x m input matrix */
    ARM_MATRIX_INSTANCE * pSrcB, /* m x l input matrix */
    ARM_MATRIX_INSTANCE * pDstC) /* n x l output matrix */
{
    FLOAT_T *pA = pSrcA->pData, *pInB = pSrcB->pData, *pOutC = pDstC->pData;
    uint16_t n = pSrcA->numRows, m = pSrcA->numCols, l = pSrcB->numCols;

    uint16_t i, j, k;
    uint16_t numRows4 = n - (n % 2);
    uint16_t numRows4space = numRows4 >> 1;

    for(i = 0; i < l; ++i)
    {
        FLOAT_T *pInA_0 = pA;
        FLOAT_T *pInA_1 = pInA_0+m*numRows4space;
        FLOAT_T *pOutC_0 = pOutC;
        FLOAT_T *pOutC_1 = pOutC_0+l*numRows4space;

        for(j = 0; j < numRows4space; ++j)
        {
            FLOAT_T sum0 = 0.0, sum1 = 0.0;
            FLOAT_T *pIn1_0 = pInA_0;
            FLOAT_T *pIn1_1 = pInA_1;
            FLOAT_T *pIn2 = pInB;

            for(k = 0; k < m; ++k)
            {
                sum0 += *pIn1_0**pIn2;
                sum1 += *pIn1_1**pIn2;

                ++pIn1_0;
                ++pIn1_1;

                pIn2 += l;
            }
            *pOutC_0 = sum0;
            *pOutC_1 = sum1;
            pOutC_0 += l;
            pOutC_1 += l;
            pInA_0 += m;
            pInA_1 += m;
        }
        for(j = 0; j < n - numRows4; ++j)
        {
            FLOAT_T sum = 0.0;
            FLOAT_T *pIn1 = pInA_1;
            FLOAT_T *pIn2 = pInB;

            for(k = 0; k < m; ++k)
            {
                sum += *pIn1++**pIn2;
                pIn2 += l;
            }
            *pOutC_1 = sum;
            pOutC_1 += l;
            pInA_1 += m;
        }
        ++pInB;
        ++pOutC;
    }
}
