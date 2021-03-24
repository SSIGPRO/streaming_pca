#include "ext_math_c.h"

void ext_mat_mult_overwriteA(
    ARM_MATRIX_INSTANCE * pSrcA, /* n x m input matrix, result will be stored here */
    ARM_MATRIX_INSTANCE * pSrcB, /* m x m square input matrix */
                FLOAT_T * buff)  /* m sized buffer */
{
    FLOAT_T *pInA = pSrcA->pData, *pB = pSrcB->pData, *pOutA = pSrcA->pData;
    FLOAT_T *p_buff = buff;
    uint16_t n = pSrcA->numRows, m = pSrcA->numCols, l = pSrcB->numCols;

    uint16_t i, j, k;
    uint16_t numCols4 = l - (l % 2);
    uint16_t numCols4space = numCols4 >> 1;

    for(i = 0; i < n; ++i)
    {
        FLOAT_T *pInB_0 = pB;
        FLOAT_T *pInB_1 = pInB_0+numCols4space;
        FLOAT_T *pOutA_0 = pOutA;
        FLOAT_T *pOutA_1 = pOutA_0+numCols4space;

        // copy row of A in buffer
        FLOAT_T *pInA_cpy = pInA, *p_buff_cpy = p_buff;
        for(j = 0; j < m; ++j)
        {
            *p_buff_cpy++ = *pInA_cpy++;
        }

        for(j = 0; j < numCols4space; ++j)
        {
            FLOAT_T sum0 = 0.0, sum1 = 0.0;
            FLOAT_T *pIn1 = p_buff;
            FLOAT_T *pIn2_0 = pInB_0;
            FLOAT_T *pIn2_1 = pInB_1;

            for(k = 0; k < m; ++k)
            {
                sum0 += *pIn1**pIn2_0;
                sum1 += *pIn1++**pIn2_1;
                pIn2_0 += l;
                pIn2_1 += l;
            }
            *pOutA_0++ = sum0;
            *pOutA_1++ = sum1;
            ++pInB_0;
            ++pInB_1;
        }
        for(j = 0; j < l - numCols4; ++j)
        {
            FLOAT_T sum = 0.0;
            FLOAT_T *pIn1 = p_buff;
            FLOAT_T *pIn2 = pInB_1;

            for(k = 0; k < m; ++k)
            {
                sum += *pIn1++**pIn2;
                pIn2 += l;
            }
            *pOutA_1++ = sum;
            ++pInB_1;
        }
        pInA += m;
        pOutA += l;
    }
}
