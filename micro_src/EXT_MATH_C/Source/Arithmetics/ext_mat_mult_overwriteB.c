#include "ext_math_c.h"

void ext_mat_mult_overwriteB(
    ARM_MATRIX_INSTANCE * pSrcA,  /* n x m input matrix */
    ARM_MATRIX_INSTANCE * pSrcB,  /* m x l input matrix, result will be stored here (NOTE works iff n >= m !!!) */
                FLOAT_T * buff)   /* m sized buffer */
{
    FLOAT_T *pA = pSrcA->pData, *pInB = pSrcB->pData, *pOutB = pSrcB->pData;
    FLOAT_T * p_buff = buff;
    uint16_t n = pSrcA->numRows, m = pSrcA->numCols, l = pSrcB->numCols;

    uint16_t i, j, k;
    uint16_t numRows4 = n - (n % 2);
    uint16_t numRows4space = numRows4 >> 1;

    for(i = 0; i < l; ++i)
    {
        FLOAT_T *pInA_0 = pA;
        FLOAT_T *pInA_1 = pInA_0+m*numRows4space;
        FLOAT_T *pOutB_0 = pOutB;
        FLOAT_T *pOutB_1 = pOutB_0+l*numRows4space;

        // copy B column in buffer
        FLOAT_T *pInB_cpy = pInB, *p_buff_cpy = p_buff;
        for(j = 0; j < m; ++j)
        {
            *p_buff_cpy++ = *pInB_cpy;
            pInB_cpy += l;
        }

        for(j = 0; j < numRows4space; ++j)
        {
            FLOAT_T sum0 = 0.0, sum1 = 0.0;
            FLOAT_T *pIn1_0 = pInA_0;
            FLOAT_T *pIn1_1 = pInA_1;
            FLOAT_T *pIn2 = p_buff;

            for(k = 0; k < m; ++k)
            {
                sum0 += *pIn1_0**pIn2;
                sum1 += *pIn1_1**pIn2;

                ++pIn1_0;
                ++pIn1_1;

                ++pIn2;
            }
            *pOutB_0 = sum0;
            *pOutB_1 = sum1;
            pOutB_0 += l;
            pOutB_1 += l;
            pInA_0 += m;
            pInA_1 += m;
        }
        for(j = 0; j < n - numRows4; ++j)
        {
            FLOAT_T sum = 0.0;
            FLOAT_T *pIn1 = pInA_1;
            FLOAT_T *pIn2 = p_buff;

            for(k = 0; k < m; ++k)
            {
                sum += *pIn1++**pIn2++;
            }
            *pOutB_1 = sum;
            pOutB_1 += l;
            pInA_1 += m;
        }
        ++pInB;
        ++pOutB;
    }
}
