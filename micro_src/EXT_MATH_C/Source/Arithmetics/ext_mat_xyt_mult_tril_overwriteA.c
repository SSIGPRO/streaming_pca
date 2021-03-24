#include "ext_math_c.h"

void ext_mat_xyt_mult_tril_overwriteA(
          ARM_MATRIX_INSTANCE * pSrcA, /* n x m input matrix, result will be
                                        * saved here (overwrite) */
    const ARM_MATRIX_INSTANCE * pSrcB, /* m x m square input matrix */
                      FLOAT_T * buff)  /* m sized buffer */
{
    FLOAT_T *pInA = pSrcA->pData, *pB = pSrcB->pData, *pOut = pSrcA->pData;
    FLOAT_T *p_buff = buff;
    uint16_t n = pSrcA->numRows, m = pSrcA->numCols;
    uint16_t i, j, k;

    for(i = 0; i < n; ++i)
    {
        FLOAT_T *pInB = pB;

        // copy A row in buffer
        FLOAT_T *pInA_cpy = pInA, *p_buff_cpy = p_buff;
        for(j = 0; j < m; ++j)
        {
            *p_buff_cpy++ = *pInA_cpy++;
        }

        for(j = 0; j < m; ++j)
        {
            FLOAT_T sum = 0.0;
            FLOAT_T *pIn1 = p_buff, *pIn2 = pInB;
            for(k = 0; k < j+1; ++k)
            {
                sum += *pIn1++**pIn2++;
            }
            *pOut++ = sum;
            pInB += m;
        }
        pInA += m;
    }
}

