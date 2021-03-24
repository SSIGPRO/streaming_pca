#include "ext_math_c.h"

void ext_mat_square_transpose_overwrite(
    ARM_MATRIX_INSTANCE * pSrcA) /* n x n square matrix, transpose will be saved here */
{
    uint16_t n = pSrcA->numRows;
    FLOAT_T *pA1 = pSrcA->pData+1, *pA2 = pSrcA->pData+n;
    uint16_t i, j;

    for(i = 1; i < n; ++i)
    {
        FLOAT_T *p1 = pA1;
        FLOAT_T *p2 = pA2;

        for(j = i; j < n; ++j)
        {
            FLOAT_T buff = *p1;
            *p1++ = *p2;
            *p2 = buff;
            p2 += n;
        }

        pA1 += n+1;
        pA2 += n+1;
    }
}
