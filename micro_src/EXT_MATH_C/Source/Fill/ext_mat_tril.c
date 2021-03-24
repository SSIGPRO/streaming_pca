#include "ext_math_c.h"

void ext_mat_tril(ARM_MATRIX_INSTANCE *pSrcA)
{
    uint16_t n = pSrcA->numRows, m = pSrcA->numCols;
    FLOAT_T *pA = pSrcA->pData;
    uint16_t i, j, tot_cols, cols;

    for(i = 0; i < n; ++i) // for each row
    {
        tot_cols = MIN(i, m);
        pA += i+1; // start from the diagonal

#if defined (ARM_MATH_LOOPUNROLL)

        cols = tot_cols >> 2U;

        for(j = 0; j < cols; ++j) // insert zeros until diagonal is reached
        {
            *pA++ = 0.0;
            *pA++ = 0.0;
            *pA++ = 0.0;
            *pA++ = 0.0;
        }

        cols = tot_cols % 0x4U; // remaining columns

#else

        cols = tot_cols;

#endif

        for(j = 0; j < cols; ++j) // insert zeros until diagonal is reached
        {
            *pA++ = 0.0;
        }
    }
}
