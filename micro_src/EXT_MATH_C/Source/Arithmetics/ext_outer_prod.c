#include "ext_math_c.h"

void ext_outer_prod(
                FLOAT_T * pSrcA,
                FLOAT_T * pSrcB,
    ARM_MATRIX_INSTANCE * pDstC)
{
    uint16_t n = pDstC->numRows, m = pDstC->numCols;
    FLOAT_T * pIn1 = pSrcA;
    FLOAT_T * pIn2 = pSrcB;
    FLOAT_T * pOut = pDstC->pData;

    uint16_t i, col;

    for(i = 0; i < n; ++i) // loop through rows
    {
#if defined (ARM_MATH_LOOPUNROLL)
        col = m >> 2;
        while(col > 0) // loop through 4 columns
        {
            *pOut++ = *pIn1 * *pIn2++;
            *pOut++ = *pIn1 * *pIn2++;
            *pOut++ = *pIn1 * *pIn2++;
            *pOut++ = *pIn1 * *pIn2++;
            col--;
        }
        col = m % 4;
#else
        col = m;
#endif
        while(col > 0) // loop through remaining columns
        {
            *pOut++ = *pIn1 * *pIn2++;
            col--;
        }
        ++pIn1;
        pIn2 = pSrcB;
    }
}
