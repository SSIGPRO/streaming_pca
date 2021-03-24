#include "ext_math_c.h"

void ext_add_outer_prod(
    ARM_MATRIX_INSTANCE * pSrcA, // matrix A
                FLOAT_T * pSrcv, // vector v
                FLOAT_T * pSrcu, // vector a
    ARM_MATRIX_INSTANCE * pDstB) // result = A + outer(v, a)
{
    uint16_t n = pSrcA->numRows, m = pSrcA->numCols;
    FLOAT_T * pIn0 = pSrcA->pData;
    FLOAT_T * pIn1 = pSrcv;
    FLOAT_T * pIn2 = pSrcu;
    FLOAT_T * pOut = pDstB->pData;

    uint16_t i, col;

    for(i = 0; i < n; ++i) // loop through rows
    {
#if defined (ARM_MATH_LOOPUNROLL)
        col = m >> 2;
        while(col > 0) // loop through 4 columns
        {
            *pOut++ = *pIn0++ + *pIn1 * *pIn2++;
            *pOut++ = *pIn0++ + *pIn1 * *pIn2++;
            *pOut++ = *pIn0++ + *pIn1 * *pIn2++;
            *pOut++ = *pIn0++ + *pIn1 * *pIn2++;
            col--;
        }
        col = m % 4;
#else
        col = m;
#endif
        while(col > 0) // loop through remaining columns
        {
            *pOut++ = *pIn0++ + *pIn1 * *pIn2++;
            col--;
        }
        ++pIn1;
        pIn2 = pSrcu;
    }
}
