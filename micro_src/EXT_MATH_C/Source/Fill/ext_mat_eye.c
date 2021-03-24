#include "ext_math_c.h"

void ext_mat_eye(ARM_MATRIX_INSTANCE * pDstA)
{
    uint16_t n = pDstA->numRows, m = pDstA->numCols;
    FLOAT_T *pA = pDstA->pData;
    uint16_t i, tot_rows, rows;

    tot_rows = MIN(n, m);

    // put ones on the diagonal

#if defined (ARM_MATH_LOOPUNROLL)

    rows = tot_rows >> 2U;

    for(i = 0; i < rows; ++i) // insert ones on the diagonal
    {
        *pA = 1.0;
        pA += m+1;

        *pA = 1.0;
        pA += m+1;

        *pA = 1.0;
        pA += m+1;

        *pA = 1.0;
        pA += m+1;
    }

    rows = tot_rows % 0x4U; // remaining columns

#else

    rows = tot_rows;

#endif

    for(i = 0; i < rows; ++i) // insert ones on the diagonal
    {
        *pA = 1.0;
        pA += m+1;
    }
}
