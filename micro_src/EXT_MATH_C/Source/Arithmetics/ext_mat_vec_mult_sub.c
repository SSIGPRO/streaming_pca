#include "ext_math_c.h"

void ext_mat_vec_mult_sub(
    const arm_matrix_instance_f32 * pSrcMat,
                    const FLOAT_T * pVec,
                          FLOAT_T * pDst)
{
    uint32_t numRows = pSrcMat->numRows;
    uint32_t numCols = pSrcMat->numCols;
    const FLOAT_T *pSrcA = pSrcMat->pData;
    const FLOAT_T *pInA1;      /* input data matrix pointer A of Q31 type */
    const FLOAT_T *pInA2;      /* input data matrix pointer A of Q31 type */
    const FLOAT_T *pInA3;      /* input data matrix pointer A of Q31 type */
    const FLOAT_T *pInA4;      /* input data matrix pointer A of Q31 type */
    const FLOAT_T *pInVec;     /* input data matrix pointer B of Q31 type */
    FLOAT_T *px;               /* Temporary output data matrix pointer */
    uint16_t i, row, colCnt; /* loop counters */
    FLOAT_T matData, matData2, vecData, vecData2;


    /* Process 4 rows at a time */
    row = numRows >> 2;
    i = 0u;
    px = pDst;

    /* The following loop performs the dot-product of each row in pSrcA with the
     * vector */
    /* row loop */
    while (row > 0) {
        /* For every row wise process, the pInVec pointer is set
         ** to the starting address of the vector */
        pInVec = pVec;

        /* Initialize accumulators */
        FLOAT_T sum1 = 0.0f;
        FLOAT_T sum2 = 0.0f;
        FLOAT_T sum3 = 0.0f;
        FLOAT_T sum4 = 0.0f;

        /* Loop unrolling: process 2 columns per iteration */
        colCnt = numCols;

        /* Initialize pointers to the starting address of the column being
         * processed */
        pInA1 = pSrcA + i;
        pInA2 = pInA1 + numCols;
        pInA3 = pInA2 + numCols;
        pInA4 = pInA3 + numCols;


        // Main loop: matrix-vector multiplication
        while (colCnt > 0u) {
            // Read 2 values from vector
            vecData = *(pInVec)++;
            /* Read 8 values from the matrix - 2 values from each of 4 rows, and
             * do multiply accumulate */
            matData = *(pInA1)++;
            sum1 += matData * vecData;
            matData = *(pInA2)++;
            sum2 += matData * vecData;
            matData = *(pInA3)++;
            sum3 += matData * vecData;
            matData = *(pInA4)++;
            sum4 += matData * vecData;

            // Decrement the loop counter
            colCnt--;
        }

        /* Saturate and subtract the result from the destination buffer */
        *px++ -= sum1;
        *px++ -= sum2;
        *px++ -= sum3;
        *px++ -= sum4;

        i = i + numCols * 4;

        /* Decrement the row loop counter */
        row--;
    }

    /* process any remaining rows */
    row = numRows & 3u;
    while (row > 0) {

        FLOAT_T sum = 0.0f;
        pInVec = pVec;
        pInA1 = pSrcA + i;

        colCnt = numCols >> 1;
        while (colCnt > 0) {
            vecData = *(pInVec)++;
            vecData2 = *(pInVec)++;
            matData = *(pInA1)++;
            matData2 = *(pInA1)++;
            sum += matData * vecData;
            sum += matData2 * vecData2;
            colCnt--;
        }
        // process remainder of row
        colCnt = numCols & 1u;


        while (colCnt > 0) {
            sum += *pInA1++ * *pInVec++;
            colCnt--;
        }

        *px++ -= sum;
        i = i + numCols;
        row--;
    }
}
