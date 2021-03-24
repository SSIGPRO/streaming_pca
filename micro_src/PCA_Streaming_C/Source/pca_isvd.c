#include "pca_streaming_c.h"

#define BETA 0.9

void pca_isvd(
                FLOAT_T * pSrcx,          /* n sized vector x (MUST BE
                                           * contiguous after pSrcUT memory) */
    ARM_MATRIX_INSTANCE * pSrcUT,         /* m x n matrix UT (U transposed) */
                FLOAT_T * pSrcs,          /* m+1 sized vector s */
                FLOAT_T * p_mp1xmp1_buff, /* m+1 x m+1 sized buffer vector */
                FLOAT_T * p_mp1_buff)     /* m+1 sized buffer vector */
{
    uint16_t n = pSrcUT->numCols, m = pSrcUT->numRows;
    uint16_t mp1 = m+1;

    FLOAT_T * pUT = pSrcUT->pData;
    FLOAT_T * px = pSrcx;
    FLOAT_T * ps = pSrcs;
    FLOAT_T * pw = p_mp1_buff;
    FLOAT_T * pr = pSrcx;
    FLOAT_T * pK = p_mp1xmp1_buff;
    FLOAT_T * pUh = p_mp1xmp1_buff;

    ARM_MATRIX_INSTANCE K, Uh, UhTnolastrow, UrT;
    ARM_MAT_INIT(&K, mp1, mp1, pK);
    ARM_MAT_INIT(&Uh, mp1, mp1, pUh);
    ARM_MAT_INIT(&UhTnolastrow, m, mp1, pUh);
    ARM_MAT_INIT(&UrT, mp1, n, pUT);

    FLOAT_T * ptempA;
    FLOAT_T * ptempB;
    FLOAT_T rnorm;
    uint16_t i;


    ext_mat_vec_mult(pSrcUT, px, pw); // w = UT@x
    ext_vec_mat_mult_sub(pw, pSrcUT, pr); // r = x - w@UT (x and r share memory)


    /* KT build *******************************************/

    ext_mat_zeros(&K); // fill with zeros

    // evaluate rnorm and normalize r
    ARM_DOT_PROD(pr, pr, n, &rnorm);
    ARM_SQRT(rnorm, &rnorm);
    ARM_SCALE(pr, 1.0/rnorm, pr, n);

    // put s on the diagonal
    ptempA = ps;
    ptempB = pK;
    for(i = 0; i < m; ++i)
    {
        *ptempB = BETA * *ptempA++;
        ptempB += mp1+1; // next row, next column
    }

    // put w on the last column
    ptempA = pw;
    ptempB = pK + m; // select last column
    for(i = 0; i < m; ++i)
    {
        *ptempB = *ptempA++;
        ptempB += mp1; // next row
    }

    *ptempB = rnorm; // put rnorm on the last element of the matrix

    /* KT build END ***************************************/

    // Note: mem(K) == mem(Uh)
    ext_mat_svd_compact(&Uh, ps, p_mp1_buff); /* SVD(K) -> Uh, s (mem(K) =
                                               * mem(Uh); VhT not needed) */

    ext_mat_square_transpose_overwrite(&Uh); // Uh = Uh.T

    ext_mat_mult_overwriteB(&UhTnolastrow, &UrT, p_mp1_buff);
    /* UT = UhT(without last row)@UrT */

    // UT overwrites m rows of UrT
}
