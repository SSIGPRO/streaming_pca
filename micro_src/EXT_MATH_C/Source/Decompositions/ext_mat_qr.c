#include "ext_math_c.h"

void ext_mat_qr(
    const ARM_MATRIX_INSTANCE * pSrcA, /* n x m input matrix A */
          ARM_MATRIX_INSTANCE * pDstQ, /* n x m output Q matrix */
          ARM_MATRIX_INSTANCE * pDstR) /* m x m output R matrix */
{
    uint16_t m = pSrcA->numCols;

    FLOAT_T *pQ = pDstQ->pData;

    ARM_MATRIX_INSTANCE ATA;
    ARM_MAT_INIT(&ATA, m, m, pQ);

    ext_mat_xtx_mult(pSrcA, &ATA); // MTM = M.T@M

    ext_mat_cholesky_t(&ATA, pDstR); // cholesky(MTM) -> L.T -> R

    ext_mat_inverse_triu_t(pDstR, pDstR); /* R -> R^-1 (save in lower triangular
                                           * part of R) */

    ext_mat_xyt_mult_tril(pSrcA, pDstR, pDstQ); /* Q = M@R^-1 (use lower
                                                 * triangular part of R) */

    ext_mat_inverse_diag(pDstR, pDstR); // Revert diagonal
    ext_mat_triu(pDstR); // Remove inverse
}
