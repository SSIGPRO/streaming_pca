#include "ext_math_c.h"

void ext_mat_qr_overwrite(
    ARM_MATRIX_INSTANCE * pSrcQ, /* n x m input matrix A (will be overwritten
                                  * by Q) */
    ARM_MATRIX_INSTANCE * pDstR, /* m x m output matrix R */
                FLOAT_T * buff)  /* m sized buffer */
{
	uint16_t m = pSrcQ->numCols;

	FLOAT_T *pR = pDstR->pData;

	ARM_MATRIX_INSTANCE ATA;
	ARM_MAT_INIT(&ATA, m, m, pR);

    ext_mat_xtx_mult(pSrcQ, &ATA); // MTM = M.T@M

    ext_mat_cholesky_t(&ATA, pDstR); // cholesky(MTM) -> L.T -> R

    ext_mat_inverse_triu_t(pDstR, pDstR); /* R -> R^-1 (save in lower triangular
                                           * part of R) */

    ext_mat_xyt_mult_tril_overwriteA(pSrcQ, pDstR, buff);
    /* Q = M@R^-1 (use lower triangular part of R) */

    ext_mat_inverse_diag(pDstR, pDstR); // Revert diagonal
    ext_mat_triu(pDstR); // Remove inverse
}
