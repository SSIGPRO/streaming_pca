#ifndef PCA_STREAMING_C
#define PCA_STREAMING_C


#include "ext_math_c.h"




/* OJA PCA STREAMING METHOD *********************/

void pca_oja(
                FLOAT_T * pSrcx,       /* n-size data vector x */
    ARM_MATRIX_INSTANCE * pSrcU,       /* n x m matrix U  */
                FLOAT_T * p_nxm_buff,  /* n x m sized vector buffer */
                FLOAT_T * p_mxm_buff,  /* m x m sized vector buffer */
                 uint32_t tau);        /* step count   */

void pca_oja_minor(
                FLOAT_T * pSrcx,       /* n-size data vector x */
    ARM_MATRIX_INSTANCE * pSrcU,       /* n x m matrix U  */
                FLOAT_T * p_nxm_buff,  /* n x m sized vector buffer */
                FLOAT_T * p_mxm_buff,  /* m x m sized vector buffer */
                 uint32_t tau);        /* step count   */

/* Low memory requirements version */

/* Note: vector x content will be destroyed
 */

void pca_oja_lowmem(
                FLOAT_T * pSrcx,       /* n-size data vector x */
    ARM_MATRIX_INSTANCE * pSrcU,       /* n x m matrix U  */
                FLOAT_T * p_mxm_buff,  /* m x m sized vector buffer */
                 uint32_t tau);        /* step count   */

void pca_oja_minor_lowmem(
                FLOAT_T * pSrcx,       /* n-size data vector x */
    ARM_MATRIX_INSTANCE * pSrcU,       /* n x m matrix U  */
                FLOAT_T * p_mxm_buff,  /* m x m sized vector buffer */
                 uint32_t tau);        /* step count   */




/* KRASULINA PCA STREAMING METHOD ***************/

/* Note: vector x content will be destroyed
 */

void pca_krasulina(
                FLOAT_T * pSrcx,       /* n-size data vector x */
    ARM_MATRIX_INSTANCE * pSrcU,       /* n x m matrix U  */
                FLOAT_T * p_nxm_buff,  /* n x m sized vector buffer */
                FLOAT_T * p_mxm_buff,  /* m x m sized vector buffer */
                 uint32_t tau);        /* step count */

void pca_krasulina_minor(
                FLOAT_T * pSrcx,       /* n-size data vector x */
    ARM_MATRIX_INSTANCE * pSrcU,       /* n x m matrix U  */
                FLOAT_T * p_nxm_buff,  /* n x m sized vector buffer */
                FLOAT_T * p_mxm_buff,  /* m x m sized vector buffer */
                 uint32_t tau);        /* step count */

/* Low memory requirements version */

void pca_krasulina_lowmem(
                FLOAT_T * pSrcx,       /* n-size data vector x */
    ARM_MATRIX_INSTANCE * pSrcU,       /* n x m matrix U  */
                FLOAT_T * p_mxm_buff,  /* m x m sized vector buffer */
                uint32_t tau);         /* step count */

void pca_krasulina_minor_lowmem(
                FLOAT_T * pSrcx,       /* n-size data vector x */
    ARM_MATRIX_INSTANCE * pSrcU,       /* n x m matrix U  */
                FLOAT_T * p_mxm_buff,  /* m x m sized vector buffer */
                 uint32_t tau);        /* step count */




/* GROUSE PCA STREAMING METHOD ******************/

/* Note: vector x content will be destroyed
 */

void pca_grouse(
                FLOAT_T * pSrcx,       /* n-size data vector x */
    ARM_MATRIX_INSTANCE * pSrcU,       /* n x m matrix U  */
                FLOAT_T * p_nxm_buff,  /* n x m sized vector buffer */
                FLOAT_T * p_m_buff);   /* m sized vector buffer */




/* PAST PCA STREAMING METHOD ********************/

/* Note: vector x content will be destroyed
 * Note: matrix PT content must be preserved from one call to another 
 *       (initialize to identity)
 */

void pca_past(
                FLOAT_T * pSrcx,       /* n-size data vector x */
    ARM_MATRIX_INSTANCE * pSrcU,       /* n x m matrix U  */
    ARM_MATRIX_INSTANCE * pSrcPT,      /* m x m matrix PT  */
                FLOAT_T * p_nxm_buff,  /* n x m sized vector buffer */
                FLOAT_T * p_m_buff);   /* m sized vector buffer */




/* ISVD PCA STREAMING METHOD ********************/

/* Note: vector s content must be preserved from one call to another
 *       (initialize to all zeros)
 * Note: vector x content will be destroyed
 * Note: the memory of vector x MUST BE after and contiguous the memory of
 *       matrix UT
 */

void pca_isvd(
                FLOAT_T * pSrcx,          /* n sized vector x (MUST BE 
                                           * contiguous after pSrcUT memory) */
    ARM_MATRIX_INSTANCE * pSrcUT,         /* m x n matrix UT (U transposed) */
                FLOAT_T * pSrcs,          /* m+1 sized vector s */
                FLOAT_T * p_mp1xmp1_buff, /* m+1 x m+1 sized buffer vector */
                FLOAT_T * p_mp1_buff);    /* m+1 sized buffer vector */



/* HFRANS PCA STREAMING METHOD *******************/

/* Note: vector x content will be destroyed
 */

void pca_hfrans_minor(
                FLOAT_T * pSrcx,       /* n-size data vector x */
    ARM_MATRIX_INSTANCE * pSrcU,       /* n x m matrix U  */
                FLOAT_T * p_n_buff,    /* n sized vector buffer */
                FLOAT_T * p_m_buff);   /* m sized vector buffer */



#endif /* PCA_STREAMING_C */
