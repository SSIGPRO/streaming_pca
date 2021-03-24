#include "ext_math_c.h"

arm_status ext_mat_cholesky_t(
  const ARM_MATRIX_INSTANCE * pSrc,
        ARM_MATRIX_INSTANCE * pDst)
{

  arm_status status;


#ifdef ARM_MATH_MATRIX_CHECK

  /* Check for matrix mismatch condition */
  if ((pSrc->numRows != pSrc->numCols) ||
      (pDst->numRows != pDst->numCols) ||
      (pSrc->numRows != pDst->numRows)   )
  {
    /* Set status as ARM_MATH_SIZE_MISMATCH */
    status = ARM_MATH_SIZE_MISMATCH;
  }
  else

#endif /* #ifdef ARM_MATH_MATRIX_CHECK */

  {
    int i,j,k;
    int n = pSrc->numRows;
    FLOAT_T invSqrtVj;
    FLOAT_T *pA,*pG;

    pA = pSrc->pData;
    pG = pDst->pData;


    for(i=0 ; i < n ; i++)
    {
       for(j=i ; j < n ; j++)
       {
          pG[i * n + j] = pA[j * n + i];

          for(k=0; k < i ; k++)
          {
             pG[i * n + j] = pG[i * n + j] - pG[k * n + i] * pG[k * n + j];
          }
       }

       if (pG[i * n + i] <= 0.0f)
       {
         return(ARM_MATH_DECOMPOSITION_FAILURE);
       }

       invSqrtVj = 1.0f/sqrtf(pG[i * n + i]);
       for(j=i ; j < n ; j++)
       {
         pG[i * n + j] = pG[i * n + j] * invSqrtVj ;
       }
    }

    status = ARM_MATH_SUCCESS;

  }


  /* Return to application */
  return (status);
}
