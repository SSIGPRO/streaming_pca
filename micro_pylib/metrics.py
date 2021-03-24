import numpy as np

eps = np.finfo(np.float64).eps
    
def norm_error(A, Ahat):
    return np.linalg.norm(A-Ahat)
    
def norm_error_rel(A, Ahat):
    err = norm_error(A, Ahat)
    if err == 0:
        err = eps
    return np.linalg.norm(A)/err
    
def norm_error_db(A, Ahat):
    rel_error = norm_error_rel(A, Ahat)
    if rel_error == 0:
        rel_error = eps
    return 20*np.log10(rel_error)
    
def print_errors(M, Mhat, name1="A", name2="B", show=False):
    # Print the matrices and the error introduced by the micro
    M = np.array(M)
    Mhat = np.array(Mhat)
    if len(M.shape) == 0:
        object = "Value"
    elif len(M.shape) == 1: 
        object = "Vector"
    elif len(M.shape) == 2: 
        object = "Matrix"
    else:
        object = "Tensor"
    if show:
        print("\n{} {}".format(object, name1))
        print(M)
        print("{} {}".format(object, name2))
        print(Mhat)
    print("{} {} vs {} error".format(object, name1, name2))
    print("    Norm error (abs) = {:.5e}".format(norm_error(M, Mhat)))
    print("    Norm error (dB)  = {:.2f} dB".format(norm_error_db(M, Mhat)))