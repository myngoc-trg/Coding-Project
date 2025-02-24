from scipy.linalg import hessenberg
import numpy as np

def count_negative_eigenvalues(T, shift):
    n = np.shape(T)[0]
    d = T[0,0] - shift
    neg_count = 1 if d<0 else 0
    for i in range(1,n):
        if d!=0:
            d = (T[i,i] - shift) - (T[i,i-1]**2)/d
        else:
            d = np.inf
        if d<0:
            neg_count += 1
    return neg_count

def bijection_algorithm(T,a,b, tol=1e-10, epsilon = 1e-5):
    eigenvalues_A = []
    n = np.shape(T)[0]
    for i in range(n):
        left_end, right_end = a, b
        while right_end - left_end > tol:
            mid_point = (left_end + right_end)/2
            neg_count = count_negative_eigenvalues(T, mid_point)
            if neg_count > i:
                right_end = mid_point
            else:
                left_end = mid_point
        eigenvalue = (left_end + right_end)/2
        if a + epsilon < eigenvalue < b - epsilon:
            eigenvalues_A.append(eigenvalue)
    return sorted(set(eigenvalues_A))

A = np.array([[1,1,0,0],
             [1,0,1,0],
             [0,1,2,1],
             [0,0,1,-1]])

eigval_A, _ = np.linalg.eig(A)
print("All Numpy eigenvalues of tridiagonal matrix A: \n", eigval_A)

a = -2
b = 1
eigenval_bijection = bijection_algorithm(A, a, b)

print(f"Eigenvalues by bijection algorithm of tridiagonal matrix A in {[a,b]}: \n", eigenval_bijection)
print("-------------------------------Testing for Hassenberg matrix -----------------------------------")

n = int(input("Enter the size of the symmetric matrix A: n = "))

random_matrix = np.random.rand(n,n)
symmetric_matrix = 1/2*(random_matrix + random_matrix.T)

hassenberg_matrix = hessenberg(symmetric_matrix)

eigval_hassenberg, _ = np.linalg.eig(hassenberg_matrix)
eigval_symmetric, _ = np.linalg.eig(symmetric_matrix)

print("Symmetric random matrix = \n", symmetric_matrix)
print("Hassenberg form of the symmetric matrix = \n", hassenberg_matrix)
print("All Numpy eigenvalues of symmetric matrix: \n", eigval_hassenberg)
print("All Numpy eigenvalues of Hassenberg matrix: \n", eigval_hassenberg)

a_H = -1
b_H = 2

print(f"Eigenvalues by bijection algorithm of Hassenberg matrix in {[a_H,b_H]}: \n", bijection_algorithm(hassenberg_matrix, a_H, b_H))





                