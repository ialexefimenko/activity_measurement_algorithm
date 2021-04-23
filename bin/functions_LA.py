'''
Функции линейной алгебры
'''
def zeros_matrix(rows, cols):

    M = []
    while len(M) < rows:
        M.append([])
        while len(M[-1]) < cols:
            M[-1].append(0.0)

    return M

def identity_matrix(n):
   
    I = zeros_matrix(n, n)
    for i in range(n):
        I[i][i] = 1.0

    return I

def copy_matrix(M):
   
    rows = len(M); cols = len(M[0])

    MC = zeros_matrix(rows, cols)


    for i in range(rows):
        for j in range(cols):
            MC[i][j] = M[i][j]

    return MC

def print_matrix(M, decimals=3):
  
    for row in M:
        print([round(x,decimals)+0 for x in row])

def transpose(M):
  
    if not isinstance(M[0],list):
        M = [M]

  
    rows = len(M); cols = len(M[0])

  
    MT = zeros_matrix(cols, rows)

  
    for i in range(rows):
        for j in range(cols):
            MT[j][i] = M[i][j]

    return MT

def matrix_addition(A, B):
  
    rowsA = len(A); colsA = len(A[0])
    rowsB = len(B); colsB = len(B[0])
    if rowsA != rowsB or colsA != colsB:
        raise ArithmeticError('Matrices are NOT the same size.')

  
    C = zeros_matrix(rowsA, colsB)

  
    for i in range(rowsA):
        for j in range(colsB):
            C[i][j] = A[i][j] + B[i][j]

    return C

def matrix_subtraction(A, B):
  
    rowsA = len(A); colsA = len(A[0])
    rowsB = len(B); colsB = len(B[0])
    if rowsA != rowsB or colsA != colsB:
        raise ArithmeticError('Matrices are NOT the same size.')

  
    C = zeros_matrix(rowsA, colsB)


    for i in range(rowsA):
        for j in range(colsB):
            C[i][j] = A[i][j] - B[i][j]

    return C

def matrix_multiply(A, B):
  
    rowsA = len(A); colsA = len(A[0])
    rowsB = len(B); colsB = len(B[0])
    if colsA != rowsB:
        raise ArithmeticError(
            'Number of A columns must equal number of B rows.')

 
    C = zeros_matrix(rowsA, colsB)
    for i in range(rowsA):
        for j in range(colsB):
            total = 0
            for ii in range(colsA):
                total += A[i][ii] * B[ii][j]
            C[i][j] = total

    return C

def multiply_matrices(list):
   
    matrix_product = list[0]

   
    for matrix in list[1:]:
        matrix_product = matrix_multiply(matrix_product, matrix)

    return matrix_product

def check_matrix_equality(A, B, tol=None):
  
    if len(A) != len(B) or len(A[0]) != len(B[0]):
        return False

  
    for i in range(len(A)):
        for j in range(len(A[0])):
            if tol == None:
                if A[i][j] != B[i][j]:
                    return False
            else:
                if round(A[i][j],tol) != round(B[i][j],tol):
                    return False

    return True

def dot_product(A, B):
   
    rowsA = len(A); colsA = len(A[0])
    rowsB = len(B); colsB = len(B[0])
    if rowsA != rowsB or colsA != colsB:
        raise ArithmeticError('Matrices are NOT the same size.')

   
    total = 0
    for i in range(rowsA):
        for j in range(colsB):
            total += A[i][j] * B[i][j]

    return total

def unitize_vector(vector):
   
    if len(vector) > 1 and len(vector[0]) > 1:
        raise ArithmeticError(
            'Vector must be a row or column vector.')

   
    rows = len(vector); cols = len(vector[0])
    mag = 0
    for row in vector:
        for value in row:
            mag += value ** 2
    mag = mag ** 0.5

  
    new = copy_matrix(vector)

  
    for i in range(rows):
        for j in range(cols):
            new[i][j] = new[i][j] / mag

    return new

def scale_matrix(scaler, M):
    
    new = copy_matrix(M)
    rows = len(new); cols = len(new[0])

    
    for i in range(rows):
        for j in range(cols):
            new[i][j] = new[i][j] * scaler

    return new 

def scale_matrix_by_max(A):
    
    
    max = 0
    for row in A:
        for col in row:
            if abs(col) > max:
                max = col

    
    new = copy_matrix(A)
    rows = len(new)
    cols = len(new[0])

    
    for i in range(rows):
        for j in range(cols):
            new[i][j] = new[i][j] / max

    return new 


def replace_nth_column_of_matrix(column_vector, M, column_num):
  
    rows = len(M); cols = len(M[0])
    
  
    if not isinstance(column_vector,list):
        column_value = column_vector
        column_vector = [] 
        for i in range(rows): 
            column_vector.append([column_value]) 


    if rows != len(column_vector):
        raise ArithmeticError('Column and Matrix rows do NOT match.')


    for i in range(rows):
        M[i][column_num] = column_vector[i][0]

    return M

def check_squareness(A):
    if len(A) != len(A[0]):
        raise ArithmeticError("Matrix must be square to inverse.")

def determinant_recursive(A, total=0):
    indices = list(range(len(A)))
    

    if len(A) == 2 and len(A[0]) == 2:
        val = A[0][0] * A[1][1] - A[1][0] * A[0][1]
        return val


    for fc in indices: 
        As = copy_matrix(A) 
        As = As[1:] 
        height = len(As)

        for i in range(height): 
            As[i] = As[i][0:fc] + As[i][fc+1:] 

        sign = (-1) ** (fc % 2) 
        sub_det = determinant_recursive(As) 
        total += sign * A[0][fc] * sub_det 

    return total

def determinant_fast(A):
    
    n = len(A)
    AM = copy_matrix(A)

   
    for fd in range(n): 
        if AM[fd][fd] == 0: 
            AM[fd][fd] = 1.0e-18 
        for i in range(fd+1,n): 
            crScaler = AM[i][fd] / AM[fd][fd] 
            for j in range(n): 
                AM[i][j] = AM[i][j] - crScaler * AM[fd][j]
    
 
    product = 1.0
    for i in range(n):
        product *= AM[i][i] 

    return product

def check_non_singular(A):
  
    det = determinant_fast(A)
    if det != 0:
        return det
    else:
        raise ArithmeticError("Singular Matrix!")

def invert_matrix(A, tol=None):
   
    check_squareness(A)
    check_non_singular(A)

  
    n = len(A)
    AM = copy_matrix(A)
    I = identity_matrix(n)
    IM = copy_matrix(I)

 
    indices = list(range(n)) 
    for fd in range(n): 
        fdScaler = 1.0 / AM[fd][fd]
    
        for j in range(n): 
            AM[fd][j] *= fdScaler
            IM[fd][j] *= fdScaler
     
        for i in indices[0:fd] + indices[fd+1:]: 
            crScaler = AM[i][fd] 
            for j in range(n): 
                AM[i][j] = AM[i][j] - crScaler * AM[fd][j]
                IM[i][j] = IM[i][j] - crScaler * IM[fd][j]
    
    return IM
    

def solve_equations(A, B, tol=None):
  
    check_squareness(A)
    check_non_singular(A)


    n = len(A)
    AM = copy_matrix(A)
    I = identity_matrix(n)
    BM = copy_matrix(B)

  
    indices = list(range(n)) 
    for fd in range(n): 
        if AM[fd][fd] == 0:
            AM[fd][fd] = 1.0e-18
        fdScaler = 1.0 / AM[fd][fd]
     
        for j in range(n): 
            AM[fd][j] *= fdScaler
        BM[fd][0] *= fdScaler
     
        for i in indices[0:fd] + indices[fd+1:]: 
            crScaler = AM[i][fd] 
            for j in range(n): 
                AM[i][j] = AM[i][j] - crScaler * AM[fd][j]
            BM[i][0] = BM[i][0] - crScaler * BM[fd][0]

    return BM
   

def least_squares(X, Y, tol=2):
   
 #    if not isinstance(X[0],list):
#         X = [X]
#     if not isinstance(type(Y[0]),list):
#         Y = [Y]
# 
#   
#     if len(X) < len(X[0]):
#         X = transpose(X)
#     if len(Y) < len(Y[0]):
#         Y = transpose(Y)

    AT = transpose(X)
    ATA = matrix_multiply(AT, X)
    ATB = matrix_multiply(AT, Y)
    coefs = solve_equations(ATA,ATB,tol=tol)
    
    return coefs