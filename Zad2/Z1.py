import copy

def eliminacjaGaussa(A, b):
    A = copy.deepcopy(A)
    b = copy.deepcopy(b)
    n = len(A)

    def znajdz_pivot(i_col):
        max_row = i_col
        max_val = abs(A[i_col][i_col])
        for i in range(i_col + 1, n):
            if abs(A[i][i_col]) > max_val:
                max_val = abs(A[i][i_col])
                max_row = i
        return max_row

    def wykonaj_eliminacje(i_col):
        for i_row in range(i_col + 1, n):
            factor = A[i_row][i_col] / A[i_col][i_col]

            for col in range(i_col, n):
                A[i_row][col] -= factor * A[i_col][col]

            b[i_row] -= factor * b[i_col]

    for i in range(n):
        pivot_row = znajdz_pivot(i)
        A[i], A[pivot_row] = A[pivot_row], A[i]
        b[i], b[pivot_row] = b[pivot_row], b[i]

        if abs(A[i][i]) < 1e-12:
            raise ValueError("Macierz jest osobliwa lub prawie osobliwa")

        wykonaj_eliminacje(i)

    x = [0.0] * n
    for i in range(n-1, -1, -1):
        suma = b[i]
        for j in range(i+1, n):
            suma -= A[i][j] * x[j]
        x[i] = suma / A[i][i]  # <- BEZ int()

    return x

