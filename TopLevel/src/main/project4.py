def jacobi(n, A, b, XO, TOL, N):
    x = XO.copy()
    x_new = [0.0 for _ in range(n)]
    for _ in range(N):
        for i in range(n):
            sigma = sum(A[i][j] * x[j] for j in range(n) if j != i)
            x_new[i] = (b[i] - sigma) / A[i][i]
        if all(abs(x_new[i] - x[i]) < TOL for i in range(n)):
            return x_new
        x = x_new.copy()
    print("Max # iterations exceeded.")
    return None

#Example numbers
n = 3
A = [[5, 1, 2], [1, 3, 1], [1, 2, 3]]
b = [2, 1, 5]
XO = [0, 0, 0]
TOL = 1e-5
N = 100

solution = jacobi(n, A, b, XO, TOL, N)
if solution:
    print("Approximate solution:", solution)
    
    
    
def gauss(n, A, b, XO, TOL, N):
    x = XO.copy()
    for iteration in range(N):
        x_new = x.copy()
        for i in range(n):
            sum1 = sum(A[i][j] * x_new[j] for j in range(i))
            sum2 = sum(A[i][j] * x[j] for j in range(i + 1, n))
            x_new[i] = (b[i] - sum1 - sum2) / A[i][i]
        #check convergence
        if all(abs(x_new[i] - x[i]) < TOL for i in range(n)):
            return x_new
        x = x_new
    print("Maximum # iterations exceeded.")
    return None

#Example numbers
n = 3
A = [[5, 1, 2], [1, 3, 1], [1, 2, 3]]
b = [2, 1, 5]
XO = [1, 1, 1]
TOL = 1e-5
N = 100

solution = gauss(n, A, b, XO, TOL, N)
if solution:
    print("Approximate solution:", solution)
else:
    print("No solution found")
    
    
def sor(n, A, b, XO, omega, TOL, N):
    x = XO[:]
    k = 1

    while k <= N:
        x_old = x[:]
        for i in range(n):
            sum1 = sum(A[i][j] * x[j] for j in range(i))
            sum2 = sum(A[i][j] * x_old[j] for j in range(i + 1, n))
            x[i] = (1 - omega) * x_old[i] + (omega / A[i][i]) * (b[i] - sum1 - sum2)

        #check convergence
        if all(abs(x[i] - x_old[i]) < TOL for i in range(n)):
            return x  

        k += 1  

        

    #If 100 iterations pass without a solution, this message is printed
    print("Maximum number of iterations exceeded")
    return None

#Example numbers
n = 3
A = [
    [5, 1, 2],
    [1, 3, 1],
    [1, 2, 3]
]
b = [2, 1, 5]
XO = [1, 1, 1] 
omega = 1.25  
TOL = 1e-5
N = 100

solution = sor(n, A, b, XO, omega, TOL, N)
if solution:
    print("Approximate solution:", solution)
else:
    print("Solution not found")
    
    
def gauss_elimination(A, b):
    n = len(A)
    M = [row[:] for row in A]  
    x = b[:]  

    for i in range(n):
        
        max_row = max(range(i, n), key=lambda r: abs(M[r][i]))
        if i != max_row:
            M[i], M[max_row] = M[max_row], M[i]
            x[i], x[max_row] = x[max_row], x[i]

        
        for j in range(i + 1, n):
            factor = M[j][i] / M[i][i]
            x[j] -= factor * x[i]
            for k in range(i, n):
                M[j][k] -= factor * M[i][k]

    #substituting back
    for i in range(n - 1, -1, -1):
        for j in range(i + 1, n):
            x[i] -= M[i][j] * x[j]
        x[i] = x[i] / M[i][i]

    return x

def residual(A, x, b):
    n = len(A)
    r = [b[i] - sum(A[i][j] * x[j] for j in range(n)) for i in range(n)]
    return r

def norm(v):
    return max(abs(vi) for vi in v)

def iterative_refinement(A, b, N, TOL, t):
    n = len(A)
    x = gauss_elimination(A, b)
    cnd = None

    for k in range(1, N + 1):
        r = residual(A, x, b)
        y = gauss_elimination(A, r)
        xx = [x[i] + y[i] for i in range(n)]

        if k == 1:
            cnd = norm(y) / norm(xx) * (10 ** t)

        if norm([xx[i] - x[i] for i in range(n)]) < TOL:
            return xx, cnd, "The procedure was successful."

        x = xx[:]  #update x for next gen

    return x, cnd, "Maximum number of iterations exceeded"

#Example numbers
A = [
    [5, 1, 2],
    [1, 3, 1],
    [1, 2, 3]
]
b = [2, 1, 5]
N = 100
TOL = 1e-5
t = 5

solution, cnd, message = iterative_refinement(A, b, N, TOL, t)
print("Approximate Solution:", solution)
print("Condition number:", cnd)
print("Status:", message)
