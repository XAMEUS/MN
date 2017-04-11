// Question 3
function [l, m]=factorisation_cholesky(D, SD)
    // arg: la diagonale et la sous-diagonale d’une matrice M symétrique définie 
    //      positive et tridiagonale.
    // return: deux vecteurs l (diagonale) et m (sous-diagonale) issus de
    //      la matrice résultante de la factorisation de cholesky.
    l=zeros(1, length(D))
    m=zeros(1, length(SD))
    l(1) = sqrt(D(1))
    for i = 1:length(D)-1
        m(i) = SD(i) / l(i);
        l(i+1) = sqrt(D(i) - m(i)^2);
    end
endfunction

// Unit Tests:
M = [2, -1, 0, 0, 0;
    -1, 2, -1, 0, 0;
    0, -1, 2, -1, 0;
    0, 0, -1, 2, -1;
    0, 0, 0, -1, 2;]
disp(M)
[l, m] = factorisation_cholesky(diag(M), diag(M, -1))
disp(l, m)
L = diag(l) + diag(m, -1)
disp(L)
A = L * L'
disp(A)
assert_checkalmostequal(A, M, 1.0D-10);

// Question 4
function Z=descente(l, m, Y)
    // arg: l la sous-diagonale et m la diagonale de L, issus de la factorisation de Cholesky d'une matrice n*n. Y vecteur de taille n.
    // return: Vecteur Z de taille n, tel que LZ=Y.
    Z = zeros(length(Y), 1)
	Z(1) = Y(1) / l(1)
    for i = 2:length(Z)
        Z(i) = (Y(i) - m(i-1) * Z(i-1)) / l(i)
    end
endfunction

Y = [1; 2; 3; 4; 5]
Z = descente(l, m, Y)
W = L * Z
assert_checkalmostequal(Y, W, 1.0D-10);

//Question 5
function X=remonte(l, m, Z)
    // arg: l la sous-diagonale et m la diagonale de L, issus de la factorisation de Cholesky d'une matrice n*n. Z vecteur de taille n.
    // return: Vecteur X de taille n, tel que T(L)X=Z.
    X = zeros(length(Z), 1)
	X(length(X)) = Z(length(X)) / l(length(X))
    for i = length(X)-1:-1:1
        X(i) = (Z(i) - m(i) * X(i+1)) / l(i)
    end
endfunction

X = remonte(l, m, Z)
disp(X)
assert_checkalmostequal(Z, L' * X, 1.0D-10);

// Question 7
n = 10  
l = 10
A = zeros(n, n)
A(1, 2) = -exp(-(l - 3*l/(n + 1))/l)
A(1, 1) = exp(-(l - l/(n + 1))/l) - A(1, 2)
for i = 2:n-1
    A(i, i-1) = A(i-1, i)
    A(i, i+1) = -exp((2 * l / (n+1) * (i+1/2)-l)/l)
    A(i, i) = -A(i, i-1) - A(i, i+1)
end
A(n, n-1) = A(n-1, n)
A(n, n) = -A(n, n-1) + exp(-(l - l*(2*n + 1)/(n + 1))/l)

B = zeros(1, n)
B(1) = exp(-3/4)

[l, m] = factorisation_cholesky(diag(A), diag(A, -1))
U = descente(l, m, B)
X = remonte(l, m, U)
plot(X')
