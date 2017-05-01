clear
stacksize(268435454)
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
        l(i+1) = sqrt(D(i+1) - m(i)^2);
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
n = 100
l = 10
function c=C(x, l)
    c = exp(-x/l)
endfunction
function c=Ci(i, n, l)
    c = i * 2 * l / (n + 1) - l
    c = C(c, l)
endfunction
function [D, SD]=gen_matriceA(n, l)
    D = zeros(n, 1)
    SD = zeros(n-1, 1)
    for i=1:n-1
        SD(i) = -Ci(i+1/2, n, l)
        D(i) = Ci(i-1/2, n, l) + Ci(i+1/2, n, l)
    end
    D(n) = Ci(n-1/2, n, l) + Ci(n+1/2, n, l)
endfunction

[A_D, A_SD] = gen_matriceA(n, l)
B = zeros(n, 1)
B(1) = Ci(1/2, n, l)
[d, m] = factorisation_cholesky(A_D, A_SD)
U = descente(d, m, B)
X = remonte(d, m, U)

//plot(X)
x = [-l:2*l/(n-1):l]
plot(x, X')
plot(x, (exp(x/l)-exp(1))/(exp(-1)-exp(1)), 'r-')

// Question 10
T = 10
n_t = 1000
delta_t = T / n_t
mu = delta_t * (n+1)**2 / (2 * l)**2
M_D = ones(n, 1) + 0.5 * mu * A_D
M_SD = 0.5 * mu * A_SD
N_D = ones(n, 1) - 0.5 * mu * A_D
N_SD = - 0.5 * mu * A_SD
U_actuel = zeros(n, 1)

scf()
//MU = Y
[d, m] = factorisation_cholesky(M_D, M_SD)
h = 20000
nb_lines = 30
for i = 1:h
    Y = zeros(n, 1)
    Y(1) = mu * B(1) + N_D(1) * U_actuel(1) + N_SD(1) * U_actuel(2)
    for i=2:n-1
        Y(i) = mu * B(i) + N_SD(i-1) * U_actuel(i-1) + N_D(i) * U_actuel(i) + N_SD(i) * U_actuel(i+1)
    end
    Y(n) = mu * B(n) + N_D(n) * U_actuel(n) + N_SD(n-1) * U_actuel(n-1)
    U = descente(d, m, Y)
    U_actuel = remonte(d, m, U)
    ix = floor(modulo(i, h/nb_lines))
    if ix == 0 then
        plot(x, U_actuel');
        e = gce()
        e.children(1).foreground=color(0, 255-255*i/(h/nb_lines)/nb_lines, 255*i/(h/nb_lines)/nb_lines);
    end
end


