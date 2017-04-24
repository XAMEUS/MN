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
n = 1000
l = 10
function c=C(x, l)
    c = exp(-x/l)
endfunction
function c=Ci(i, n, l)
    c = i * 2 * l / (n + 1) - l
    c = C(c, l)
endfunction
function A=gen_matrice(n, l)
    A = zeros(n, n)
    A(1, 2) = -Ci(1+1/2, n, l)
    A(1, 1) = Ci(1-1/2, n, l) - A(1, 2)
    for i = 2:n-1
        A(i, i-1) = A(i-1, i)
        A(i, i+1) = -Ci(i+1/2, n, l)
        A(i, i) = -A(i, i-1) - A(i, i+1)
    end
    A(n, n-1) = A(n-1, n)
    A(n, n) = -A(n, n-1) + Ci(n+1/2, n, l)
endfunction

A = gen_matrice(n, l)
B = zeros(n, 1)
B(1) = Ci(1/2, n, l)
[d, m] = factorisation_cholesky(diag(A), diag(A, -1))
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
M = eye(n, n) + 0.5 * mu * A
N = eye(n, n) - 0.5 * mu * A
U_actuel = zeros(n, 1)

scf()
//MU = Y
[d, m] = factorisation_cholesky(diag(M), diag(M, -1))
h = 1//20000
nb_lines = 30
for i = 1:h
    Y = N * U_actuel + mu * B
    U = descente(d, m, Y)
    U_actuel = remonte(d, m, U)
    ix = floor(modulo(i, h/nb_lines))
    if ix == 0 then
        plot(x, U_actuel');
        e = gce()
        e.children(1).foreground=color(0, 255-255*i/(h/nb_lines)/nb_lines, 255*i/(h/nb_lines)/nb_lines);
    end
end

// Question 11
// Attention: On redéfinit des constantes et une fonction dans cette partie!
funcprot(0)
a = 0.8
l = 10
T = 60
n = 2000
delta_x = 2*l / n
n_t = 3000
delta_t = T / n_t
t_inter = 2 * T / 3
t_fin = T
mu = delta_t * (n+1)**2 / (2 * l)**2
F_cible = [-0.1, -0.18]
function c=C(x, x_d)
    c = 1 - a * exp(-(x - x_d)**2 / 4)
endfunction
function c=Ci(i, n, l, x_d)
    c = i * 2 * l / (n + 1) - l
    c = C(c, x_d)
endfunction
function A=gen_matrice(n, l, x_d)
    A = zeros(n, n)
    A(1, 2) = -Ci(1+1/2, n, l, x_d)
    A(1, 1) = Ci(1-1/2, n, l, x_d) - A(1, 2)
    for i = 2:n-1
        A(i, i-1) = A(i-1, i)
        A(i, i+1) = -Ci(i+1/2, n, l, x_d)
        A(i, i) = -A(i, i-1) - A(i, i+1)
    end
    A(n, n-1) = A(n-1, n)
    A(n, n) = -A(n, n-1) + Ci(n+1/2, n, l, x_d)
endfunction

function [f_inter, f_fin] = flux(x_d)
    A = gen_matrice(n, l, x_d)
    M = eye(n, n) + 0.5 * mu * A
    N = eye(n, n) - 0.5 * mu * A
    [d, m] = factorisation_cholesky(diag(M), diag(M, -1))
    U_actuel = zeros(n, 1)
    for t = 0:n_t
        Y = N * U_actuel
        // B plus constante.
        Y(1) = Y(1) + mu * Ci(1/2, n, l, x_d) * (delta_t/T)**2 * (t**2 + (t+1)**2) / 2
        U_actuel = remonte(d, m, descente(d, m, Y))
        if t == int(2 * n_t / 3) then
            f_inter = Ci(1/2, n, l, x_d) * (U_actuel(1) - (t * delta_t / T)**2) / delta_x - (delta_x / T)**2 * t
        end
    end
    f_fin = Ci(1/2, n, l, x_d) * (U_actuel(1) - (t * delta_t / T)**2) / delta_x - (delta_x / T)**2 * t
endfunction

//Question 12
function norme=J(x_d)
    num = flux(x_d) - F_cible
    norme = (num(1)**2 + num(2)**2) / (F_cible(1)**2 + F_cible(2)**2)
endfunction
//TODO tracer la courbe

//Question 13
function res=dichotomie(funct, epsilon, a, b) // A tester!
    while(b - a >= epsilon)
        disp(a, b)
        Ja = funct(a + (b-a) / 4)
        Jb = funct(a + (b-a) / 2)
        Jc = funct(a + (b-a) * 3 / 4)
        if Ja <= Jb then
            b = a + (b-a) / 2
        elseif Jb <= Jc then
            swpa = a + (b-a) / 4
            swpb = a + (b-a) * 3 / 4
            a = swpa
            b = swpb
        else
            a = a + (b-a) / 2
        end
    end
    res = b - a
endfunction

//res_dich = dichotomie(J, 1e-5, -l, l)
