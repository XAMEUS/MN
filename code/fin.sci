clear
stacksize(268435454)
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
function Z=descente(l, m, Y)
    // arg: l la sous-diagonale et m la diagonale de L, issus de la factorisation de Cholesky d'une matrice n*n. Y vecteur de taille n.
    // return: Vecteur Z de taille n, tel que LZ=Y.
    Z = zeros(length(Y), 1)
	Z(1) = Y(1) / l(1)
    for i = 2:length(Z)
        Z(i) = (Y(i) - m(i-1) * Z(i-1)) / l(i)
    end
endfunction
function X=remonte(l, m, Z)
    // arg: l la sous-diagonale et m la diagonale de L, issus de la factorisation de Cholesky d'une matrice n*n. Z vecteur de taille n.
    // return: Vecteur X de taille n, tel que T(L)X=Z. 
    X = zeros(length(Z), 1)
	X(length(X)) = Z(length(X)) / l(length(X))
    for i = length(X)-1:-1:1
        X(i) = (Z(i) - m(i) * X(i+1)) / l(i)
    end
endfunction

// Question 11
a = 0.8
l = 10
T = 60
n = 2000
delta_x = 2*l / (n + 1)
n_t = 3000
delta_t = T / n_t
t_inter = 2 * T / 3
t_fin = T
mu = delta_t * (n+1)**2 / (2 * l)**2
F_cible = [-0.1, -0.18]
function c=C(x, x_d)
    c = 1 - a * exp(-(x - x_d)**2 / 4)
endfunction
function c=Ci(i, x_d)
    c = C(i * delta_x - l, x_d)
endfunction
function [D, SD]=gen_matriceA(x_d)
    D = zeros(n, 1)
    SD = zeros(n-1, 1)
    for i=1:n-1
        SD(i) = -Ci(i+1/2, x_d)
        D(i) = Ci(i-1/2, x_d) + Ci(i+1/2, x_d)
    end
    D(n) = Ci(n-1/2, x_d) + Ci(n+1/2, x_d)
endfunction

function [f_inter, f_fin] = flux(x_d)
    [A_D, A_SD] = gen_matriceA(x_d)
    M_D = ones(n, 1) + 0.5 * mu * A_D
    M_SD = 0.5 * mu * A_SD
    N_D = ones(n, 1) - 0.5 * mu * A_D
    N_SD = - 0.5 * mu * A_SD
    [d, m] = factorisation_cholesky(M_D, M_SD)
    U_actuel = zeros(n, 1)
    for t = 0:n_t
        Y = zeros(n, 1)
        Y(1) = N_D(1) * U_actuel(1) + N_SD(1) * U_actuel(2)
        for i=2:n-1
            Y(i) = N_SD(i-1) * U_actuel(i-1) + N_D(i) * U_actuel(i) + N_SD(i) * U_actuel(i+1)
        end
        Y(n) = N_D(n) * U_actuel(n) + N_SD(n-1) * U_actuel(n-1)
        // B plus constante.
        Y(1) = Y(1) + mu * Ci(1/2, x_d) * (t**2 + (t+1)**2) / 2 / (n_t**2)
        U_actuel = remonte(d, m, descente(d, m, Y))
        if t == int(2 * n_t / 3) then
            f_inter = (Ci(1/2, x_d) * (U_actuel(1) - (t / n_t)**2) / delta_x) - delta_x * t / T / n_t
        end
    end
    f_fin = (Ci(1/2, x_d) * (U_actuel(1) - (t / n_t)**2) / delta_x) - delta_x * t / T / n_t
endfunction

//Question 12
function norme=J(x_d)
    num = flux(x_d) - F_cible
    norme = (num(1)**2 + num(2)**2) / (F_cible(1)**2 + F_cible(2)**2)
endfunction

scf()
aaa = -6
bbb = 3
ppp = 5
xxx = aaa:abs(aaa-bbb)/(ppp-1):bbb
yyy = eye(1, ppp)
for i=1:ppp
    disp(xxx, yyy)
    yyy(1, i) = J(xxx(1, i));
end
//plot(xxx, yyy)

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
