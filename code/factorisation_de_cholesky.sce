// Question 3
function [l, m]=factorisation_cholesky(D, SD)
    // arg: la diagonale et la sous-diagonale d’une matrice M symétrique définie 
    //      positive et tridiagonale.
    // return: deux vecteurs l (sous-diagonale) et m (diagonale) issus de
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

