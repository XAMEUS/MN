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
x = [-l:2*l/(n-1):l]

