clear
stacksize(268435454)
exec("q3a5.sci")
exec("ressources_q7aq10.sci")

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
    for j=2:n-1
        Y(j) = mu * B(j) + N_SD(j-1) * U_actuel(j-1) + N_D(j) * U_actuel(j) + N_SD(j) * U_actuel(j+1)
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
