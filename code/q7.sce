clear
stacksize(268435454)
exec("q3a5.sci")
exec("ressources_q7aq10.sci")

// Question 7
[d, m] = factorisation_cholesky(A_D, A_SD)
U = descente(d, m, B)
X = remonte(d, m, U)

//plot(X)
plot(x, X')
plot(x, (exp(x/l)-exp(1))/(exp(-1)-exp(1)), 'r-')

