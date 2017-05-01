clear
stacksize(268435454)
exec("q3a5.sci")
exec("ressources_q11aq14.sci")

//Question 12
scf()
aaa = -6
bbb = 3
ppp = 50
xxx = aaa:abs(aaa-bbb)/(ppp-1):bbb
yyy = zeros(1, ppp)
for i=1:ppp
    disp(i)
    yyy(1, i) = J(xxx(1, i));
end
plot(xxx, yyy)

