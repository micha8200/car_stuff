function a = acc_func(v, P, m, g_max)

a = zeros(size(v));
vc = P/m/(9.81*g_max);
ok = v>0 & v<=vc;
a(ok) = P/m./vc;
ok = v>0 & v>=vc;
a(ok) = P/m./v(ok);

end