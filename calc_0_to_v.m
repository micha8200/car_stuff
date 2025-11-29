function te = calc_0_to_v(v, p2m, u)
gg = 9.81;
vc = p2m/u/gg;
tc = p2m/(u*gg).^2;

t1 = tc;
t2 = 0.5*(v.^2-vc.^2)/p2m;

te = t1 + t2;
end