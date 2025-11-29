function [v, P2m, u, tc] = v_t_car_acc(t, P2m, u)
gg = 9.81;
vc = P2m/u/gg;
tc = P2m/(u*gg).^2;
vlo     = t<tc;
v       = zeros(size(t));
v(vlo)  = u*gg*t(vlo);
v(~vlo) = sqrt(vc.^2+2*P2m*(t(~vlo)-tc));
end