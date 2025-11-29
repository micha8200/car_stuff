function v0m = extrap_interp_start_t(t, v)
for i=1:3
    dv              = v(2+i) - v(2);
    dt              = t(2+i) - t(2);
    dt1             = t(2) - t(1);
    a0              = dv/dt;
    if a0>0
        v0m  = v(2) - a0*dt1;
        if v0m<0
            return
        end
    end
end
v0m = inf;
end