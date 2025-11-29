function ploy_t_v_a_curves(t, v, tq, vq, aq, cc, PtoM, g_max_friction)
% t[s]
% v[m/s]
% a[m/s/s]
figure()
if nargin<6
    cc = lines(1);
end
vc = PtoM/(9.81*g_max_friction);
t_max = max(t);

s(1) = subplot(2, 2, 1);
hold on
plot(t, v*3.6, 'o', 'Color',cc);
plot(tq, vq*3.6, '.-', 'Color',cc);
t1=0:0.1:10; 
plot(t1, v_t_car_acc(t1, PtoM, g_max_friction)*3.6); 

xlim([0 t_max])
ylim([0 150])
ylabel('velocity[km/h]')
grid on

s(2) = subplot(2, 2, 3);
hold on
% plot(dataAcc.sT, dataAcc.a/10-10);
plot(tq, aq/9.81, '.-', 'Color',cc);
ylim([0 1.1])
xlim([0 t_max])
ylabel('acceleration[g]')
xlabel('time[sec]')
grid on
legend show

s(3)  = subplot(2,2, [2 4]);
hold on
ok = vq>15/3.6 & aq > 0 & vq<115/3.6; % ;
plot(vq(ok)*3.6, aq(ok)/9.81, '.-', 'Color',cc);
ok = vq>0 & aq > 0 & vq<115/3.6 & vq>=vc ;% ;
x = [vc vq(ok)]*3.6;
y = PtoM*1./([vc vq(ok)])/9.81;
plot(x, y, ':', 'DisplayName','P', 'Color',cc);
ok = vq>0 & aq > 0 & vq<115/3.6 & vq<=vc ;% ;
vcs = vc * ones(size(vq));
x = [vq(ok) vc]*3.6;
y = PtoM*1./([vcs(ok) vc])/9.81;
plot(x, y, ':', 'DisplayName','P', 'Color',cc);


xlim([0 150])
ylim([0 1])
linkaxes(s(1:2), 'x')
end