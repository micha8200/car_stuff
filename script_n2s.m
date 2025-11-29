buf0 = '123456789ABCD';
val = 123456;
[buf, ixe] = n2s(buf0, 8, val, 0, 'aaa', 1);
disp(buf)
[buf, ixe] = n2s(buf0, 8, val, 0);
disp(buf)