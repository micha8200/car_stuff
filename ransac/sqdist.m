function d = sqdist(a)
N = size(a, 1);
d =zeros(N*(N-1)/2, 1);
k = 1;
for i=1:N
    for j=i+1:N
        d(k) = sum((a(i,:)-a(j,:)).^2);
        k = k+1;
    end
end
end