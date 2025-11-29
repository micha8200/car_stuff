function d = sqdist(a, iseuclid)
N = size(a, 1);
d = zeros(N*(N-1)/2, 1);
k = 1;
for i=1:N
    for j=i+1:N
        if iseuclid
            d(k) = sum((a(i,:)-a(j,:)).^2);
        else
            d(k) = sum(abs(a(i,:)-a(j,:)));
        end
        k = k+1;
    end
end
end