function d = hamming(x1,x2)
    N = size(x1,1);
    d = zeros(N,1);
    for n=1:N
        d(n) = hamming_uint8(x1(n),x2(n));
    end
    d = sum(d);
end

