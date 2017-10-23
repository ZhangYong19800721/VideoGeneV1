function d = hamming_uint8(x1,x2)
    b1 = dec2bin(x1,8);
    b2 = dec2bin(x2,8);
    d = pdist2(b1,b2,'hamming') * 8;
end

