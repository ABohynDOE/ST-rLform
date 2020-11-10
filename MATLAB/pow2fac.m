function powers = pow2fac(x,rpad)
    i = 1;
    index = 1;
    powers = zeros(1,ceil(log2(x)));
    while i <= x
        if bitand(x,i)
            powers(1,index) = 1;
        end
        i = bitshift(i,1);
        index = index +1;
    end
    if nargin ==1
        rpad = length(powers);
    end
    if length(powers)<rpad 
        powers = [powers zeros(1,rpad-length(powers))];
    end
end