function p = PWF(N, i, d)
%   This function calculates the Present Worth Factor for N inflating costs
%   with inflation rate i and discount rate d.  The total present worth is
%   the annual cost multiplied by the Present Worth Factor.

if i==d
    p = N/(1+i);
else
    p = (1/(d-i)) * (1 - ((1+i)/(1+d))^N);
end

end
