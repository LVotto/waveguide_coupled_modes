function s = int_term(i, j, k0, ba, bb, coeffs_a, coeffs_b, ...
                      nsa, nsb, ta, tb, lo, hi)
    ga = sqrt(ba^2 - k0^2 * nsa(i)^2);
    gb = sqrt(bb^2 - k0^2 * nsb(j)^2);
    a = coeffs_a(1, i); b = coeffs_a(2, i);
    c = coeffs_b(1, j); d = coeffs_b(2, j);
    if lo >= hi
        s = 0;
        return
    end

    syms x;
    s = 0;
    if ga + gb ~= 0
        s = s + a * c / (ga + gb) * exp((ga + gb) * x - ga * ta - gb * tb);
        s = s - b * d / (ga + gb) * exp(-(ga+ gb) * x + ga * ta + gb * tb);
    else
        s = s + a * c * exp(-ga * ta - gb * tb) * x;
        s = s + b * d * exp( ga * ta + gb * tb) * x;
    end
    
    if ga - gb ~= 0
        s = s + a * d / (ga - gb) * exp((ga - gb) * x - ga * ta + gb * tb);
        s = s - b * c / (ga - gb) * exp(-(ga- gb) * x + ga * ta - gb * tb);
    else
        s = s + a * d * exp(-ga * ta + gb * tb) * x;
        s = s + b * c * exp( ga * ta - gb * tb) * x;
    end
    
    s_hi = 0;
    if abs(hi) ~= Inf
        s_hi = subs(s, hi);
    end
    s_lo = 0;
    if abs(lo) ~= Inf
        s_lo = subs(s, lo);
    end
%      s_hi = subs(s, hi); s_lo = subs(s, lo);
     s = double(s_hi - s_lo);
end