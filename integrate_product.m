function I = integrate_product(lo, hi, k0, ba, bb, ws, wsa, wsb,...
                               ns, nsa, nsb, coeffs_a, coeffs_b,...
                               t0a, t0b, delta)
    if ~exist('delta', 'var')
        delta = false;
    end
    
    I = 0;
    if lo >= hi
        return
    end
    
    permt0 = 8.8541878128E-12;      % Free space permittivity
    ts = zeros(1, length(ws));
    for i = 1:length(ts)
        ts(i) = sum(ws(1:i));
    end
    tsa = zeros(1, length(wsa));
    for i = 1:length(tsa)
        tsa(i) = sum(wsa(1:i)) + t0a;
    end
    tsb = zeros(1, length(wsb));
    for i = 1:length(tsb)
        tsb(i) = sum(wsb(1:i)) + t0b;
    end
    
    i = 1;
    while lo >= ts(i)
        i = i + 1;
    end
    bounds = [lo, ts(i:end)];
    i = length(bounds);
    while hi <= bounds(i)
        i = i - 1;
    end
    bounds = [bounds(1:i), hi];
    for i = 1:length(bounds) - 1
        tlo = bounds(i); thi = bounds(i + 1);
        x = (tlo + thi) / 2;
        ii = get_ind(x, ws);
        ia = get_ind(x, wsa, t0a);
        ib = get_ind(x, wsb, t0b);
%         ta = tsa(ia)
%         tb = tsb(ib)
        term = int_term(ia, ib, k0, ba, bb, coeffs_a, coeffs_b, ...
                        nsa, nsb, tsa(ia), tsb(ib), tlo, thi);
        if ((delta) && (term ~= 0))
            term = term * permt0 * (ns(ii) ^ 2 - nsb(ib) ^ 2);
        end
        I = I + term;
    end
end