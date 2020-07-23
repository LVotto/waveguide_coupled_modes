function s = pterm(j, k0, beta, a, b, ns, ws, pol)
    gamma = sqrt(beta^2 - k0^2 * ns(j)^2);
    c0 = physconst('LightSpeed');
    omega = k0 * c0;
    permt0 = 8.8541878128E-12;      % Free space permittivity
%     permb0 = 1 / c0 ^ 2 / permt0;
    permb0 = 1.25663706212E-6;
    
    s = 0;
    if j == 1
        s = s + abs(a) ^ 2 / 2 / real(gamma);
    elseif j == length(ws)
        s = s + abs(b) ^ 2 / 2 / real(gamma);   
    else
        if real(gamma) ~= 0
            s = s + abs(a) ^ 2 * (1 - exp(-2 * real(gamma) * ws(j))) ...
                  / 2 / real(gamma);
            s = s - abs(b) ^ 2 * (1 - exp( 2 * real(gamma) * ws(j))) ...
                  / 2 / real(gamma);
        else
            s = s + (abs(a) ^ 2 + abs(b) ^ 2) * ws(j);
        end
        if imag(gamma) ~= 0
            s = s + a * conj(b) * (1 - exp(-2 * 1j * imag(gamma) * ws(j))) ...
                  / 2 / 1j / imag(gamma);
            s = s - b * conj(a) * (1 - exp( 2 * 1j * imag(gamma) * ws(j))) ...
                  / 2 / 1j / imag(gamma);
        else
            s = s + (a * conj(b) + b * conj(a)) * ws(j);
        end
    end
    
    if all(pol == 'te')
        mul = beta / 2 / omega / permb0;
    else
        mul = beta / 2 / omega / permt0 / ns(j) ^ 2;
    end
    s = mul * s;
end