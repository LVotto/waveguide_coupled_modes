function coeffs = tmt_coeffs(beta, k0, ws, ns, pol, normalize)
    syms b;
    if ~exist('normalize', 'var')
        normalize = false;
    end
    if ~exist('pol', 'var')
        pol = 'tm';
    end
    coeffs = zeros(2, length(ns), 'sym');
    coeffs(2, end) = b;
    for i = length(ns):-1:2
        coeffs(:, i - 1) = tmt_matrix(i, beta, k0, ws, ns, pol) ... 
                         * coeffs(:, i);
    end
    
    
    % For normalization:
    if normalize
        S2 = zeros(1, length(ns), 'sym');
        for i = 1:length(ns)
            S2(i) = pterm(i, k0, beta, coeffs(1, i), coeffs(2, i),  ...
                          ns, ws, pol);
        end
        B = 1 / sqrt(sum(subs(S2, b, 1)));
    else
        B = 1;
    end
    
    coeffs = subs(coeffs, b, B);
    coeffs(2, 1) = 0; coeffs(1, end) = 0;
    coeffs = double(coeffs);
end

