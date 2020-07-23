function betas = tmt_betas(k0, ns, ws, pol)
    syms b;
    tm = eye(2, 2);
    for j = 2:length(ns)
        tm = tm * tmt_matrix(j, b, k0, ws, ns, pol);
    end
    tosolve = tm(2, 2);
    bmin = k0 * min(ns); bmax = k0 * max(ns);
    bs = linspace(bmin, bmax, 10);
    betas = zeros(1, length(bs));
    for i = 1:length(bs) - 1
        bb = vpasolve(tosolve, [bs(i) bs(i + 1)]);
        if ~isempty(bb)
            betas(i) = bb;
        end
    end
    betas = unique(betas);
    betas = betas(~isnan(betas));
    betas = sort(betas);
    betas = betas(betas~=0);
end