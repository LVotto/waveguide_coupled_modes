function ind = get_ind(x, ws, x0)
    if ~exist('x0', 'var') 
        x0 = 0; 
    end
    ind = 1;
    for i = 1:length(ws) - 1
        t = sum(ws(1:i)) + x0;
        ind = ind + heaviside(x - t);
    end
end