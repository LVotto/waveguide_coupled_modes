function n = nfromx(x, ws, ns, x0)
    ind = get_ind(x, ws, x0);
    n = ns(ind);
end