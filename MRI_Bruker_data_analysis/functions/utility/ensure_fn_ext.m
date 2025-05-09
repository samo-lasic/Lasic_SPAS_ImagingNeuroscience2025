function fn = ensure_fn_ext(fn, ext)
% ensure filename extension

if iscell(fn)
    ind = ~contains(fn, ext);
    fn(ind) = strcat(fn(ind),ext);
elseif ~contains(fn, ext) % one string
    fn = strcat(fn,ext);
end

end