function res = CachedFactorial(n)
% This function does speed things up... trust me
persistent cache;
if(isempty(cache))
    cache = zeros(100,1);
end
res = cache(n+1);
if(res)
    return;
else
    res = factorial(n);
    cache(n+1) = res;
end
end
