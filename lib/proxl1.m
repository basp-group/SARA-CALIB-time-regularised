function p = proxl1(x, gamma)

% proximal operator of gamma*l1(x)
[Ny, Nx] = size(x);
if min(Ny, Nx)>1
    x = x(:);
end
p = sign(x) .* max(abs(x) - gamma, 0);
if min(Ny, Nx)>1
    p = reshape(p, Ny, Nx);
end
end
% % % p = y + sign(x-y).*max(abs(x-y)-gamma,0);