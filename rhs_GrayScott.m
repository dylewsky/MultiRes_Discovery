function rhs = rhs_GrayScott(uv, del_2, n, ru, rv, f, k)

u = uv(1:n);
v = uv(n+1:end);

ut = ru * del_2 * u - u .* v .* v + f * (ones(size(u)) - u);
vt = rv * del_2 * v + u .* v .* v - (f + k) * v;

rhs = [ut; vt];