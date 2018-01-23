function Z = zernike(i, r, theta)
load('zernike_index');
n = zernike_index(i,1);
m = zernike_index(i,2);

if m == 0
    Z = sqrt(n+1) * zrf(n,0,r);
else
    if mod(i,2) == 0
        Z = sqrt(2*(n+1)) * zrf(n,m,r) .* cos(m*theta);
    else
        Z = sqrt(2*(n+1)) * zrf(n,m,r) .* sin(m*theta);
    end
end
return

%zernike radial function
function R = zrf(n, m, r)
R = 0;
for s = 0:((n-m)/2)
    num = (-1)^s * gamma(n-s+1);
    denom = gamma(s+1) * gamma((n+m)/2 - s+1)...
        * gamma((n-m)/2 - s+1);
    R = R + num/denom * r.^(n- 2*s);
end