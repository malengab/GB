function up = xp6vec(~,u,Y1,Y2,Y3)
% c(x,y) = Y1*exp(-Y2*x^2 -Y3*y^2)
n = length(u)/8;

x = u(1:n);
y = u(n+1:2*n);
px = u(2*n+1:3*n);
py = u(3*n+1:4*n);
m00 = u(4*n+1:5*n);
m10 = u(5*n+1:6*n);
m11 = u(6*n+1:7*n);
a0 = u(7*n+1:end);

pn = sqrt(px.^2+py.^2);
ppx = px./pn;
ppy = py./pn;

g = Y1*exp(-Y2*x.^2 - Y3*y.^2);
c0 = 1 - g;
cx = 2*g.*x*Y2;
cy = 2*g.*y*Y3;

cxx =  -(-2*Y2 + 4*Y2^2*x.^2).*g;
cxy =  -4*Y2*Y3*x.*y.*g;
cyy=  -(-2*Y3 + 4*Y3^2*y.^2).*g;

% mp = M*pp;
mp1 = m00.*ppx + m10.*ppy;
mp2 = m10.*ppx + m11.*ppy;

% Elements of D, M*B and M*C*M

d00 = pn.*cxx;
mb00 = mp1.*cx;
mcm00 = (m00.*m00+m10.*m10-mp1.*mp1).*c0./pn;

d11 = pn.*cyy;
mb11 = mp2.*cy;
mcm11 = (m11.*m11+m10.*m10-mp2.*mp2).*c0./pn;

d10 = pn.*cxy;
mb10 = mp1.*cy;
mcm10 = (m00.*m10+m11.*m10-mp1.*mp2).*c0./pn;

mb01 = mp2.*cx;

mm1 = d00+2*mb00+mcm00;
mm2 = d10+mb10+mb01+mcm10;
mm3 = d11+2*mb11+mcm11;

a0p =  a0./(2*pn).*(-c0.*(ppx.*mp1+ppy.*mp2) + c0.*(m00+m11) - cx.*px-cy.*py); 

up = [-c0.*ppx; ...  % x
    -c0.*ppy; ...    % y
    cx.*pn; ...      % px
    cy.*pn; ...      % py
    mm1; ...         % m00
    mm2; ...         % m10
    mm3; ...         % m11
    a0p];            % a0
end