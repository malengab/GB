function [d24,xi,xix,g] = statpts(C1,C2,M1,M2,M4)
% stationary points
d24 = M2./M4;
xi = -C2/2./M4;
xix = (-C1.*M4+M2.*C2)./(2*M1.*M4-2*M2.^2);
xiyx = xi - d24.*xix; % xi_y in xi_x
% xiyx = xiy(xix);    
g = C1.*xix+C2.*xiyx+M1.*xix.^2+M4.*xiyx.^2+2*M2.*xix.*xiyx;