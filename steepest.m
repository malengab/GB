function I1 = steepest(f,omega,n,sid4,d24,xi,aa,xix)
% steepest descent example with stationary point
% usign generalized Laguerre polynomials
% 
% % stationary points
% d24 = d2./d4;
% xi = -c2/2./d4;
% xiy = @(x) xi - d24.*x;
% xix = (-c1.*d4+d2.*c2)./(2*d1.*d4-2*d2.^2);
% xiyx = xiy(xix);    % xi_y in xi_x
% % plot(imag(xix),imag(xiyx),'.');
% % plot3(real(xix),imag(xix),real(xiyx),'.')
% % figure
% % plot3(real(xix),imag(xix),imag(xiyx),'.')
% % oscillator evaluated at the stationary point
% g = c1.*xix+c2.*xiyx+d1.*xix.^2+d4.*xiyx.^2+2*d2.*xix.*xiyx; 
% gb = g + b;

% %% auxiliary variables
% aa = d1-d2.^2./d4;
% bb = c1-c2.*d2./d4;
% cc = bb.^2 + 4*(g + c2.^2/4./d4).*aa;

%% Gauss-Laguerre quadrature's weigths and nodes
% n = 1;     % nr of points
alpha = -1/2;  % power of generalized G-L
[pj,w] = GaussLaguerre(n,alpha);    % points + weights
pj = pj/omega;

%%
I1 = zeros(size(sid4));
for ii = 1:n
    for jj = 1:n
        p = pj(ii); q = pj(jj);
        % paths in x
%         uin = sqrt_horse(cc + 4i*aa*p,-1i./aa);     
%         S = sqrt_horse(p.*aa,-aa);   
%         UX2 = (-bb + uin)/2./aa;      % steepest descent path nr 1
%         UX1 = (-bb - uin)/2./aa;      % steepest descent path nr 2
        S = aa.*sqrt(p);
        UX1 = xix + S;      % steepest descent path nr 1
        UX2 = xix - S;      % steepest descent path nr 2

        % get the directions of the paths right
% %         dir1 = imag(uin) > 0; %
        dir1 = real(UX1)>real(xix); % just for the check
%         UX1 = dir1.*ux1 + ~dir1.*ux2;
%         UX2 = dir1.*ux2 + ~dir1.*ux1;
        if any(~dir1)
            error('shout1')
        end

        % derivatives of paths
        UP1 = aa./sqrt(p)/2; %aa./S/2;  % its derivative
%         UP2 = -UP1;    % path nr 2 der
%         UP1 = dir1.*up1 + ~dir1.*up2;
%         UP2 = dir1.*up2 + ~dir1.*up1;
        
        % auxiliary % xiy = @(x) xi - d24.*x;
        % xi = -C2/2./M4 and d24 = M2./M4
        xiyu1 = xi - d24.*UX1; %xiy(UX1);   % xi_y in u_x1 
        xiyu2 = xi - d24.*UX2; %xiy(UX2);   
%         au1 = c2 + 2*d2.*UX1;
%         au2 = c2 + 2*d2.*UX2;
        
        %paths in y   
%       sd4 = sqrt(1i./d4);  % newly sid4
      hin = sid4.*sqrt(q);
%     hin = sqrt_horse(1i*q./d4,-1i./d4);
%     hin2 = sqrt_horse((au2+ 2*d4.*xiyu2).^2 + 4i*q*d4,-1i./d4);
%     Sy = sqrt_horse(1i*d24.^2*p./aa + 1i*q./d4,-1i./d4);
%     Syn = sqrt_horse(4*d4.^2.*(xiyu1-xiyx).^2 + 4i*q.*d4,-1i./d4);

    VY11 = xiyu1 + hin;
    VY12 = xiyu2 + hin;
    VY21 = xiyu1 - hin;
    VY22 = xiyu2 - hin;
%       
    % path directions in y
% %     dir2 = imag(hin1) > 0;
% %     dir3 = imag(hin2) > 0;   %real(vy11)>real(xiyu1);
    dir2 = real(VY11)>real(xiyu1);
    dir3 = real(VY12)>real(xiyu2);
    if any(~dir2)
        error('shout2')
    elseif any(~dir3)
        error('shout3')
    end
%     VY11 = dir2.*vy11 + ~dir2.*vy21;
%     VY12 = dir3.*vy12 + ~dir3.*vy22;
%     VY21 = dir2.*vy21 + ~dir2.*vy11;
%     VY22 = dir3.*vy22 + ~dir3.*vy12;
    
  % derivative of paths
    VP11 = sid4./sqrt(q)/2; %1i./hin/2./d4;
%     VP12 = VP11;
%     VP21 = -VP11;
%     VP22 = -VP11;
 
%     % reshuffle them 
%         VP11 = dir2.*vp11 + ~dir2.*vp21;
%         VP21 = dir2.*vp21 + ~dir2.*vp11;
%         VP12 = dir3.*vp12 + ~dir3.*vp22;
%         VP22 = dir3.*vp22 + ~dir3.*vp12;
% %         
    % integral contributions
    UV = UP1.*VP11;
%         F121 = f(UX1,VY11);  
% %         F122 = f(UX2,VY12).*UP2.*VP12;  % VP12 = VP11, UP2 = -UP1
%         F122 = -f(UX2,VY12);
% %         F211 = f(UX1,VY21).*UP1.*VP21;  % VP21 = VP22, VP22 = -VP11
%         F211 = -f(UX1,VY21);
% %         F212 = f(UX2,VY22).*UP2.*VP22;  % UP2 = -UP1, VP22 = -VP11
%         F212 = f(UX2,VY22);
        
        F = f([UX1 UX2 UX1 UX2],[VY11 VY12 VY21 VY22]);
        F12 = F(:,1)+F(:,2)+F(:,3) + F(:,4);
    % integral increment    
    incr = omega^(-2-2*alpha).*F12.*UV*(p*q)^(-alpha)*w(ii)*w(jj);
%             (F121-F122-F211+F212).*UV*(p*q)^(-alpha)*w(ii)*w(jj);
        I1 = I1 + incr;
    end
end
% I1 = exp(1i*omega*gb).*I1;