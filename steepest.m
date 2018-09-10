function I1 = steepest(f,omega,n,sid4,d24,xi,aa,xix)
% steepest descent example with stationary point
% usign generalized Laguerre polynomials
% 

%% Gauss-Laguerre quadrature's weigths and nodes
% n = 1;     % nr of points
alpha = -1/2;  % power of generalized G-L
[pj,w] = GaussLaguerre(n,alpha);    % points + weights
pj = pj/omega;  % rescale

%%
I1 = zeros(size(sid4));
% loop through all p's and q's
for ii = 1:n
    for jj = 1:n
        p = pj(ii); q = pj(jj);
        S = aa.*sqrt(p);
        UX1 = xix + S;      % steepest descent path nr 1
        UX2 = xix - S;      % steepest descent path nr 2

        % get the directions of the paths right
        dir1 = real(UX1)>real(xix); % just for the check
        if any(~dir1)
            error('shout1')
        end

        % derivatives of paths
        UP1 = aa./sqrt(p)/2; %aa./S/2;  % its derivative
        
        % auxiliary
        xiyu1 = xi - d24.*UX1; 
        xiyu2 = xi - d24.*UX2;
        
        %paths in y   
      hin = sid4.*sqrt(q);
    VY11 = xiyu1 + hin;
    VY12 = xiyu2 + hin;
    VY21 = xiyu1 - hin;
    VY22 = xiyu2 - hin;
%       
    % path directions in y: JUST TO CHECK
    dir2 = real(VY11)>real(xiyu1);
    dir3 = real(VY12)>real(xiyu2);
    if any(~dir2)
        error('shout2')
    elseif any(~dir3)
        error('shout3')
    end
    
  % derivative of paths
    VP11 = sid4./sqrt(q)/2;
 
    % integral contributions
    UV = UP1.*VP11;
        
        F = f([UX1 UX2 UX1 UX2],[VY11 VY12 VY21 VY22]); % push all at the same time
        if ~isempty(F)
            F12 = F(:,1)+F(:,2)+F(:,3) + F(:,4);
          % integral increment    
            incr = omega^(-2-2*alpha).*F12.*UV*(p*q)^(-alpha)*w(ii)*w(jj);
            I1 = I1 + incr;
        end
    end
end