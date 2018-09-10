function [vend1,c1,tt] = GBsol(epsilon,xx,yy,phi01,u1,x0,Rpsi,tol)
 % save transition points
    xnew1 = u1(1);
    ynew1 = u1(2);
    pnew1 = u1(3);
    pnew2 = u1(4);
    m1 = u1(5);
    m2 = u1(6);
    m3 = u1(7);
    m4 = u1(8);
    A1 = u1(9);
%% RAY 1    
% tol = 1e-10;             % smaller contribution is omitted
if (tol<abs(A1));       % small GB amplitude does not contribute
    % GB radius
    m1i = imag(m1); m2i = imag(m2); m4i = imag(m4);
    sqrtm = sqrt((m1i-m4i).^2+4*m2i.^2);
    lam1 = (m1i+m4i+sqrtm)/2;
    lam2 = (m1i+m4i-sqrtm)/2;
    lam = min(lam1,lam2);
   
    K = log(tol./abs(A1));  % circle where the beam value is tol
    R = sqrt(-2*K*epsilon./lam);    % radius of each GB
%     in = K<=0;   % discard those that are smaller than that
% discard also those that do not make it to the support
% R2 = rad(imag([m1,m2;m3,m4]),A1,epsilon,tol); % wrong
    % centre grid around the ray y-x(t)
       
if (norm([xnew1-x0,ynew1])< Rpsi + R )        % if GB makes it to the support
    % local evaluation
        xxs1 = xx - xnew1;
        yys1 = yy - ynew1;
% circle around the termination point with radius ~sqrt(epsilon)
tic
c1 = (xxs1.^2 + yys1.^2<=R^2);      % find all such grid points
% tt = toc;

% c1 = (sqrt(xxs1.^2+yys1.^2)<Inf);
xt = xxs1(c1);
yt = yys1(c1);
tt = toc;
% length(nonzeros(c1))/numel(xx)

% tic
 % solution at time T
    phi1 = phi01 + xt*pnew1 + yt*pnew2 +...
        1/2*(m1*xt.^2 + (m2 + m3)*xt.*yt + m4*yt.^2);
    vend1 = A1*exp(1i*phi1/epsilon);
% tt = toc;
else
    vend1 = [];
    c1 = [];
    tt = 0;
end
else
    vend1 = [];
    c1 = [];
    tt = 0;
end
