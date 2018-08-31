% plots the stationary points, discards some combinations either 
% based dependent on the distance or the radius
close all
kdata = 10*2.^(0:8);% out of memory for 0:7
time1 = []; 
tot =[];
time2 = []; 
Qdata = []; 
sumdata = [];
for k = kdata;
epsilon = 1/k;      % 1/wave number  
% Y1 = .8; Y2 = .75; Y3 = 1;  % former random variables
Y1 = 0; Y2 = 0; Y3 = 0;
% Y1 = 0.1; Y2 = 0.2; Y3 = 0.3;

%% GB preliminaries
% sx = -.5;            % initial bump position
ppw = 10;           % points per wavelength
T = 1;              % final time 
Lx = 4; Ly = 3;        % length of the domain
% n = 1;               %nr of GL points

% initial amplitude
d = 10;
% % ONE PULSE
% a = @(x,y) exp(-d*((x-sx).^2 + y.^2));
% phi = @(x,y) -x;
% px = @(x,y) -1*ones(size(x));   % p_x
% py = @(x,y) 0*ones(size(x));    % p_y
% m00 = @(x,y) (0 + 1i)*ones(size(x)); % p_xx + 1i
% m10 = @(x,y) 0*ones(size(x)); % p_xy
% m11 = @(x,y) (0 + 1i)*ones(size(x)); % p_yy + 1i

% % two pulses
sx1 = -1; sx2 = 1;
a = @(x,y) exp(-d*((x-sx1).^2 + y.^2)) + exp(-d*((x-sx2).^2 + y.^2));
phi = @(x,y) abs(x);
px = @(x,y) sign(x).*ones(size(x));   % p_x
py = @(x,y) 0*ones(size(x));    % p_y
m00 = @(x,y) (0 + 1i)*ones(size(x)); % p_xx + 1i
m10 = @(x,y) 0*ones(size(x)); % p_xy
m11 = @(x,y) (0 + 1i)*ones(size(x)); % p_yy + 1i

% test bump function 
x0 = 0;
psif = @(x,y) exp(-5*((x-x0).^2+y.^2)); %(x.^2 + y.^2<e^2).*exp(-(x.^2+y.^2)./(e^2-x.^2-y.^2));
Rpsi = sqrt(-1/5*log(1e-12));   % 'radius' of the test func

% cutoff
% % choice = 3;
% % iff choice = 1
% c = 10;
% dist = sqrt(epsilon)*c; % distance at which beams interact
% % if choice = 2
tol = 1e-7;  % radius based
% % if choice == 3 bound the imag Psi > bound
% bound = 1;

tol2 = 1e-11; % how precise should NSDM be
nmax = 256; % max number of GL points
nmin = 1;

%% %%% NO NEED TO CHANGE ANYTHING UNDER THIS POINT %%%%%%%%%%%%%

[xx,yy,h,ssx,ssy,cgap] = GBgrid(ppw,Lx,Ly,epsilon);      % GB grid points
% ssx = [0;0.1;0;0.1]; ssy = [0;0.1;0.1;0];
% ssx = 0; ssy = 0;
nnn = length(ssx(:));
         
    %% Compute Gaussian beams
    
    U0=[ssx; ...            % x
        ssy; ...            % y
        px(ssx,ssy); ...    %-ones(nnn,1); ...      % px
        py(ssx,ssy); ...    %zeros(nnn,1); ...      % py
        m00(ssx,ssy); ...   %1i*ones(nnn,1); ...    % m00
        m10(ssx,ssy); ...   %zeros(nnn,1); ...      % m10
        m11(ssx,ssy); ...   1i*ones(nnn,1); ...    % m11
        a(ssx,ssy)];        % A  
    
     options = odeset('RelTol',1e-12);
     [~,U1] = ode45(@xp6vec, [0,T], U0, options, Y1,Y2,Y3);
     
 %% initiate auxiliary variables
         u1 = reshape(U1(end,:),nnn,8);
        phi01 = phi(ssx,ssy);               % phi0
        xnew1 = u1(:,1);
        ynew1 = u1(:,2);
        pnew1 = u1(:,3);
        pnew2 = u1(:,4);
        m1 = u1(:,5);
        m2 = u1(:,6);
        m3 = u1(:,6);
        m4 = u1(:,7);
        A1 = u1(:,8);

 %%     % plot initial data
% v0 = zeros(size(xx));
% for ind = 1:nnn;
% v00 = a(ssx(ind),ssy(ind))*exp(1i*k*(phi(ssx(ind),ssy(ind))+...
%     (xx-ssx(ind))*px(ssx(ind),ssy(ind))+(yy-ssy(ind))*py(ssx(ind),ssy(ind))+...
%     1/2*m00(ssx(ind),ssy(ind))*(xx-ssx(ind)).^2 + ...
%     m10(ssx(ind),ssy(ind)).*(xx-ssx(ind)).*(yy-ssy(ind))+...
%     + 1/2*m11(ssx(ind),ssy(ind)).*(yy-ssy(ind)).^2));
% v0 = v0 + v00;
% end
% v0 = 1/(2*pi*epsilon)*cgap^2*v0;
% figure;
% surf(xx,yy,abs(v0));
% shading interp
% colorbar
% view([0 90])
% title('t = 0','FontSize',16)
% print('-dpdf','sol0.pdf')

%% compare one GB results
% tic
% % slow sol to compare
% vv = zeros(size(xx));
% for ind = 1:length(ssx);
% v = A1(ind)*exp(1i*k*(phi01(ind) + (xx-xnew1(ind))*pnew1(ind) + (yy-ynew1(ind))*pnew2(ind)+...
%     1/2*m1(ind)*(xx-xnew1(ind)).^2 + m2(ind).*(xx-xnew1(ind)).*(yy-ynew1(ind))+...
%     + 1/2*m4(ind).*(yy-ynew1(ind)).^2));
% vv = vv + v;
% end
% v1 = 1/(2*pi*epsilon)*cgap^2*vv;
% figure;
% surf(xx,yy,abs(v1));
% shading interp
% colorbar
% view([0 90])
% title('T = 3','FontSize',16)
% print('-dpdf','sol3.pdf')

% Q  = sum(sum(abs(v1).^2.*psif(xx,yy)))*h^2;
% t2 = toc;
% Qdata = [Qdata Q];
% time2 = [time2 t2];

% tic
%% eigenvalues
m1i = imag(m1); m2i = imag(m2); m4i = imag(m4);
sqrtm = sqrt((m1i-m4i).^2+4*m2i.^2);
lam1 = (m1i+m4i+sqrtm)/2;
lam2 = (m1i+m4i-sqrtm)/2;
lam = min(lam1,lam2);
   
K = log(tol./abs(A1));  % circle where the beam value is tol

R = sqrt(-2*K/k./lam);    % radius of each GB
in = K<=0;   % discard those that are smaller than that
% discard also those that do not make it to the support
in2 = sqrt((x0-xnew1).^2+ynew1.^2)<Rpsi+R;
in = in & in2;

R = R(in);
xnew1 = xnew1(in); ynew1 = ynew1(in); pnew1 = pnew1(in);
pnew2 = pnew2(in); m1 = m1(in); m2 = m2(in); m4 = m4(in);
phi01 = phi01(in); A1 = A1(in);

% auxiliary
b1 = phi01-xnew1.*pnew1-ynew1.*pnew2+...
    1/2*(m1.*xnew1.^2+2*m2.*xnew1.*ynew1+m4.*ynew1.^2); 
c11 = pnew1-(m1.*xnew1+m2.*ynew1);
c12 = pnew2-(m2.*xnew1+m4.*ynew1);

%% batch size
bs = 1e3;
ne = length(xnew1);     % nr of elements
nr = ceil(ne/bs);   % max nr of runs
summ = 0;

% x1 =  []; y1 = [];  % statpts all
% x2 = []; y2 = [];   % statpts only contributing
% us2 = []; usA = []; g1 = []; g2 = [];

us = []; donesum = 0; notdonesum = 0; totalsum = 0;
xdone = []; ydone = []; ndone = [];
t1 = 0;
% to avoid big matrices, we combine the points by batches
for kk = 0:(nr-1);
    for ll = 0:(nr-1);  % all combinations w kk
        if kk == nr-1   % if at last step take the leftovers
             ind = bs*kk+1:ne; 
        else         % take a given batch
            ind = (bs*kk+1):(bs*(kk+1));
        end
        if ll == nr-1   % leftovers
            ind2 = bs*ll+1:ne;
        else
            ind2 =(bs*ll+1):(bs*(ll+1));
        end
    [ba,bb] = meshgrid(b1(ind),b1(ind2));
    B = ba-conj(bb);
    [c1a,c1b] = meshgrid(c11(ind),c11(ind2));
    C1 = c1a-conj(c1b);
    [c2a,c2b] = meshgrid(c12(ind),c12(ind2));
    C2 = c2a-conj(c2b);
    [m11a,m11b] = meshgrid(1/2*m1(ind),1/2*m1(ind2));
    [m12a,m12b] = meshgrid(1/2*m2(ind),1/2*m2(ind2));
    % m13 = m12;
    [m14a,m14b] = meshgrid(1/2*m4(ind),1/2*m4(ind2));
    M1 = m11a-conj(m11b);
    M2 = m12a-conj(m12b);
    % M3 = M2;
    M4 = m14a-conj(m14b);
    [A1a,A1b] = meshgrid(A1(ind),A1(ind2));
    A = A1a.*conj(A1b);
% all stationary points

[d24,xi,xix,g] = statpts(C1,C2,M1,M2,M4);
sid4 = sqrt(1i./M4);

% gb = g + B;

% x1 = [x1; xix(:)]; y1 = [y1; xiyx(:)];  % save
% x1 = xix(:); y1 = xiyx(:);
% g values
% g = C1.*xix+C2.*xiyx+M1.*xix.^2+M4.*xiyx.^2+2*M2.*xix.*xiyx+ B;
%%
% if any(imag(g)<-1e-10);
%     warning('Negative stat points!!!')
% end
%%
% g1 = [g1;g(:)];

%% choose only some
%   % which combinations of points are relevant
    [XX1,XX2] = meshgrid(xnew1(ind),xnew1(ind2));
    [YY1,YY2] = meshgrid(ynew1(ind),ynew1(ind2));
%     if choice == 1
%         CC = sqrt((XX2-XX1).^2+(YY2-YY1).^2)<dist;
% % %     figure
% % %     spy(CC)
% %     disp(sum(sum(CC))/nnn^2);
%     elseif choice == 2 
%   %alternatively based on size of GB
% m1i = imag(M1); m2i = imag(M2); m4i = imag(M4);
% sqrtm = sqrt((m1i-m4i).^2+4*m2i.^2);
% lam1 = (m1i+m4i+sqrtm)/2;
% lam2 = (m1i+m4i-sqrtm)/2;
% lam = min(lam1,lam2);
% K = log(tol./abs(A));  % circle where the beam value is tol
% % in = K<=max(K);   % discard those that are smaller than that
% R = sqrt(-2*K/k./lam);    % radius of each GB
    [RR1,RR2] = meshgrid(R(ind),R(ind2));
    CC = sqrt((XX2-XX1).^2+(YY2-YY1).^2)<=RR1+RR2;

% CC = RR1<=max(max(RR1));
% %     figure
% %     spy(CC2)   
% %     CC = CC2;
%     elseif choice == 3
%         CC = imag(g)<bound;           
%     end
%     disp(sum(sum(CC))/nnn^2);

%%
B = sparse(B(CC));
C1 = sparse(C1(CC));
C2 = sparse(C2(CC));
M1 = sparse(M1(CC));
M2 = sparse(M2(CC));
M4 = sparse(M4(CC));
A = sparse(A(CC));
% 
% % on the top of it, discard those whose amplitude is too small to
% % contribute
% CCC = abs(A)>1e-9;
% B = sparse(B(CCC));
% C1 = sparse(C1(CCC));
% C2 = sparse(C2(CCC));
% M1 = sparse(M1(CCC));
% M2 = sparse(M2(CCC));
% M4 = sparse(M4(CCC));
% A = sparse(A(CCC));

% only the stationary points that make it to the integral
d24 = d24(CC);%M2./M4;
xi = xi(CC);%-C2/2./M4;
sid4 = sid4(CC);
% xiy = @(x) xi - d24.*x;
xix = xix(CC);%(-C1.*M4+M2.*C2)./(2*M1.*M4-2*M2.^2);
% xiyx = xiyx(CC); %xiy(xix);    % xi_y in xi_x
% x2 = [x2; xix(:)]; y2 = [y2; xiyx(:)];  % save
% x2 = xix(:); y2 = xiyx(:);
% g values
g = g(CC);%C1.*xix+C2.*xiyx+M1.*xix.^2+M4.*xiyx.^2+2*M2.*xix.*xiyx+ B;
% gb = gb(CC);
Aexpg = A.*exp(1i*k*(g+B));

%% auxiliary variables
aa = sqrt(1i./(M1-M2.^2./M4));    % q := 1i./aa_old
bb = C1-C2.*M2./M4;
% cc = bb.^2 + 4*(g + C2.^2/4./M4).*aa; % iszero
% g2 = [g2; gg(:)];
% us = zeros(size(B));
notdone = 1; n = nmin; uspold = zeros(size(B)); Ipartold = uspold;
totalsum = totalsum + length(B); uspart = []; 
xdonepart = []; ydonepart = []; ndonepart = [];
% for nn = 1:5;
tic
while ~isempty(notdone) && n <= nmax;
    usp = steepest(psif,k,n,sid4,d24,xi,aa,xix);
%     us = [us; usp];
    Ipart = Aexpg.*usp; % partial integral
    done = abs(Ipart-Ipartold)<tol2;
    donesum = donesum + sum(done);
%     uspart = [uspart; Ipart(done)];% A(done).*usp(done)];
%     xdonepart = [xdonepart; xix(done)]; % save the stat points
%     ydonepart = [ydonepart; xiyx(done)];    % coordinates
%     ndonepart = [ndonepart; n*ones(sum(done),1)]; % how many p point to resolve
    notdone = ~(done);% + abs(Ipart)>1e20 + isnan(Ipart));   % discard diverging
%     B = B(notdone); 
%     C1 = C1(notdone); C2 = C2(notdone);
%     M1 = M1(notdone); 
%     M2 = M2(notdone); 
    M4 = M4(notdone);
    d24 = d24(notdone); xi = xi(notdone); xix = xix(notdone);
    sid4 = sid4(notdone);
% %     xiyx =xiyx(notdone); 
% %     g = g(notdone); 
    Aexpg = Aexpg(notdone);
% %     A = A(notdone);
    aa = aa(notdone); bb = bb(notdone); %cc = cc(notdone);
    n = 2*n;
%     uspold = usp(notdone);
    Ipartold = Ipart(notdone);
    summ = summ + 1/(2*pi/k)^2*sum(sum(Ipart(done)))*cgap^4;
end
    notdonesum = notdonesum + sum(notdone);
%     us = [us; uspart]; 
%     ndone = [ndone; ndonepart];
%     xdone = [xdone; xdonepart]; ydone = [ydone; ydonepart];
%     us2 = [us2; us(:)];
%     us = A.*us;
%     usA = [usA; us(:)];
tt = toc;
        t1 = t1+tt;
    end
end
%     summ = summ + 1/(2*pi/k)^2*sum(sum(us))*cgap^4;
% t1 = toc;
% issue a warning if not all integrals resolved
if donesum ~= totalsum
warning([num2str(notdonesum),' not resolved, ', ...
    num2str(totalsum-donesum-notdonesum),' discarded out of ',...
    num2str(totalsum)])
end
% disp(totalsum)
tot = [tot totalsum];
time1 = [time1 t1];
sumdata = [sumdata, summ];
end
%%
% load('Qdata.mat'); load('time2.mat');

% loglog(1./kdata,time1,'*-',...
%     1./kdata,1e-2*kdata.^1,'--',...
%     1./kdata,1e-2*kdata.^2,'-.',...
%     1./kdata,1e-3*kdata.^3,':')
% legend('NSDM','\epsilon^{-1}',...
%     '\epsilon^{-2}','\epsilon^{-3}')
% title('Computation cost NSDM','FontSize',16)
% xlabel('\epsilon'); ylabel('tic-toc time');
% % print('-dpdf','time_NSDM.pdf')
% 
% figure
% loglog(1./kdata,tot,'*-',1./kdata,1e3*kdata,'--',1./kdata,1e2*kdata.^2,'-.')
% legend('# of integrals','\epsilon^{-1}','\epsilon^{-2}')
% xlabel('\epsilon'); ylabel('# of integrals')
% title('Number of integrals','FontSize',16)
% % print('-dpdf','nrintegrals2.pdf')


loglog(...%1./kdata,time2(1:length(kdata)),'o-',...
    1./kdata,time1,'*-',...
    1./kdata,1e-2*kdata.^1,'--',...
    1./kdata,1e-2*kdata.^2,'-.',...
    1./kdata,1e-3*kdata.^3,':')
legend(...%'trap',...
    'NSDM','\epsilon^{-1}',...
    '\epsilon^{-2}','\epsilon^{-3}')
title('Computation cost trapezoidal rule vs NSDM','FontSize',16)
xlabel('\epsilon'); ylabel('tic-toc time');
% print('-dpdf','time3.pdf')

% figure
% loglog(1./kdata,abs((sumdata-Qdata(1:length(kdata)))./sumdata),'o-')
% xlabel('\epsilon'); ylabel('relative error')
% title('Relative error of trapezodial rule to NSDM','FontSize',16)
% % print('-dpdf','error11.pdf')

figure
loglog(1./kdata,tot,'*-',1./kdata,1e3*kdata,'--',1./kdata,1e2*kdata.^2,'-.')
legend('# of integrals','\epsilon^{-1}','\epsilon^{-2}')
xlabel('\epsilon'); ylabel('# of integrals')
title('Number of integrals','FontSize',16)
% print('-dpdf','nrintegrals3.pdf')

% disp(totalsum)
% 
%%
% disp([Q, summ]);
% disp(['Relative error ',num2str( abs((Q-summ)/summ))]);
