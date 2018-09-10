% QOI computed using the NSDM
close all
kdata = 10*2.^(0:5);    % wave numbers % out of memory for 0:8

time2 = []; 
Qdata = []; 
time12 = []; time22 = []; time23 = [];
tot2 =[];
Qdata3 = [];
sumdata2 = [];

tol = 1e-5;  % radius based tolerance
tol3 = 1e-11; % test function radius
cdata = [];

tol2 = 1e-11; % how precise should NSDM be
nmax = 256; % max number of GL points
nmin = 1;   % min number of GL points

% test function 
x0 = 0;
psif = @(x,y) exp(-5*((x-x0).^2+y.^2)); 
Rpsi = sqrt(-1/5*log(tol3));   % 'radius' of the test func


for k = kdata;
epsilon = 1/k;      % 1/wave number  
Y1 = .8; Y2 = .75; Y3 = 1;  % former random variables, example 1
% Y1 = 0; Y2 = 0; Y3 = 0;   % example 2
% Y1 = 0.1; Y2 = 0.2; Y3 = 0.3; % example 3

%% GB preliminaries
ppw = 10;               % points per wavelength
T = 3;                  % final time 
Lx = 4; Ly = 4;         % length of the domain

% initial amplitude
d = 10;

% % ONE PULSE
% sx = -.5;             % initial bump position
% a = @(x,y) exp(-d*((x-sx).^2 + y.^2));
% phi = @(x,y) -x;
% px = @(x,y) -1*ones(size(x));   % p_x
% py = @(x,y) 0*ones(size(x));    % p_y
% m00 = @(x,y) (0 + 1i)*ones(size(x)); % p_xx + 1i
% m10 = @(x,y) 0*ones(size(x)); % p_xy
% m11 = @(x,y) (0 + 1i)*ones(size(x)); % p_yy + 1i

% % two pulses
sx1 = -1; sx2 = 1;      % centers
a = @(x,y) exp(-d*((x-sx1).^2 + y.^2)) + exp(-d*((x-sx2).^2 + y.^2));
phi = @(x,y) abs(x);    % phi_0
px = @(x,y) sign(x).*ones(size(x));   % p_x
py = @(x,y) 0*ones(size(x));    % p_y
m00 = @(x,y) (0 + 1i)*ones(size(x)); % p_xx + 1i
m10 = @(x,y) 0*ones(size(x)); % p_xy
m11 = @(x,y) (0 + 1i)*ones(size(x)); % p_yy + 1i


%% %%% NO NEED TO CHANGE ANYTHING UNDER THIS POINT %%%%%%%%%%%%%

[xx,yy,h,ssx,ssy,cgap] = GBgrid(ppw,Lx,Ly,epsilon);      % GB grid points
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

%% uncomment if you want a slow trapezoidal rule reference solution
% tic        
% Q = GBtrap(xx,yy,nnn,ssx,ssy,a,phi,k,m00,m10,m11,px,py,...
%     xnew1,ynew1,A1,phi01,m1,m2,m4,pnew1,pnew2,cgap,psif,h,epsilon);
% Qdata = [Qdata Q];
% t2 = toc;
% time2 = [time2 t2];


%% old GB solution from Paper 2
      % % initiate auxiliary variables
%       uu = zeros(size(xx));
%         xxx = xx(:); xxx = repmat(xxx,1,nnn);
%         xxs1 = xxx - repmat(xnew1,;
%         yys1 = yy - ynew1;
% sumc = 0;
% t3 = 0;
%       % Loop over GB solutions 
% for ii = 1:length(ssx);
%          uuu = U1(end,ii:nnn:end); uuu=[uuu(1:6) uuu(6) uuu(7:8)];  % GB parameters      
%          s1 = ssx(ii);
%          s2 = ssy(ii);
%          phi0 = phi(s1,s2);               % phi0
% %          tic
%         [vend1,c1,tt] = GBsol(epsilon,xx,yy,phi0,uuu,x0,Rpsi,tol);
%         sumc = sumc + length(nonzeros(c1));
%         % Sum solution values inside QOI support
%         uu(c1) = uu(c1) + vend1;
% %         tt = toc;
% 
% t3 = t3 + tt;
% end
%         cdata = [cdata sumc];
%       %% Compute Q by trapezoidal rule
%       uu = 1/(2*pi*epsilon)*cgap^2*uu;
%       Q = h^2*sum(sum(abs(uu).^2.*psif(xx-x0,yy)));   % quantity of interest
% %       t3 = toc;
%       Qdata3 = [Qdata3 Q];
%       time3 = [time3 t3];

%% eigenvalues
tic   % tic, pair 1

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
summ2 = 0;

us2 = []; donesum2 = 0; notdonesum2 = 0; totalsum2 = 0;
xdone2 = []; ydone2 = []; ndone2 = [];
t22 = 0; t23 = 0;
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
 %% choose only some

%   % which combinations of points are relevant
    [XX1,XX2] = meshgrid(xnew1(ind),xnew1(ind2));
    [YY1,YY2] = meshgrid(ynew1(ind),ynew1(ind2));
tStart = tic;   % tic, pair 2
    [RR1,RR2] = meshgrid(R(ind),R(ind2));
    CC = (XX2-XX1).^2+(YY2-YY1).^2<=(RR1+RR2).^2;  % combinations of z and z' that overlap
    [j2,j1] = find(CC); % find nonzero indices
tt = toc(tStart);   % toc, pair 2

    j2 = j2 + ll*bs; j1 = j1 + kk*bs;   % shift the index w the current batch

    ba2 = b1(j1); bb2 = b1(j2);
    c1a2 = c11(j1); c1b2 = c11(j2);
    c2a2 = c12(j1); c2b2 = c12(j2);
    m11a2 = 1/2*m1(j1); m11b2 = 1/2*m1(j2);
    m12a2 = 1/2*m2(j1); m12b2 = 1/2*m2(j2);
    m14a2 = 1/2*m4(j1); m14b2 = 1/2*m4(j2);
    A1a2 = A1(j1); A1b2 = A1(j2);
 
     B2 = ba2-conj(bb2);
    C12 = c1a2-conj(c1b2);
    C22 = c2a2-conj(c2b2);    
    M12 = m11a2-conj(m11b2);
    M22 = m12a2-conj(m12b2);
    M42 = m14a2-conj(m14b2);
    A2 = A1a2.*conj(A1b2);   
    
% stationary points and auxiliary variables
[d242,xi2,xix2,g2] = statpts(C12,C22,M12,M22,M42);
sid42 = sqrt(1i./M42);

%% only the stationary points that make it to the integral
Aexpg2 = A2.*exp(1i*k*(g2+B2));
aa2 = sqrt(1i./(M12-M22.^2./M42));    % q := 1i./aa_old
bb2 = C12-C22.*M22./M42;

notdone2 = 1; n = nmin; uspold2 = zeros(size(B2)); Ipartold2 = uspold2;
totalsum2 = totalsum2 + length(B2); uspart2 = []; 
xdonepart2 = []; ydonepart2 = []; ndonepart2 = [];

tStart2 = tic;  % tic, pair 3
while ~isempty(notdone2) && n <= nmax;
    usp2 = steepest(psif,k,n,sid42,d242,xi2,aa2,xix2);
    Ipart2 = Aexpg2.*usp2; % partial integral
    done2 = abs(Ipart2-Ipartold2)<tol2;
    donesum2 = donesum2 + sum(done2);
    notdone2 = ~(done2);    
    M42 = M42(notdone2);
    d242 = d242(notdone2); xi2 = xi2(notdone2); xix2 = xix2(notdone2);
    sid42 = sid42(notdone2);
    Aexpg2 = Aexpg2(notdone2);
    aa2 = aa2(notdone2); bb2 = bb2(notdone2); %cc = cc(notdone);
   
    n = 2*n;    % number of GL points increases 
    Ipartold2 = Ipart2(notdone2);   % push only not finished integrals
    summ2 = summ2 + 1/(2*pi/k)^2*sum(sum(Ipart2(done2)))*cgap^4;
end
tt2 = toc(tStart2); % toc, pair 3
    notdonesum2 = notdonesum2 + sum(notdone2);

        t22 = t22+tt;
        t23 = t23+tt2;
    end
end
t12 = toc;  % toc, pair 2

% compute the number of intergals
if donesum2 ~= totalsum2
warning([num2str(notdonesum2),' not resolved, ', ...
    num2str(totalsum2-donesum2-notdonesum2),' discarded out of ',...
    num2str(totalsum2)])
end
% disp(totalsum)
tot2 = [tot2 totalsum2];
time12 = [time12 t12];
time22 = [time22 t22];
time23 = [time23 t23];
sumdata2 = [sumdata2, summ2];
end
%%
% reference solution and time saved
load('Qdata.mat'); load('time2.mat'); % caustic case
% load('Qdata2.mat'); load('time22.mat'); % non-caustic case

ln = length(time2);

loglog(1./kdata(1:ln), time2, 'o-',... %1./kdata,time2(1:length(kdata)),'o-',...
    1./kdata,time12,'*-',...%     1./kdata,time3,'s-',... 1./kdata,time22,'s-',...1./kdata,time23,'<-',...1./kdata,time12-time22,'>-',...
    1./kdata,1e-2*kdata.^1,'--',...
    1./kdata,1e-2*kdata.^2,'-.',...
    1./kdata,1e-3*kdata.^3,':')
lgd = legend('trap','NSDM',...%'sort','SD alone',...%'old',
    '\epsilon^{-1}','\epsilon^{-2}','\epsilon^{-3}');
set(lgd,'FontSize',16);
title('Computation cost trapezoidal rule vs NSDM','FontSize',16)
xlabel('\epsilon','FontSize',16); ylabel('tic-toc time','FontSize',16);
axis([1/kdata(end) 1/kdata(1) 1e-1 1e5])
% print('-dpdf','time2.pdf')

figure
loglog(1./kdata(1:ln),abs(sumdata2(1:ln)-Qdata),'o-')
% loglog(1./kdata,abs((sumdata2-Qdata(1:length(kdata)))./sumdata2),'o-',...
%     1./kdata, abs((Qdata3-Qdata(1:length(kdata)))./Qdata3))
xlabel('\epsilon','FontSize',16); ylabel('relative error','FontSize',16)
title('Relative error of trapezodial rule to NSDM','FontSize',16)
% print('-dpdf','error2.pdf')

figure
loglog(1./kdata,tot2,'*-',1./kdata,1e3*kdata,'--',1./kdata,1e2*kdata.^2,'-.')
lgd = legend('# of integrals','\epsilon^{-1}','\epsilon^{-2}');
set(lgd,'FontSize',16);
xlabel('\epsilon','FontSize',16); ylabel('# of integrals','FontSize',16)
title('Number of integrals','FontSize',16)
% print('-dpdf','nrintegrals2.pdf')
%