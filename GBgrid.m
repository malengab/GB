function [xx,yy,h,ssx,ssy,cgap] = GBgrid(ppw,Lx,Ly,eps)
% number of spatial steps
N = ceil(ppw*Lx/(2*pi*eps));   % number of nodes in x       
% Ny = ceil(ppw*Ly/(2*pi*eps));   % number of nodes in y       
h = Lx/(N-1);                  % spatial step length
x = -Lx/2:h:Lx/2;               % domain
y = -Ly/2:h:Ly/2;
[xx,yy] = meshgrid(x,y);

%% GB starting points
% coarser grid with ~sqrt(epsilon)-spacing
cgap = 0.75*sqrt(eps);            % gap between GB:s
xg = -Lx/2:cgap:Lx/2;
yg = -Ly/2:cgap:Ly/2;
[xxg,yyg] = meshgrid(xg,yg);

ssx = xxg(:);      % and concatenate them into a vector
ssy = yyg(:);
end
