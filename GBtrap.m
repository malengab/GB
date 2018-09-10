function Q = GBtrap(xx,yy,nnn,ssx,ssy,a,phi,k,m00,m10,m11,px,py,...
    xnew1,ynew1,A1,phi01,m1,m2,m4,pnew1,pnew2,cgap,psif,h,epsilon)
%% trapezoidal rule for GB solution     
% plot initial data
v0 = zeros(size(xx));
for ind = 1:nnn;
v00 = a(ssx(ind),ssy(ind))*exp(1i*k*(phi(ssx(ind),ssy(ind))+...
    (xx-ssx(ind))*px(ssx(ind),ssy(ind))+(yy-ssy(ind))*py(ssx(ind),ssy(ind))+...
    1/2*m00(ssx(ind),ssy(ind))*(xx-ssx(ind)).^2 + ...
    m10(ssx(ind),ssy(ind)).*(xx-ssx(ind)).*(yy-ssy(ind))+...
    + 1/2*m11(ssx(ind),ssy(ind)).*(yy-ssy(ind)).^2));
v0 = v0 + v00;
end
v0 = 1/(2*pi*epsilon)*cgap^2*v0;
figure;
surf(xx,yy,abs(v0));
shading interp
colorbar
view([0 90])
title('T = 0','FontSize',16)
axis([-2 2 -2 2])
xlabel('x','FontSize',16); ylabel('y','FontSize',16)
% print('-dpdf','sol02.pdf')
print('-dpdf','sol04.pdf')

%% compare one GB results
% tic
% slow sol to compare
vv = zeros(size(xx));
for ind = 1:length(ssx);
v = A1(ind)*exp(1i*k*(phi01(ind) + (xx-xnew1(ind))*pnew1(ind) + (yy-ynew1(ind))*pnew2(ind)+...
    1/2*m1(ind)*(xx-xnew1(ind)).^2 + m2(ind).*(xx-xnew1(ind)).*(yy-ynew1(ind))+...
    + 1/2*m4(ind).*(yy-ynew1(ind)).^2));
vv = vv + v;
end
v1 = 1/(2*pi*epsilon)*cgap^2*vv;
figure;
surf(xx,yy,abs(v1));
shading interp
colorbar
view([0 90])
title('T = 1','FontSize',16)
xlabel('x','FontSize',16); ylabel('y','FontSize',16)
axis([-2 2 -2 2])
% print('-r1200','-dpdf','sol2.pdf')
print('-dpdf','sol4.pdf')
Q  = sum(sum(abs(v1).^2.*psif(xx,yy)))*h^2;
% t2 = toc;
% Qdata = [Qdata Q];
% time2 = [time2 t2];