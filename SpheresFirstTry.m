clear all
close all

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated2d.gif';


L=7;%Length
me=9.10938356*10^-31;%Mass of electron
DT=.05;%DeltaT
NT=500;%Iterations
NTOUT=25;%Not used???
NG=8;%Number of Grid Points
N=500;%Number of Particles
WP=1;
QM=-1;
V0=0.2;
VXT=0.0;
VYT=0.0;
XP1=1;
YP1=1;
VX1=0.0;
VY1=0.0;
mode=1;
Q=WP^2/(QM*N/L);
    
dx=L/NG;%Delta X
BoundIter=1;
acc=1

%sphere1
sxc1=L/4;
syc1=L/4;
sr1=L/8;
vsx1=0;
vsy1=0;


sx1=[sxc1-sr1:2*sr1/100:sxc1+sr1];
sy1=(sr1.^2-(sx1-sxc1).^2).^(1/2)+syc1;
nsy1=-(sr1.^2-(sx1-sxc1).^2).^(1/2)+syc1;

%sphere2
sxc2=3*L/4;
syc2=3*L/4;
sr2=L/8;
vsx2=0;
vsy2=0;

sx2=[sxc2-sr2:2*sr2/100:sxc2+sr2];
sy2=(sr2.^2-(sx2-sxc2).^2).^(1/2)+syc2;
nsy2=-(sr2.^2-(sx2-sxc2).^2).^(1/2)+syc2;


rho_back=Q*N/(L^2-pi*(sr1^2+sr2^2))+.5*[0 0 0 0 0 0 0 0
    0 0 -1 0 0 0 0 0
    0 -1 -4 -1 0 0 0 0
    0 0 -1 0 0 0 0 0
    0 0 0 0 0 1 0 0
    0 0 0 0 1 4 1 0
    0 0 0 0 0 1 0 0
    0 0 0 0 0 0 0 0]

%color
c=linspace(1,10,N);

% initial positions
%xp=linspace(0,L-L/N,N)';
%xp=linspace(0,L/2,N)';
xp=L*rand(N,1);
%yp=linspace(0,L-L/N,N)';
%yp=L*rand(N,1)/2;
yp=L*rand(N,1);
vxp=VXT*randn(N,1);
vyp=VYT*randn(N,1);
pm=[1:N]';
pm=1-2*mod(pm,2);
vxp=vxp;%+pm.*V0;
vyp=vyp;%+pm.*V0;
% Perturbation
vxp=vxp+VX1*sin(2*pi*xp/L*mode);
vyp=vyp+VY1*sin(2*pi*yp/L*mode);
p=1:N;p=[p p;p p];
v=acc*ones(NG^2,1);
v1=ones(NG^2,1);
Poisson=-4*diag(v1)+diag(v(1:NG^2-1),1)+diag(v(1:NG^2-1),-1)+diag(v(1:NG^2-NG),-NG)+diag(v(1:NG^2-NG),NG);


i=1;

for i=1:(NG-1)
    Poisson=Poisson-full(sparse(NG*i,NG*i+1,acc,NG^2,NG^2))-full(sparse(NG*i+1,NG*i,acc,NG^2,NG^2));
end
InvPoi=inv(Poisson);


% Main computational cycle
for it=1:NT
% update xp

scatter(xp,yp,20,'b','filled')
hold on
sx1=[sxc1-sr1:2*sr1/100:sxc1+sr1];
sy1=(sr1.^2-(sx1-sxc1).^2).^(1/2)+syc1;
nsy1=-(sr1.^2-(sx1-sxc1).^2).^(1/2)+syc1;
sx2=[sxc2-sr2:2*sr2/100:sxc2+sr2];
sy2=(sr2.^2-(sx2-sxc2).^2).^(1/2)+syc2;
nsy2=-(sr2.^2-(sx2-sxc2).^2).^(1/2)+syc2;
plot(sx1,sy1,'r')
plot(sx1,nsy1,'r')
plot(sx2,sy2,'r')
plot(sx2,nsy2,'r')
axis([0 L 0 L])
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if it == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
      hold off
xp=xp+vxp*DT;
yp=yp+vyp*DT;
sxc1=sxc1+vsx1*DT;
syc1=syc1+vsy1*DT;
sxc2=sxc2+vsx2*DT;
syc2=syc2+vsy2*DT;
% apply bc on the particle positions
n=1;
for n=1:BoundIter
out=(xp<0); 
xp(out)=-xp(out);
vxp(out)=-vxp(out);
out=(xp>=L);
xp(out)=-xp(out)+2*L;
vxp(out)=-vxp(out);
out=(yp<0); 
yp(out)=-yp(out);
vyp(out)=-vyp(out);
out=(yp>=L);
yp(out)=-yp(out)+2*L;
vyp(out)=-vyp(out);
if sxc1<0
    sxc1=sxc1+L;
elseif sxc1>L
    sxc1=sxc1-L
end
if sxc2<0
    sxc2=sxc2+L;
elseif sxc1>L
    sxc2=sxc2-L
end
end
n=1;
for n=1:BoundIter
out=(((xp-sxc1).^2+(yp-syc1).^2)<sr1^2);
dfc=(((xp-sxc1).^2+(yp-syc1).^2)).^.5;
vxp(out)=vxp(out)-2*(vxp(out).*(xp(out)-sxc1)+vyp(out).*(yp(out)-syc1)).*(xp(out)-sxc1)./(dfc(out).^2);
vyp(out)=vyp(out)-2*(vxp(out).*(xp(out)-sxc1)+vyp(out).*(yp(out)-syc1)).*(yp(out)-syc1)./(dfc(out).^2);
xp(out)=xp(out)+2*(sr1-dfc(out)).*(xp(out)-sxc1)./dfc(out);
yp(out)=yp(out)+2*(sr1-dfc(out)).*(yp(out)-syc1)./dfc(out);
end
for n=1:BoundIter
out=(((xp-sxc2).^2+(yp-syc2).^2)<sr2^2);
dfc=(((xp-sxc2).^2+(yp-syc2).^2)).^.5;
vxp(out)=vxp(out)-2*(vxp(out).*(xp(out)-sxc2)+vyp(out).*(yp(out)-syc2)).*(xp(out)-sxc2)./(dfc(out).^2);
vyp(out)=vyp(out)-2*(vxp(out).*(xp(out)-sxc2)+vyp(out).*(yp(out)-syc2)).*(yp(out)-syc2)./(dfc(out).^2);
xp(out)=xp(out)+2*(sr2-dfc(out)).*(xp(out)-sxc2)./dfc(out);
yp(out)=yp(out)+2*(sr2-dfc(out)).*(yp(out)-syc2)./dfc(out);
end


% projection p->g
g1=cat(2,ceil(xp/dx-.5),ceil(yp/dx-.5));%put on grid points
g=[g1;g1+1];%put on grid points
%fraz1=1-abs(xp/dx-g1+.5);%weighting
%fraz=[fraz1;1-fraz1]%weighting
fraz1=((1-abs(xp/dx-g1(1:N,1)+.5))'*diag((1-abs(yp/dx-g1(1:N,2)+.5))'))';
fraz2=((abs(xp/dx-g1(1:N,1)+.5))'*diag((1-abs(yp/dx-g1(1:N,2)+.5))'))';
fraz3=((1-abs(xp/dx-g1(1:N,1)+.5))'*diag((abs(yp/dx-g1(1:N,2)+.5))'))';
fraz4=((abs(xp/dx-g1(1:N,1)+.5))'*diag((abs(yp/dx-g1(1:N,2)+.5))'))';
fraz=[fraz1 fraz2;fraz3 fraz4];

% apply bc on the projection
n=1;
for n=1:BoundIter
out=(g<1);
g(out)=g(out)+NG;
out=(g>NG);
g(out)=g(out)-NG;
n=n+1;
end

n=1;
for n=1:N
    if n==1
        weight=full(sparse([g(1,1) g(1,1) g(1+N,1) g(1+N,1)],[g(1,2) g(1+N,2) g(1,2) g(1+N,2)],[fraz(1,1) fraz(1+N,1) fraz(1,2) fraz(1+N,2)] ,NG,NG));%matrix stating weights of each particle on each grid point
    else
        weight=cat(3,weight,full(sparse([g(n,1) g(n,1) g(n+N,1) g(n+N,1)],[g(n,2) g(n+N,2) g(n,2) g(n+N,2)],[fraz(n,1) fraz(n+N,1) fraz(n,2) fraz(n+N,2)] ,NG,NG)));%matrix stating weights of each particle on each grid point
    end
    n=n+1;
end


rho=-((Q/dx^2)*sum(weight,3))'+rho_back;%rho of each grid point
n=1;
for n=1:NG
    if n==1
        rho1=rho(1:NG,1);
    else
        rho1=cat(1,rho1,rho(1:NG,n));
    end
    n=n+1;
end
test=sum(rho1)

% computing fields
Phi=InvPoi*(-rho1(1:NG^2)*dx^2);
n=1;
for n=1:NG
    if n==1
        Phi1=Phi(1:NG);
    else
        Phi1=cat(2,Phi1,Phi((n-1)*NG+1:n*NG));
    end
    n=n+1;
end
%Electric Field in x direciton
n=1;
for n=1:NG^2
    if n==1
        Ex=[(Phi(1+NG))/(2*dx)];
    else
        in1=n+NG;
        in2=n-NG;
        if in1>NG^2
            Ex=[Ex;(-Phi(in2))/(2*dx)];
        elseif in2<1
            Ex=[Ex;(Phi(in1))/(2*dx)];
        else
            Ex=[Ex;(Phi(in1)-Phi(in2))/(2*dx)];
        end
    end
    n=n+1;
end
n=1;
for n=1:NG
    if n==1
        Ex1=Ex(1:NG);
    else
        Ex1=cat(2,Ex1,Ex((n-1)*NG+1:n*NG));
    end
    n=n+1;
end


%Electric Field in Y Direction
n=1;
for n=1:NG^2
    if n==1
        Ey=[(Phi(2)-Phi(NG))/(2*dx)];
    else
        in1=n+1;
        in2=n-1;
        if mod(in1,NG)==1
            Ey=[Ey;(-Phi(in2))/(2*dx)];
        elseif mod(in2,NG)==0
            Ey=[Ey;(Phi(in1))/(2*dx)];
        else
            Ey=[Ey;(Phi(in1)-Phi(in2))/(2*dx)];
        end
    end
    n=n+1;
end
n=1;
for n=1:NG
    if n==1
        Ey1=Ey(1:NG);
    else
        Ey1=cat(2,Ey1,Ey((n-1)*NG+1:n*NG));
    end
    n=n+1;
end



EEnergy=sum(Ex1.^2+Ey1.^2,'all')

% projection q->p and update of vp
n=1;
for n=1:N
    vxp(n)=vxp(n)+sum(weight(:,:,n).*Ex1','all')*QM*DT;
    n=n+1;
end
for n=1:N
    vyp(n)=vyp(n)+sum(weight(:,:,n).*Ey1','all')*QM*DT;
    n=n+1;
end
end