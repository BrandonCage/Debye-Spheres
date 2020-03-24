clear all
close all

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated2d.gif';

totEn=0;
dE=0;
iterations=0;

L=150;%Length
me=9.10938356*10^-31;%Mass of electron
mp=1.6726219*10^-27;%mass of proton
e=1.60217662*10^-19;%electron charge
e0=8.854*10^-12;%Epsilon Naught
k=1.38064852*10^-23;%Boltzman
DT=.01;%DeltaT
NT=1500;%Iterations
NTOUT=25;%Not used???
NG=10;%Number of Grid Points
N=1000;%Number of Particles
WP=1;
QM=N*e;
V0=0.2;
VXT=0.0;
VYT=0.0;
XP1=1;
YP1=1;
VX1=0.0;
VY1=0.0;
mode=1;
Q=N*e/L^2%WP^2/(QM*N/L);
rho_back=0%Q*N/L^2;
dx=L/NG;%Delta X
BoundIter=10;
acc=1;



%color
col=[linspace(1,1,N/2) linspace(10,10,N/2)];

% initial loading for the 2 Stream instability
%xp=linspace(0,L-L/N,N)';
%xp=linspace(0,L/2,N)';
xp=L*rand(N,1);
%yp=linspace(0,L-L/N,N)';
%yp=L*rand(N,1)/2;
yp=L*rand(N,1);

Con=1;
Temp=10^-31;
randvtot=(-2*log((1-rand(N,1)))/Con).^.5;
randang=2*pi*rand(N,1);
vxp=0%randvtot.*cos(randang);
vyp=0%randvtot.*sin(randang);
pm=[1:N]';
pm=1-2*mod(pm,2);
vxp=vxp;%+pm.*V0;
vyp=vyp;%+pm.*V0;
% Perturbation
vxp=vxp+VX1*sin(2*pi*xp/L*mode);
vyp=vyp+VY1*sin(2*pi*yp/L*mode);

n=0;
for n=1:N
    if n==1
        massmat=[me];
    elseif n<N/2+1
        massmat=[massmat me];
    else
        massmat=[massmat mp];
    end
end

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
scatter(xp,yp,20,col,'filled') 
axis([0 L 0 L])




%hold off;

%Cone=me/k/Temp;
%Conp=mp/k/Temp;
%xtemp=[0:.04:2];
%ytemp=N*.04*Cone*xtemp.*exp(-Cone*xtemp.^2/2);
%vtot=[(vxp(1:N).^2+vyp(1:N).^2).^.5;4]
%histogram([vtot;.4],100)
%axis([0 2 0 N/5])
%hold on;
%plot(xtemp,ytemp)
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
xp=xp+vxp*DT;
yp=yp+vyp*DT;
% apply bc on the particle positions
n=1;
for n=1:BoundIter
out=(xp<0);
j=find(out);
if isempty(j)==0
    weight(:,:,j)=[];
end
xp(out)=[];
yp(out)=[];
vxp(out)=[];
vyp(out)=[];
massmat(out)=[];
col(out)=[];
N=N-sum(out,'all');

out=(xp>=L);
j=find(out);
if isempty(j)==0
    weight(:,:,j)=[];
end
xp(out)=[];
yp(out)=[];
vxp(out)=[];
vyp(out)=[];
massmat(out)=[];
col(out)=[];
N=N-sum(out,'all');

out=(yp<0); 
j=find(out);
if isempty(j)==0
    weight(:,:,j)=[];
end
xp(out)=[];
yp(out)=[];
vxp(out)=[];
vyp(out)=[];
massmat(out)=[];
col(out)=[];
N=N-sum(out,'all');

out=(yp>=L);
j=find(out);
if isempty(j)==0
    weight(:,:,j)=[];
end
xp(out)=[];
yp(out)=[];
vxp(out)=[];
vyp(out)=[];
massmat(out)=[];
col(out)=[];
N=N-sum(out,'all');
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


rho=-((e)*sum(weight(:,:,1:N/2),3))'+((e)*sum(weight(:,:,N/2+1:N),3))'+rho_back;%rho of each grid point
n=1;
for n=1:NG
    if n==1
        rho1=rho(1:NG,1);
    else
        rho1=cat(1,rho1,rho(1:NG,n));
    end
    n=n+1;
end
%test=sum(rho1,'all')

% computing fields
Phi=InvPoi*(-rho1(1:NG^2)*dx^2)/e0;
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
        Ex=[(Phi(1+NG)-Phi(NG^2-NG+1))/(2*dx)];
    else
        in1=n+NG;
        in2=n-NG;
        if in1>NG^2
            in1=in1-NG^2;
        end
        if in2<1
            in2=in2+NG^2;
        end
        Ex=[Ex;(Phi(in1)-Phi(in2))/(2*dx)];
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
            in1=in1-NG;
        end
        if mod(in2,NG)==0
            in2=in2+NG;
        end
        Ey=[Ey;(Phi(in1)-Phi(in2))/(2*dx)];
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

% projection q->p and update of vp
n=1;

for n=1:N
    if massmat(n)==me
        vxp(n)=vxp(n)+sum(weight(:,:,n).*Ex1','all')*e*DT/me;
    else
        vxp(n)=vxp(n)-sum(weight(:,:,n).*Ex1','all')*e*DT/mp;
    end
    n=n+1;
end
n=1;
for n=1:N
    if massmat(n)==me
        vyp(n)=vyp(n)+sum(weight(:,:,n).*Ey1','all')*e*DT/me;
    else
        vyp(n)=vyp(n)-sum(weight(:,:,n).*Ey1','all')*e*DT/mp;
    end
    n=n+1;
end

%Momentum test
xmom=massmat*vxp
ymom=massmat*vyp

%Energy test

dE=[dE;.5*massmat*(vxp.^2+vyp.^2)+.5*e0*sum(Ex1.^2+Ey1.^2,'all')-totEn];
totEn=.5*massmat*(vxp.^2+vyp.^2)+.5*e0*sum(Ex1.^2+Ey1.^2,'all')
Temp=massmat*(vxp.^2+vyp.^2)/N/2/k
iterations=iterations+1;
end