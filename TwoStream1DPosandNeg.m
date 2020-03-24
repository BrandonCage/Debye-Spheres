clearvars 
close all

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated2d.gif';

Ex1=zeros(10,10)

dE=0;
dKE=0;
dUE=0;
UdUE=0;
LdUE=0;
dxmom=0;
dymom=0;
iterations=0;

L=2*pi;%Length
e0=1;%Epsilon Naught
k=1.68637205*10^-10;
DT=.005*2*pi;%DeltaT
NT=400;%Iterations
NTOUT=25;%Not used???
NG=16;%Number of Grid Points
nden=400;
N=1000;%Number of Particles
me=nden/N;%Mass of electron
mp=100*nden/N;%mass of proton
e=nden/N;%electron charge
WP=1;
QM=N*e;
V0=0.2;
VXT=0.0;
VYT=0.0;
XP1=60;
YP1=0;
VX1=0;
VY1=0;
mode=1;
Q=N*e/L^2%WP^2/(QM*N/L);
dx=L/NG;%Delta X
rho_back=0;
BoundIter=10;
acc=.999999;




%color
c=[linspace(1,1,N/2) linspace(10,10,N/2)];

% initial loading for the 2 Stream instability
%xp=linspace(0,L-L/N,N)';
xp=[linspace(0,L-L*2/N,N/2)';linspace(0,L-L*2/N,N/2)'];
%xp=L*rand(N,1);
%yp=linspace(0,L-L/N,N)';
%yp=L*rand(N,1)/2;
%yp=L*rand(N,1);
%xp=[0;L/2]
%yp=[.05+L/2;.05+L/2]
yp=(L/2-L/2/NG)*ones(N,1);


Temp=10000000000000;
Con=me/k/Temp;
randvtot=(-2*log((1-rand(N,1)))/Con).^.5;
randang=2*pi*rand(N,1);
vxp=randvtot.*cos(randang);
vyp=randvtot.*sin(randang);
pm=[1:N]';
pm=1-2*mod(pm,2);
vxp=[ones(N/2,1);-ones(N/2,1)];%+pm.*V0;
vyp=zeros(N,1);%+pm.*V0;
% Perturbation
vxp=vxp+VX1*sin(2*pi*xp/L*mode);
vyp=vyp+VY1*sin(2*pi*yp/L*mode);
xp=xp-XP1*(L/N)*sin(2*pi*xp/L*mode);


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
Poisson=Poisson+diag(v(1:NG),NG^2-NG)+diag(v(1:NG),-NG^2+NG);

i=1;

for i=1:NG
    Poisson=Poisson+full(sparse(NG*i,NG*i-NG+1,acc,NG^2,NG^2))+full(sparse(NG*i-NG+1,NG*i,acc,NG^2,NG^2));
end
for i=1:(NG-1)
    Poisson=Poisson-full(sparse(NG*i,NG*i+1,acc,NG^2,NG^2))-full(sparse(NG*i+1,NG*i,acc,NG^2,NG^2));
end
InvPoi=inv(Poisson);

% Main computational cycle
for it=1:NT-1
% update xp
scatter(xp,vxp,20,c,'filled') 
axis([0 L -3 3])
%plot([Ex1(5,1:10) Ex1(5,1:10)])



%hold off;

%Cone=me/k/Temp;
%Conp=mp/k/Temp;
%xtemp=[0:.04:2];
%ytemp=N*.04*Cone*xtemp.*exp(-Cone*xtemp.^2/2);
%vtot=[(vxp(1:N).^2+vyp(1:N).^2).^.5;4];
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
      if it~=1
        xp=xp+vxp*DT;
      end
yp=yp;%+vyp*DT;
% apply bc on the particle positions
n=1;
for n=1:BoundIter
out=(xp<0); 
xp(out)=xp(out)+L;
out=(xp>=L);
xp(out)=xp(out)-L;
out=(yp<0); 
yp(out)=yp(out)+L;
out=(yp>=L);
yp(out)=yp(out)-L;
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

testmat=weight(1:NG,2,1:N)
rho=(-((e)*sum(weight(:,:,1:N/2),3))'+((e)*sum(weight(:,:,N/2+1:N),3))')/dx^2+rho_back;%rho of each grid point
plot(rho(2,1:NG))
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
for n=1:N/2
    vxp(n)=vxp(n)+sum(weight(:,:,n).*Ex1','all')*e*DT/me;
    n=n+1;
end
for n=N/2+1:N
    vxp(n)=vxp(n)-sum(weight(:,:,n).*Ex1','all')*e*DT/mp;
    n=n+1;
end
for n=1:N/2
    vyp(n)=vyp(n)+sum(weight(:,:,n).*Ey1','all')*e*DT/me;
    n=n+1;
end
for n=N/2+1:N
    vyp(n)=vyp(n)-sum(weight(:,:,n).*Ey1','all')*e*DT/mp;
    n=n+1;
end

%Momentum test
xmom=sum(massmat*vxp,'all')
if xmom==[0]
    dxmom=xmom;
else
dxmom=[dxmom;sum(massmat*vxp,'all')];
end
ymom=sum(massmat*vyp,'all')
if ymom==[0]
    dymom=ymom;
else
dymom=[dymom;ymom];
end

%Energy test

totEn=.5*sum(massmat*(vxp.^2+vyp.^2),'all')+.5*dx^2*e0*sum(Ex1.^2+Ey1.^2,'all')
if dE==[0]
    dE=totEn;
else
dE=[dE;.5*sum(massmat*(vxp.^2+vyp.^2),'all')+.5*dx^2*e0*sum(Ex1.^2+Ey1.^2,'all')];
end

totUEn=.5*dx^2*e0*sum(Ex1.^2+Ey1.^2,'all');
if dUE==[0]
    dUE=totUEn;
else
dUE=[dUE;totUEn];
end



UtotUEn=0;
for n=1:NG^2
    in1=n;
    in2=n+1;
    in3=n+NG;
    in4=n+NG+1;
    if mod(n,NG)==0
        in2=in2-NG;
        in4=in4-NG;
    end
    if in3>NG^2
        in3=in3-NG;
        in4=in4-NG;
    end
    findmax=max([Ex(in1)^2+Ey(in1)^2;Ex(in2)^2+Ey(in2)^2;Ex(in3)^2+Ey(in3)^2;Ex(in4)^2+Ey(in4)^2]);
    UtotUEn=UtotUEn+findmax;
end
UtotUEn=.5*e0*dx^2*UtotUEn;
if UdUE==[0]
    UdUE=UtotUEn;
else
UdUE=[UdUE;UtotUEn];
end

LtotUEn=0;
for n=1:NG^2
    in1=n;
    in2=n+1;
    in3=n+NG;
    in4=n+NG+1;
    if mod(n,NG)==0
        in2=in2-NG;
        in4=in4-NG;
    end
    if in3>NG^2
        in3=in3-NG;
        in4=in4-NG;
    end
    findmin=min([Ex(in1)^2+Ey(in1)^2;Ex(in2)^2+Ey(in2)^2;Ex(in3)^2+Ey(in3)^2;Ex(in4)^2+Ey(in4)^2]);
    LtotUEn=LtotUEn+findmin;
end
LtotUEn=.5*e0*dx^2*LtotUEn;
if LdUE==[0]
    LdUE=LtotUEn;
else
LdUE=[LdUE;LtotUEn];
end


totKEn=.5*sum(massmat*(vxp.^2+vyp.^2),'all');
if dKE==[0]
    dKE=totKEn;
else
    dKE=[dKE;totKEn];
end

it
totUEn/totKEn;

Temp=(sum(massmat*(vxp.^2+vyp.^2),'all')/N)/2/k;
LD =(e0*totKEn/(4*pi*(nden/L^2)*e^2))^.5;
if .5*LD<dx
    "Error, Mesh not fine enough or not hot enough"
    dx
    LD
end
PP=(pi*me*L^2/(N*e^2))^.5;
if .25*PP<DT
    "Error, large time step"
    PP 
end

%CFL Condition
C=DT/dx;
if C>1
    "May be unstable, C>1"
    C
end

end
figure
plot(dE)
hold on
%plot(dKE+UdUE)
%plot(dKE+LdUE)
plot(dUE)
plot(dKE)
%plot(UdUE)
%plot(LdUE)
legend("Energy","Potential","Kinetic")
hold off

figure
hold on
plot(dxmom(2:NT))
%plot(dymom(2:NT))
hold off

figure
plot(dE)