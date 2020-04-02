clearvars 
close all

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated2d.gif';


dE=0;
dKE=0;
dUE=0;
dxmom=0;
dymom=0;
iterations=0;

L=1*2.284*10^10;%Length
e0=1;%Epsilon Naught
k=1.68637205*10^-10;
DT=.001*2.284*10^10;%DeltaT
NT=200;%Iterations
NTOUT=25;%Not used???
NG=15;%Number of Grid Points
nden=10000;

N=600;%Number of Particles
me=nden/N;%Mass of electron
mp=nden/N;%mass of proton
e=nden/N;%electron charge
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
Q=N*e/L^2;%WP^2/(QM*N/L);
rho_back=0;%Q*N/L^2;
dx=L/NG;%Delta X
BoundIter=1;
acc=.9999999;





% initial loading for the 2 Stream instability
%xp=linspace(0,L-L/N,N)';
%xp=linspace(0,L/2,N)';
xp=L*rand(N,1);
%yp=linspace(0,L-L/N,N)';
%yp=L*rand(N,1)/2;
yp=L*rand(N,1);
%xp=[0;L/2]
%yp=[.05+L/2;.05+L/2]

Temp=10000000000000;
Con=me/k/Temp;
randvtot=(-2*log((1-rand(N,1)))/Con).^.5;
randang=2*pi*rand(N,1);
vxp=randvtot.*cos(randang);
vyp=randvtot.*sin(randang);
pm=[1:N]';
pm=1-2*mod(pm,2);
vxp=vxp;%+pm.*V0;
vyp=vyp;%+pm.*V0;
% Perturbation
vxp=vxp+VX1*sin(2*pi*xp/L*mode);
vyp=vyp+VY1*sin(2*pi*yp/L*mode);

massmat=mp*ones(1,N);
for n=1:N
    if n<N/2+1
        massmat(n)=me;
    end
end

p=1:N;p=[p p;p p];
v=acc*ones(NG^2,1);
v1=ones(NG^2,1);
Poisson=-4*diag(v1)+diag(v(1:NG^2-1),1)+diag(v(1:NG^2-1),-1)+diag(v(1:NG^2-NG),-NG)+diag(v(1:NG^2-NG),NG);
Poisson=Poisson+diag(v(1:NG),NG^2-NG)+diag(v(1:NG),-NG^2+NG);


for i=1:NG
    Poisson=Poisson+full(sparse(NG*i,NG*i-NG+1,acc,NG^2,NG^2))+full(sparse(NG*i-NG+1,NG*i,acc,NG^2,NG^2));
end
for i=1:(NG-1)
    Poisson=Poisson-full(sparse(NG*i,NG*i+1,acc,NG^2,NG^2))-full(sparse(NG*i+1,NG*i,acc,NG^2,NG^2));
end

Ex=zeros(NG^2,1);
Ey=zeros(NG^2,1);
xmom=zeros(NT,1);
ymom=zeros(NT,1);
totEn=zeros(NT,1);
totUEn=zeros(NT,1);
totKEn=zeros(NT,1);
UdUE=zeros(NT,1);
LdUE=zeros(NT,1);


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

In1=(((xp-sxc1).^2+(yp-syc1).^2)<sr1^2);
In2=(((xp-sxc2).^2+(yp-syc2).^2)<sr2^2);

%color
c1=[linspace(1,1,N/2)];
c1(In1(1:N/2))=3;
c1(In2(1:N/2))=3;
c2=[linspace(10,10,N/2)];
c2(In1(N/2+1:N))=8;
c2(In2(N/2+1:N))=8;
c=[c1 c2];
% Main computational cycle
for it=1:NT
% update xp
scatter(xp,yp,20,c,'filled') 

hold on
sx1=[sxc1-sr1:2*sr1/100:sxc1+sr1];
sy1=(sr1.^2-(sx1-sxc1).^2).^(1/2)+syc1;
nsy1=-(sr1.^2-(sx1-sxc1).^2).^(1/2)+syc1;
sx2=[sxc2-sr2:2*sr2/100:sxc2+sr2];
sy2=real((sr2.^2-(sx2-sxc2).^2).^(1/2)+syc2);
nsy2=real(-(sr2.^2-(sx2-sxc2).^2).^(1/2)+syc2);
plot(sx1,sy1,'r')
plot(sx1,nsy1,'r')
plot(sx2,sy2,'r')
plot(sx2,nsy2,'r')
axis([0 L 0 L])
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
hold off
xp=xp+vxp*DT;
yp=yp+vyp*DT;
% apply bc on the particle positions
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
%Apply Sphere Bounce
out=(((xp-sxc1).^2+(yp-syc1).^2)<sr1^2);
out=xor(out,In1);
dfc=(((xp-sxc1).^2+(yp-syc1).^2)).^.5;
vxptemp(out)=vxp(out)-2*(vxp(out).*(xp(out)-sxc1)+vyp(out).*(yp(out)-syc1)).*(xp(out)-sxc1)./(dfc(out).^2);
vyptemp(out)=vyp(out)-2*(vxp(out).*(xp(out)-sxc1)+vyp(out).*(yp(out)-syc1)).*(yp(out)-syc1)./(dfc(out).^2);
vxp(out)=vxptemp(out);
vyp(out)=vyptemp(out);
xp(out)=xp(out)+2*(sr1-dfc(out)).*(xp(out)-sxc1)./dfc(out);
yp(out)=yp(out)+2*(sr1-dfc(out)).*(yp(out)-syc1)./dfc(out);

out=(((xp-sxc2).^2+(yp-syc2).^2)<sr2^2);
out=xor(out,In2);
dfc=(((xp-sxc2).^2+(yp-syc2).^2)).^.5;
vxptemp(out)=vxp(out)-2*(vxp(out).*(xp(out)-sxc2)+vyp(out).*(yp(out)-syc2)).*(xp(out)-sxc2)./(dfc(out).^2);
vyptemp(out)=vyp(out)-2*(vxp(out).*(xp(out)-sxc2)+vyp(out).*(yp(out)-syc2)).*(yp(out)-syc2)./(dfc(out).^2);
vxp(out)=vxptemp(out);
vyp(out)=vyptemp(out);
xp(out)=xp(out)+2*(sr2-dfc(out)).*(xp(out)-sxc2)./dfc(out);
yp(out)=yp(out)+2*(sr2-dfc(out)).*(yp(out)-syc2)./dfc(out);

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
for n=1:BoundIter
out=(g<1);
g(out)=g(out)+NG;
out=(g>NG);
g(out)=g(out)-NG;
end

for n=1:N
    if n==1
        weight=sparse([index2(g(1,1),g(1,2),NG) index2(g(1,1),g(N+1,2),NG) index2(g(N+1,1),g(1,2),NG) index2(g(N+1,1),g(N+1,2),NG)],[1 1 1 1],[fraz1(1) fraz2(1) fraz3(1) fraz4(1)] ,NG^2,N);%matrix stating weights of each particle on each grid point
    else
        weight=weight+sparse([index2(g(n,1),g(n,2),NG) index2(g(n,1),g(N+n,2),NG) index2(g(N+n,1),g(n,2),NG) index2(g(N+n,1),g(N+n,2),NG)],[n n n n],[fraz1(n) fraz2(n) fraz3(n) fraz4(n)],NG^2,N);%matrix stating weights of each particle on each grid point
    end
end


rho=(-((e)*sum(weight(:,1:N/2),2))+((e)*sum(weight(:,N/2+1:N),2)))/dx^2+rho_back;%rho of each grid point
%test=sum(rho1,'all')

% computing fields
Phi=Poisson\(-rho(1:NG^2)*dx^2)/e0;
%Electric Field in x direciton
for n=1:NG^2
        in1=n+NG;
        in2=n-NG;
        if in1>NG^2
            in1=in1-NG^2;
        end
        if in2<1
            in2=in2+NG^2;
        end
        Ex(n)=(Phi(in1)-Phi(in2))/(2*dx);
end


%Electric Field in Y Direction
for n=1:NG^2
        in1=n+1;
        in2=n-1;
        if mod(in1,NG)==1
            in1=in1-NG;
        end
        if mod(in2,NG)==0
            in2=in2+NG;
        end
        Ey(n)=(Phi(in1)-Phi(in2))/(2*dx);
end


% projection q->p and update of vp
for n=1:N/2
    vxp(n)=vxp(n)+sum(weight(:,n).*Ex','all')*e*DT/me;
end
for n=N/2+1:N
    vxp(n)=vxp(n)-sum(weight(:,n).*Ex','all')*e*DT/mp;
end
for n=1:N/2
    vyp(n)=vyp(n)+sum(weight(:,n).*Ey','all')*e*DT/me;
end
for n=N/2+1:N
    vyp(n)=vyp(n)-sum(weight(:,n).*Ey','all')*e*DT/mp;
end

%Momentum test
xmom(it)=sum(massmat*vxp);
ymom(it)=sum(massmat*vyp);

%Energy test

totEn(it)=.5*sum(massmat*(vxp.^2+vyp.^2),'all')+.5*dx^2*e0*sum(Ex.^2+Ey.^2,'all');

totUEn(it)=.5*dx^2*e0*sum(Ex.^2+Ey.^2,'all');


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
end
UdUE(it)=.5*e0*dx^2*findmax;

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
end
LdUE(it)=.5*e0*dx^2*findmin;



totKEn(it)=.5*sum(massmat*(vxp.^2+vyp.^2),'all');


it


Temp=(sum(massmat*(vxp.^2+vyp.^2),'all')/N)/2/k;
LD =(e0*totKEn(it)/(4*pi*(nden/L^2)*e^2))^.5;
if .5*LD<dx
    "Error, Mesh not fine enough or not hot enough"
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
hold on
plot(totEn)
%plot(dKE+UdUE)
%plot(dKE+LdUE)
plot(totUEn)
plot(totKEn)
%plot(UdUE)
%plot(LdUE)
legend("Energy","Potential","Kinetic")
hold off

figure
hold on
plot(xmom)
plot(ymom)
hold off

figure
plot(totEn)