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

L=1;%Length
e0=1;%Epsilon Naught
k=1.68637205*10^-10;
DT=.05;%DeltaT
NT=200;%Iterations
NTOUT=25;%Not used???
NG=15;%Number of Grid Points
nden=5;%Total charge of electrons in simulation


N=1000;%Number of Particles
me=nden/N;%Mass of electron
mp=nden/N;%mass of proton
e=nden/N;%electron charge

Nin=N;
QM=N*e;


VX1=0.0;%Velocity Petubation Magnitude
VY1=0.0;%Velocity Petubation Magnitude
mode=1;%wavelengths of sin
rho_back=0;%Background charge density
rho_sphere1=0;
rho_sphere2=0;
dx=L/NG;%Delta X
acc=.9999999;%Number close to 1 as to make the Poisson Matrix nonsingular

Nin=N;

% initial loading for Particle Positions
%(Commented out some other configurations)
%xp=linspace(0,L-L/N,N)';
%xp=linspace(0,L/2,N)';
xp=L*rand(N,1);
%yp=linspace(0,L-L/N,N)';
%yp=L*rand(N,1)/2;
yp=L*rand(N,1);
%xp=[0;L/2]
%yp=[.05+L/2;.05+L/2]


Temp=1000000;%Initial Temperature
Con=me/k/Temp;%Constant in Maxwellian Distribution
randvtot=(-2*log((1-rand(N,1)))/Con).^.5;%Initialize the speeds into a Maxwellian Distribution
randang=2*pi*rand(N,1);%Choose random direction for velocities
vxp=randvtot.*cos(randang);%put random speed and random direction together
vyp=randvtot.*sin(randang);


% Perturbation (Utilized in 2 Stream)
vxp=vxp+VX1*sin(2*pi*xp/L*mode);
vyp=vyp+VY1*sin(2*pi*yp/L*mode);


%Mass Matrix
massmat=mp*ones(1,N);%Make all into ions
parfor n=1:N
    if n<N/2+1
        massmat(n)=me;%Change the first half into electrons
    end
end
%Charge Matrix
emat=e*ones(1,N);%Make all into ions
for n=1:N
    if n<N/2+1
        emat(n)=-e;%Change the first half into electrons
    end
end

%Build Poisson Matrix
v=acc*ones(NG^2,1);%Vector of values close to one
v1=ones(NG^2,1);%Vector of values equal to one
%Put -4 on the diagonal and 1 on surrounding values
Poisson=-4*diag(v1)+diag(v(1:NG^2-1),1)+diag(v(1:NG^2-1),-1)+diag(v(1:NG^2-NG),-NG)+diag(v(1:NG^2-NG),NG);
Poisson=Poisson+diag(v(1:NG),NG^2-NG)+diag(v(1:NG),-NG^2+NG);

%Put ones where there are some missing
for i=1:NG
    Poisson=Poisson+full(sparse(NG*i,NG*i-NG+1,acc,NG^2,NG^2))+full(sparse(NG*i-NG+1,NG*i,acc,NG^2,NG^2));
end
%Remove ones from places where they shouldn't be
for i=1:(NG-1)
    Poisson=Poisson-full(sparse(NG*i,NG*i+1,acc,NG^2,NG^2))-full(sparse(NG*i+1,NG*i,acc,NG^2,NG^2));
end

%Preallocate memory
Ex=zeros(NG^2,1);
Ey=zeros(NG^2,1);
xmom=zeros(NT,1);
ymom=zeros(NT,1);
totEn=zeros(NT,1);
totUEn=zeros(NT,1);
totKEn=zeros(NT,1);
UdUE=zeros(NT,1);
LdUE=zeros(NT,1);
lostenergy=0;


bins=10;%Bins for charge around the circle
s1cdx=linspace(2*pi/bins,2*pi,bins);%Show Angle
s1cdy=zeros(bins,1);%Initialize at 0 charge
s2cdx=linspace(2*pi/bins,2*pi,bins);%Show Angle
s2cdy=zeros(bins,1);%Initialize at 0 charge

%sphere1
sxc1=L/2;%Sphere 1 x center
syc1=L/2;%Sphere 1 y center
sr1=L/4;%Sphere 1 radius
vsx1=0;%X Velocity of sphere 1
vsy1=0;%Y velocity of sphere 1
ms1=1;%Mass of sphere


sx1=[sxc1-sr1:2*sr1/100:sxc1+sr1];%Data to plot the sphere on x
sy1=(sr1.^2-(sx1-sxc1).^2).^(1/2)+syc1;%Data to plot sphere on y
nsy1=-(sr1.^2-(sx1-sxc1).^2).^(1/2)+syc1;%Data to plot negative sphere on y

%sphere2
sxc2=8*L/10;%Sphere 2 x center
syc2=L/2;%Sphere 2 y center
sr2=L/20000;%Sphere 2 radius
vsx2=0;%X Velocity of sphere 2
vsy2=0;%Y velocity of sphere 
ms2=1;%Mass of sphere

sx2=[sxc2-sr2:2*sr2/100:sxc2+sr2];%Data to plot the sphere on x
sy2=(sr2.^2-(sx2-sxc2).^2).^(1/2)+syc2;%Data to plot sphere on y
nsy2=-(sr2.^2-(sx2-sxc2).^2).^(1/2)+syc2;%Data to plot negative sphere on y

In1=(((xp-sxc1).^2+(yp-syc1).^2)<sr1^2);%Particles in Sphere 1
In2=(((xp-sxc2).^2+(yp-syc2).^2)<sr2^2);%Particles in sphere 2

emat(In1)=0;%Eliminate charge in sphere
emat(In2)=0;%Elimnate charge in sphere


%color
c1=[linspace(1,1,N/2)];%Ions
c2=[linspace(10,10,N/2)];%Electrons
c=[c1 c2];%&PUT it all together

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
axis([L/NG L-L/NG L/NG L-L/NG])

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

%Particle Delete if in a bordering grid space or out of bounds
out=xp<L/NG;
xp(out)=[];
yp(out)=[];
vxp(out)=[];
vyp(out)=[];
emat(out)=[];
massmat(out)=[];
c(out)=[];
N=N-sum(out,'all');

out=xp>L-L/NG;
xp(out)=[];
yp(out)=[];
vxp(out)=[];
vyp(out)=[];
emat(out)=[];
massmat(out)=[];
c(out)=[];
N=N-sum(out,'all');

out=yp<L/NG;
xp(out)=[];
yp(out)=[];
vxp(out)=[];
vyp(out)=[];
emat(out)=[];
massmat(out)=[];
c(out)=[];
N=N-sum(out,'all');

out=yp>L-L/NG;
xp(out)=[];
yp(out)=[];
vxp(out)=[];
vyp(out)=[];
emat(out)=[];
massmat(out)=[];
c(out)=[];
N=N-sum(out,'all');

%Add Fresh Plasma
NGnew=4*NG-4;%Border Grid spaces getting replaced
Nnew=floor(NGnew*Nin/NG^2+.5);%Put in a number of particles to keep N about the same
randmat=randi([1 4],Nnew,1);%Which side will they appear?
in1=sum(randmat(:) == 1);%Put on Bottom
xp=[xp;(L-L/NG)*rand(in1,1)];
yp=[yp;(L/NG)*rand(in1,1)];
in1=sum(randmat(:) == 2);%Put on right
xp=[xp;L-(L/NG)*rand(in1,1)];
yp=[yp;(L-L/NG)*rand(in1,1)];
in1=sum(randmat(:) == 3);%Put on top
xp=[xp;L/NG+(L-L/NG)*rand(in1,1)];
yp=[yp;L-(L/NG)*rand(in1,1)];
in1=sum(randmat(:) == 4);%Put on left
xp=[xp;(L/NG)*rand(in1,1)];
yp=[yp;L/NG+(L-L/NG)*rand(in1,1)];

randvtot=(-2*log((1-rand(Nnew,1)))/Con).^.5;%Maxwellian
randang=2*pi*rand(Nnew,1);%Random direction
vxp=[vxp;randvtot.*cos(randang)];
vyp=[vyp;randvtot.*sin(randang)];

randmat=randi([0 1],Nnew,1);%Ion or Electron?

emat=[emat 2*e*randmat'-e];%Set new charges
massmat=[massmat (mp-me)*randmat'+me];%Set new masses
c=[c 9*randmat'+1];%Set correct colors

N=N+Nnew;%Add new particles to N


% Put particles on grid points
g1=cat(2,ceil(xp/dx-.5),ceil(yp/dx-.5));%put onto the lower grid points
g=[g1;g1+1];%Put onto lower and higher grid points
fraz1=((1-abs(xp/dx-g1(1:N,1)+.5))'*diag((1-abs(yp/dx-g1(1:N,2)+.5))'))';%lower left weight
fraz2=((abs(xp/dx-g1(1:N,1)+.5))'*diag((1-abs(yp/dx-g1(1:N,2)+.5))'))';%lower right weight
fraz3=((1-abs(xp/dx-g1(1:N,1)+.5))'*diag((abs(yp/dx-g1(1:N,2)+.5))'))';%upper left weight
fraz4=((abs(xp/dx-g1(1:N,1)+.5))'*diag((abs(yp/dx-g1(1:N,2)+.5))'))';%upper right weight
fraz=[fraz1 fraz2;fraz3 fraz4];%Put it all together now yall

% If it's on a nonexistant grip point, shift it to the other sides
out=(g<1);
g(out)=g(out)+NG;
out=(g>NG);
g(out)=g(out)-NG;

%Put the weights into a 3d matrix, where each row in the third
    %dimension is a new particle. Adding all allements on each row should
    %equal 1 in 4 adjacent indices
for n=1:N
    if n==1
        weight=sparse([index2(g(1,1),g(1,2),NG) index2(g(1,1),g(N+1,2),NG) index2(g(N+1,1),g(1,2),NG) index2(g(N+1,1),g(N+1,2),NG)],[1 1 1 1],[fraz1(1) fraz2(1) fraz3(1) fraz4(1)] ,NG^2,N);%matrix stating weights of each particle on each grid point
    else
        weight=weight+sparse([index2(g(n,1),g(n,2),NG) index2(g(n,1),g(N+n,2),NG) index2(g(N+n,1),g(n,2),NG) index2(g(N+n,1),g(N+n,2),NG)],[n n n n],[fraz1(n) fraz2(n) fraz3(n) fraz4(n)],NG^2,N);%matrix stating weights of each particle on each grid point
    end
end


%Compute the charge density
rho=sum(emat.*weight(:,1:N),2)/dx^2+rho_sphere1+rho_sphere2;%rho of each grid point



% compute Potential with Poissons equation
Phi=Poisson\(-rho*dx^2)/e0;

%Particle Delete
out=(((xp-sxc1).^2+(yp-syc1).^2)<sr1^2);%If in sphere
j=find(out);%Which are out
if isempty(j)==0%Add charge to sphere
    rho_sphere1=rho_sphere1+sum(emat(j).*weight(:,j),2)/dx^2;
    weight(:,j)=[];
end
%Determine where the charge hit and add it to the ghraph
angle1=(sign(-xp(out)+sxc1)+1)*.5*pi+(atan((yp(out)-syc1)./(xp(out)-sxc1)));
binnum=floor(angle1*bins/(2*pi)+bins/4+.5);
binnum(find(~binnum))=binnum(find(~binnum))+10;
s1cdy=s1cdy+sparse(binnum,1,emat(out),bins,1);
lostenergy=lostenergy+.5*massmat(out)*(vxp(out).^2+vyp(out).^2);%Kinetic energy of particles deleted
%Delete variables
xp(out)=[];
yp(out)=[];
vxp(out)=[];
vyp(out)=[];
emat(out)=[];
massmat(out)=[];
c(out)=[];
N=N-sum(out,'all');%Fix N

out=(((xp-sxc2).^2+(yp-syc2).^2)<sr2^2);%If in sphere
j=find(out);%Which are out
if isempty(j)==0%Add charge to sphere
    rho_sphere2=rho_sphere2+sum(emat(j).*weight(:,j),2)/dx^2;
    weight(:,j)=[];
end
%Determine where the charge hit and add it to the ghraph

angle2=(sign(-xp(out)+sxc2)+1)*.5*pi+(atan((yp(out)-syc2)./(xp(out)-sxc2)));
binnum=floor(angle2*bins/(2*pi)+bins/4+.5);
binnum(find(~binnum))=binnum(find(~binnum))+10;
s2cdy=s2cdy+sparse(binnum,1,emat(out),bins,1);
lostenergy=lostenergy+.5*massmat(out)*(vxp(out).^2+vyp(out).^2);
%Delete Variable
xp(out)=[];
yp(out)=[];
vxp(out)=[];
vyp(out)=[];
emat(out)=[];
massmat(out)=[];
c(out)=[];
N=N-sum(out,'all');%Fix N




%Electric Field in x direciton
parfor n=1:NG^2
        in1=n+NG;%Index of right node
        in2=n-NG;%Index of left node
        if in1>NG^2
            in1=in1-NG^2;%Fix if on right boundary
        end
        if in2<1
            in2=in2+NG^2;%Fix if on left boundary
        end
        Ex(n)=(Phi(in1)-Phi(in2))/(2*dx);%Electric field in x
end


%Electric Field in Y Direction
parfor n=1:NG^2
        in1=n+1;%Index of above node
        in2=n-1;%index of lower lode
        if mod(in1,NG)==1
            in1=in1-NG;%Fix if on the top
        end
        if mod(in2,NG)==0
            in2=in2+NG;%fix if on the bottom
        end
        Ey(n)=(Phi(in1)-Phi(in2))/(2*dx);%Electric field in y
end



% projection q->p and update of vp
    vxp=vxp-weight(:,1:N)'*Ex.*emat'*DT./massmat';
    vyp=vyp-weight(:,1:N)'*Ey.*emat'*DT./massmat';
    

%Momentum test
xmom(it)=sum(massmat*vxp);
ymom(it)=sum(massmat*vyp);

%Energy test

totEn(it)=.5*sum(massmat*(vxp.^2+vyp.^2),'all')+.5*dx^2*e0*sum(Ex.^2+Ey.^2,'all');

totEnWLost(it)=totEn(it)+lostenergy;

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

LeftSphereForceLateral(it)=-rho_sphere1'*Ex;
LeftSphereForceLongitude=-rho_sphere1'*Ey;
%RightSphereForceLateral(it)=-rho_sphere2'*Ex;
%RightSphereForceLongitude=-rho_sphere2'*Ey;

Temp=(sum(massmat*(vxp.^2+vyp.^2),'all')/N)/2/k;
LD =(e0*totKEn(it)/(4*pi*(nden/L^2)*e^2))^.5;
if .5*LD<dx
    "Error, Mesh not fine enough or not hot enough"
    LD
    dx
end
PP=(pi*me*L^2/(N*e^2))^.5;
if .25*PP<DT
    "Error, large time step"
    PP
    DT
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
plot(totEnWLost)
%plot(UdUE)
%plot(LdUE)
legend("Energy","Potential","Kinetic","Energy Including Heat")
hold off

figure
hold on
plot(xmom)
plot(ymom)
hold off

figure
plot(LeftSphereForceLateral)
legend("Force on Sphere1")
xlabel("time")
ylabel("Force (Positive is Attractive)")