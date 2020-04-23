clearvars 
close all

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated2d.gif';



L=1*2.284*10^10;%Length
e0=1;%Epsilon Naught
k=1.68637205*10^-10;%Boltzman
DT=.01*2.284*10^10;%DeltaT
NT=100;%Iterations
NG=10;%Number of Grid Points
nden=5;%Total charge of electrons in simulation

N=600;%Number of Particles
me=nden/N;%Mass of electron
mp=nden/N;%mass of proton
e=nden/N;%electron charge


VX1=0.0;%Velocity Petubation Magnitude
VY1=0.0;%Velocity Petubation Magnitude
mode=1;%wavelengths of sin
rho_back=0;%Background charge density
dx=L/NG;%Delta X
acc=.9999999;%Number close to 1 as to make the Poisson Matrix nonsingular

dE=0;%Initialization of Energy Graph Value
dKE=0;
dUE=0;
dxmom=0;
dymom=0;


%color
c=[linspace(1,1,N/2) linspace(10,10,N/2)];%Color matrix distinguishing electrons and ions

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


Temp=10000000;%Initial Temperature
Con=me/k/Temp;%Constant in Maxwellian Distribution
randvtot=(-2*log((1-rand(N,1)))/Con).^.5;%Initialize the speeds into a Maxwellian Distribution
randang=2*pi*rand(N,1);%Choose random direction for velocities
vxp=randvtot.*cos(randang);%put random speed and random direction together
vyp=randvtot.*sin(randang);

% Perturbation (Utilized in 2 Stream)
vxp=vxp+VX1*sin(2*pi*xp/L*mode);
vyp=vyp+VY1*sin(2*pi*yp/L*mode);

%Mass Matrix
massmat=mp*ones(1,600);%Make all into ions
parfor n=1:N
    if n<N/2+1
        massmat(n)=me;%Change the first half into electrons
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
UdUE=zeros(NT,1)
LdUE=zeros(NT,1);


% Main computational cycle
for it=1:NT
% update xp
scatter(xp,yp,20,c,'filled') 
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
xp=xp+vxp*DT;%Change particle position
yp=yp+vyp*DT;%Change particle position

% If position moves out of the square, loop it around to the other side
out=(xp<0); 
xp(out)=xp(out)+L;
out=(xp>=L);
xp(out)=xp(out)-L;
out=(yp<0); 
yp(out)=yp(out)+L;
out=(yp>=L);
yp(out)=yp(out)-L;


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

for n=1:N
    %Put the weights into a 3d matrix, where each row in the third
    %dimension is a new particle. Adding all allements on each row should
    %equal 1 in 4 adjacent indices
    if n==1
        weight=full(sparse([g(1,1) g(1,1) g(1+N,1) g(1+N,1)],[g(1,2) g(1+N,2) g(1,2) g(1+N,2)],[fraz(1,1) fraz(1+N,1) fraz(1,2) fraz(1+N,2)] ,NG,NG));%matrix stating weights of each particle on each grid point
    else
        weight=cat(3,weight,full(sparse([g(n,1) g(n,1) g(n+N,1) g(n+N,1)],[g(n,2) g(n+N,2) g(n,2) g(n+N,2)],[fraz(n,1) fraz(n+N,1) fraz(n,2) fraz(n+N,2)] ,NG,NG)));%matrix stating weights of each particle on each grid point
    end
end


%Compute the charge density
rho=(-((e)*sum(weight(:,:,1:N/2),3))'+((e)*sum(weight(:,:,N/2+1:N),3))')/dx^2+rho_back;%rho of each grid point
for n=1:NG
    if n==1
        rho1=rho(1:NG,1);
    else
        rho1=cat(1,rho1,rho(1:NG,n));
    end
end

% compute Potential with Poissons equation
Phi=Poisson\(-rho1(1:NG^2)*dx^2)/e0;
%Put Potential into a square Matrix, for visualation purposes
for n=1:NG
    if n==1
        Phi1=Phi(1:NG);
    else
        Phi1=cat(2,Phi1,Phi((n-1)*NG+1:n*NG));
    end
end

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

%Assemble into square matrix, for visualization
for n=1:NG
    if n==1
        Ex1=Ex(1:NG);
    else
        Ex1=cat(2,Ex1,Ex((n-1)*NG+1:n*NG));
    end
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

%Put into square matrix for visualization
for n=1:NG
    if n==1
        Ey1=Ey(1:NG);
    else
        Ey1=cat(2,Ey1,Ey((n-1)*NG+1:n*NG));
    end
end

% Update velocity
parfor n=1:N/2%For x velocity of electrons
    vxp(n)=vxp(n)+sum(weight(:,:,n).*Ex1','all')*e*DT/me;
end
parfor n=N/2+1:N%For x velocity of ions
    vxp(n)=vxp(n)-sum(weight(:,:,n).*Ex1','all')*e*DT/mp;
end
parfor n=1:N/2%For y velocity of electrons
    vyp(n)=vyp(n)+sum(weight(:,:,n).*Ey1','all')*e*DT/me;
end
parfor n=N/2+1:N%For y velocity of ions
    vyp(n)=vyp(n)-sum(weight(:,:,n).*Ey1','all')*e*DT/mp;
end

%Momentum test (Check to see if momentum is constant
xmom(it)=sum(massmat*vxp);
ymom(it)=sum(massmat*vyp);

%Energy test (Check if Total energy is constant)

totEn(it)=.5*sum(massmat*(vxp.^2+vyp.^2),'all')+.5*dx^2*e0*sum(Ex1.^2+Ey1.^2,'all');

%Potential energy
totUEn(it)=.5*dx^2*e0*sum(Ex1.^2+Ey1.^2,'all');

%Kinetic Energy
totKEn(it)=.5*sum(massmat*(vxp.^2+vyp.^2),'all');

%Display Iteration
it


Temp=(sum(massmat*(vxp.^2+vyp.^2),'all')/N)/2/k;%Calculate the temperature
LD =(e0*totKEn(it)/(4*pi*(nden/L^2)*e^2))^.5;%Calculate the debye length
if .5*LD<dx%Compare to delta x, if too small, display error
    "Error, Mesh not fine enough or not hot enough"
    LD
    dx
end
PP=(pi*me*L^2/(N*e^2))^.5;%Compute plasma period 
if .25*PP<DT%If plasma period is small display error
    "Error, large time step"
    PP 
end

%CFL Condition
C=DT/dx;
if C>1%If time step is larger than dx
    "May be unstable, C>1"
    C
end

end
figure%Plot Energy
plot(totEn)
hold on
%plot(dKE+UdUE)
%plot(dKE+LdUE)
plot(totUEn)
plot(totKEn)
%plot(UdUE)
%plot(LdUE)
legend("Energy","Potential","Kinetic")
hold off

figure%Plot momentums
plot(xmom)
legend("x Momentum")
figure
plot(ymom)
legend("y Momentum")

figure%Plot Energy
plot(totEn)
legend("Energy")
