clear all
close all

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated2.gif';

dE=[0];
dKE=[0];
dUE=[0];

L=2*pi;%Length
me=9.10938356*10^-31;%Mass of electron
DT=.5;%DeltaT
NT=50;%Iterations
NTOUT=25;%Not used???
NG=32;%Number of Grid Points
N=100;%Particles
WP=1;
QM=-1;
V0=0.0;
VT=0.0;
XP1=1;
V1=0.0;
mode=1;
Q=WP^2/(QM*N/L);
rho_back=-Q*N/L;
dx=L/NG;%Delta X

%color
c=linspace(1,10,N);

% initial loading for the 2 Stream instability
%xp=linspace(0,L-L/N,N)';
xp=L*rand(N,1)
yp=linspace(0,0,N)';
vp=VT*randn(N,1);
pm=[1:N]';pm=1-2*mod(pm,2);
vp=vp+pm.*V0;
% Perturbation
vp=vp+V1*sin(2*pi*xp/L*mode);
xp=xp+XP1*(L/N)*sin(2*pi*xp/L*mode);
p=1:N;p=[p p];
un=ones(NG-1,1);
Poisson=spdiags([un -2*un un],[-1 0 1],NG-1,NG-1)
% Main computational cycle
for it=1:NT
% update xp
scatter(xp,yp,20,c,'filled') 
axis([0 2*pi -1 1])
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
xp=xp+vp*DT;
% apply bc on the particle positions
out=(xp<0); 
xp(out)=xp(out)+L;
out=(xp>=L);
xp(out)=xp(out)-L;
% projection p->g
g1=ceil(xp/dx-.5);%put on grid points
g=[g1;g1+1];%put on grid points
fraz1=1-abs(xp/dx-g1+.5);%weighting
fraz=[fraz1;1-fraz1];%weighting

% apply bc on the projection
out=(g<1);
g(out)=g(out)+NG;
out=(g>NG);
g(out)=g(out)-NG;

mat=sparse(p,g,fraz,N,NG);%matrix stating weights of each particle on each grid point
rho=full((Q/dx)*sum(mat))'+rho_back%rho of each grid point

% computing fields
Phi=Poisson\(-rho(1:NG-1)*dx^2)
Phi=[Phi;0];
Eg=([Phi(NG); Phi(1:NG-1)]-[Phi(2:NG);Phi(1)])/(2*dx);
% projection q->p and update of vp
vp=vp+mat*QM*Eg*DT;
it


%Energy test

totEn=.5*sum((vp.^2),'all')+.5*dx*sum(Eg.^2,'all');
if dE==[0]
    dE=totEn;
else
    dE=[dE;totEn];
end
totKEn=.5*sum((vp.^2),'all');
if dKE==[0]
    dKE=totKEn;
else
    dKE=[dKE;totKEn];
end
totUEn=.5*dx*sum(Eg.^2,'all');
if dUE==[0]
    dUE=totUEn;
else
    dUE=[dUE;totUEn];
end
end
figure
plot(dE)
hold on
plot(dKE)
plot(dUE)

