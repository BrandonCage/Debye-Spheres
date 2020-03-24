clear all
close all

time_start = clock();

% Simulation Parameters
L=2*pi/3.0600; % Simulation box
DT=0.1;		   % time step (in wp^-1)
NT=1000;        % number of time steps
NG=64;		   % number of cells
dx=L/NG;       % grid spacing
N=1000;		   % number of particles
WP=1;          % plasma frequency
QM=-1;         % charge to mass ratio (= 1 for electrons)
V0=0.2;        % drift velocity
VT=0.0;        % thermal velocity
Q=WP^2/(QM*N/L); % charge of the computational particle
% Protons (ions) do not move during the simulation. This approximation is valid when the phenomena
% under study develops over time scale
% The plasma is neutral and ion contribute with their charge to neutralize the system
rho_back=-Q*N/L;

% Simulate the 2 stream instability: simulate two electron beam moving in opposite directions
xp=linspace(0,L-L/N,N)';
vp=VT*randn(N,1);
pm=[1:N]';
pm=1-2*mod(pm,2);
vp=vp+pm.*V0;




% Use a perturbation to initiate the instability
XP1=0.1; 
V1=0.0;
mode=1;
xp=xp+XP1*(L/N)*sin(2*pi*xp/L*mode);

p=1:N;p=[p p];
un=ones(NG-1,1);
Poisson=spdiags([un -2*un un],[-1 0 1],NG-1,NG-1);

% Output vectors
histTotEnergy = []; histKinEnergy = []; histPotEnergy = []; histMomentum = []; histThermalVelocity = []; 
histSpectrum = []; histEfield = []; histPhi    = [];
histPARTx = []; histPARTv = [];


h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated2d.gif';
c=[];
for n=1:N/2
    c=[c;1;10];
end

% start the computational cycle
for it=1:NT
    
scatter(xp,vp,20,c,'filled') 
axis([0 L -1 1])
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
	% calculate new particle position
	xp=xp+vp*DT;  
	% apply periodic boundary conditions
	out=(xp<0); xp(out)=xp(out)+L;
	out=(xp>=L);xp(out)=xp(out)-L;

	% interpolation particle to grid
	% to calculate the charge density on the grid
	% charge density is caculated on the grid cell centers
	g1=floor(xp/dx-.5)+1;
	g=[g1;g1+1];
	fraz1=1-abs(xp/dx-g1+.5);
	fraz=[fraz1;1-fraz1];	
	out=(g<1);g(out)=g(out)+NG;
	out=(g>NG);g(out)=g(out)-NG;
	mat=full(sparse(p,g,fraz,N,NG));
	rho=full((Q/dx)*sum(mat))'+rho_back; %'
	
	size(rho)
	% Field solver
	% The electrostatic potential Phi is calculated on the grid cell centers
	% by solving numercially \nabla^2 Phi = - rho
	size(Poisson);
	size(rho);
	Phi=Poisson\(-rho(1:NG-1)*dx^2);
	Phi=[Phi;0];
	size(Phi);

% The electric field is calculated on the grid cell centers
	% by central finite difference of Phi; E = - \nabla Phi
	Eg=([Phi(NG); Phi(1:NG-1)]-[Phi(2:NG);Phi(1)])/(2*dx);
	size(Eg);
	% interpolation particle to grid, using mat calculated 
    % in the interpolation grid to particle
	% and velocity update
	vp=vp+mat*QM*Eg*DT;
	
	Ekin   = 0.5*abs(Q)*sum(vp.^2);
	Efield = 0.5*sum(Eg.^2)*dx;
	Etot   =  Ekin + Efield ;
	histMomentum  = [histMomentum  sum(abs(Q)*vp)];
	histKinEnergy = [histKinEnergy Ekin];
	histPotEnergy   = [histPotEnergy Efield];
	histTotEnergy = [histTotEnergy Etot];
	histThermalVelocity = [histThermalVelocity cov(vp(1:2:end))]; % Take only the even particle (just one beam)
	histEfield = [histEfield Eg];
	histPhi    = [histPhi Phi];
	% take the Fourier transform of the electric field for spectral analysis
	NFFT = 2^nextpow2(length(Eg)); % Next power of 2 from length of Eg
	Y = fft(Eg,NFFT)/length(Eg);
	histSpectrum = [histSpectrum 2*abs(Y(1:NFFT/2+1))];
	if (mod(it,NT/4)==0) % save particle phase space every 1/4 of the simulation
			 histPARTx = [histPARTx xp(1:3:end)]; % take only every 3 particles
			 histPARTv = [histPARTv vp(1:3:end)]; % take only every 3 particles
	end
	
end
% timing
time_end = clock();
time_elapsed = etime(time_end,time_start)


time = linspace(0,NT*DT,NT);
k = 2*pi*(1/(2*dx))*linspace(0,1,NFFT/2+1);
space = linspace(0,L,NG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Examples of the plots you can make   %
%remove the % to use it               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
semilogy(time,histTotEnergy,time,histKinEnergy,time,histPotEnergy); xlabel('time'); ylabel('Energies'); legend('Tot En','Kin En','Pot En' );
%print("energy1D.png")
%replot
figure(2)
plot(time,histMomentum); xlabel('time'); ylabel('momentum'); % Momentum history
%figure(3)
%plot(time,histThermalVelocity); % Thermal velocity history
%figure(4)
%surf(time,space,-histPhi); xlabel('time'); ylabel('space'); zlabel('Phi'); % Potential history
%figure(5)
%hist(xp,NG); xlabel('x'); title('Number of computational electron per cell'); %electron space distribution at the last time
%figure(6)
%hist(histPARTx(:,2),NG); xlabel('x'); title('Number of computational electron per cell');
%figure(7)
%hist(vp,100); xlabel('v'); title('Electron distribution function'); % distribution function
%figure(8)
%hist(histPARTv(:,2),100); xlabel('v'); title('Electron distribution function'); % distribution function at half simulation
%figure(9)
%pcolor(time,space,histEfield); shading interp; xlabel('time'); ylabel('space'); colorbar;
%figure(10)
%pcolor(k,time,histSpectrum'); ylabel('time'); xlabel('k');
%figure(11)
%bar(k,histSpectrum(:,end)); xlabel('x'); % spectrum at the last cycle
%figure(12)
%plot(xp,vp,'.'); xlabel('xp'); ylabel('vp');     % phase space at last cycle'
%figure(13)
%plot(histPARTx(:,2),histPARTv(:,2),'.'); xlabel('xp'); ylabel('vp'); % phase space at half simulation







