/*
 * TestMat3.cpp
 *
 *  Created on: Nov 5, 2020
 *      Author: shadowcage72
 */






#include <iostream>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

int Matrix()
{
  MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  VectorXd v(2);
  v<<1,-1;
  return 0;
}

int index2(double y,double z,int NG)
{
	if(z>NG)
	{
		z-=NG;
	}
	if(z<1)
	{
		z+=NG;
	}
	if(y>NG)
	{
		y-=NG;
	}
	if(y<1)
	{
		y+=NG;
	}
	return (y-1)*NG+z-1;
}








int main()
{


	double e0=1;//Epsilon Naught
	double k=.000000000168637205;//Boltzman
	double meterconv=0.00000000000003511282;//Convert from our units to meters
	double L=9/meterconv;//Length
	double DT=L/150;//DeltaT
	int NT=500000;//Iterations
	int NG=80;//Number of Grid Points
	double NG1=NG;
	double ndensity=2000000000000*meterconv*meterconv*meterconv;
	double totalparticles=.001*2*ndensity*L*L;//Total number of electrons in simulation


	//At least 25 per cell is ideal

	int N=400000;//Number of Particles
	double me=totalparticles/N;//Mass of electron
	double mp=100*totalparticles/N;//mass of proton
	double e=totalparticles/N;//electron charge

	double Nin=N;
	int Nin1=N;
	double dx=L/NG;//Delta X
	double acc=1;//Number close to 1 as to make the Poisson Matrix nonsingular

	double dE=0;//Initialization of Energy Graph Value
	double dKE=0;
	double dUE=0;
	double dxmom=0;
	double dymom=0;
	double lostenergy=0;
	double esphere1=0;
	double esphere2=0;



	// initial loading for Particle Positions
	VectorXd xp(2*N);
	VectorXd yp(2*N);
	int n;
	for(n=0;n<N;n++)
	{
	xp(n)=L*((double) rand() / (RAND_MAX));
	yp(n)=L*((double) rand() / (RAND_MAX));
	xp(n+N)=0;
	yp(n+N)=0;
	};

	//Mass Matrix
	VectorXd massmat(2*N);
	for(n=0;n<N/2;n++)
		{
	massmat(n)=me;//Make all into ions
	massmat(n+N/2)=mp;
	massmat(n+N)=0;
	massmat(n+N+N/2)=0;
		}

	//Charge Matrix
	VectorXd emat(2*N);
		for(n=0;n<N/2;n++)
			{
		emat(n)=-e;//Make all into ions
		emat(n+N/2)=e;
		emat(n+N)=0;
		emat(n+N+N/2)=0;
			}

	double Tempin=60*11700;//Initial Temperature
	double Temp=Tempin;//Temperature that will change in time

	VectorXd vxp(2*N);
	VectorXd vyp(2*N);
	double pi=3.141592653589793;
	double vrand;
	double anglerand;
	for(n=0;n<N;n++)
	{
		vrand=sqrt((-2*log((1-((double) rand() / (RAND_MAX)))))/(massmat(n)/me/k/Temp));
		anglerand=2*pi*((double) rand() / (RAND_MAX));
		vxp(n)=vrand*cos(anglerand);
		vxp(n+N)=0;
		vyp(n)=vrand*sin(anglerand);//Make all into ions
		vyp(n+N)=0;
	}

	MatrixXd Poisson(NG*NG,NG*NG);
	for(int n=0;n<NG*NG;n++)
	{
		for(int m=0;m<NG*NG;m++)
		{
			Poisson(n,m)=0;
		}
	}
	for(n=0;n<NG*NG;n++)
	{
		Poisson(n,n)=-4;
		if(n>0)
		{
			Poisson(n,n-1)=1;
		}
		if(n<NG*NG-1)
		{
			Poisson(n,n+1)=1;
		}
		if(n>NG-1)
		{
			Poisson(n,n-NG)=1;
		}
		if(n<NG*NG-NG)
		{
			Poisson(n,n+NG)=1;
		}
	}

	for(n=1;n<NG;n++)
	{
		Poisson(NG*n-1,NG*n)=0;
		Poisson(NG*n,NG*n-1)=0;
	}


MatrixXd invPoi(NG*NG,NG*NG);
invPoi=Poisson.inverse();

//Preallocate memory
	VectorXd Ex(NG*NG);
	VectorXd Ey(NG*NG);
	VectorXd xmom(NT);
	VectorXd ymom(NT,1);
	VectorXd ForceXSphere1(NT,1);
	VectorXd ForceYSphere1(NT,1);
	VectorXd ForceXSphere2(NT,1);
	VectorXd ForceYSphere2(NT,1);
	VectorXd totEn(NT,1);
	VectorXd totUEn(NT,1);
	VectorXd totKEn(NT,1);
	VectorXd UdUE(NT,1);
	VectorXd LdUE(NT,1);
	VectorXd Volt1(NT,1);
	VectorXd Volt2(NT,1);
	VectorXd Charge2(NT,1);
	VectorXd Charge1(NT,1);
	MatrixXd rho1allit(NG,NG);
	MatrixXd Phi1allit(NG,NG);
	for(n=0;n<NG;n++)
	{
		for(int m=0;m<NG;m++)
		{
			rho1allit(n,m)=0;
			Phi1allit(n,m)=0;
		}
	}
	VectorXd TotalMomentumAbsorbed1X(NT,1);
	VectorXd TotalMomentumAbsorbed1Y(NT,1);
	VectorXd TotalMomentumAbsorbed2X(NT,1);
	VectorXd TotalMomentumAbsorbed2Y(NT,1);
	VectorXd Temp1(NT,1);
	VectorXd TotalCharge1(NT,1);
	VectorXd totEnWLost(NT,1);
	VectorXd ForceYmean(NT,1);
	MatrixXd rho1(NG,NG);
	MatrixXd Phi1(NG,NG);
	MatrixXd Ex1(NG,NG);
	MatrixXd Ey1(NG,NG);


	double LD =sqrt(k*Temp/ndensity);


	//sphere1
	double sxc1=.5*L;//Sphere 1 x center
	double syc1=.5*L;//Sphere 1 y center
	double sr1=.08*L;//Sphere 1 radius


	//sphere2
	double sxc2=.8654*L;//Sphere 2 x center
	double syc2=.55504*L;//Sphere 2 y center
	double sr2=.000000000001*L;//Sphere 2 radius

	for(n=0;n<N;n++)
	{
		if(((xp(n)-sxc1)*(xp(n)-sxc1)+(yp(n)-syc1)*(yp(n)-syc1))<sr1*sr1)
		{
			emat(n)=0;
		}
		if(((xp(n)-sxc2)*(xp(n)-sxc2)+(yp(n)-syc2)*(yp(n)-syc2))<sr2*sr2)
		{
			emat(n)=0;
		}
	}

	//Capacity Matrix

	//Put a particle on every grid point
	int i=0;
	int j=NG-1;
	VectorXd xpcap(NG*NG);
	VectorXd ypcap(NG*NG);
	for(int n=0;n<NG*NG;n++)
	{
		xpcap(n)=dx*j;
		ypcap(n)=dx*i;
		if(i==NG-1)
		{
			i=0;
			j--;
		}
		else
		{
			i++;
		}
	}

	//Now a particle is on each point




	VectorXd ematc1(NG*NG);
	VectorXd ematc2(NG*NG);
	VectorXd ematc3(NG*NG);
	int lengthin1=0;
	int lengthin2=0;
	int lengthin3=0;
	for(n=0;n<NG*NG;n++)
	{
		ematc3(n)=0;
		if(((xpcap(n)-sxc1)*(xpcap(n)-sxc1)+(ypcap(n)-syc1)*(ypcap(n)-syc1))<sr1*sr1)
		{
			ematc1(n)=1;
			ematc3(n)=1;
			lengthin1++;
			lengthin3++;
		}
		else
		{
			ematc1(n)=0;
		}

		if(((xpcap(n)-sxc2)*(xpcap(n)-sxc2)+(ypcap(n)-syc2)*(ypcap(n)-syc2))<sr2*sr2)
		{
			ematc2(n)=1;
			ematc3(n)=1;
			lengthin2++;
			lengthin3++;
		}
		else
		{
			ematc2(n)=0;
		}
	}
	cout <<"\nlengthin2="<<lengthin2;

	VectorXd rhotemp1(NG*NG);
	VectorXd Phitemp1(NG*NG);
	MatrixXd Bprime1(lengthin1,lengthin1);
	MatrixXd Bprime2(lengthin2,lengthin2);
	MatrixXd Bprime3(lengthin3,lengthin3);
	for(n=0;n<NG*NG;n++)
	{
		rhotemp1(n)=0;
	}
	i=0;
	j=0;
	for(n=0;n<NG*NG;n++)
	{
		if(ematc1(n)==1)
		{
			rhotemp1(n)=1/dx/dx;
			Phitemp1=-invPoi*rhotemp1*dx*dx/e0;
			for (int m=0;m<NG*NG;m++)
			{
				if(ematc1(m)==1)
				{
					Bprime1(i,j)=Phitemp1(m);
					i++;
				}
			}
			i=0;
			j++;
			rhotemp1(n)=0;
		}
	}
	MatrixXd Cap1=Bprime1.inverse();
	cout <<"\nCheckpoint2";
	i=0;
	j=0;
	for(n=0;n<NG*NG;n++)
	{
		if(ematc2(n)==1)
		{
			rhotemp1(n)=1/dx/dx;
			Phitemp1=-invPoi*rhotemp1*dx*dx/e0;
			for (int m=0;m<NG*NG;m++)
			{
				if(ematc2(m)==1)
				{
					Bprime2(i,j)=Phitemp1(m);
					i++;
				}
			}
			i=0;
			j++;
			rhotemp1(n)=0;
		}
	}
	MatrixXd Cap2=Bprime2.inverse();

	i=0;
	j=0;
	for(n=0;n<NG*NG;n++)
	{
		if(ematc3(n)==1)
		{
			rhotemp1(n)=1/dx/dx;
			Phitemp1=-invPoi*rhotemp1*dx*dx/e0;
			for (int m=0;m<NG*NG;m++)
			{
				if(ematc3(m)==1)
				{
					Bprime3(i,j)=Phitemp1(m);
					i++;
				}
			}
			i=0;
			j++;
			rhotemp1(n)=0;
		}
	}
	MatrixXd Cap3=Bprime3.inverse();

	cout <<"\nCheckpoint3";



	//Circular Matrix representing 1 particle around Sphere 1
		int NCirc=2;//Number of particles for circle
		double NCirc1=NCirc;
		VectorXd xpc(NCirc);
		VectorXd ypc(NCirc);
		for (n=0;n<NCirc;n++)
		{
			xpc(n)=sxc1;//+sr1*sin(2*pi/NCirc*n);
			ypc(n)=syc1;//+sr1*cos(2*pi/NCirc*n);
		}
		double charge=1/NCirc1;

		MatrixXd gc(2*NCirc,2);
		for (n=0;n<NCirc;n++)
		{
			gc(n,0)=floor(xpc(n)/dx+.5);
			gc(n,1)=floor(ypc(n)/dx+.5);
			gc(NCirc+n,0)=floor(xpc(n)/dx+.5)+1;
			gc(NCirc+n,1)=floor(ypc(n)/dx+.5)+1;
		}


		MatrixXd frazc(2*NCirc,2);
		for (n=0;n<NCirc;n++)
		{
			frazc(n,0)=(1-abs(xpc(n)/dx-gc(n,0)))*(1-abs(ypc(n)/dx-gc(n,1)));
			frazc(n,1)=(abs(xpc(n)/dx-gc(n,0)))*(1-abs(ypc(n)/dx-gc(n,1)));
			frazc(NCirc+n,0)=(1-abs(xpc(n)/dx-gc(n,0)))*(abs(ypc(n)/dx-gc(n,1)));
			frazc(NCirc+n,1)=(abs(xpc(n)/dx-gc(n,0)))*(abs(ypc(n)/dx-gc(n,1)));
		}
		cout <<"\nCheckpoint4";

		//Make the 1 particle sphere matrix
		VectorXd rho_sphere1(NG*NG);
		for (n=0;n<NG*NG;n++)
		{
			rho_sphere1(n)=0;
		}
		for(n=0;n<NCirc;n++)
		{
			rho_sphere1(index2(gc(n,0),gc(n,1),NG))+=charge*frazc(n,0)/dx/dx;
			rho_sphere1(index2(gc(n,0),gc(NCirc+n,1),NG))+=charge*frazc(n,1)/dx/dx;
			rho_sphere1(index2(gc(NCirc+n,0),gc(n,1),NG))+=charge*frazc(NCirc+n,0)/dx/dx;
			rho_sphere1(index2(gc(NCirc+n,0),gc(NCirc+n,1),NG))+=charge*frazc(NCirc+n,1)/dx/dx;
		}



		//Circular Matrix representing 1 particle around Sphere 2
			for (n=0;n<NCirc;n++)
			{
				xpc(n)=sxc2;//+sr2*sin(2*pi/NCirc*n);
				ypc(n)=syc2;//+sr2*cos(2*pi/NCirc*n);
			}

			for (n=0;n<NCirc;n++)
			{
				gc(n,0)=floor(xpc(n)/dx+.5);
				gc(n,1)=floor(ypc(n)/dx+.5);
				gc(NCirc+n,0)=floor(xpc(n)/dx+.5)+1;
				gc(NCirc+n,1)=floor(ypc(n)/dx+.5)+1;
			}


			for (n=0;n<NCirc;n++)
			{
				frazc(n,0)=(1-abs(xpc(n)/dx-gc(n,0)))*(1-abs(ypc(n)/dx-gc(n,1)));
				frazc(n,1)=(abs(xpc(n)/dx-gc(n,0)))*(1-abs(ypc(n)/dx-gc(n,1)));
				frazc(NCirc+n,0)=(1-abs(xpc(n)/dx-gc(n,0)))*(abs(ypc(n)/dx-gc(n,1)));
				frazc(NCirc+n,1)=(abs(xpc(n)/dx-gc(n,0)))*(abs(ypc(n)/dx-gc(n,1)));
			}

			//Make the 1 particle sphere matrix
					VectorXd rho_sphere2(NG*NG);
					for (n=0;n<NG*NG;n++)
					{
						rho_sphere2(n)=0;
					}
					for(n=0;n<NCirc;n++)
					{
						rho_sphere2(index2(gc(n,0),gc(n,1),NG))+=charge*frazc(n,0)/dx/dx;
						rho_sphere2(index2(gc(n,0),gc(NCirc+n,1),NG))+=charge*frazc(n,1)/dx/dx;
						rho_sphere2(index2(gc(NCirc+n,0),gc(n,1),NG))+=charge*frazc(NCirc+n,0)/dx/dx;
						rho_sphere2(index2(gc(NCirc+n,0),gc(NCirc+n,1),NG))+=charge*frazc(NCirc+n,1)/dx/dx;
					}

cout <<"\nCheckpoint5";

//
	//color
	VectorXd c(2*N);
	for (n=0;n<N/2;n++)
		{
		c(n)=1;
		c(n+N/2)=10;
		c(n+N)=0;
		c(n+N+N/2)=0;
		}


	ofstream myfile1;
		myfile1.open("4Volt1.txt");

		ofstream myfile2;
		myfile2.open("4Volt2.txt");

		ofstream myfile3;
		myfile3.open("4Charge1.txt");

		ofstream myfile4;
		myfile4.open("4Charge2.txt");


		ofstream myfile5;
		myfile5.open("4rho1allit.txt");
		myfile5 << rho1allit;
		myfile5.close();

		ofstream myfile6;
		myfile6.open("4Phi1allit.txt");
		myfile6 << Phi1allit;
		myfile6.close();

		ofstream myfile7;
		myfile7.open("4Temp1.txt");

		ofstream myfile8;
		myfile8.open("4velocitysquared.txt");
		myfile8 << vxp.array().square()+vyp.array().square();
		myfile8.close();

		ofstream myfile9;
		myfile9.open("4totUEn.txt");

		ofstream myfile10;
		myfile10.open("4totKEn.txt");

		ofstream myfile11;
		myfile11.open("4totEn.txt");

		ofstream myfile12;
		myfile12.open("4Phi1.txt");
		myfile12 << Phi1;
		myfile12.close();

		ofstream myfile13;
		myfile13.open("4rho1.txt");
		myfile13 << rho1;
		myfile13.close();

		ofstream myfile14;
		myfile14.open("4emat.txt");
		myfile14 << emat;
		myfile14.close();
		cout <<"\nCheckpoint6";
		ofstream myfile15;
		myfile15.open("4xp.txt");
		myfile15 << xp;
		myfile15.close();
		cout <<"\nCheckpoint7";
		ofstream myfile16;
		myfile16.open("4yp.txt");
		myfile16 << yp;
		myfile16.close();
		cout <<"\nCheckpoint8";
		ofstream myfile17;
		myfile17.open("4vxp.txt");
		myfile17 << vxp;
		myfile17.close();
		cout <<"\nCheckpoint9";
		ofstream myfile18;
		myfile18.open("4vyp.txt");
		myfile18 << vyp;
		myfile18.close();
		cout <<"\nCheckpoint10";
		ofstream myfile19;
		myfile19.open("4invPoi88.txt");
		myfile19 << invPoi;
		myfile19.close();
		cout <<"\nCheckpoint11";
		ofstream myfile20;
		myfile20.open("4Cap1.txt");
		myfile20 << Cap1;
		myfile20.close();
		cout <<"\nCheckpoint12";
		ofstream myfile21;
		myfile21.open("4Cap2.txt");
		myfile21 << Cap2;
		myfile21.close();
		cout <<"\nCheckpoint13";
		ofstream myfile22;
		myfile22.open("4Cap3.txt");
		myfile22 << Cap3;
		myfile22.close();
		cout <<"\nCheckpoint14";
	// Main computational cycle
	for(int it=0;it<NT;it++)
	{

	// update xp
	for(n=0;n<N;n++)
	{
		xp(n)+=vxp(n)*DT;
		yp(n)+=vyp(n)*DT;
	}



	//Capacity Matrix

	//Particle Delete if in a bordering grid space or out of bounds
	for(n=0;n<N;n++)
	{
		if(xp(n)<L/NG1||xp(n)>L-L/NG1||yp(n)<L/NG1||yp(n)>L-L/NG1)
		{
			xp(n)=0;
			yp(n)=0;
			vxp(n)=0;
			vyp(n)=0;
			emat(n)=0;
			massmat(n)=0;
			c(n)=0;

		}
	}


	n=0;
	while(n<N)
	{
		if(xp(n)==0)
		{
			if((N-n-1)>0)
			{
				while(xp(N)==0)
				{
					N=N-1;
				}
				if(n<N)
				{
					xp(n)=xp(N);
					yp(n)=yp(N);
					vxp(n)=vxp(N);
					vyp(n)=vyp(N);
					emat(n)=emat(N);
					massmat(n)=massmat(N);
					c(n)=c(N);
					xp(N)=0;
					yp(N)=0;
					vxp(N)=0;
					vyp(N)=0;
					emat(N)=0;
					massmat(N)=0;
					c(N)=0;
					N=N-1;
				}
				else
				{
					N++;
				}
			}
			else
			{
				N=N-1;
			}
		}
		n++;
	}
	while(xp(N)!=0)
	{
		N++;
	}

	cout  << "\n N="<<N;




	//Add Fresh Plasma
	double NGnew=4*NG-4;//Border Grid spaces getting replaced
	double Nnew=floor(NGnew*Nin/NG1/NG1+.5);//Put in a number of particles to keep N about the same
	double randmat;//Which side will they appear?
	for(n=0;n<Nnew;n++)
	{
		randmat=rand()%4;
		if(randmat==0)
		{
			xp(N+n)=(L-L/NG)*((double) rand() / (RAND_MAX));
			yp(N+n)=(L/NG)*((double) rand() / (RAND_MAX));
		}
		if(randmat==1)
		{
			xp(N+n)=L-(L/NG)*((double) rand() / (RAND_MAX));
			yp(N+n)=(L-L/NG)*((double) rand() / (RAND_MAX));
		}
		if(randmat==2)
		{
			xp(N+n)=L/NG+(L-L/NG)*((double) rand() / (RAND_MAX));
			yp(N+n)=L-(L/NG)*((double) rand() / (RAND_MAX));
		}
		if(randmat==3)
		{
			xp(N+n)=L/NG+(L-L/NG)*((double) rand() / (RAND_MAX));
			yp(N+n)=L-(L/NG)*((double) rand() / (RAND_MAX));
		}
	}

	for(n=0;n<Nnew;n++)
	{
		randmat=rand()%2;

		massmat(N+n)=(mp-me)*randmat+me;//Set new masses
		vrand=sqrt((-2*log((1-((double) rand() / (RAND_MAX)))))/(massmat(n+N)/me/k/Temp));
		anglerand=2*pi*((double) rand() / (RAND_MAX));//Random direction
		vxp(N+n)=vrand*cos(anglerand);
		vyp(N+n)=vrand*sin(anglerand);
		emat(N+n)=2*e*randmat-e;//Set new charges
		c(N+n)=9*randmat+1;//Set correct colors
	}
	N=N+Nnew;//Add new particles to N




	MatrixXd g(4*Nin1,2);
	for (n=0;n<N;n++)
	{
		g(n,0)=floor(xp(n)/dx+.5);
		g(n,1)=floor(yp(n)/dx+.5);
		g(N+n,0)=floor(xp(n)/dx+.5)+1;
		g(N+n,1)=floor(yp(n)/dx+.5)+1;
	}

	MatrixXd fraz(4*Nin1,2);
	for (n=0;n<N;n++)
	{
		fraz(n,0)=(1-abs(xp(n)/dx-g(n,0)))*(1-abs(yp(n)/dx-g(n,1)));
		fraz(n,1)=(abs(xp(n)/dx-g(n,0)))*(1-abs(yp(n)/dx-g(n,1)));
		fraz(N+n,0)=(1-abs(xp(n)/dx-g(n,0)))*(abs(yp(n)/dx-g(n,1)));
		fraz(N+n,1)=(abs(xp(n)/dx-g(n,0)))*(abs(yp(n)/dx-g(n,1)));
	}

	    //Put the weights into a 3d matrix, where each row in the third
	    //dimension is a new particle. Adding all allements on each row should
	    //equal 1 in 4 adjacent indices

	//Compute the charge density
	VectorXd rho(NG*NG);
	for (n=0;n<NG*NG;n++)
	{
		rho(n)=0;
	}
	for(n=0;n<N;n++)
	{
		rho(index2(g(n,0),g(n,1),NG))+=emat(n)*fraz(n,0)/dx/dx;
		rho(index2(g(n,0),g(N+n,1),NG))+=emat(n)*fraz(n,1)/dx/dx;
		rho(index2(g(N+n,0),g(n,1),NG))+=emat(n)*fraz(N+n,0)/dx/dx;
		rho(index2(g(N+n,0),g(N+n,1),NG))+=emat(n)*fraz(N+n,1)/dx/dx;
	}

	esphere2=0;
	rho=rho+esphere1*rho_sphere1+esphere2*rho_sphere2;




	VectorXd Phi(NG*NG);
	double PhiC1;
	Phi=-invPoi*rho*dx*dx/e0;



	VectorXd Phiin1(lengthin1);
	VectorXd Phitemp1(lengthin1);
	int m=0;
	for(n=0;n<NG*NG;n++)
	{
		if(ematc1(n)==1)
		{
			Phiin1(m)=Phi(n);
			m++;
		}
	}

	Phitemp1=Cap1*Phiin1;
	PhiC1=Phitemp1.sum()/Cap1.sum();

	double PhiC2;
	VectorXd Phiin2(lengthin2);
	VectorXd Phitemp2(lengthin2);
	m=0;
	for(n=0;n<NG*NG;n++)
	{
		if(ematc2(n)==1)
		{
			Phiin2(m)=Phi(n);
			m++;
		}
	}

	Phitemp2=Cap2*Phiin2;
	PhiC2=Phitemp2.sum()/Cap2.sum();

	VectorXd PhiC1C2(lengthin3);
	m=0;
	i=0;
	j=0;
	for (n=0;n<NG*NG;n++)
	{
		if(ematc1(n)==1)
		{
			PhiC1C2(m)=PhiC1-Phi(n);
			i++;
			m++;
		}
		if(ematc2(n)==1)
		{
			PhiC1C2(m)=PhiC2-Phi(n);
			j++;
			m++;
		}
	}

	VectorXd drho3=Cap3.transpose()*PhiC1C2/dx/dx;
	m=0;
	for(n=0;n<NG*NG;n++)
	{
		if(ematc3(n)==1)
		{
			rho(n)+=drho3(m);
			m++;
		}
	}


	Volt1(it)=PhiC1;
	Volt2(it)=PhiC2;
	Charge1(it)=esphere1;
	Charge2(it)=esphere2;



	MatrixXd rho1(NG,NG);
	for(n=0;n<NG;n++)
	{
	    for(m=0;m<NG;m++)
	    {
	    	rho1(n,m)=rho(n*NG+m);
	    }
	}

	if(it>19999)
	{
		rho1allit=(it-20000)*rho1allit/(it-19999)+rho1/(it-19999);
	}

	// compute Potential with Poissons equation
	Phi=-invPoi*rho*dx*dx/e0;//New Voltage

	for(n=0;n<NG;n++)
		{
			Phi(n)=0;
			Phi(n+NG)=0;
			Phi(n+2*NG)=0;
			Phi(NG*NG-n-1)=0;
			Phi(NG*NG-n-NG-1)=0;
			Phi(NG*NG-n-2*NG-1)=0;

			Phi(n*NG)=0;
			Phi(n*NG+1)=0;
			Phi(n*NG+2)=0;
			Phi(n*NG+NG-1)=0;
			Phi(n*NG+NG-2)=0;
			Phi(n*NG+NG-3)=0;
		}


	MatrixXd Phi1(NG,NG);
	for(n=0;n<NG;n++)
	{
	    for(m=0;m<NG;m++)
	    {
	    	Phi1(n,m)=Phi(n*NG+m);
	    }
	}


if(it>19999)
{
	Phi1allit=(it-20000)*Phi1allit/(it-19999)+Phi1/(it-19999);
}
	for(n=0;n<N;n++)
		{
		if((((xp(n)-sxc1)*(xp(n)-sxc1)+(yp(n)-syc1)*(yp(n)-syc1))<sr1*sr1)||(((xp(n)-sxc2)*(xp(n)-sxc2)+(yp(n)-syc2)*(yp(n)-syc2))<sr2*sr2))
			{
				if((((xp(n)-sxc1)*(xp(n)-sxc1)+(yp(n)-syc1)*(yp(n)-syc1))<sr1*sr1))
				{
					esphere1+=emat(n);
					TotalMomentumAbsorbed1X(it)+=massmat(n)*vxp(n);
					TotalMomentumAbsorbed1Y(it)+=massmat(n)*vyp(n);
				}
				if((((xp(n)-sxc2)*(xp(n)-sxc2)+(yp(n)-syc2)*(yp(n)-syc2))<sr2*sr2))
				{
					esphere2+=emat(n);
					TotalMomentumAbsorbed2X(it)+=massmat(n)*vxp(n);
					TotalMomentumAbsorbed2Y(it)+=massmat(n)*vyp(n);
				}
				xp(n)=0;
				yp(n)=0;
				vxp(n)=0;
				vyp(n)=0;
				emat(n)=0;
				massmat(n)=0;
				c(n)=0;
			}
		}

	n=0;
		while(n<N)
		{
			if(xp(n)==0)
			{
				if((N-n-1)>0)
				{
					while(xp(N)==0)
					{
						N=N-1;
					}
					if(n<N)
					{
						xp(n)=xp(N);
						yp(n)=yp(N);
						vxp(n)=vxp(N);
						vyp(n)=vyp(N);
						emat(n)=emat(N);
						massmat(n)=massmat(N);
						c(n)=c(N);
						xp(N)=0;
						yp(N)=0;
						vxp(N)=0;
						vyp(N)=0;
						emat(N)=0;
						massmat(N)=0;
						c(N)=0;
						N=N-1;
					}
					else
					{
						N++;
					}
				}
				else
				{
					N=N-1;
				}
			}
			n++;
		}

		while(xp(N)!=0)
		{
			N++;
		}

	//Electric Field in x direction
	int in1;
	int in2;
	for(n=0;n<NG*NG;n++)
	{
	        in1=n+NG;//Index of right node
	        in2=n-NG;//Index of left node
	        if(in1>NG*NG-1)
	        {
	            Ex(n)=(-Phi(in2))/(2*dx);//Electric field in x, zero out edge
	        }
	        else if (in2<0)
			{
	            Ex(n)=(Phi(in1))/(2*dx);//Electric field in x, zero out edge
			}
	        else
	        {
	        	Ex(n)=(Phi(in1)-Phi(in2))/(2*dx);//Electric field in x
	        }
	}

		for(n=0;n<NG;n++)
		{
		    for(m=0;m<NG;m++)
		    {
		    	Ex1(n,m)=Ex(n*NG+m);
		    }
		}

	//Electric Field in Y Direction
	for (n=0;n<NG*NG;n++)
	{
	        in1=n+1;//Index of above node
	        in2=n-1;//index of lower node
	        if (in1%NG==0)
	        {
	            Ey(n)=(-Phi(in2))/(2*dx);//Electric field in y, zere out edge
	        }
	        else if ((in2+NG)%NG==NG-1)
	        {
	            Ey(n)=(Phi(in1))/(2*dx);//Electric field in y, zero out edge
	        }
	        else
	        {
	        	Ey(n)=(Phi(in1)-Phi(in2))/(2*dx);//Electric field in y
	        }
	}

			for(n=0;n<NG;n++)
			{
			    for(m=0;m<NG;m++)
			    {
			    	Ey1(n,m)=Ey(n*NG+m);
			    }
			}

	    //Update Velocity
	for(n=0;n<N;n++)
	{
		vxp(n)+=-Ex(index2(g(n,0),g(n,1),NG))*fraz(n,0)*emat(n)*DT/massmat(n);
		vxp(n)+=-Ex(index2(g(n,0),g(N+n,1),NG))*fraz(n,1)*emat(n)*DT/massmat(n);
		vxp(n)+=-Ex(index2(g(N+n,0),g(n,1),NG))*fraz(n+N,0)*emat(n)*DT/massmat(n);
		vxp(n)+=-Ex(index2(g(N+n,0),g(N+n,1),NG))*fraz(n+N,1)*emat(n)*DT/massmat(n);


		vyp(n)+=-Ey(index2(g(n,0),g(n,1),NG))*fraz(n,0)*emat(n)*DT/massmat(n);
		vyp(n)+=-Ey(index2(g(n,0),g(N+n,1),NG))*fraz(n,1)*emat(n)*DT/massmat(n);
		vyp(n)+=-Ey(index2(g(N+n,0),g(n,1),NG))*fraz(n+N,0)*emat(n)*DT/massmat(n);
		vyp(n)+=-Ey(index2(g(N+n,0),g(N+n,1),NG))*fraz(n+N,1)*emat(n)*DT/massmat(n);
	}
	cout << "\nit="<<it;
	cout << "\nPhic1="<<Volt1(it)*510999;
	cout << "\nPhic2="<<Volt2(it)*510999;

	double Temp=0;
	for(n=0;n<N;n++)
	{
		Temp+=massmat(n)/me*(vxp(n)*vxp(n)+vyp(n)*vyp(n))/N/2/k;
	}

	totUEn(it)=0;
	for(n=0;n<NG;n++)
	{
		totUEn(it)+=.5*dx*dx*e0*(Ex(n)*Ex(n)+Ey(n)*Ey(n));
	}

	totKEn(it)=0;

	for(n=0;n<N;n++)
	{
		totKEn(it)+=.5*massmat(n)*(vxp(n)*vxp(n)+vyp(n)*vyp(n));
	}

	totEn(it)=totUEn(it)+totKEn(it);

	cout << "\nTotal Energy: " << totEn(it);


	cout << "\nTemp=" << Temp/11700;

	cout << "\nCharge1=" << esphere1/e;
	cout << "\nCharge2=" << esphere2/e;

	Temp1(it)=Temp;

	LD =sqrt(k*Temp/ndensity);
			if(.5*LD<dx)//Compare to delta x, if too small, display error
			{
			    cout<<"\nError, Mesh not fine enough or not hot enough\nLD="<<LD<<"\ndx="<<dx;
			}
			double PP=sqrt(pi*me*L*L/(N*e*e));//Compute plasma period
			if(.25*PP<DT)//If plasma period is small display error
			{
				cout<<"/nError, large time step\nPP="<<PP<<"\nDT="<<DT;
			}

			//CFL Condition
			double C=DT/dx;
			if(C>1)//If time step is larger than dx
			{
			    cout <<"\nMay be unstable, C>1.\nC="<<C;
			}



	if(it%10000==9999)
	{

		for(n=0;n<10000;n++)
		{
			myfile1 <<Volt1(it-9999+n)<<"\n";

			myfile2 << Volt2(it-9999+n)<<"\n";

			myfile3 << Charge1(it-9999+n)<<"\n";

			myfile4 << Charge2(it-9999+n)<<"\n";

			myfile7 << Temp1(it-9999+n)<<"\n";

			myfile9 << totUEn(it-9999+n)<<"\n";

			myfile10 << totKEn(it-9999+n)<<"\n";

			myfile11 << totEn(it-9999+n)<<"\n";
		}

		myfile5.open("4rho1allit.txt");
		myfile5 << rho1allit;
		myfile5.close();

		myfile6.open("4Phi1allit.txt");
		myfile6 << Phi1allit;
		myfile6.close();

		myfile8.open("4velocitysquared.txt");
		myfile8 << vxp.array().square()+vyp.array().square();
		myfile8.close();

		myfile12.open("4Phi1.txt");
		myfile12 << Phi1;
		myfile12.close();

		myfile13.open("4rho1.txt");
		myfile13 << rho1;
		myfile13.close();

		myfile14.open("4emat.txt");
		myfile14 << emat;
		myfile14.close();

		myfile15.open("4xp.txt");
		myfile15 << xp;
		myfile15.close();

		myfile16.open("4yp.txt");
		myfile16 << yp;
		myfile16.close();

		myfile17.open("4vxp.txt");
		myfile17 << vxp;
		myfile17.close();

		myfile18.open("4vyp.txt");
		myfile18 << vyp;
		myfile18.close();
	}


	}


	myfile1.close();
	myfile2.close();
	myfile3.close();
	myfile4.close();
	myfile7.close();
	myfile9.close();
	myfile10.close();
	myfile11.close();








return 0;
}

//
//	 ForceXSphere1(it)=dx^2*(drho1tot+esphere1*rho_sphere1)'*Ex;
//	 ForceYSphere1(it)=dx^2*(drho1tot+esphere1*rho_sphere1)'*Ey;
//	 ForceXSphere2(it)=dx^2*(drho2tot+esphere2*rho_sphere2)'*Ex;
//	 ForceYSphere2(it)=dx^2*(drho2tot+esphere2*rho_sphere2)'*Ey;
//
//	 Force1=ForceYSphere1(it)
//	 Force2=ForceYSphere2(it)
//
//	 vsx1=vsx1;//-ForceXSphere1*DT/ms1;//Change sphere 1 x velocity
//	 vsy1=vsy1;//-ForceYSphere1*DT/ms1;//Change sphere 1 y velocity
//	 vsx2=vsx2;//-ForceXSphere2*DT/ms2;//Change sphere 2 x velocity
//	 vsy2=vsy2;//-ForceYSphere2*DT/ms2;//Change sphere 2 y velocity
//
//
//	 Charge1(it)=esphere1;
//	 Charge2(it)=esphere2;
//	 Temp1(it)=Temp;
//	 TotalCharge1(it)=sum(emat);
//
//	//Energy test (Check if Total energy is constant)
//
//	totEn(it)=.5*sum(massmat*(vxp.^2+vyp.^2))+.5*dx^2*e0*sum(Ex.^2+Ey.^2);
//
//
//	totEnWLost(it)=totEn(it)+lostenergy;
//
//	totUEn(it)=.5*dx^2*e0*sum(Ex.^2+Ey.^2);
//
//
//	Kinetic Energy
//	totKEn(it)=.5*sum(massmat*(vxp.^2+vyp.^2));
//
//	//Display Iteration
//	it
//
//
//	Temp=(sum(massmat/me*(vxp.^2+vyp.^2))/N)/2/k;//Calculate the temperature
//	Temperature=Temp/11700
//	LD =(k*Temp/ndensity)^.5;
//	if .5*LD<dx//Compare to delta x, if too small, display error
//	    "Error, Mesh not fine enough or not hot enough"
//	    LD
//	    dx
//	end
//	PP=(pi*me*L^2/(N*e^2))^.5;//Compute plasma period
//	if .25*PP<DT//If plasma period is small display error
//	    "Error, large time step"
//	    PP
//	    DT
//	end
//
//	//CFL Condition
//	C=DT/dx;
//	if C>1//If time step is larger than dx
//	    "May be unstable, C>1"
//	    C
//	end
//
//	//PhiC2*510999
//	end
//
//
//}

