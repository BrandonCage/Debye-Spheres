/*
 * TestMat3.cpp
 *
 *  Created on: Nov 5, 2020
 *      Author: shadowcage72
 */






#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;


int index2(double y,double z,int NG)
{
	if(z>NG-1)
	{
		z-=NG;
	}
	if(z<0)
	{
		z+=NG;
	}
	if(y>NG-1)
	{
		y-=NG;
	}
	if(y<0)
	{
		y+=NG;
	}
	return y*NG+z;
}








int main()
{


	double starttime = clock();
	double e0=1;//Epsilon Naught
	double k=.000000000168637205;//Boltzman
	double meterconv=0.00000000000003511282;//Convert from our units to meters
	double L=9/meterconv;//Length
	double DT=L/100;//L/60;//DeltaT
	int NT=500000;//Iterations
	int NG=52;//Number of Grid Points
	double NG1=NG;
	double ndensity=5000000000*meterconv*meterconv*meterconv;
	double totalparticles=2*ndensity*L*L;//Total number of electrons in simulation


	//At least 25 per cell is ideal

	int N=100000;//Number of Particles
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
		while(isinf(vrand))
		{
			vrand=sqrt((-2*log((1-((double) rand() / (RAND_MAX)))))/(massmat(n)/me/k/Temp));
		}
		anglerand=2*pi*((double) rand() / (RAND_MAX));
		vxp(n)=vrand*cos(anglerand);
		vxp(n+N)=0;
		vyp(n)=vrand*sin(anglerand);//Make all into ions
		vyp(n+N)=0;
	}
MatrixXd invPoi(NG*NG,NG*NG);
	for(int n=0;n<NG*NG;n++)
	{
		for(int m=0;m<NG*NG;m++)
		{
			invPoi(n,m)=0;
		}
	}
	for(n=0;n<NG*NG;n++)
	{
		invPoi(n,n)=-4;
		if(n>0)
		{
			invPoi(n,n-1)=1;
		}
		if(n<NG*NG-1)
		{
			invPoi(n,n+1)=1;
		}
		if(n>NG-1)
		{
			invPoi(n,n-NG)=1;
		}
		if(n<NG*NG-NG)
		{
			invPoi(n,n+NG)=1;
		}
	}
	for(n=1;n<NG;n++)
	{
		invPoi(NG*n-1,NG*n)=0;
		invPoi(NG*n,NG*n-1)=0;
	}

invPoi=invPoi.inverse();
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
	double sr1=.06*L;//Sphere 1 radius


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
		xpcap(n)=dx*(j+.5);
		ypcap(n)=dx*(i+.5);
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



	//Circular Matrix representing 1 particle around Sphere 1
		int NCirc=64;//Number of particles for circle
		double NCirc1=NCirc;
		VectorXd xpc(NCirc);
		VectorXd ypc(NCirc);
		for (n=0;n<NCirc;n++)
		{
			xpc(n)=sxc1+sr1*sin(2*pi/NCirc*n);
			ypc(n)=syc1+sr1*cos(2*pi/NCirc*n);
		}
		double charge=1/NCirc1;

		MatrixXd gc(2*NCirc,2);
		for (n=0;n<NCirc;n++)
		{
			gc(n,0)=floor(xpc(n)/dx-.5);
			gc(n,1)=floor(ypc(n)/dx-.5);
			gc(NCirc+n,0)=(floor(xpc(n)/dx-.5)+1);
			gc(NCirc+n,1)=(floor(ypc(n)/dx-.5)+1);
		}


		MatrixXd frazc(2*NCirc,2);
		for (n=0;n<NCirc;n++)
		{
			frazc(n,0)=(1-(xpc(n)/dx-.5-gc(n,0)))*(1-(ypc(n)/dx-.5-gc(n,1)));
			frazc(n,1)=(xpc(n)/dx-.5-gc(n,0))*(1-(ypc(n)/dx-.5-gc(n,1)));
			frazc(NCirc+n,0)=(1-(xpc(n)/dx-.5-gc(n,0)))*(ypc(n)/dx-.5-gc(n,1));
			frazc(NCirc+n,1)=(xpc(n)/dx-.5-gc(n,0))*(ypc(n)/dx-.5-gc(n,1));
		}


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
				xpc(n)=sxc2+sr2*sin(2*pi/NCirc*n);
				ypc(n)=syc2+sr2*cos(2*pi/NCirc*n);
			}

			for (n=0;n<NCirc;n++)
			{
				gc(n,0)=floor(xpc(n)/dx-.5);
				gc(n,1)=floor(ypc(n)/dx-.5);
				gc(NCirc+n,0)=(floor(xpc(n)/dx-.5)+1);
				gc(NCirc+n,1)=(floor(ypc(n)/dx-.5)+1);
			}


			for (n=0;n<NCirc;n++)
			{
				frazc(n,0)=(1-(xpc(n)/dx-.5-gc(n,0)))*(1-(ypc(n)/dx-.5-gc(n,1)));
				frazc(n,1)=(xpc(n)/dx-.5-gc(n,0))*(1-(ypc(n)/dx-.5-gc(n,1)));
				frazc(NCirc+n,0)=(1-(xpc(n)/dx-.5-gc(n,0)))*(ypc(n)/dx-.5-gc(n,1));
				frazc(NCirc+n,1)=(xpc(n)/dx-.5-gc(n,0))*(ypc(n)/dx-.5-gc(n,1));
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



cout <<"lengthin1="<<lengthin1;
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
		myfile1.open("6Volt1.txt");

		ofstream myfile2;
		myfile2.open("6Volt2.txt");

		ofstream myfile3;
		myfile3.open("6Charge1.txt");

		ofstream myfile4;
		myfile4.open("6Charge2.txt");


		ofstream myfile5;
		myfile5.open("6rho1allit.txt");
		myfile5 << rho1allit;
		myfile5.close();

		ofstream myfile6;
		myfile6.open("6Phi1allit.txt");
		myfile6 << Phi1allit;
		myfile6.close();

		ofstream myfile7;
		myfile7.open("6Temp1.txt");

		ofstream myfile8;
		myfile8.open("6velocitysquared.txt");
		myfile8 << vxp.array().square()+vyp.array().square();
		myfile8.close();

		ofstream myfile9;
		myfile9.open("6totUEn.txt");

		ofstream myfile10;
		myfile10.open("6totKEn.txt");

		ofstream myfile11;
		myfile11.open("6totEn.txt");

		ofstream myfile12;
		myfile12.open("6Phi1.txt");
		myfile12 << Phi1;
		myfile12.close();

		ofstream myfile13;
		myfile13.open("6rho1.txt");
		myfile13 << rho1;
		myfile13.close();

		ofstream myfile14;
		myfile14.open("6emat.txt");
		myfile14 << emat;
		myfile14.close();
		cout <<"\nCheckpoint6";
		ofstream myfile15;
		myfile15.open("6xp.txt");
		myfile15 << xp;
		myfile15.close();
		cout <<"\nCheckpoint7";
		ofstream myfile16;
		myfile16.open("6yp.txt");
		myfile16 << yp;
		myfile16.close();
		cout <<"\nCheckpoint8";
		ofstream myfile17;
		myfile17.open("6vxp.txt");
		myfile17 << vxp;
		myfile17.close();
		cout <<"\nCheckpoint9";
		ofstream myfile18;
		myfile18.open("6vyp.txt");
		myfile18 << vyp;
		myfile18.close();
		cout <<"\nCheckpoint10";
		ofstream myfile19;
		myfile19.open("6invPoi88.txt");
		//myfile19 << invPoi;
		myfile19.close();
		cout <<"\nCheckpoint11";
		ofstream myfile20;
		myfile20.open("6Cap1.txt");
		myfile20 << Cap1;
		myfile20.close();
		cout <<"\nCheckpoint12";
		ofstream myfile21;
		myfile21.open("6Cap2.txt");
		myfile21 << Cap2;
		myfile21.close();
		cout <<"\nCheckpoint13";
		ofstream myfile22;
		myfile22.open("6Cap3.txt");
		myfile22 << Cap3;
		myfile22.close();
		cout <<"\nCheckpoint14";
		ofstream myfile23;
		myfile23.open("6TotalMomentumAbsorbed1X.txt");
		ofstream myfile24;
		myfile24.open("6TotalMomentumAbsorbed1Y.txt");
		ofstream myfile25;
		myfile25.open("TotalMomentumAbsorbed2X.txt");
		ofstream myfile26;
		myfile26.open("6TotalMomentumAbsorbed2Y.txt");
	// Main computational cycle
	for(int it=0;it<NT;it++)
	{
	// update xp
#pragma omp parallel sections num_threads(4)
		{
#pragma omp section
			{
				for(n=0;n<N/4;n++)
					{
						xp(n)+=vxp(n)*DT/2;
						yp(n)+=vyp(n)*DT/2;
					}
			}
#pragma omp section
			{
				for(n=N/4;n<N/2;n++)
					{
						xp(n)+=vxp(n)*DT/2;
						yp(n)+=vyp(n)*DT/2;
					}
			}
#pragma omp section
			{
				for(n=N/2;n<N/4*3;n++)
					{
						xp(n)+=vxp(n)*DT/2;
						yp(n)+=vyp(n)*DT/2;
					}
			}
#pragma omp section
			{
				for(n=N/4*3;n<N;n++)
					{
						xp(n)+=vxp(n)*DT/2;
						yp(n)+=vyp(n)*DT/2;
					}
			}
		}


	//Capacity Matrix

	//Particle Delete if in a bordering grid space or out of bounds
#pragma omp parallel sections num_threads(4)
		{
#pragma omp section
			{
				for(n=0;n<N/4;n++)
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
			}
#pragma omp section
			{
				for(n=N/4;n<N/2;n++)
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
			}
#pragma omp section
			{
				for(n=N/2;n<N/4*3;n++)
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
			}
#pragma omp section
			{
				for(n=N/4*3;n<N;n++)
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
#pragma omp parallel sections num_threads(4)
		{
#pragma omp section
			{
				for(n=0;n<Nnew/4;n++)
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
							xp(N+n)=(L/NG)*((double) rand() / (RAND_MAX));
							yp(N+n)=L/NG+(L-L/NG)*((double) rand() / (RAND_MAX));
						}
						randmat=rand()%2;
						massmat(N+n)=(mp-me)*randmat+me;//Set new masses
						vrand=sqrt((-2*log((1-((double) rand() / (RAND_MAX)))))/(massmat(n+N)/me/k/Temp));
						while(isinf(vrand))
						{
							vrand=sqrt((-2*log((1-((double) rand() / (RAND_MAX)))))/(massmat(n)/me/k/Temp));
						}
						anglerand=2*pi*((double) rand() / (RAND_MAX));//Random direction
						vxp(N+n)=vrand*cos(anglerand);
						vyp(N+n)=vrand*sin(anglerand);
						emat(N+n)=2*e*randmat-e;//Set new charges
						c(N+n)=9*randmat+1;//Set correct colors
					}
			}
#pragma omp section
			{
				for(n=Nnew/4;n<Nnew/2;n++)
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
							xp(N+n)=(L/NG)*((double) rand() / (RAND_MAX));
							yp(N+n)=L/NG+(L-L/NG)*((double) rand() / (RAND_MAX));
						}
						randmat=rand()%2;
						massmat(N+n)=(mp-me)*randmat+me;//Set new masses
						vrand=sqrt((-2*log((1-((double) rand() / (RAND_MAX)))))/(massmat(n+N)/me/k/Temp));
						while(isinf(vrand))
								{
									vrand=sqrt((-2*log((1-((double) rand() / (RAND_MAX)))))/(massmat(n)/me/k/Temp));
								}
						anglerand=2*pi*((double) rand() / (RAND_MAX));//Random direction
						vxp(N+n)=vrand*cos(anglerand);
						vyp(N+n)=vrand*sin(anglerand);
						emat(N+n)=2*e*randmat-e;//Set new charges
						c(N+n)=9*randmat+1;//Set correct colors
					}
			}
#pragma omp section
			{
				for(n=Nnew/2;n<Nnew/4*3;n++)
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
							xp(N+n)=(L/NG)*((double) rand() / (RAND_MAX));
							yp(N+n)=L/NG+(L-L/NG)*((double) rand() / (RAND_MAX));
						}
						randmat=rand()%2;
						massmat(N+n)=(mp-me)*randmat+me;//Set new masses
						vrand=sqrt((-2*log((1-((double) rand() / (RAND_MAX)))))/(massmat(n+N)/me/k/Temp));
						while(isinf(vrand))
								{
									vrand=sqrt((-2*log((1-((double) rand() / (RAND_MAX)))))/(massmat(n)/me/k/Temp));
								}
						anglerand=2*pi*((double) rand() / (RAND_MAX));//Random direction
						vxp(N+n)=vrand*cos(anglerand);
						vyp(N+n)=vrand*sin(anglerand);
						emat(N+n)=2*e*randmat-e;//Set new charges
						c(N+n)=9*randmat+1;//Set correct colors
					}
			}
#pragma omp section
			{
				for(n=Nnew/4*3;n<Nnew;n++)
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
							xp(N+n)=(L/NG)*((double) rand() / (RAND_MAX));
							yp(N+n)=L/NG+(L-L/NG)*((double) rand() / (RAND_MAX));
						}
						randmat=rand()%2;
						massmat(N+n)=(mp-me)*randmat+me;//Set new masses
						vrand=sqrt((-2*log((1-((double) rand() / (RAND_MAX)))))/(massmat(n+N)/me/k/Temp));
						while(isinf(vrand))
								{
									vrand=sqrt((-2*log((1-((double) rand() / (RAND_MAX)))))/(massmat(n)/me/k/Temp));
								}
						anglerand=2*pi*((double) rand() / (RAND_MAX));//Random direction
						vxp(N+n)=vrand*cos(anglerand);
						vyp(N+n)=vrand*sin(anglerand);
						emat(N+n)=2*e*randmat-e;//Set new charges
						c(N+n)=9*randmat+1;//Set correct colors
					}
			}
		}
	N=N+Nnew;//Add new particles to N




	MatrixXd g(4*Nin1,2);
#pragma omp parallel sections num_threads(4)
				{
		#pragma omp section
					{
						for(n=0;n<N/4;n++)
						{
							g(n,0)=floor(xp(n)/dx-.5);
							g(n,1)=floor(yp(n)/dx-.5);
							g(N+n,0)=(floor(xp(n)/dx-.5)+1);
							g(N+n,1)=(floor(yp(n)/dx-.5)+1);
						}
					}
		#pragma omp section
					{
						for(n=N/4;n<N/2;n++)
						{
							g(n,0)=floor(xp(n)/dx-.5);
							g(n,1)=floor(yp(n)/dx-.5);
							g(N+n,0)=(floor(xp(n)/dx-.5)+1);
							g(N+n,1)=(floor(yp(n)/dx-.5)+1);
						}
					}
		#pragma omp section
					{
						for(n=N/2;n<N/4*3;n++)
						{
							g(n,0)=floor(xp(n)/dx-.5);
							g(n,1)=floor(yp(n)/dx-.5);
							g(N+n,0)=(floor(xp(n)/dx-.5)+1);
							g(N+n,1)=(floor(yp(n)/dx-.5)+1);
						}
					}
		#pragma omp section
					{
						for(n=N/4*3;n<N;n++)
						{
							g(n,0)=floor(xp(n)/dx-.5);
							g(n,1)=floor(yp(n)/dx-.5);
							g(N+n,0)=(floor(xp(n)/dx-.5)+1);
							g(N+n,1)=(floor(yp(n)/dx-.5)+1);
						}
					}
				}
	MatrixXd fraz(4*Nin1,2);
#pragma omp parallel sections num_threads(4)
		{
#pragma omp section
			{
				for(n=0;n<N/4;n++)
				{
					fraz(n,0)=(1-(xp(n)/dx-.5-g(n,0)))*(1-(yp(n)/dx-.5-g(n,1)));
					fraz(n,1)=(xp(n)/dx-.5-g(n,0))*(1-(yp(n)/dx-.5-g(n,1)));
					fraz(N+n,0)=(1-(xp(n)/dx-.5-g(n,0)))*(yp(n)/dx-.5-g(n,1));
					fraz(N+n,1)=(xp(n)/dx-.5-g(n,0))*(yp(n)/dx-.5-g(n,1));
				}
			}
#pragma omp section
			{
				for(n=N/4;n<N/2;n++)
				{
					fraz(n,0)=(1-(xp(n)/dx-.5-g(n,0)))*(1-(yp(n)/dx-.5-g(n,1)));
					fraz(n,1)=(xp(n)/dx-.5-g(n,0))*(1-(yp(n)/dx-.5-g(n,1)));
					fraz(N+n,0)=(1-(xp(n)/dx-.5-g(n,0)))*(yp(n)/dx-.5-g(n,1));
					fraz(N+n,1)=(xp(n)/dx-.5-g(n,0))*(yp(n)/dx-.5-g(n,1));
				}
			}
#pragma omp section
			{
				for(n=N/2;n<N/4*3;n++)
				{
					fraz(n,0)=(1-(xp(n)/dx-.5-g(n,0)))*(1-(yp(n)/dx-.5-g(n,1)));
					fraz(n,1)=(xp(n)/dx-.5-g(n,0))*(1-(yp(n)/dx-.5-g(n,1)));
					fraz(N+n,0)=(1-(xp(n)/dx-.5-g(n,0)))*(yp(n)/dx-.5-g(n,1));
					fraz(N+n,1)=(xp(n)/dx-.5-g(n,0))*(yp(n)/dx-.5-g(n,1));
				}
			}
#pragma omp section
			{
				for(n=N/4*3;n<N;n++)
				{
					fraz(n,0)=(1-(xp(n)/dx-.5-g(n,0)))*(1-(yp(n)/dx-.5-g(n,1)));
					fraz(n,1)=(xp(n)/dx-.5-g(n,0))*(1-(yp(n)/dx-.5-g(n,1)));
					fraz(N+n,0)=(1-(xp(n)/dx-.5-g(n,0)))*(yp(n)/dx-.5-g(n,1));
					fraz(N+n,1)=(xp(n)/dx-.5-g(n,0))*(yp(n)/dx-.5-g(n,1));
				}
			}
		}
	    //Put the weights into a 3d matrix, where each row in the third
	    //dimension is a new particle. Adding all allements on each row should
	    //equal 1 in 4 adjacent indices

	//Compute the charge density
	VectorXd rho(NG*NG);
#pragma omp parallel sections num_threads(4)
		{
#pragma omp section
			{
				for(n=0;n<NG*NG/4;n++)
				{
						rho(n)=0;
					}
			}
#pragma omp section
			{
				for(n=NG*NG/4;n<NG*NG/2;n++)
				{
						rho(n)=0;
					}
			}
#pragma omp section
			{
				for(n=NG*NG/2;n<NG*NG/4*3;n++)
				{
						rho(n)=0;
					}
			}
#pragma omp section
			{
				for(n=NG*NG/4*3;n<NG*NG;n++)
				{
						rho(n)=0;
					}
			}
		}

#pragma omp parallel sections num_threads(4)
		{
#pragma omp section
			{
				for(n=0;n<N/4;n++)
				{
						rho(index2(g(n,0),g(n,1),NG))+=emat(n)*fraz(n,0)/dx/dx;
						rho(index2(g(n,0),g(N+n,1),NG))+=emat(n)*fraz(n,1)/dx/dx;
						rho(index2(g(N+n,0),g(n,1),NG))+=emat(n)*fraz(N+n,0)/dx/dx;
						rho(index2(g(N+n,0),g(N+n,1),NG))+=emat(n)*fraz(N+n,1)/dx/dx;
					}
			}
#pragma omp section
			{
				for(n=N/4;n<N/2;n++)
				{
						rho(index2(g(n,0),g(n,1),NG))+=emat(n)*fraz(n,0)/dx/dx;
						rho(index2(g(n,0),g(N+n,1),NG))+=emat(n)*fraz(n,1)/dx/dx;
						rho(index2(g(N+n,0),g(n,1),NG))+=emat(n)*fraz(N+n,0)/dx/dx;
						rho(index2(g(N+n,0),g(N+n,1),NG))+=emat(n)*fraz(N+n,1)/dx/dx;
					}
			}
#pragma omp section
			{
				for(n=N/2;n<N/4*3;n++)
				{
						rho(index2(g(n,0),g(n,1),NG))+=emat(n)*fraz(n,0)/dx/dx;
						rho(index2(g(n,0),g(N+n,1),NG))+=emat(n)*fraz(n,1)/dx/dx;
						rho(index2(g(N+n,0),g(n,1),NG))+=emat(n)*fraz(N+n,0)/dx/dx;
						rho(index2(g(N+n,0),g(N+n,1),NG))+=emat(n)*fraz(N+n,1)/dx/dx;
					}
			}
#pragma omp section
			{
				for(n=N/4*3;n<N;n++)
				{
						rho(index2(g(n,0),g(n,1),NG))+=emat(n)*fraz(n,0)/dx/dx;
						rho(index2(g(n,0),g(N+n,1),NG))+=emat(n)*fraz(n,1)/dx/dx;
						rho(index2(g(N+n,0),g(n,1),NG))+=emat(n)*fraz(N+n,0)/dx/dx;
						rho(index2(g(N+n,0),g(N+n,1),NG))+=emat(n)*fraz(N+n,1)/dx/dx;
					}
			}
		}

	esphere2=0;
	rho=rho+esphere1*rho_sphere1+esphere2*rho_sphere2;


	VectorXd Phiseg1(NG*NG/4);
	VectorXd Phiseg2(NG*NG/4);
	VectorXd Phiseg3(NG*NG/4);
	VectorXd Phiseg4(NG*NG/4);

	double PhiC1;
#pragma omp parallel sections num_threads(4)
		{
#pragma omp section
			{
				Phiseg1=-invPoi.block(0,0,NG*NG/4,NG*NG)*rho*dx*dx/e0;
			}
#pragma omp section
			{
				Phiseg2=-invPoi.block(NG*NG/4,0,NG*NG/4,NG*NG)*rho*dx*dx/e0;
			}
#pragma omp section
			{
				Phiseg3=-invPoi.block(NG*NG/2,0,NG*NG/4,NG*NG)*rho*dx*dx/e0;
			}
#pragma omp section
			{
				Phiseg4=-invPoi.block(3*NG*NG/4,0,NG*NG/4,NG*NG)*rho*dx*dx/e0;
			}
		}



	VectorXd Phiin1(lengthin1);
	VectorXd Phitemp1(lengthin1);
	int m=0;
	for(n=0;n<NG*NG;n++)
	{
		if(ematc1(n)==1)
		{
			if(n<NG*NG/4)
			{
				Phiin1(m)=Phiseg1(n);
			}
			else if(n<NG*NG/2)
			{
				Phiin1(m)=Phiseg2(n-NG*NG/4);
			}
			else if(n<3*NG*NG/4)
			{
				Phiin1(m)=Phiseg3(n-NG*NG/2);
			}
			else
			{
				Phiin1(m)=Phiseg4(n-3*NG*NG/4);
			}
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
			if(n<NG*NG/4)
			{
				Phiin2(m)=Phiseg1(n);
			}
			else if(n<NG*NG/2)
			{
				Phiin2(m)=Phiseg2(n-NG*NG/4);
			}
			else if(n<3*NG*NG/4)
			{
				Phiin2(m)=Phiseg3(n-NG*NG/2);
			}
			else
			{
				Phiin2(m)=Phiseg4(n-3*NG*NG/4);
			}
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
			if(n<NG*NG/4)
			{
				PhiC1C2(m)=PhiC1-Phiseg1(n);
			}
			else if(n<NG*NG/2)
			{
				PhiC1C2(m)=PhiC1-Phiseg2(n-NG*NG/4);
			}
			else if(n<3*NG*NG/4)
			{
				PhiC1C2(m)=PhiC1-Phiseg3(n-NG*NG/2);
			}
			else
			{
				PhiC1C2(m)=PhiC1-Phiseg4(n-3*NG*NG/4);
			}
			i++;
			m++;
		}
		if(ematc2(n)==1)
		{
			if(n<NG*NG/4)
			{
				PhiC1C2(m)=PhiC2-Phiseg1(n);
			}
			else if(n<NG*NG/2)
			{
				PhiC1C2(m)=PhiC2-Phiseg2(n-NG*NG/4);
			}
			else if(n<3*NG*NG/4)
			{
				PhiC1C2(m)=PhiC2-Phiseg3(n-NG*NG/2);
			}
			else
			{
				PhiC1C2(m)=PhiC2-Phiseg4(n-3*NG*NG/4);
			}
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
#pragma omp parallel sections num_threads(4)
		{
#pragma omp section
			{
				Phiseg1=-invPoi.block(0,0,NG*NG/4,NG*NG)*rho*dx*dx/e0;
			}
#pragma omp section
			{
				Phiseg2=-invPoi.block(NG*NG/4,0,NG*NG/4,NG*NG)*rho*dx*dx/e0;
			}
#pragma omp section
			{
				Phiseg3=-invPoi.block(NG*NG/2,0,NG*NG/4,NG*NG)*rho*dx*dx/e0;
			}
#pragma omp section
			{
				Phiseg4=-invPoi.block(3*NG*NG/4,0,NG*NG/4,NG*NG)*rho*dx*dx/e0;
			}
		}



	MatrixXd Phi1(NG,NG);
	for(n=0;n<NG;n++)
	{
	    for(m=0;m<NG;m++)
	    {
			if(n*NG+m<NG*NG/4)
			{
		    	Phi1(n,m)=Phiseg1(n*NG+m);
			}
			else if(n*NG+m<NG*NG/2)
			{
		    	Phi1(n,m)=Phiseg2(n*NG+m-NG*NG/4);
			}
			else if(n*NG+m<3*NG*NG/4)
			{
		    	Phi1(n,m)=Phiseg3(n*NG+m-NG*NG/2);
			}
			else
			{
		    	Phi1(n,m)=Phiseg4(n*NG+m-3*NG*NG/4);
			}
	    }
	}


if(it>19999)
{
	Phi1allit=(it-20000)*Phi1allit/(it-19999)+Phi1/(it-19999);
}

#pragma omp parallel sections num_threads(4)
		{
#pragma omp section
			{
				for(n=0;n<N/4;n++)
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
			}
#pragma omp section
			{
				for(n=N/4;n<N/2;n++)
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
			}
#pragma omp section
			{
				for(n=N/2;n<N/4*3;n++)
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
			}
#pragma omp section
			{
				for(n=N/4*3;n<N;n++)
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
	double Phiind1;
	double Phiind2;
#pragma omp parallel sections num_threads(4)
				{
		#pragma omp section
					{
						for(n=0;n<NG*NG/4;n++)
						{
							        in1=n+NG;//Index of right node
							        in2=n-NG;//Index of left node
							        Phiind1=0;
							        Phiind2=0;
							        if(in1<NG*NG)
							        {
							        	if(in1<NG*NG/4)
							        				{
							        					Phiind1=Phiseg1(in1);
							        				}
							        				else if(in1<NG*NG/2)
							        				{
							        					Phiind1=Phiseg2(in1-NG*NG/4);
							        				}
							        				else if(in1<3*NG*NG/4)
							        				{
							        					Phiind1=Phiseg3(in1-NG*NG/2);
							        				}
							        				else
							        				{
							        					Phiind1=Phiseg4(in1-3*NG*NG/4);
							        				}
							        }
							        if(in2>=0)
							        {
							        	if(in2<NG*NG/4)
							        				{
							        					Phiind2=Phiseg1(in2);
							        				}
							        				else if(in2<NG*NG/2)
							        				{
							        					Phiind2=Phiseg2(in2-NG*NG/4);
							        				}
							        				else if(in2<3*NG*NG/4)
							        				{
							        					Phiind2=Phiseg3(in2-NG*NG/2);
							        				}
							        				else
							        				{
							        					Phiind2=Phiseg4(in2-3*NG*NG/4);
							        				}
							        }

							        	Ex(n)=(Phiind1-Phiind2)/(2*dx);//Electric field in x
							}
					}
		#pragma omp section
					{
						for(n=NG*NG/4;n<NG*NG/2;n++)
						{
							        in1=n+NG;//Index of right node
							        in2=n-NG;//Index of left node
							        Phiind1=0;
							        Phiind2=0;
							        if(in1<NG*NG)
							        {
							        	if(in1<NG*NG/4)
							        				{
							        					Phiind1=Phiseg1(in1);
							        				}
							        				else if(in1<NG*NG/2)
							        				{
							        					Phiind1=Phiseg2(in1-NG*NG/4);
							        				}
							        				else if(in1<3*NG*NG/4)
							        				{
							        					Phiind1=Phiseg3(in1-NG*NG/2);
							        				}
							        				else
							        				{
							        					Phiind1=Phiseg4(in1-3*NG*NG/4);
							        				}
							        }
							        if(in2>=0)
							        {
							        	if(in2<NG*NG/4)
							        				{
							        					Phiind2=Phiseg1(in2);
							        				}
							        				else if(in2<NG*NG/2)
							        				{
							        					Phiind2=Phiseg2(in2-NG*NG/4);
							        				}
							        				else if(in2<3*NG*NG/4)
							        				{
							        					Phiind2=Phiseg3(in2-NG*NG/2);
							        				}
							        				else
							        				{
							        					Phiind2=Phiseg4(in2-3*NG*NG/4);
							        				}
							        }

							        	Ex(n)=(Phiind1-Phiind2)/(2*dx);//Electric field in x
							}
					}
		#pragma omp section
					{
						for(n=NG*NG/2;n<3*NG*NG/4;n++)
						{
							        in1=n+NG;//Index of right node
							        in2=n-NG;//Index of left node
							        Phiind1=0;
							        Phiind2=0;
							        if(in1<NG*NG)
							        {
							        	if(in1<NG*NG/4)
							        				{
							        					Phiind1=Phiseg1(in1);
							        				}
							        				else if(in1<NG*NG/2)
							        				{
							        					Phiind1=Phiseg2(in1-NG*NG/4);
							        				}
							        				else if(in1<3*NG*NG/4)
							        				{
							        					Phiind1=Phiseg3(in1-NG*NG/2);
							        				}
							        				else
							        				{
							        					Phiind1=Phiseg4(in1-3*NG*NG/4);
							        				}
							        }
							        if(in2>=0)
							        {
							        	if(in2<NG*NG/4)
							        				{
							        					Phiind2=Phiseg1(in2);
							        				}
							        				else if(in2<NG*NG/2)
							        				{
							        					Phiind2=Phiseg2(in2-NG*NG/4);
							        				}
							        				else if(in2<3*NG*NG/4)
							        				{
							        					Phiind2=Phiseg3(in2-NG*NG/2);
							        				}
							        				else
							        				{
							        					Phiind2=Phiseg4(in2-3*NG*NG/4);
							        				}
							        }

							        	Ex(n)=(Phiind1-Phiind2)/(2*dx);//Electric field in x
							}
					}
		#pragma omp section
					{
						for(n=3*NG*NG/4;n<NG*NG;n++)
						{
							        in1=n+NG;//Index of right node
							        in2=n-NG;//Index of left node
							        Phiind1=0;
							        Phiind2=0;
							        if(in1<NG*NG)
							        {
							        	if(in1<NG*NG/4)
							        				{
							        					Phiind1=Phiseg1(in1);
							        				}
							        				else if(in1<NG*NG/2)
							        				{
							        					Phiind1=Phiseg2(in1-NG*NG/4);
							        				}
							        				else if(in1<3*NG*NG/4)
							        				{
							        					Phiind1=Phiseg3(in1-NG*NG/2);
							        				}
							        				else
							        				{
							        					Phiind1=Phiseg4(in1-3*NG*NG/4);
							        				}
							        }
							        if(in2>=0)
							        {
							        	if(in2<NG*NG/4)
							        				{
							        					Phiind2=Phiseg1(in2);
							        				}
							        				else if(in2<NG*NG/2)
							        				{
							        					Phiind2=Phiseg2(in2-NG*NG/4);
							        				}
							        				else if(in2<3*NG*NG/4)
							        				{
							        					Phiind2=Phiseg3(in2-NG*NG/2);
							        				}
							        				else
							        				{
							        					Phiind2=Phiseg4(in2-3*NG*NG/4);
							        				}
							        }

							        	Ex(n)=(Phiind1-Phiind2)/(2*dx);//Electric field in x
							}
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
#pragma omp parallel sections num_threads(4)
		{
#pragma omp section
			{
				for(n=0;n<NG*NG/4;n++)
				{
			        in1=n+1;//Index of above node
			        in2=n-1;//index of lower node
					        Phiind1=0;
					        Phiind2=0;
					        if(in1%NG!=0)
					        {
					        	if(in1<NG*NG/4)
					        				{
					        					Phiind1=Phiseg1(in1);
					        				}
					        				else if(in1<NG*NG/2)
					        				{
					        					Phiind1=Phiseg2(in1-NG*NG/4);
					        				}
					        				else if(in1<3*NG*NG/4)
					        				{
					        					Phiind1=Phiseg3(in1-NG*NG/2);
					        				}
					        				else
					        				{
					        					Phiind1=Phiseg4(in1-3*NG*NG/4);
					        				}
					        }
					        if((in2+NG)%NG!=NG-1)
					        {
					        	if(in2<NG*NG/4)
					        				{
					        					Phiind2=Phiseg1(in2);
					        				}
					        				else if(in2<NG*NG/2)
					        				{
					        					Phiind2=Phiseg2(in2-NG*NG/4);
					        				}
					        				else if(in2<3*NG*NG/4)
					        				{
					        					Phiind2=Phiseg3(in2-NG*NG/2);
					        				}
					        				else
					        				{
					        					Phiind2=Phiseg4(in2-3*NG*NG/4);
					        				}
					        }

					        	Ey(n)=(Phiind1-Phiind2)/(2*dx);//Electric field in x
					}
			}
#pragma omp section
			{
				for(n=NG*NG/4;n<NG*NG/2;n++)
				{
			        in1=n+1;//Index of above node
			        in2=n-1;//index of lower node
					        Phiind1=0;
					        Phiind2=0;
					        if(in1%NG!=0)
					        {
					        	if(in1<NG*NG/4)
					        				{
					        					Phiind1=Phiseg1(in1);
					        				}
					        				else if(in1<NG*NG/2)
					        				{
					        					Phiind1=Phiseg2(in1-NG*NG/4);
					        				}
					        				else if(in1<3*NG*NG/4)
					        				{
					        					Phiind1=Phiseg3(in1-NG*NG/2);
					        				}
					        				else
					        				{
					        					Phiind1=Phiseg4(in1-3*NG*NG/4);
					        				}
					        }
					        if((in2+NG)%NG!=NG-1)
					        {
					        	if(in2<NG*NG/4)
					        				{
					        					Phiind2=Phiseg1(in2);
					        				}
					        				else if(in2<NG*NG/2)
					        				{
					        					Phiind2=Phiseg2(in2-NG*NG/4);
					        				}
					        				else if(in2<3*NG*NG/4)
					        				{
					        					Phiind2=Phiseg3(in2-NG*NG/2);
					        				}
					        				else
					        				{
					        					Phiind2=Phiseg4(in2-3*NG*NG/4);
					        				}
					        }

					        	Ey(n)=(Phiind1-Phiind2)/(2*dx);//Electric field in x
					}
			}
#pragma omp section
			{
				for(n=NG*NG/2;n<NG*NG/4*3;n++)
				{
			        in1=n+1;//Index of above node
			        in2=n-1;//index of lower node
					        Phiind1=0;
					        Phiind2=0;
					        if(in1%NG!=0)
					        {
					        	if(in1<NG*NG/4)
					        				{
					        					Phiind1=Phiseg1(in1);
					        				}
					        				else if(in1<NG*NG/2)
					        				{
					        					Phiind1=Phiseg2(in1-NG*NG/4);
					        				}
					        				else if(in1<3*NG*NG/4)
					        				{
					        					Phiind1=Phiseg3(in1-NG*NG/2);
					        				}
					        				else
					        				{
					        					Phiind1=Phiseg4(in1-3*NG*NG/4);
					        				}
					        }
					        if((in2+NG)%NG!=NG-1)
					        {
					        	if(in2<NG*NG/4)
					        				{
					        					Phiind2=Phiseg1(in2);
					        				}
					        				else if(in2<NG*NG/2)
					        				{
					        					Phiind2=Phiseg2(in2-NG*NG/4);
					        				}
					        				else if(in2<3*NG*NG/4)
					        				{
					        					Phiind2=Phiseg3(in2-NG*NG/2);
					        				}
					        				else
					        				{
					        					Phiind2=Phiseg4(in2-3*NG*NG/4);
					        				}
					        }

					        	Ey(n)=(Phiind1-Phiind2)/(2*dx);//Electric field in x
					}
			}
#pragma omp section
			{
				for(n=NG*NG/4*3;n<NG*NG;n++)
				{
			        in1=n+1;//Index of above node
			        in2=n-1;//index of lower node
					        Phiind1=0;
					        Phiind2=0;
					        if(in1%NG!=0)
					        {
					        	if(in1<NG*NG/4)
					        				{
					        					Phiind1=Phiseg1(in1);
					        				}
					        				else if(in1<NG*NG/2)
					        				{
					        					Phiind1=Phiseg2(in1-NG*NG/4);
					        				}
					        				else if(in1<3*NG*NG/4)
					        				{
					        					Phiind1=Phiseg3(in1-NG*NG/2);
					        				}
					        				else
					        				{
					        					Phiind1=Phiseg4(in1-3*NG*NG/4);
					        				}
					        }
					        if((in2+NG)%NG!=NG-1)
					        {
					        	if(in2<NG*NG/4)
					        				{
					        					Phiind2=Phiseg1(in2);
					        				}
					        				else if(in2<NG*NG/2)
					        				{
					        					Phiind2=Phiseg2(in2-NG*NG/4);
					        				}
					        				else if(in2<3*NG*NG/4)
					        				{
					        					Phiind2=Phiseg3(in2-NG*NG/2);
					        				}
					        				else
					        				{
					        					Phiind2=Phiseg4(in2-3*NG*NG/4);
					        				}
					        }

					        	Ey(n)=(Phiind1-Phiind2)/(2*dx);//Electric field in x
					}
			}
		}

			for(n=0;n<NG;n++)
			{
			    for(m=0;m<NG;m++)
			    {
			    	Ey1(n,m)=Ey(n*NG+m);
			    }
			}

			//Update Position for leap frog
#pragma omp parallel sections num_threads(4)
		{
#pragma omp section
			{
				for(n=0;n<N/4;n++)
					{
						xp(n)+=vxp(n)*DT/2;
						yp(n)+=vyp(n)*DT/2;
					}
			}
#pragma omp section
			{
				for(n=N/4;n<N/2;n++)
					{
						xp(n)+=vxp(n)*DT/2;
						yp(n)+=vyp(n)*DT/2;
					}
			}
#pragma omp section
			{
				for(n=N/2;n<N/4*3;n++)
					{
						xp(n)+=vxp(n)*DT/2;
						yp(n)+=vyp(n)*DT/2;
					}
			}
#pragma omp section
			{
				for(n=N/4*3;n<N;n++)
					{
						xp(n)+=vxp(n)*DT/2;
						yp(n)+=vyp(n)*DT/2;
					}
			}
		}

			 //Update Velocity
					#pragma omp parallel sections num_threads(4)
							{
					#pragma omp section
								{
									for(n=0;n<N/4;n++)
									{
											vxp(n)+=-Ex(index2(g(n,0),g(n,1),NG))*fraz(n,0)*emat(n)*DT/massmat(n);
											vxp(n)+=-Ex(index2(g(N+n,0),g(n,1),NG))*fraz(n,1)*emat(n)*DT/massmat(n);
											vxp(n)+=-Ex(index2(g(n,0),g(N+n,1),NG))*fraz(n+N,0)*emat(n)*DT/massmat(n);
											vxp(n)+=-Ex(index2(g(N+n,0),g(N+n,1),NG))*fraz(n+N,1)*emat(n)*DT/massmat(n);


											vyp(n)+=-Ey(index2(g(n,0),g(n,1),NG))*fraz(n,0)*emat(n)*DT/massmat(n);
											vyp(n)+=-Ey(index2(g(N+n,0),g(n,1),NG))*fraz(n,1)*emat(n)*DT/massmat(n);
											vyp(n)+=-Ey(index2(g(n,0),g(N+n,1),NG))*fraz(n+N,0)*emat(n)*DT/massmat(n);
											vyp(n)+=-Ey(index2(g(N+n,0),g(N+n,1),NG))*fraz(n+N,1)*emat(n)*DT/massmat(n);
										}
								}
					#pragma omp section
								{
									for(n=N/4;n<N/2;n++)
									{
											vxp(n)+=-Ex(index2(g(n,0),g(n,1),NG))*fraz(n,0)*emat(n)*DT/massmat(n);
											vxp(n)+=-Ex(index2(g(N+n,0),g(n,1),NG))*fraz(n,1)*emat(n)*DT/massmat(n);
											vxp(n)+=-Ex(index2(g(n,0),g(N+n,1),NG))*fraz(n+N,0)*emat(n)*DT/massmat(n);
											vxp(n)+=-Ex(index2(g(N+n,0),g(N+n,1),NG))*fraz(n+N,1)*emat(n)*DT/massmat(n);


											vyp(n)+=-Ey(index2(g(n,0),g(n,1),NG))*fraz(n,0)*emat(n)*DT/massmat(n);
											vyp(n)+=-Ey(index2(g(N+n,0),g(n,1),NG))*fraz(n,1)*emat(n)*DT/massmat(n);
											vyp(n)+=-Ey(index2(g(n,0),g(N+n,1),NG))*fraz(n+N,0)*emat(n)*DT/massmat(n);
											vyp(n)+=-Ey(index2(g(N+n,0),g(N+n,1),NG))*fraz(n+N,1)*emat(n)*DT/massmat(n);
										}
								}
					#pragma omp section
								{
									for(n=N/2;n<3*N/4;n++)
									{
											vxp(n)+=-Ex(index2(g(n,0),g(n,1),NG))*fraz(n,0)*emat(n)*DT/massmat(n);
											vxp(n)+=-Ex(index2(g(N+n,0),g(n,1),NG))*fraz(n,1)*emat(n)*DT/massmat(n);
											vxp(n)+=-Ex(index2(g(n,0),g(N+n,1),NG))*fraz(n+N,0)*emat(n)*DT/massmat(n);
											vxp(n)+=-Ex(index2(g(N+n,0),g(N+n,1),NG))*fraz(n+N,1)*emat(n)*DT/massmat(n);


											vyp(n)+=-Ey(index2(g(n,0),g(n,1),NG))*fraz(n,0)*emat(n)*DT/massmat(n);
											vyp(n)+=-Ey(index2(g(N+n,0),g(n,1),NG))*fraz(n,1)*emat(n)*DT/massmat(n);
											vyp(n)+=-Ey(index2(g(n,0),g(N+n,1),NG))*fraz(n+N,0)*emat(n)*DT/massmat(n);
											vyp(n)+=-Ey(index2(g(N+n,0),g(N+n,1),NG))*fraz(n+N,1)*emat(n)*DT/massmat(n);
										}
								}
					#pragma omp section
								{
									for(n=3*N/4;n<N;n++)
									{
											vxp(n)+=-Ex(index2(g(n,0),g(n,1),NG))*fraz(n,0)*emat(n)*DT/massmat(n);
											vxp(n)+=-Ex(index2(g(N+n,0),g(n,1),NG))*fraz(n,1)*emat(n)*DT/massmat(n);
											vxp(n)+=-Ex(index2(g(n,0),g(N+n,1),NG))*fraz(n+N,0)*emat(n)*DT/massmat(n);
											vxp(n)+=-Ex(index2(g(N+n,0),g(N+n,1),NG))*fraz(n+N,1)*emat(n)*DT/massmat(n);


											vyp(n)+=-Ey(index2(g(n,0),g(n,1),NG))*fraz(n,0)*emat(n)*DT/massmat(n);
											vyp(n)+=-Ey(index2(g(N+n,0),g(n,1),NG))*fraz(n,1)*emat(n)*DT/massmat(n);
											vyp(n)+=-Ey(index2(g(n,0),g(N+n,1),NG))*fraz(n+N,0)*emat(n)*DT/massmat(n);
											vyp(n)+=-Ey(index2(g(N+n,0),g(N+n,1),NG))*fraz(n+N,1)*emat(n)*DT/massmat(n);
										}
								}
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






	if(it%100000==99999)
	{

		for(n=0;n<100000;n++)
		{
			myfile1 <<Volt1(it-99999+n)<<"\n";

			myfile2 << Volt2(it-99999+n)<<"\n";

			myfile3 << Charge1(it-99999+n)<<"\n";

			myfile4 << Charge2(it-99999+n)<<"\n";

			myfile7 << Temp1(it-99999+n)<<"\n";

			myfile9 << totUEn(it-99999+n)<<"\n";

			myfile10 << totKEn(it-99999+n)<<"\n";

			myfile11 << totEn(it-99999+n)<<"\n";

			myfile23 << TotalMomentumAbsorbed1X(it-99999+n)<<"\n";

			myfile24 << TotalMomentumAbsorbed1Y(it-99999+n)<<"\n";

			myfile25 << TotalMomentumAbsorbed2X(it-99999+n)<<"\n";

			myfile26 << TotalMomentumAbsorbed2Y(it-99999+n)<<"\n";
		}

		myfile5.open("6rho1allit.txt");
		myfile5 << rho1allit;
		myfile5.close();

		myfile6.open("6Phi1allit.txt");
		myfile6 << Phi1allit;
		myfile6.close();

		myfile8.open("6velocitysquared.txt");
		myfile8 << vxp.array().square()+vyp.array().square();
		myfile8.close();

		myfile12.open("6Phi1.txt");
		myfile12 << Phi1;
		myfile12.close();

		myfile13.open("6rho1.txt");
		myfile13 << rho1;
		myfile13.close();

		myfile14.open("6emat.txt");
		myfile14 << emat;
		myfile14.close();

		myfile15.open("6xp.txt");
		myfile15 << xp;
		myfile15.close();

		myfile16.open("6yp.txt");
		myfile16 << yp;
		myfile16.close();

		myfile17.open("6vxp.txt");
		myfile17 << vxp;
		myfile17.close();

		myfile18.open("6vyp.txt");
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
	myfile23.close();
	myfile24.close();
	myfile25.close();
	myfile26.close();




	cout<< "\nTime="<<(clock ()-starttime)/CLOCKS_PER_SEC;



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
//
//
//
//}

