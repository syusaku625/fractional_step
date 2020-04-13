#include<iostream>
#include<vector>
#include<cmath>
#include<fstream>
using namespace std;
#define xn 31
#define yn 31
#define Re 100
#define lm 20000

inline void update(double dx, double dy, double dt, vector<double> &u, vector<double> &v, vector<double> &p, double &divv, double uwall)
{
	int i,j;
	
	//u//
	for(i=2;i<xn;i++){
		for(j=1;j<yn;j++){
			u[i*(yn+1)+j]=u[i*(yn+1)+j]+dt*(p[(i-1)*(yn+1)+j]-p[i*(yn+1)+j])/dx;
		}
	}
	//v//
	for(i=1;i<xn;i++){
		for(j=2;j<yn;j++){
			v[i*(yn+2)+j]=v[i*(yn+2)+j]+dt*(p[i*(yn+1)+j-1]-p[i*(yn+1)+j])/dy;
		}
	}
	divv=0.0;
	for(i=1;i<xn;i++){
		for(j=1;j<yn;j++){
			divv+=fabs((u[(i+1)*(yn+1)+j]-u[i*(yn+1)+j])/dx+(v[i*(yn+2)+j+1]-v[i*(yn+2)+j])/dy);
		}
	}
}


inline void poi(int &km, double dx, double dy, double dt, double &err, vector<double> &u, vector<double> &v, vector<double> &p)
{
	int i,j,k;
	double C1=dy*dy/(2.0*(dx*dx+dy*dy));
	double C2=dx*dx/(2.0*(dx*dx+dy*dy));
	double C3=dx*dx*dy*dy/(2.0*(dx*dx+dy*dy))/dt;
	double in;
	
	//Poisson equation//
	for(k=1;k<=km;k++){
		err=0.;
		//Neumann BC//
		for(j=0;j<yn+1;j++){
			p[0*(yn+1)+j]=p[1*(yn+1)+j];
			p[xn*(yn+1)+j]=p[(xn-1)*(yn+1)+j];
		}
		for(i=0;i<xn+1;i++){
			p[i*(yn+1)+0]=p[i*(yn+1)+1];
			p[i*(yn+1)+yn]=p[i*(yn+1)+yn-1];
		}
		//iteration//
		for(i=1;i<xn;i++){
			for(j=1;j<yn;j++){
				in=p[i*(yn+1)+j];
				p[i*(yn+1)+j]=C1*(p[(i+1)*(yn+1)+j]+p[(i-1)*(yn+1)+j])+C2*(p[i*(yn+1)+j+1]+p[i*(yn+1)+j-1])-C3*((u[(i+1)*(yn+1)+j]-u[i*(yn+1)+j])/dx+(v[i*(yn+2)+j+1]-v[i*(yn+2)+j])/dy);
				err+=(p[i*(yn+1)+j]-in)*(p[i*(yn+1)+j]-in);
			}
		}
		if(err<=0.001) break;
	}
}


inline void vel(vector<double> &u, vector<double> &v, double dx, double dy, double dt, double uwall)
{
    double uad,vad;
	double udif,vdif;
	double umid,vmid;
    
    for(int i=0; i<yn+1; i++){
        //BC for left wall
        u[(yn+1)+i]=0.0;
        u[i]=u[2*(yn+1)+i];
        v[i]=-v[(yn+1)+i];
        //BC for right wall
        u[xn*(yn+1)+i]=0.0;
        u[(xn+1)*(yn+1)+i]=u[(xn-1)*(yn+1)+i];
        v[xn*(yn+2)+i]=-v[(xn-1)*(yn+2)+i];
    }
    //corner 
    v[yn+1]=-v[(yn+2)+yn+1];
    v[xn*(yn+2)+yn+1]=-v[(xn-1)*(yn+2)+yn+1];

    for(int i=0; i<xn+1; i++){
        //BC for bottom
        v[i*(yn+2)+1]=0.0;
        v[i*(yn+2)]=v[i*(yn+2)+2];
        u[i*(yn+1)]=-u[i*(yn+1)+1];
        //BC for top
        v[i*(yn+2)+yn]=0.0;
        v[i*(yn+2)+yn+1]=v[i*(yn+2)+yn-1];
        u[i*(yn+1)+yn]=2.0*uwall-u[i*(yn+1)+yn-1];
    }
    //corner
    u[(xn+1)*(yn+1)]=-u[xn*(yn+1)];
	u[(xn+1)*(yn+1)+yn]=-u[xn*(yn+1)+yn];

    //u
    for(int i=2;i<xn;i++){
		for(int j=1;j<yn;j++){
			vmid=(v[i*(yn+2)+j]+v[i*(yn+2)+j+1]+v[(i-1)*(yn+2)+j+1]+v[(i-1)*(yn+2)+j])/4.0;
			uad=u[i*(yn+1)+j]*(u[(i+1)*(yn+1)+j]-u[(i-1)*(yn+1)+j])/2.0/dx+vmid*(u[i*(yn+1)+j+1]-u[i*(yn+1)+j-1])/2.0/dy;
			udif=(u[(i+1)*(yn+1)+j]-2.0*u[i*(yn+1)+j]+u[(i-1)*(yn+1)+j])/dx/dx+(u[i*(yn+1)+j+1]-2.0*u[i*(yn+1)+j]+u[i*(yn+1)+j-1])/dy/dy;
			u[i*(yn+1)+j]=u[i*(yn+1)+j]+dt*(-uad+1.0/Re*udif);
		}
	}

    //v//
	for(int i=1;i<xn;i++){
		for(int j=2;j<yn;j++){
			umid=(u[i*(yn+1)+j]+u[(i+1)*(yn+1)+j]+u[(i+1)*(yn+1)+j-1]+u[i*(yn+1)+j-1])/4.0;
			vad=umid*(v[(i+1)*(yn+2)+j]-v[(i-1)*(yn+2)+j])/2.0/dx+v[i*(yn+2)+j]*(v[i*(yn+2)+j+1]-v[i*(yn+2)+j-1])/2.0/dy;
			vdif=(v[(i+1)*(yn+2)+j]-2.0*v[i*(yn+2)+j]+v[(i-1)*(yn+2)+j])/dx/dx+(v[i*(yn+2)+j+1]-2.0*v[i*(yn+2)+j]+v[i*(yn+2)+j-1])/dy/dy;
			v[i*(yn+2)+j]=v[i*(yn+2)+j]+dt*(-vad+1.0/Re*vdif);
		}
	}
}

int main()
{
    vector<double> u((xn+2)*(yn+1)), v((xn+1)*(yn+2)), p((xn+1)*(yn+1));
    double lx=1.0, ly=1.0;
    double dx=lx/(xn-1);
    double dy=ly/(yn-1);
    double dt=0.001;
    double uwall=1.0;
    double err;
    int km=100;
    double divv;
    ofstream fk,ff;
	fk.open("velocity.vtk");
	ff.open("pressure.vtk");
    //initialization
    for(int i=0; i<xn+1; i++){
        for(int j=0; j<yn+1; j++){
            p[i*(xn+1)+j]=0.0;           
        }
    }
    for(int i=0; i<xn+2; i++){
        for(int j=0; j<yn+1; j++){
            u[i*(yn+1)+j]=0.0;
        }
    }
    for(int i=0; i<xn+1; i++){
        for(int j=0; j<yn+2; j++){
            v[i*(yn+1)+j]=0.0;
        }
    }
    //time step
    for(int l=0; l<=lm; l++){
        for(l=1;l<=lm;l++){	
            vel(u, v, dx, dy, dt, uwall);
            poi(km, dx, dy, dt, err, u, v, p);
            update(dx, dy, dt, u, v, p, divv, uwall);
            if(l%1000==0) cout<<l<<" "<<err<<" "<<divv<<endl;
        }
    }
    //output
    fk << "# vtk DataFile Version 3.0" << endl;
    fk << "vtk output" << endl;
    fk << "ASCII" << endl;
    fk << "DATASET UNSTRUCTURED_GRID" << endl;
    fk << "POINTS " << xn*yn << " double" << endl;
    for(int i=0; i<xn; i++){
	    for(int j=0; j<yn; j++){
            fk<< i*dx << " " << j*dy << " " << 0 << endl;
	    }
	}
    fk << "POINT_DATA " << xn*yn << endl;
    fk << "VECTORS velocity[m/s] double" << endl; 
	for(int i=0; i<xn; i++){
	    for(int j=0; j<yn; j++){
	        fk << u[i*(yn+1)+(yn+2)+j] << " " << v[i*(yn+2)+(yn+2)+1+j] << " " << 0 << endl;
	    }
	}

    ff<< "# vtk DataFile Version 3.0" << endl;
    ff<< "vtk output" << endl;
    ff<< "ASCII" << endl;
    ff<< "DATASET STRUCTURED_POINTS" << endl;
    


    for(int i=1; i<xn; i++){
		for(int j=1; j<yn; j++){
			//fk<<double((i-0.5)*dx)<<" "<<double((j-0.5)*dy)<<" "<<(u[i*(yn+1)+j]+u[(i+1)*(yn+1)+j])/2.0<<" "<<(v[i*(yn+2)+j]+v[i*(yn+2)+j+1])/2.0<<endl;
        	ff<<double((i-0.5)*dx)<<" "<<double((j-0.5)*dy)<<" "<<p[i*(yn+1)+j]<<endl;
		}
	}
}