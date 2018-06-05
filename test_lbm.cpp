#include "./solve/solve.h"
#include "./lbm/Domain.h"
#include <math.h>


void Setup (LBM::Domain &lbm_dom, void * UD)
{
    // myUserData &dat = (*static_cast<myUserData *> (UD));    
    size_t nx = lbm_dom.nx;
    size_t ny = lbm_dom.ny;
    size_t nz = lbm_dom.nz;
    size_t Nneigh = lbm_dom.Nneigh;
    size_t dx = lbm_dom.dx;
    size_t dt = lbm_dom.dt;
    double cswl = sqrt(3)/3.0;
    double csw = 1497.0;
    double ratiov = cswl/csw;//Ul/U
    double ratiot = dt/1e-5;//tl/t
    double ratio = (41.0*7.0)/0.41;
    ratiov = ratio/ratiot;
    // std::cout<<"v = "<<1*ratiov<<std::endl;
    Vec3_t vb(0.001,0.0,0.0);
		for(size_t j =0; j<ny; j++)
		{
			if(lbm_dom.IsSolid[0][j][0]) continue;
			double L  = ny - 2;                       
        	double yp = j;                      
        	double vvx = 1.5*1*ratiov*4/(L*L)*(L*yp - yp*yp);
			vb = vvx, 0.0, 0.0;			
			double *f = lbm_dom.F[0][j][0];		
			for(size_t k=0; k<Nneigh; k++)
			{	
				double VdotC = dot(vb,lbm_dom.C[k]);
                double VdotV = dot(vb,vb); 
				double rho = lbm_dom.Rho[0][j][0];
                double Feq   = lbm_dom.W[k]*rho*(1.0 + 3.0*VdotC/lbm_dom.Cs + 4.5*VdotC*VdotC/(lbm_dom.Cs*lbm_dom.Cs) - 1.5*VdotV/(lbm_dom.Cs*lbm_dom.Cs));
				VdotC = dot(lbm_dom.Vel[1][j][0],lbm_dom.C[k]);
                VdotV = dot(lbm_dom.Vel[1][j][0],lbm_dom.Vel[1][j][0]);
				rho = lbm_dom.Rho[1][j][0];
                double Feq1   = lbm_dom.W[k]*rho*(1.0 + 3.0*VdotC/lbm_dom.Cs + 4.5*VdotC*VdotC/(lbm_dom.Cs*lbm_dom.Cs) - 1.5*VdotV/(lbm_dom.Cs*lbm_dom.Cs));
				f[k] = Feq + lbm_dom.F[1][j][0][k] - Feq1;
            }
            
			lbm_dom.Vel   [0][j][0] = Vec3_t(0.0,0.0,0.0);
	        lbm_dom.Rho   [0][j][0] = 0.0;
    	    for (size_t k=0;k<Nneigh;k++)
    	    {
    	    	lbm_dom.Rho[0][j][0] +=  lbm_dom.F[0][j][0][k];
    	        lbm_dom.Vel[0][j][0] +=  lbm_dom.F[0][j][0][k]*lbm_dom.C[k];
    	    }
            lbm_dom.Vel[0][j][0] /= lbm_dom.Rho[0][j][0];
            
		}

        for(size_t j =0; j<ny; j++)
        {
        if(lbm_dom.IsSolid[nx-1][j][0]) continue;        
        double *f = lbm_dom.F[nx-1][j][0];
        for(size_t k=0; k<Nneigh; k++)
        {	
            f[k] = lbm_dom.F[nx-2][j][0][k];
            }
    
        lbm_dom.Vel   [nx-1][j][0] = Vec3_t(0.0,0.0,0.0);
            lbm_dom.Rho   [nx-1][j][0] = 0.0;
        for (size_t k=0;k<Nneigh;k++)
        {
            lbm_dom.Rho[nx-1][j][0] +=  lbm_dom.F[nx-1][j][0][k];
            lbm_dom.Vel[nx-1][j][0] +=  lbm_dom.F[nx-1][j][0][k]*lbm_dom.C[k];
        }
            lbm_dom.Vel[nx-1][j][0] /= lbm_dom.Rho[nx-1][j][0];
        
        }
}

void Setup1 (LBM::Domain &lbm_dom, void * UD)
{
    size_t nx = lbm_dom.nx;
    size_t ny = lbm_dom.ny;
    size_t nz = lbm_dom.nz;
    size_t Nneigh = lbm_dom.Nneigh;
    size_t dx = lbm_dom.dx;
    size_t dt = lbm_dom.dt;
    size_t Nproc = lbm_dom.Nproc;
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(size_t iy=0; iy<ny; ++iy)
    {
        if(!lbm_dom.IsSolid[0][iy][0])
        {
            for(size_t k=0; k<Nneigh; ++k)
            {
                lbm_dom.F[0][iy][0][k] = lbm_dom.Feq(k,1.05,lbm_dom.Vel[1][iy][0]) + lbm_dom.F[1][iy][0][k] - lbm_dom.Feq(k,lbm_dom.Rho[1][iy][0],lbm_dom.Vel[1][iy][0]);
            }
        }
        if(!lbm_dom.IsSolid[nx-1][iy][0])
        {
            for(size_t k=0; k<Nneigh; ++k)
            {
                lbm_dom.F[nx-1][iy][0][k] = lbm_dom.Feq(k,1.0,lbm_dom.Vel[nx-2][iy][0]) + lbm_dom.F[nx-2][iy][0][k] - lbm_dom.Feq(k,lbm_dom.Rho[nx-2][iy][0],lbm_dom.Vel[nx-2][iy][0]);
            }
        }
    }   
}

int main (int argc, char **argv) try
{    
    LBM::Domain lbm_dom;
    
    lbm_dom.Nproc = 1;
    
    
    
    size_t nx = 200;
    size_t ny = 50;//0.41 m
    size_t nz = 1;
    std::cout<<"nx = "<<nx<<std::endl;
    std::cout<<"ny = "<<ny<<std::endl;
    double dx = 1.0;
    double dt = 1.0;

    double nu = 0.5/3.0;
    std::cout<<"nu = "<<nu<<std::endl;
    lbm_dom.AllocateDomainMem(nx,ny,nz,nu,dx,dt);
    lbm_dom.Tau = 1;
    lbm_dom.Sc = 0.0;
	//initial
    double rho = 1.0;
    lbm_dom.Rho0 = rho;
    Vec3_t v0(0.0,0.0,0.0);
    for(size_t ix=0; ix<nx; ix++)
    for(size_t iy=0; iy<ny; iy++)
    for(size_t iz=0; iz<nz; iz++)
    {
        lbm_dom.Rho[ix][iy][iz] = rho;
        lbm_dom.Vel[ix][iy][iz] = v0;
		lbm_dom.Check[ix][iy][iz] = 0.0;
        // lbm_dom.BForce[ix][iy][iz] = 0.000141*3,0.0,0.0;
        // size_t s = 0;
        for(size_t k=0; k<lbm_dom.Nneigh; k++)
        {
            lbm_dom.F[ix][iy][iz][k] = lbm_dom.Feq(k,rho,v0);
        }
    }
    
	
    for (size_t ix=0;ix<nx;ix++)
	for (size_t iz=0;iz<nz;iz++)
	{
		lbm_dom.IsSolid[ix][0][iz] = true;
		lbm_dom.IsSolid[ix][ny-1][iz] = true;
    }
    

    
    

    double tt = 2.5e5;
    Solve( lbm_dom, tt, 250,Setup1,NULL,"my_test1", true);    

}MECHSYS_CATCH