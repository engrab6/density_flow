#include "./solve/solve.h"


void Setup3(LBM::Domain &fluid_dom, SCALAR::Domain &con_dom, SCALAR::Domain &temp_dom, void *UD)
{
    size_t Nproc = fluid_dom.Nproc;
    size_t nx = fluid_dom.nx;
    size_t ny = fluid_dom.ny;
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif    
    for(size_t ix=0; ix<fluid_dom.nx; ++ix)
    for(size_t iz=0; iz<fluid_dom.nz; ++iz)
    {
            // std::cout<<fluid_dom.IsSolid[0][iy][iz]<<std::endl;        
        size_t index = ny-2;
        if(!temp_dom.IsSolid[ix][index][iz])
        {
            // temp_dom.F[index][iy][iz][0] = 0.5 - temp_dom.F[index][iy][iz][2]; 
            temp_dom.F[ix][index][iz][0] = 0.0 - temp_dom.F[ix][index][iz][2]-temp_dom.F[ix][index][iz][1]-temp_dom.F[ix][index][iz][3]; 
            con_dom.F[ix][index][iz][0] = 1.0 - con_dom.F[ix][index][iz][2]-con_dom.F[ix][index][iz][1]-con_dom.F[ix][index][iz][3]; 
            // temp_dom.Rho   [index][iy][iz] = 0.0;
            // for (size_t k=0;k<temp_dom.Nneigh;k++)
            // {
            //     temp_dom.Rho[index][iy][iz] +=  temp_dom.F[index][iy][iz][k];
            // }
            // std::cout<<index<<" "<<iy<<" "<<temp_dom.Rho[index][iy][iz]<<std::endl;
        }
        
    	// index = nx-2;
        // if(!temp_dom.IsSolid[index][iy][iz])
		// {
        //     temp_dom.F[index][iy][iz][2] =  0.0 - temp_dom.F[index][iy][iz][0]-temp_dom.F[index][iy][iz][1]-temp_dom.F[index][iy][iz][3];
        //     // temp_dom.Rho   [index][iy][iz] = 0.0;
        //     // for (size_t k=0;k<temp_dom.Nneigh;k++)
        //     // {
        //     //     temp_dom.Rho[index][iy][iz] +=  temp_dom.F[nx-1][iy][iz][k];
        //     // }
        // }
        index = 1;
        if(!temp_dom.IsSolid[ix][index][iz])
		{
            temp_dom.F[ix][index][iz][0] = 1.0 - temp_dom.F[ix][index][iz][2]-temp_dom.F[ix][index][iz][1]-temp_dom.F[ix][index][iz][3]; 
            con_dom.F[ix][index][iz][0] = 0.0 - con_dom.F[ix][index][iz][2]-con_dom.F[ix][index][iz][1]-con_dom.F[ix][index][iz][3]; 
            
            
        }
    }
    //adiadtic
    // #ifdef USE_OMP
    // #pragma omp parallel for schedule(static) num_threads(Nproc)
    // #endif
    // for(size_t ix=0; ix<nx; ++ix)
    // {
    //     for (size_t k=0; k<temp_dom.Nneigh; ++k)
    //     {
    //         temp_dom.F[ix][1][0][k] = temp_dom.F[ix][2][0][k];
    //         temp_dom.F[ix][ny-2][0][k] = temp_dom.F[ix][ny-3][0][k];
            
    //     }
    // }
}


int main () try
{
    LBM::Domain fluid_dom;
    SCALAR::Domain temp_dom;
    SCALAR::Domain con_dom;
    int Nproc = 6;
    fluid_dom.Nproc = Nproc;
    temp_dom.Nproc = Nproc;
    con_dom.Nproc = Nproc;
    // myUserData my_dat;
    // UserData = &my_dat;

    size_t nx = 102;
    size_t ny = 102;
    size_t nz = 1;
    double dx = 1;
    double dt = 1;
    double mu = 0.01;
    double Pr = 7;//mu/alpha
    double Ra = 1e6;//g*beta*dT*M^3/(alpha*mu)
    double N = 1.22; //betac*dC/betat*dT
    double Le = 100; //alpha/D
    double alpha = mu/Pr;
    double D = alpha/Le;
    fluid_dom.AllocateDomainMem(nx,ny,nz,mu,dx,dt);
    temp_dom.AllocateDomainMem(nx,ny,nz,alpha,dx,dt);
    con_dom.AllocateDomainMem(nx,ny,nz,D,dx,dt);
    fluid_dom.Sc = 0.0;
    temp_dom.beta = 2;//volumn thermal expansion coefficient
    con_dom.beta = N*temp_dom.beta;
    double M = (double) ny;
    fluid_dom.g = 0.0,-Ra*alpha*mu/(M*M*M*1.0*temp_dom.beta),0.0;
    std::cout<<"g= "<<fluid_dom.g(1)<<std::endl;
    std::cout<<fluid_dom.g(1)*temp_dom.beta<<std::endl;
    //initial
    
    double rho = 1.0;
    double temp = 0.0;
    double con = 0.0;
    Vec3_t v0(0.0,0.0,0.0);
    for(size_t ix=0; ix<nx; ix++)
    for(size_t iy=0; iy<ny; iy++)
    for(size_t iz=0; iz<nz; iz++)
    {
        fluid_dom.Rho[ix][iy][iz] = rho;
        fluid_dom.Vel[ix][iy][iz] = v0;
        temp_dom.Rho[ix][iy][iz] = temp;
        temp_dom.Vel[ix][iy][iz] = v0;
        con_dom.Rho[ix][iy][iz] = con;
        con_dom.Vel[ix][iy][iz] = v0;
        fluid_dom.Check[ix][iy][iz] = 0.0;
        for(size_t k=0; k<fluid_dom.Nneigh; ++k)
        {
            fluid_dom.F[ix][iy][iz][k] = fluid_dom.Feq(k,rho,v0);            
        }
        for(size_t k=0; k<temp_dom.Nneigh; ++k)
        {
            temp_dom.F[ix][iy][iz][k] = temp_dom.Feq(k,temp,v0);            
        }
        for(size_t k=0; k<con_dom.Nneigh; ++k)
        {
            con_dom.F[ix][iy][iz][k] = con_dom.Feq(k,con,v0);            
        }
    }
    //boundary condition
    
    for(size_t ix=0; ix<nx; ix++)
    for(size_t iz=0; iz<nz; iz++)
    {
        fluid_dom.IsSolid[ix][0][iz] = true;
        fluid_dom.IsSolid[ix][ny-1][iz] = true;    
    }
        // fluid_dom.IsSolid[1][1][0] = true;
        // fluid_dom.IsSolid[nx-2][1][0] = true;
        // fluid_dom.IsSolid[1][ny-2][0] = true;
        // fluid_dom.IsSolid[nx-2][1][0] = true;
    
   
    for(size_t ix=1; ix<nx-1; ix++)
    for(size_t iz=0; iz<nz; iz++)
    {
        fluid_dom.IsBoundary[ix][1][iz] = true;
        fluid_dom.IsBoundary[ix][ny-2][iz] = true;

        temp_dom.IsBoundary[ix][1][iz] = true;
        temp_dom.IsBoundary[ix][ny-2][iz] = true;
        con_dom.IsBoundary[ix][1][iz] = true;
        con_dom.IsBoundary[ix][ny-2][iz] = true;
    }
        // fluid_dom.IsBoundary[1][1][0] = false;
        // fluid_dom.IsBoundary[nx-2][1][0] = false;
        // fluid_dom.IsBoundary[1][ny-2][0] = false;
        // fluid_dom.IsBoundary[nx-2][1][0] = false;
    for(size_t ix=0; ix<nx; ix++)
    {
        temp_dom.IsSolid[ix][0][0] = true;
        temp_dom.IsSolid[ix][ny-1][0] = true;
        con_dom.IsSolid[ix][0][0] = true;
        con_dom.IsSolid[ix][ny-1][0] = true;
    }

    double Ttotal = 2.5e6;
    double tt = Ttotal/1e3;
    // Ttotal = 100;
    // tt = 1;

    Solve(fluid_dom,con_dom,temp_dom,Ttotal,tt,NULL,Setup3,NULL,"test_sta",true);
    // Solve(fluid_dom, Ttotal,tt,NULL,NULL,"test1",true);

}MECHSYS_CATCH