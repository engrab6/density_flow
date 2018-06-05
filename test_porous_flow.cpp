#include "./solve/solve.h"

struct myUserData
{
    std::ofstream oss_ss;
};

void Report(LBM::Domain &fluid_dom, SCALAR::Domain &con_dom, SCALAR::Domain &temp_dom, void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    size_t nx = fluid_dom.nx;
    size_t ny = fluid_dom.ny;
    size_t nz = fluid_dom.nz;
    double dx = fluid_dom.dx;
    if(fluid_dom.Time <1e-6)
    {
        String fs;
        fs.Printf("%s.out","temperature1");
        dat.oss_ss.open(fs.CStr(),std::ios::out);
        dat.oss_ss<<Util::_10_6<<"Time"<<Util::_8s<<"NuH"<<Util::_8s<<"NuC"<<Util::_8s<<"Nua"<<Util::_8s<<"Temperature\n";
    }else{
        size_t index = 1;
        double S = 0;
        for(size_t iy = 0; iy<ny; ++iy)
        {
            double sita = temp_dom.Rho[index][iy][0] - temp_dom.Rho[index+1][iy][0];
            S += sita;
        }
        index = nx-3;
        double S1 = 0;
        for(size_t iy = 0; iy<ny; ++iy)
        {
            double sita = temp_dom.Rho[index][iy][0] - temp_dom.Rho[index+1][iy][0];
            S1 += sita;
        }
        dat.oss_ss<<Util::_10_6<<fluid_dom.Time<<Util::_8s<<S<<Util::_8s<<S1<<Util::_8s<<0.5*(S+S1)<<Util::_8s<<temp_dom.Rho[49][49][0]<<std::endl;
    }
}



void Setup2(LBM::Domain &fluid_dom, SCALAR::Domain &con_dom, SCALAR::Domain &temp_dom, void *UD)
{
    size_t Nproc = fluid_dom.Nproc;
    size_t nx = fluid_dom.nx;
    size_t ny = fluid_dom.ny;

    if(fluid_dom.Time > fluid_dom.Ttotal-4.0*fluid_dom.tout)
    {
        isF = true;
    }

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif    
    for(size_t ix=0; ix<fluid_dom.nx; ++ix)
    for(size_t iz=0; iz<fluid_dom.nz; ++iz)
    {
            // std::cout<<fluid_dom.IsSolid[0][iy][iz]<<std::endl;        
        size_t index = 1;
        double *ff = fluid_dom.F[ix][index][iz];
        double *ff1 = fluid_dom.F[ix][index+1][iz];
        for(size_t k=0; k< fluid_dom.Nneigh; ++k)
        {
            Vec3_t vel = fluid_dom.Vel[ix][index+1][iz];
            double rho = fluid_dom.Rho[ix][index+1][iz];
            double rhob = 1.0;
            ff[k] = fluid_dom.Feq(k,rhob,vel) +  ff1[k] -fluid_dom.Feq(k,rho,vel);
            
        }
        index = ny-1;
        // for(size_t k=0; k< fluid_dom.Nneigh; ++k)
        // {
        //     Vec3_t vel = fluid_dom.Vel[ix][index-1][iz];
        //     double rho = fluid_dom.Rho[ix][index-1][iz];
        //     double rho1 = fluid_dom.Rho[ix][index][iz];
        //     Vec3_t velb(0.1,0.0,0.0);
        //     fluid_dom.F[ix][index][iz][k] = fluid_dom.Feq(k,rho,velb) +  fluid_dom.F[ix][index-1][iz][k] -fluid_dom.Feq(k,rho,vel);
        // }
        // fluid_dom.Vel   [ix][index][iz] = Vec3_t(0.0,0.0,0.0);
        // fluid_dom.Rho   [ix][index][iz] = 0.0;
        // for (size_t k=0;k<fluid_dom.Nneigh;k++)
        // {
        //     fluid_dom.Rho[ix][index][iz] +=  fluid_dom.F[ix][index][iz][k];
        //     fluid_dom.Vel[ix][index][iz] +=  fluid_dom.F[ix][index][iz][k]*fluid_dom.C[k];
        // }
        // fluid_dom.Vel[ix][index][iz] *= fluid_dom.Cs/fluid_dom.Rho[ix][index][iz];
        //Zhou He
        double ux = 0.01;
        double uy = 0.0;
        double *f = fluid_dom.F[ix][index][iz];
        
        double rho = 1.0/(1.0-uy)*(f[0]+f[1]+f[3]+2.0*(f[2]+f[5]+f[6]));
        fluid_dom.Vel[ix][index][iz] = ux,uy,0.0;
        fluid_dom.Rho[ix][index][iz] = rho;
        f[4] = f[2]- 2.0/3.0 * fluid_dom.Rho[ix][index][iz]*uy;
        f[8] = f[6] + 0.5*(f[3] - f[1]) + 0.5*rho*ux - 1.0/6.0*rho*uy;
        f[7] = f[5] + 0.5*(f[1] - f[3]) - 0.5*rho*ux - 1.0/6.0*rho*uy;

    }
        
}

void Setup3(LBM::Domain &fluid_dom, SCALAR::Domain &con_dom, SCALAR::Domain &temp_dom, void *UD)
{
    size_t Nproc = fluid_dom.Nproc;
    size_t nx = fluid_dom.nx;
    size_t ny = fluid_dom.ny;
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif    
    for(size_t ix=0; ix<fluid_dom.nx; ++ix)
    {
        
        size_t index = ny-1;
        double *ft = temp_dom.F[ix][index][0];
        double *fc = con_dom.F[ix][index][0];
        double *fc1 = con_dom.F[ix][index-1][0];
        double *ft1 = temp_dom.F[ix][index-1][0];
        Vec3_t vel1 = fluid_dom.Vel[ix][index-1][0];
        double temp1 = temp_dom.Rho[ix][index-1][0];
        double con1 = con_dom.Rho[ix][index-1][0];
        // ft[3] = 1.0 - ft[0] - ft[1] -ft[2];
        // fc[3] = 1.0 - fc[0] - fc[1] -fc[2];
        for(size_t k=0; k<temp_dom.Nneigh; ++k)
        {
            ft[k] = temp_dom.Feq(k,1.0,vel1) + ft1[k] - temp_dom.Feq(k,temp1,vel1);
            fc[k] = con_dom.Feq(k,1.0,vel1) + fc1[k] - con_dom.Feq(k,con1,vel1);  
        }
        // ft[3] = 0.5 - ft[1];
        // fc[3] = 0.5 - fc[1];


        index = 0;
        Vec3_t vel = fluid_dom.Vel[ix][index][0];
        fc1 = con_dom.F[ix][index+1][0];
        ft1 = temp_dom.F[ix][index+1][0];
        
        if(vel(1)>0)
        {
            
            ft = temp_dom.F[ix][index][0];
            fc = con_dom.F[ix][index][0];
            ft[0] = 0.0 - ft[3] - ft[1] -ft[2];
            fc[0] = 0.0 - fc[3] - fc[1] -fc[2];
        }else{
            for(size_t k=0; k<temp_dom.Nneigh; ++k)
            {
                ft[k] = ft1[k];
                fc[k] = fc1[k];
            }
        }
        
    }
}

int main () try
{
    LBM::Domain fluid_dom;
    SCALAR::Domain temp_dom;
    SCALAR::Domain con_dom;
    int Nproc = 7;
    fluid_dom.Nproc = Nproc;
    temp_dom.Nproc = Nproc;
    con_dom.Nproc = Nproc;
    myUserData my_dat;
    UserData = &my_dat;

    size_t nx = 100-4;
    size_t ny = 102;
    size_t nz = 1;
    double dx = 1;
    double dt = 1;
    double mu = 0.01;
    double Pr = 0.7;//mu/alpha
    double Ra = 1e5;//g*beta*dT*M^3/(alpha*mu)
    double N = -1.05; //betac*dC/betat*dT
    double Le = 2; //alpha/D
    double alpha = mu/Pr;
    double D = alpha/Le;
    fluid_dom.AllocateDomainMem(nx,ny,nz,mu,dx,dt);
    temp_dom.AllocateDomainMem(nx,ny,nz,alpha,dx,dt);
    con_dom.AllocateDomainMem(nx,ny,nz,D,dx,dt);
    fluid_dom.Sc = 0.0;
    temp_dom.beta = 2;//volumn thermal expansion coefficient
    con_dom.beta = N*temp_dom.beta;
    double M = (double) 350;
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
    size_t d = 10;
    size_t porex = 6;
    size_t porey = 6;
    size_t depth = 50;
    
    for(size_t iy=porey/2+4; iy<depth; iy +=d+porey)
    {
        for(size_t iiy=0; iiy<d; ++iiy)
        {
            for(size_t ix=porex/2; ix< nx-porex/2; ix += d+porex)
            {
                for(size_t iix=0; iix< d; ++iix)
                {
                    fluid_dom.IsSolid[ix+iix][iy+iiy][0] = true;
                    temp_dom.IsSolid[ix+iix][iy+iiy][0] = true;
                    con_dom.IsSolid[ix+iix][iy+iiy][0] = true;

                }
            }
        }   
    }
    for(size_t ix=0; ix<nx; ++ix)
    {
        fluid_dom.IsSolid[ix][0][0] = true;
    }
    for(size_t ix=0; ix<nx; ++ix)
    for(size_t iy=0; iy<ny; ++iy)
    {
        if(fluid_dom.IsSolid[ix][iy][0])
        {
            for(size_t k=0; k<fluid_dom.Nneigh; ++k)
            {
                size_t nix = (size_t)((int)ix + (int)fluid_dom.C[k](0) + (int)nx)%nx;
                size_t niy = (size_t)((int)iy + (int)fluid_dom.C[k](1) + (int)ny)%ny;
                if(!fluid_dom.IsSolid[nix][niy][0])
                {
                    fluid_dom.IsBoundary[nix][niy][0] = true;
                }
            }
        }
    }
    for(size_t ix=0; ix<nx; ++ix)
    {
        fluid_dom.IsBoundary[ix][ny-1][0] = false;
        fluid_dom.IsBoundary[ix][1][0] = false;
    }
    

    double Ttotal = 2.5e5;
    double tt = Ttotal/1e3;
    fluid_dom.Ttotal = Ttotal;
    fluid_dom.tout = tt;
    // Ttotal = 100;
    // tt = 1;
    Solve(fluid_dom,con_dom,temp_dom,Ttotal,tt,Setup2,NULL,NULL,"test_porous_flow",true);
    // Solve(fluid_dom, Ttotal,tt,NULL,NULL,"test1",true);

}MECHSYS_CATCH