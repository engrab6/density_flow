#ifndef MECHSYS_SOLVE_H
#define MECHSYS_SOLVE_H

// #include "../scalar/Domain.h"
#include "../lbm/Domain.h"
#include "../scalar/Domain.h"


void *UserData;
typedef void (*ptDFun_t) (LBM::Domain &fluid_dom, SCALAR::Domain &con_dom, SCALAR::Domain &temp_dom, void *UserData);
typedef void (*ptDFun_t1) (LBM::Domain &fluid_dom, void *UserData);
bool isF;

void WriteXDMF(char const * FileKey, LBM::Domain &fluid_dom, SCALAR::Domain &con_dom, SCALAR::Domain &temp_dom)
{
    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    size_t nx = fluid_dom.nx;
    size_t ny = fluid_dom.ny;
    size_t nz = fluid_dom.nz;

    size_t  Nx = nx;
    size_t  Ny = ny;
    size_t  Nz = nz;
    size_t Step = 1;
    size_t Nl = 1;
    // std::cout<<Nx<<" "<<Ny<<" "<<Nz<<std::endl;
    for (size_t j=0;j<Nl;j++)
    {
        // Creating data sets
        double * Density   = new double[  Nx*Ny*Nz];
        double * Concentration   = new double[  Nx*Ny*Nz];
        double * Temperature   = new double[  Nx*Ny*Nz];
        double * Ga     = new double[  Nx*Ny*Nz];
        double * Vvec      = new double[3*Nx*Ny*Nz];
        double * BFvec      = new double[3*Nx*Ny*Nz];
        double * Ff     = new double [fluid_dom.Nneigh*Nx*Ny*Nz];
        size_t i=0;
        for (size_t m=0;m<nz;m+=Step)
        for (size_t l=0;l<ny;l+=Step)
        for (size_t n=0;n<nx;n+=Step)
        {
            double rho    = 0.0;
            double temp   = 0.0;
            double con    = 0.0;
            double gamma  = 0.0;
            Vec3_t vel    = Vec3_t(0.0,0.0,0.0);
            Vec3_t velp    = Vec3_t(0.0,0.0,0.0);
            Vec3_t accp    = Vec3_t(0.0,0.0,0.0);
            Vec3_t BF     = Vec3_t(0.0,0.0,0.0);
            double ttemp = 0.0;
            for (size_t ni=0;ni<Step;ni++)
            for (size_t li=0;li<Step;li++)
            for (size_t mi=0;mi<Step;mi++)
            {
                rho    += fluid_dom.Rho    [n+ni][l+li][m+mi];
                temp   += temp_dom.Rho    [n+ni][l+li][m+mi];
                con   += con_dom.Rho    [n+ni][l+li][m+mi];
                ttemp    = fluid_dom.IsSolid[n+ni][l+li][m+mi] ? 2.0: 0.0;
                // gamma  += std::max(fluid_dom.Gamma[n+ni][l+li][m+mi],ttemp);
                gamma = ttemp;
                ttemp    = fluid_dom.IsBoundary[n+ni][l+li][m+mi] ? 1.0: 0.0;
                gamma += ttemp;
                vel    += fluid_dom.Vel    [n+ni][l+li][m+mi];
                BF    += fluid_dom.BForce    [n+ni][l+li][m+mi];
                velp    += fluid_dom.VelP    [n+ni][l+li][m+mi];
                accp    += fluid_dom.AccP    [n+ni][l+li][m+mi];
            }
            rho  /= Step*Step*Step;
            temp  /= Step*Step*Step;
            con  /= Step*Step*Step;
            gamma/= Step*Step*Step;
            vel  /= Step*Step*Step;
            velp  /= Step*Step*Step;
            accp  /= Step*Step*Step;
            BF   /= Step*Step*Step;
            Ga   [i]  = (double) gamma;
            Density [i]  = (double) rho;            
            Concentration [i]  = (double) con;            
            Temperature [i]  = (double) temp;            
            Vvec[3*i  ]  = (double) vel(0);//*(1.0-Ga[i]);
            Vvec[3*i+1]  = (double) vel(1);//*(1.0-Ga[i]);
            Vvec[3*i+2]  = (double) vel(2);//*(1.0-Ga[i]);
            BFvec[3*i ]   = (double) BF(0);
            BFvec[3*i+1]   = (double) BF(1);
            BFvec[3*i+2]   = (double) BF(2);
            if(isF)
            {
                for (size_t k=0; k<fluid_dom.Nneigh; k++)
                {
                    Ff[fluid_dom.Nneigh*i + k] = (double) fluid_dom.F[n][l][m][k];
                }
            }
            
            i++;
        }
        
        //Writing data to h5 file
        hsize_t dims[1];
        dims[0] = Nx*Ny*Nz;
        String dsname;
        dsname.Printf("Density_%d",j);
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Density );
        dims[0] = Nx*Ny*Nz;
        dsname.Printf("Temperature_%d",j);
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Temperature );
        dims[0] = Nx*Ny*Nz;
        dsname.Printf("Concentration_%d",j);
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Concentration );
        if (j==0)
        {
            dsname.Printf("Gamma");
            H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Ga   );
        }
        dims[0] = 3*Nx*Ny*Nz;
        dsname.Printf("Velocity_%d",j);
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Vvec    );
        dims[0] = 3*Nx*Ny*Nz;
        dsname.Printf("BForce_%d",j);
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,BFvec    );
        if(isF)
        {
            dims[0] = fluid_dom.Nneigh*Nx*Ny*Nz;
            dsname.Printf("F_%d",j);
            H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Ff    );
        }
        
        dims[0] = 1;
        int N[1];
        N[0] = Nx;
        dsname.Printf("Nx");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,N);
        dims[0] = 1;
        N[0] = Ny;
        dsname.Printf("Ny");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,N);
        dims[0] = 1;
        N[0] = Nz;
        dsname.Printf("Nz");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,N);

        delete [] Density ;
        delete [] Concentration ;
        delete [] Temperature ;
        delete [] Ga   ;
        delete [] Vvec    ;
        delete [] BFvec  ;
        if(isF)
        {
            delete []Ff;
        }
    }


    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

    // Writing xmf file
    std::ostringstream oss;


    if (nz==1)
    {
        oss << "<?xml version=\"1.0\" ?>\n";
        oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        oss << "<Xdmf Version=\"2.0\">\n";
        oss << " <Domain>\n";
        oss << "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n";
        oss << "     <Topology TopologyType=\"2DCoRectMesh\" Dimensions=\"" << ny << " " << nx << "\"/>\n";
        oss << "     <Geometry GeometryType=\"ORIGIN_DXDY\">\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2\"> 0.0 0.0\n";
        oss << "       </DataItem>\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2\"> 1.0 1.0\n";
        oss << "       </DataItem>\n";
        oss << "     </Geometry>\n";
        for (size_t j=0;j<Nl;j++)
        {
        oss << "     <Attribute Name=\"Density_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << nx << " " << ny << " " << nz << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Density_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Concentration_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << nx << " " << ny << " " << nz << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Concentration_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Temperature_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << nx << " " << ny << " " << nz << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Temperature_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Velocity_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << nx << " " << ny << " " << nz << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Velocity_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"BForce_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << nx << " " << ny << " " << nz << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/BForce_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        }
        oss << "     <Attribute Name=\"Gamma\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << nx << " " << ny << " " << nz << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Gamma\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Overlap\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << nx << " " << ny << " " << nz << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Overlap\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "   </Grid>\n";
        oss << " </Domain>\n";
        oss << "</Xdmf>\n";
    }
    
    else
    {
        std::cout<<"To Be Continue";
    }
    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();

}

void Imprint1(LBM::Domain &fluid_dom, SCALAR::Domain &con_dom, SCALAR::Domain &temp_dom)
{
    size_t Nproc = fluid_dom.Nproc;

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif    
    for(size_t ix=0; ix<fluid_dom.nx; ++ix)
    for(size_t iy=0; iy<fluid_dom.ny; ++iy)
    for(size_t iz=0; iz<fluid_dom.nz; ++iz)
    {
        double dT = (temp_dom.Rho[ix][iy][iz]-0.0);
        double dC = (con_dom.Rho[ix][iy][iz]-0.0);
        if(dT<0) dT = 0.0;
        if(dC<0) dC = 0.0;
        fluid_dom.BForce[ix][iy][iz] = (-temp_dom.beta*dT+con_dom.beta*dC)*fluid_dom.Rho[ix][iy][iz]*fluid_dom.g;
    }
    
}

void Imprint2(LBM::Domain &fluid_dom, SCALAR::Domain &con_dom, SCALAR::Domain &temp_dom)
{
    size_t Nproc = fluid_dom.Nproc;

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif    
    for(size_t ix=0; ix<fluid_dom.nx; ++ix)
    for(size_t iy=0; iy<fluid_dom.ny; ++iy)
    for(size_t iz=0; iz<fluid_dom.nz; ++iz)
    {
        temp_dom.Vel[ix][iy][iz] = fluid_dom.Vel[ix][iy][iz];
        con_dom.Vel[ix][iy][iz] = fluid_dom.Vel[ix][iy][iz];
    }
    
}


void Solve(LBM::Domain &fluid_dom, SCALAR::Domain &con_dom, SCALAR::Domain &temp_dom, double Tf, double dtOut,ptDFun_t ptSetup1,ptDFun_t ptSetup2, ptDFun_t ptReport, char const* TheFileKey, bool RenderVideo)
{
    String FileKey;
    FileKey.Printf("%s",TheFileKey);
    bool Finished = false;
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1   , TERM_RST);
    std::cout<<"Tau of fluid "<<fluid_dom.Tau<<std::endl;
    std::cout<<"Tau of temp "<<temp_dom.Tau<<std::endl;
    std::cout<<"Tau of con "<<con_dom.Tau<<std::endl;
    isF = false;
    

    double Time = 0.0;
    size_t idx_out = 0;
    double tout = Time;
    double dt = fluid_dom.dt;
    while (Time < Tf)
    {
        if (Time >= tout)
        {
            if (TheFileKey!=NULL)
            {
                String fn;
				String fn1;
                fn.Printf    ("%s_%04d", TheFileKey, idx_out);
                
                if ( RenderVideo) 
                {
                    #ifdef USE_HDF5
                    // fluid_dom.WriteXDMF  (fn.CStr());
                    WriteXDMF  (fn.CStr(),fluid_dom, con_dom, temp_dom);
                    #endif
                }
            }
    
            tout += dtOut;
            idx_out++;
            printf(" Accomplished                                  = %3.2f %% \n", Time);
            if (ptReport!=NULL) (*ptReport) (fluid_dom, con_dom, temp_dom, UserData);            
        }
        
        
        Imprint1(fluid_dom, con_dom, temp_dom);
        
        fluid_dom.CollideMRT();
        fluid_dom.Stream();
        fluid_dom.BounceBack();
        if(ptSetup1 != NULL) (*ptSetup1) (fluid_dom, con_dom, temp_dom, UserData); 
        
        fluid_dom.CalcProps();
        
        

        Imprint2(fluid_dom, con_dom, temp_dom);        
                
        temp_dom.Collide();
        temp_dom.Stream();
        // temp_dom.BounceBack();
        con_dom.Collide();
        con_dom.Stream();
        // con_dom.BounceBack();
        if(ptSetup2 != NULL) (*ptSetup2) (fluid_dom, con_dom, temp_dom, UserData);        
        temp_dom.CalcProps();
        con_dom.CalcProps();

        Time += dt;
        fluid_dom.Time = Time;
    }
    Finished = true;
    if (ptReport!=NULL) (*ptReport) (fluid_dom, con_dom, temp_dom, UserData);
    printf("%s  Final CPU time       = %s\n",TERM_CLR2, TERM_RST);

}

void Solve(LBM::Domain &fluid_dom, double Tf, double dtOut,ptDFun_t1 ptSetup1, ptDFun_t1 ptReport, char const* TheFileKey, bool RenderVideo)
{
    String FileKey;
    FileKey.Printf("%s",TheFileKey);
    bool Finished = false;
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1   , TERM_RST);
    std::cout<<"Tau of fluid "<<fluid_dom.Tau<<std::endl;


    

    double Time = 0.0;
    size_t idx_out = 0;
    double tout = Time;
    double dt = fluid_dom.dt;
    while (Time < Tf)
    {
        if (Time >= tout)
        {
            if (TheFileKey!=NULL)
            {
                String fn;
				String fn1;
                fn.Printf    ("%s_%04d", TheFileKey, idx_out);
                
                if ( RenderVideo) 
                {
                    #ifdef USE_HDF5
                    fluid_dom.WriteXDMF  (fn.CStr());
                    #endif
                }
            }
    
            tout += dtOut;
            idx_out++;
            printf(" Accomplished                                  = %3.2f %% \n", Time);
            if (ptReport!=NULL) (*ptReport) (fluid_dom, UserData);            
        }
        

        // Imprint1(fluid_dom, con_dom, temp_dom);
        
        fluid_dom.Collide();
        fluid_dom.Stream();
        // fluid_dom.BounceBack();
        
        if(ptSetup1 != NULL) (*ptSetup1) (fluid_dom, UserData); 
        fluid_dom.CalcProps();
        
        // Imprint2(fluid_dom, con_dom, temp_dom);        
        
        // temp_dom.Collide();
        // temp_dom.Stream();
        
        // con_dom.Collide();
        // con_dom.Stream();
        
        // if(ptSetup2 != NULL) (*ptSetup2) (fluid_dom, con_dom, temp_dom, UserData);        

        Time += dt;
    }
    Finished = true;
    // if (ptReport!=NULL) (*ptReport) (fluid_dom, con_dom, temp_dom, UserData);
    printf("%s  Final CPU time       = %s\n",TERM_CLR2, TERM_RST);

}




#endif
