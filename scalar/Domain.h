#ifndef MECHSYS_SCALAR_DOMAIN_H
#define MECHSYS_SCALAR_DOMAIN_H


#include <mechsys/dem/domain.h>
#include <mechsys/lbm/Lattice.h>
#include <mechsys/lbm/Interacton.h>
#include <mechsys/dem/domain.h>
#include <mechsys/lbm/Dem.h>

// STD
#include <map>
#include <vector>
#include <utility>
#include <set>
#include <ctime>
#include <ratio>
#include <chrono>
namespace SCALAR{

class Domain
{
public:
    static const double   WEIGHTSD2Q4   [ 4] ;
    static const Vec3_t   LVELOCD2Q4    [ 4] ;
    static const size_t   OPPOSITED2Q4  [ 4] ;
    Domain();

    double ****F;
    double ****Ftemp;
    double ****Omeis;
    double ***Rho;
    double Rho0;
    Vec3_t ***Vel;
    Vec3_t ***BForce;
    bool ***IsSolid;
    bool ***IsBoundary;
    double ***Gamma;
    size_t Ncells;
    size_t Nneigh;
    double Sc;
    const double *W;
    const Vec3_t *C;
    const size_t *Op;
    size_t nx ;
    size_t ny ;
    size_t nz;
    size_t Nproc;
    double nu;
    double beta;
    double Tau;
    double dx;
    double dt;
    double Cs;
    double Alpha;
    double ratio; //Ll/L
    double ratiot; //tl/t
    double ratiom; //ml/m
    double ratiov;
    double Time;
    Vec3_t g;
    Array <LBM::Disk     *>                            Disks;         ///< Array of Disks for 2D calculation

	#ifdef USE_OMP
    omp_lock_t      lck;             ///< to protect variables in multithreading
	#endif
    //function
    double Feq(size_t k, double rho, Vec3_t &vel);
    void WriteXDMF(char const * FileKey);
    void Collide();
    void Stream();
    void BounceBack();
    void CalcProps();
    void AllocateDomainMem(size_t nx, size_t ny, size_t nz, double nu0, double dx0, double dt0);
};
inline Domain::Domain ()
{
    Alpha = 1.0;
    Nproc = 1;
	omp_init_lock(&lck);
}
inline void Domain::AllocateDomainMem(size_t nx0, size_t ny0, size_t nz0,double nu0, double dx0, double dt0)
{
    nx = nx0;
    ny = ny0;
    nz = nz0;
    dx = dx0;
    dt = dt0;
    Cs = dx/dt;
    Nneigh = 4;
	Ncells = nx*ny*nz;
	W      = WEIGHTSD2Q4;
    C      = LVELOCD2Q4;
    Op     = OPPOSITED2Q4;
	Sc = 0.17;
    Alpha = 1.0;
    beta = 0.0;
    Time = 0.0;
	F = new double***[nx];
	Ftemp = new double***[nx];
	Rho = new double **[nx];
	Vel = new Vec3_t **[nx];
    BForce = new Vec3_t **[nx];
    IsSolid = new bool **[nx];
    IsBoundary = new bool **[nx];
    Gamma = new double **[nx];
	nu = nu0;
    // Tau = 3.0*nu*dt/(dx*dx)+0.5;
    Tau = 2.0*nu+0.5;
	for (size_t i=0; i<nx; i++)
	{
		F[i] = new double **[ny];
		Ftemp[i] = new double **[ny];
        Rho[i] = new double *[ny];
		Vel[i] = new Vec3_t *[ny];
        BForce[i] = new Vec3_t *[ny];
        Gamma[i] = new double *[ny];
		IsSolid[i] = new bool *[ny];
        IsBoundary[i] = new bool *[ny];
        for (size_t j=0; j<ny; j++)
		{
			F[i][j] = new double*[nz];
			Ftemp[i][j] = new double*[nz];
            Rho[i][j] = new double [nz];
			Vel[i][j] = new Vec3_t [nz];
            BForce[i][j] = new Vec3_t [nz];
			Gamma[i][j] = new double [nz];
            IsSolid[i][j] = new bool[nz];
            IsBoundary[i][j] = new bool[nz];
            for (size_t k=0; k<nz; k++)
			{
				F[i][j][k] = new double[Nneigh];
				Ftemp[i][j][k] = new double[Nneigh];
                IsSolid[i][j][k] = false;
                IsBoundary[i][j][k] = false;
                BForce[i][j][k] = 0.0, 0.0, 0.0;
				for (size_t nn=0; nn<Nneigh; nn++)
				{
					F[i][j][k][nn] = 0.0;
					Ftemp[i][j][k][nn] = 0.0;
                }
			}
		}
	}
	
}

inline void Domain::WriteXDMF(char const * FileKey)
{
    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    size_t  Nx = nx;
    size_t  Ny = ny;
    size_t  Nz = nz;
    size_t Step = 1;
    size_t Nl = 1;
    std::cout<<Nx<<" "<<Ny<<" "<<Nz<<std::endl;
    for (size_t j=0;j<Nl;j++)
    {
        // Creating data sets
        double * Density   = new double[  Nx*Ny*Nz];
        double * Ga     = new double[  Nx*Ny*Nz];
        double * Vvec      = new double[3*Nx*Ny*Nz];
        double * BFvec      = new double[3*Nx*Ny*Nz];
        double * Ff         = new double[Nneigh*Nx*Ny*Nz];
        size_t i=0;
        for (size_t m=0;m<nz;m+=Step)
        for (size_t l=0;l<ny;l+=Step)
        for (size_t n=0;n<nx;n+=Step)
        {
            double rho    = 0.0;
            double gamma  = 0.0;
            Vec3_t vel    = Vec3_t(0.0,0.0,0.0);
            Vec3_t BF     = Vec3_t(0.0,0.0,0.0);
            double temp = 0.0;
            for (size_t ni=0;ni<Step;ni++)
            for (size_t li=0;li<Step;li++)
            for (size_t mi=0;mi<Step;mi++)
            {
                rho    += Rho    [n+ni][l+li][m+mi];
                temp    = IsSolid[n+ni][l+li][m+mi] ? 1.0: 0.0;
                gamma  += std::max(Gamma[n+ni][l+li][m+mi],temp);
                vel    += Vel    [n+ni][l+li][m+mi];
                BF    += BForce    [n+ni][l+li][m+mi];
            }
            rho  /= Step*Step*Step;
            gamma/= Step*Step*Step;
            vel  /= Step*Step*Step;
            BF   /= Step*Step*Step;
            Ga   [i]  = (double) gamma;
            Density [i]  = (double) rho;            
            Vvec[3*i  ]  = (double) vel(0);//*(1.0-Ga[i]);
            Vvec[3*i+1]  = (double) vel(1);//*(1.0-Ga[i]);
            Vvec[3*i+2]  = (double) vel(2);//*(1.0-Ga[i]);
            BFvec[3*i ]   = (double) BF(0);
            BFvec[3*i+1]   = (double) BF(1);
            BFvec[3*i+2]   = (double) BF(2);
            for (size_t k=0; k<Nneigh; k++)
            {
                Ff[Nneigh*i + k] = (double) F[n][l][m][k];
            }
            i++;
        }
        
        //Writing data to h5 file
        hsize_t dims[1];
        dims[0] = Nx*Ny*Nz;
        String dsname;
        dsname.Printf("Density_%d",j);
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Density );
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
        dims[0] = Nneigh*Nx*Ny*Nz;
        dsname.Printf("F_%d",j);
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Ff    );
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
        delete [] Ga   ;
        delete [] Vvec    ;
        delete [] BFvec  ;
        delete [] Ff;
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

inline double Domain::Feq(size_t k, double rho, Vec3_t &vel)
{
    // double VdotV = dot(vel,vel);   
    double VdotC = dot(vel,C[k]);
    // double Feq = W[k]*rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
    double Feq = W[k]*rho*(1.0+2.0*VdotC/(Cs*Cs));
    return Feq;
}

inline void Domain::Collide()
{
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
	for(size_t ix=0; ix<nx;ix++)
	for(size_t iy=0; iy<ny;iy++)
	for(size_t iz=0; iz<nz;iz++)
	{
		if(!IsSolid[ix][iy][iz])
		{
    			double ws = 1.0/Tau;
                double rho = Rho[ix][iy][iz];
                Vec3_t vel = Vel[ix][iy][iz];
                // double Bn = (Gamma[ix][iy][iz]*(tau-0.5))/((1.0-Gamma[ix][iy][iz])+(tau-0.5));
                for (size_t k=0; k<Nneigh; k++)
                {
                    // std::cout<<ix<<" "<<iy<<" "<<iz<<" "<<k<<std::endl;
                    // std::cout<<BForce[0][1][0]<<std::endl;
                    // double ForceTerm = dt*3.0*W[k]*dot(BForce[ix][iy][iz],C[k])/(Cs*Cs);                    
                	// Ftemp[ix][iy][iz][k] = F[ix][iy][iz][k] -(NonEq[k]/tau)+ForceTerm;
                    Ftemp[ix][iy][iz][k] = (1.0-ws)*F[ix][iy][iz][k]+ws*Feq(k,rho,vel);
                    // if(Ftemp[ix][iy][iz][k]<0.0) Ftemp[ix][iy][iz][k] = 0.0;
                     
                }
                
    	}else{
    		for(size_t k=0; k<Nneigh; k++)
    		{
    			Ftemp[ix][iy][iz][k] = F[ix][iy][iz][Op[k]];
    		}
    	}
    }
    double **** tmp = F;
    F = Ftemp;
    Ftemp = tmp;
}

inline void Domain::Stream()
{
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
	for(size_t ix=0; ix<nx; ix++)
	for(size_t iy=0; iy<ny; iy++)
	for(size_t iz=0; iz<nz; iz++)
	for(size_t k=0; k<Nneigh; k++)
	{
		size_t nix = (size_t)((int)ix + (int)C[k](0) + (int)nx)%nx;
        size_t niy = (size_t)((int)iy + (int)C[k](1) + (int)ny)%ny;
        size_t niz = (size_t)((int)iz + (int)C[k](2) + (int)nz)%nz;
        Ftemp[nix][niy][niz][k] = F[ix][iy][iz][k];
	}

    double **** tmp = F;
    F = Ftemp;
    Ftemp = tmp;
}

inline void Domain::BounceBack()
{
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
	for(size_t ix=0; ix<nx; ix++)
	for(size_t iy=0; iy<ny; iy++)
	for(size_t iz=0; iz<nz; iz++)
    {
        if(IsBoundary[ix][iy][iz])
        {
            double Temp[Nneigh];
            for(size_t k=0; k<Nneigh; k++)
            {
                Temp[k] = F[ix][iy][iz][k];

            }
            for(size_t k=0; k<Nneigh; k++)
    	    {
				size_t nix = (size_t)((int)ix + (int)C[Op[k]](0) + (int)nx)%nx;
                size_t niy = (size_t)((int)iy + (int)C[Op[k]](1) + (int)ny)%ny;
                size_t niz = (size_t)((int)iz + (int)C[Op[k]](2) + (int)nz)%nz;
                if(IsSolid[nix][niy][niz])
                {
    				F[ix][iy][iz][k] = Temp[Op[k]];                            
                }          
            }
        }
	    
    }
}

// inline void Domain::Adiabatic()
// {
//     #ifdef USE_OMP
//     #pragma omp parallel for schedule(static) num_threads(Nproc)
//     #endif
// 	for(size_t ix=0; ix<nx; ix++)
// 	for(size_t iy=0; iy<ny; iy++)
// 	for(size_t iz=0; iz<nz; iz++)
//     {
//         if(IsSolid[ix][iy][iz])
//         {
//             for(size_t k=0; k<Nneigh; k++)
//     	    {
// 				size_t nix = (size_t)((int)ix + (int)C[k](0) + (int)nx)%nx;
//                 size_t niy = (size_t)((int)iy + (int)C[k](1) + (int)ny)%ny;
//                 size_t niz = (size_t)((int)iz + (int)C[k](2) + (int)nz)%nz;
//                 if(!IsSolid[nix][niy][niz])
//                 {
//     				F[ix][iy][iz][k] = F[nix][niy][niz][k];                      
//                 }          
//             }
//         }
	    
//     }
// }

inline void Domain::CalcProps()
{    
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
	for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        Vel   [ix][iy][iz] = Vec3_t(0.0,0.0,0.0);
        Rho   [ix][iy][iz] = 0.0;
        for (size_t k=0;k<Nneigh;k++)
        {
            Rho[ix][iy][iz] +=  F[ix][iy][iz][k];
            Vel[ix][iy][iz] +=  F[ix][iy][iz][k]*C[k];
        }
        Vel[ix][iy][iz] *= Cs/Rho[ix][iy][iz];
    }
}


/*inline void Domain::LoadResults (char const * FileKey)
{
    String fn(FileKey);
    fn.append(".hdf5");
    if (!Util::FileExists(fn)) throw new Fatal("File <%s> not found",fn.CStr());
    printf("\n%s--- Loading file %s --------------------------------------------%s\n",TERM_CLR1,fn.CStr(),TERM_RST);
    hid_t file_id;
    file_id = H5Fopen(fn.CStr(), H5F_ACC_RDONLY, H5P_DEFAULT);
    int data[1];
    H5LTread_dataset_int(file_id,"NP",data);
    size_t NP = data[0];
}
*/
// const double   Domain::WEIGHTSD2Q9   [ 9] = { 4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};
// const Vec3_t   Domain::LVELOCD2Q9    [ 9] = { {0,0,0}, {1,0,0}, {0,1,0}, {-1,0,0}, {0,-1,0}, {1,1,0}, {-1,1,0}, {-1,-1,0}, {1,-1,0} };
// const size_t   Domain::OPPOSITED2Q9  [ 9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };
// const double   Domain::MD2Q9 [ 9][ 9] =  { { 1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0},
// {-4.0, -1.0, -1.0, -1.0, -1.0,  2.0,  2.0,  2.0,  2.0},
// { 4.0, -2.0, -2.0, -2.0, -2.0,  1.0,  1.0,  1.0,  1.0},
// { 0.0,  1.0,  0.0, -1.0,  0.0,  1.0, -1.0, -1.0,  1.0},
// { 0.0, -2.0,  0.0,  2.0,  0.0,  1.0, -1.0, -1.0,  1.0},
// { 0.0,  0.0,  1.0,  0.0, -1.0,  1.0,  1.0, -1.0, -1.0},
// { 0.0,  0.0, -2.0,  0.0,  2.0,  1.0,  1.0, -1.0, -1.0},
// { 0.0,  1.0, -1.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0},
// { 0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0}};
const double Domain::WEIGHTSD2Q4 [4] = {0.25,0.25,0.25,0.25};
const Vec3_t Domain::LVELOCD2Q4  [4] = {{1,0,0},{0,1,0},{-1,0,0},{0,-1,0}};
const size_t Domain::OPPOSITED2Q4[4] = {2,3,0,1};

}; // namespace LBM

#endif // MECHSYS_LBM_DOMAIN_H
