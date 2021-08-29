#include<iostream>
#include<vector>
#include<string>
#include<cmath>
#include<fstream>
#include <iostream>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/SparseQR>
#include <H5Cpp.h>


typedef Eigen::Triplet<double> T;

using namespace std;
using namespace H5;

void initialize(vector<vector<double>> &u, vector<vector<double>> &v, vector<vector<double>> &p, vector<vector<double>> &u_tilde, vector<vector<double>> &v_tilde)
{
    for(int i=0; i<u.size(); i++){
        for(int j=0; j<v.size(); j++){
            u[i][j]=0.0;
            v[i][j]=0.0;
            u_tilde[i][j]=0.0;
            v_tilde[i][j]=0.0;
            p[i][j]=0.0;
        }
    }
}


void boundary_initialize(vector<vector<double>> &u, vector<vector<double>> &v, vector<vector<double>> &u_tilde, vector<vector<double>> &v_tilde)
{
    for(int i=0;i<u.size(); i++){
        u[i][0]=0.0; u[i][u.size()-1]=0.0; u[u.size()-1][i]=0.0;
        v[i][0]=0.0; v[0][i]=0.0; v[u.size()-1][i]=0.0; v[i][u.size()-1]=0.0;
        u_tilde[i][0]=0.0; u_tilde[i][u.size()-1]=0.0; u_tilde[u.size()-1][i]=0.0;
        v_tilde[i][0]=0.0; v_tilde[0][i]=0.0; v_tilde[u.size()-1][i]=0.0; v_tilde[i][u.size()-1]=0.0;

    }
    for(int i=0;i<u.size(); i++){
        u[0][i]=10.0;
        u_tilde[0][i]=10.0;
    }
}

void prediction_velocity_u(vector<vector<double>> u, vector<vector<double>> v, vector<vector<double>> &u_tilde, double Re, double dx, double dy, double dt){
    for(int i=1; i<u_tilde.size()-1; i++){
        for(int j=1; j<u_tilde[i].size()-1; j++){
            //calc_x_diff
            double advection_u=((u[i][j+1]-u[i][j-1])/(2.0*dx))*u[i][j];
            //calc_y_diff
            double advection_v=((u[i+1][j]-u[i-1][j])/(2.0*dy))*v[i][j];
            //calc_x_diff_diffusion
            double diffusion_x=(u[i][j+1]+u[i][j-1]-2.0*u[i][j])/(dx*dx);
            //calc_y_diff_diffusion
            double diffusion_y=(u[i+1][j]+u[i-1][j]-2.0*u[i][j])/(dy*dy);

            double right_hand=(diffusion_x+diffusion_y)*(1.0/Re);
            double left_hand=advection_u+advection_v;

            u_tilde[i][j]=u[i][j]+dt*(right_hand-left_hand);
        }
    }
}

void prediction_velocity_v(vector<vector<double>> u, vector<vector<double>> v, vector<vector<double>> &v_tilde, double Re, double dx, double dy, double dt){
    for(int i=1; i<v.size()-1; i++){
        for(int j=1; j<v[i].size()-1; j++){
            //calc_x_diff
            double advection_u=((v[i][j+1]-v[i][j-1])/(2.0*dx))*u[i][j];
            //calc_y_diff
            double advection_v=((v[i+1][j]-v[i-1][j])/(2.0*dy))*v[i][j];
            //calc_x_diff_diffusion
            double diffusion_x=(v[i][j+1]+v[i][j-1]-2.0*v[i][j])/(dx*dx);
            //calc_y_diff_diffusion
            double diffusion_y=(v[i+1][j]+v[i-1][j]-2.0*v[i][j])/(dy*dy);

            double right_hand=(diffusion_x+diffusion_y)*(1.0/Re);
            double left_hand=advection_u+advection_v;

            v_tilde[i][j]=v[i][j]+dt*(right_hand-left_hand);
        }
    }
}

void pressure(vector<vector<double>> &p, vector<vector<double>> u_tilde, vector<vector<double>> v_tilde, double y_point_num, double x_point_num, double dx, double dy, double dt)
{
    vector<vector<double>> A(y_point_num*x_point_num, vector<double>(x_point_num*y_point_num));
    vector<double> b(y_point_num*x_point_num);
    for(int i=0; i<A.size(); i++){
        for(int j=0; j<A[i].size(); j++){
            A[i][j]=0.0;
        }
    }
    int row=0;
    for(int i=0; i<y_point_num; i++){
        for(int j=0; j<x_point_num; j++){
            if(i==0){
                A[row][y_point_num*i+j]=1.0;
                A[row][y_point_num*(i+1)+j]=-1.0;
                row++;
                continue;
            } 
            if(i==y_point_num-1){
                A[row][y_point_num*i+j]=1.0;
                A[row][y_point_num*(i-1)+j]=-1.0;
                row++;
                continue;
            }
            if(j==0){
                A[row][y_point_num*i+j]=1.0;
                A[row][y_point_num*i+j+1]=-1.0;
                row++;
                continue;
            } 
            if(j==y_point_num-1){
                A[row][y_point_num*i+j]=1.0;
                A[row][y_point_num*i+j-1]=-1.0;
                row++;
                continue;
            }
            A[row][y_point_num*i+j+1]+=(1.0/(dx*dx));
            A[row][y_point_num*i+j-1]+=(1.0/(dx*dx));
            A[row][y_point_num*i+j]+=(-2.0/(dx*dx));
            A[row][y_point_num*(i+1)+j]+=(1.0/(dy*dy));
            A[row][y_point_num*(i-1)+j]+=(1.0/(dy*dy));
            A[row][y_point_num*i+j]+=(-2.0/(dy*dy));
            b[row]+=((u_tilde[i][j+1]-u_tilde[i][j-1])/(2.0*dx)+(v_tilde[i+1][j]-v_tilde[i-1][j])/(2.0*dy))/dt;
            row++;
        }
    }

    std::vector<T> tripletVec;
    for(int i=0;i<A.size(); i++){
        for(int j=0; j<A[i].size(); j++){
            if(fabs(A[i][j])>0.00001) tripletVec.push_back(T(i,j,A[i][j]));
        }
    }
    Eigen::SparseMatrix<double> M(x_point_num*y_point_num,x_point_num*y_point_num);
    M.setFromTriplets(tripletVec.begin(), tripletVec.end());

    Eigen::SparseMatrix<double> b_t(x_point_num*y_point_num,1);
    std::vector<T> tripletVecb;
    for(int i=0; i<b.size(); i++){
        if(fabs(b[i])>0.00001) tripletVecb.push_back( T(i,0,b[i]) );
    }
    b_t.setFromTriplets(tripletVecb.begin(), tripletVecb.end());

    Eigen::VectorXd x; 
    Eigen::SparseQR< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;  
    solver.compute(M);

    if( solver.info() != Eigen::Success ) {
        std::cerr << "decomposition failed" << std::endl;
    }
    x = solver.solve(b_t);
    if( solver.info() != Eigen::Success ) {
        std::cerr << "solving failed" << std::endl;
    }
    
    for(int i=0; i<y_point_num; i++){
        for(int j=0; j<x_point_num; j++){
            p[i][j]=x(i*x_point_num+j);
        }
    }

}

void update_u(vector<vector<double>> &u, vector<vector<double>> u_tilde, vector<vector<double>> p, double dt, double dx, double x_point_num, double y_point_num)
{
    for(int i=1; i<y_point_num-1; i++){
        for(int j=1; j<x_point_num-1; j++){
            u[i][j]=u_tilde[i][j]-((p[i][j+1]-p[i][j-1])/(2.0*dx))*dt;
        }
    }
}

void update_v(vector<vector<double>> &v, vector<vector<double>> v_tilde, vector<vector<double>> p, double dt, double dy, double x_point_num, double y_point_num)
{
    for(int i=1; i<y_point_num-1; i++){
        for(int j=1; j<x_point_num-1; j++){
            v[i][j]=v_tilde[i][j]-((p[i+1][j]-p[i-1][j])/(2.0*dy))*dt;
        }
    }
}

void export_vtk(string filename, string filename2, vector<vector<double>> u, vector<vector<double>> v, vector<vector<double>> p, double dx, double dy, double x_point_num, double y_point_num)
{
    ofstream ofs(filename);
    ofstream of(filename2);
    int point_count=0;
    for(int i=0; i<x_point_num; i++){
        for(int j=0; j<y_point_num; j++){
            if((fabs(u[j][i])!=0.0 && fabs(v[j][i])!=0.0)){
                point_count++;
            }
        }
    }

    of << "# vtk DataFile Version 3.0" << endl;
    of << "vtk output" << endl;
    of << "ASCII" << endl;
    of << "DATASET UNSTRUCTURED_GRID" << endl;
    of << "POINTS " << point_count << " double" << endl;
    for(int i=0; i<y_point_num; i++){
	    for(int j=0; j<x_point_num; j++){
            if((fabs(u[j][i])!=0.0 && fabs(v[j][i])!=0.0)){
                of << j*dx << " " << i*dy << " " << 0 << endl;
            }
	    }
	}
    of << "POINT_DATA " << point_count << endl;
    of << "VECTORS velocity[m/s] double" << endl; 
	for(int i=0; i<y_point_num; i++){
	    for(int j=0; j<x_point_num; j++){
            if((fabs(u[j][i])!=0.0 && fabs(v[j][i])!=0.0)){
                of << u[i][j] << " " << v[i][j] << " " << 0 << endl;
            }
	    }
	}

    ofs << "# vtk DataFile Version 3.0" << endl;
    ofs << "vtk output" << endl;
    ofs << "ASCII" << endl;
    ofs << "DATASET STRUCTURED_POINTS" << endl;
    ofs << "DIMENSIONS " << x_point_num << " " << y_point_num << " " << 1 << endl;
    ofs << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
    ofs << "SPACING " << dx << " " << dy << " " << 1.0 << endl;
    ofs << "POINT_DATA " << x_point_num*y_point_num << endl;
    ofs <<"SCALARS " << "pressure " << "double" << endl;
    ofs <<"LOOKUP_TABLE default" << endl;
    for(int i=0; i<y_point_num; i++){
        for(int j=0; j<x_point_num; j++){
            ofs << p[i][j] << endl;
        }
    }
    ofs <<"SCALARS " << "velocity_magnitude " << "double" << endl;
    ofs <<"LOOKUP_TABLE default" << endl;
    for(int i=0; i<y_point_num; i++){
        for(int j=0; j<x_point_num; j++){
            ofs << sqrt(u[i][j]*u[i][j]+v[i][j]*v[i][j]) << endl;
        }
    }
}

void write_HDF5_main(H5File file, string Gr, string dataName, string name, int num_all_points, double x_tmp[])
{
    Group group = file.openGroup(Gr.c_str());
    dataName = Gr + "/" + name;
    H5std_string DATASET_NAME(dataName.c_str());
    hsize_t dim[1] = {num_all_points}; // dataset dimensions
    H5::DataSpace dataspace(1, dim);
    H5::IntType datatype(H5::PredType::NATIVE_DOUBLE);
    datatype.setOrder(H5T_ORDER_LE);
    H5::DataSet dataset = file.createDataSet(DATASET_NAME, datatype, dataspace);
    dataset.write(x_tmp, H5::PredType::NATIVE_DOUBLE);
}

void export_initial_HDF5_1D(string H5fileName, string name, int num_all_points, vector<double> x, bool create)
{
    double *x_tmp = x.data();
    std::string dataName,Gr = "/"+to_string(0);
    H5std_string FILE_NAME(H5fileName.c_str());    
    if(create){
        H5File file(FILE_NAME, H5F_ACC_TRUNC);
        file.createGroup(Gr.c_str());
        write_HDF5_main(file, Gr, dataName, name, num_all_points, x_tmp);
    }
    else{
        H5File file(FILE_NAME, H5F_ACC_RDWR);
        write_HDF5_main(file, Gr, dataName, name, num_all_points, x_tmp);
    }
}

void export_initial_HDF5_2D(string H5fileName, string name, int num_all_points, vector<vector<double>> x, bool create)
{
    double x_tmp[num_all_points];
    for(int i=0; i<x.size(); i++){
        for(int j=0; j<x[i].size(); j++){
            x_tmp[i*x.size()+j]=x[i][j];
        }
    }
    std::string dataName,Gr = "/"+to_string(0);
    H5std_string FILE_NAME(H5fileName.c_str());    
    if(create){
        H5File file(FILE_NAME, H5F_ACC_TRUNC);
        file.createGroup(Gr.c_str());
        write_HDF5_main(file, Gr, dataName, name, num_all_points, x_tmp);
    }
    else{
        H5File file(FILE_NAME, H5F_ACC_RDWR);
        write_HDF5_main(file, Gr, dataName, name, num_all_points, x_tmp);
    }   
}

void export_HDF5_1D(string H5fileName, string name, int num_all_points, vector<double> x, int loop, bool create)
{
    double *x_tmp = x.data();
    H5std_string FILE_NAME(H5fileName.c_str());    
    H5File file(FILE_NAME, H5F_ACC_RDWR);
    std::string dataName,Gr = "/"+to_string(loop);
    if(create){
        file.createGroup(Gr.c_str());
    }
    write_HDF5_main(file, Gr, dataName, name, num_all_points, x_tmp);
}

void export_HDF5_2D(string H5fileName, string name, int num_all_points, vector<vector<double>> x, int loop, bool create)
{
    double x_tmp[num_all_points];
    for(int i=0; i<x.size(); i++){
        for(int j=0; j<x[i].size(); j++){
            x_tmp[i*x.size()+j]=x[i][j];
        }
    }
    H5std_string FILE_NAME(H5fileName.c_str());    
    H5File file(FILE_NAME, H5F_ACC_RDWR);
    std::string dataName,Gr = "/"+to_string(loop);
    if(create){
        file.createGroup(Gr.c_str());
    }
    write_HDF5_main(file, Gr, dataName, name, num_all_points, x_tmp);
}

int main()
{
    int x_point_num=11, y_point_num=11;
    int num_all_points=x_point_num*y_point_num;
    double dx=0.1; double dy=0.1;
    double dt=0.001;
    int time_step=10000;
    double Re=10;
    vector<vector<double>> u(y_point_num, vector<double>(x_point_num));
    vector<vector<double>> v(y_point_num, vector<double>(x_point_num));
    vector<vector<double>> u_tilde(y_point_num, vector<double>(x_point_num));
    vector<vector<double>> v_tilde(y_point_num, vector<double>(x_point_num));
    vector<vector<double>> p(y_point_num, vector<double>(x_point_num));
    vector<double> x(num_all_points);
    vector<double> y(num_all_points);
    bool True=1;
    bool False=0;
    string name;
    for(int i=0; i<y_point_num; i++){
	    for(int j=0; j<x_point_num; j++){
            x[i*x_point_num+j]=j*dx;
            y[i*x_point_num+j]=i*dy;
	    }
	}
    initialize(u,v,p,u_tilde,v_tilde);
    boundary_initialize(u, v, u_tilde, v_tilde);
    string H5fileName="data.h5";
    name="x";
    export_initial_HDF5_1D(H5fileName,name,num_all_points,x,True);
    name="y";
    export_initial_HDF5_1D(H5fileName,name,num_all_points,y,False);
    name="u";
    export_initial_HDF5_2D(H5fileName,name,num_all_points,u,False);
    name="v";
    export_initial_HDF5_2D(H5fileName,name,num_all_points,v,False);
    name="p";
    export_initial_HDF5_2D(H5fileName,name,num_all_points,p,False);
    for(int i=1; i<=10; i++){
        boundary_initialize(u, v, u_tilde, v_tilde);
        prediction_velocity_u(u, v, u_tilde, Re, dx, dy, dt);
        prediction_velocity_v(u, v, v_tilde, Re, dx, dy, dt);
        pressure(p, u_tilde, v_tilde, y_point_num, x_point_num, dx, dy, dt);
        update_u(u, u_tilde, p, dt, dx, x_point_num, y_point_num);
        update_v(v, v_tilde, p, dt, dy, x_point_num, y_point_num);
        cout << "iter" << " " << ":" << " " << i << endl;
        name="x";
        export_HDF5_1D(H5fileName,name,num_all_points,x,i,True);
        name="y";
        export_HDF5_1D(H5fileName,name,num_all_points,y,i,False);
        name="u";
        export_HDF5_2D(H5fileName,name,num_all_points,u,i,False);
        name="v";
        export_HDF5_2D(H5fileName,name,num_all_points,v,i,False);
        name="p";
        export_HDF5_2D(H5fileName,name,num_all_points,p,i,False);
    }
}