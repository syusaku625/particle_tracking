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

using namespace std;
using namespace H5;

//int main()
//{
//    string dataName;
//    string H5file="data.h5";
//    for(int i=0;i<=10; i++){
//        dataName = to_string(i)+"/u";
//        H5::H5File file(H5file,H5F_ACC_RDONLY);
//        H5::DataSet dataset = file.openDataSet(dataName);
//        H5::DataSpace dataspace = dataset.getSpace();  
//        hsize_t dims_out[1];
//        int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);
//        int nx = (unsigned long)(dims_out[0]);
//        double data[nx];
//        dataset.read(data, H5::PredType::NATIVE_DOUBLE);
//        vector<vector<double>> u_t(11, vector<double>(11));
//        for(int j=0; j<11; j++){
//            for(int k=0; k<11; k++){
//                u_t[j][k]=data[j*11+k];
//            }
//        }
//    }
//}