//#include<iostream>
//#include<vector>
//#include<fstream>
//#include<string>
//#include<sstream>
//#include<ctime>       
//#include<cstdlib>
//
//using namespace std;
//
//void export_particle(string filename, vector<double> particle_x, vector<double> particle_y)
//{
//    ofstream ofs(filename);
//    ofs << "# vtk DataFile Version 3.0" << endl;
//    ofs << "Data.vtk" << endl;
//    ofs << "ASCII" << endl;
//    ofs << "DATASET UNSTRUCTURED_GRID" << endl;
//
//    ofs << "POINTS " << particle_x.size() << " float" << endl; 
//    for(int i=0; i<particle_x.size(); i++){
//        ofs << particle_x[i] << " " << particle_y[i] << " " << 0.0 << endl;
//    }
//    ofs << "CELL_TYPES " << particle_x.size() << endl;
//    for(int i=0; i<particle_x.size(); i++){
//        ofs << 1 << endl;
//    }
//    ofs << "POINT_DATA " << particle_x.size() << endl;
//    ofs << "SCALARS radius float" << endl;
//    ofs << "LOOKUP_TABLE default" << endl;
//    for(int i=0; i<particle_x.size(); i++){
//        ofs << 0.1 << endl;
//    }
//}
//
//int main()
//{
//    double dx = 0.1;
//    double dy = 0.1;
//    int y_point_num=31;
//    int x_point_num=31;
//    double dt=0.001;
//    int time_step=10000;
//    vector<double> x;
//    vector<double> y;
//    vector<vector<double>> u(time_step);
//    vector<vector<double>> v(time_step);
//    for(int i=0; i<time_step; i++){
//        string filename="velocity_vector/vel"+to_string(i)+".vtk";
//        ifstream ifs(filename);
//        if(!ifs){
//            cout << "file can't open" << endl;
//            exit(1);
//        }
//        string str;
//        for(int j=0; j<5; j++) getline(ifs, str);
//        for(int j=0; j<841; j++){
//            getline(ifs,str);
//            istringstream ss(str);
//            string s;
//            int count=0;
//            while (ss >> s) {
//                if(count==0) x.push_back(stod(s));
//                if(count==1) y.push_back(stod(s));
//                count++;
//            }
//        }
//        for(int j=0; j<2; j++){
//            getline(ifs, str);
//        }
//        for(int j=0; j<841; j++){
//            getline(ifs,str);
//            istringstream ss(str);
//            string s;
//            int count=0;
//            while (ss >> s) {
//                if(count==0) u[i].push_back(stod(s));
//                if(count==1) v[i].push_back(stod(s));
//                count++;
//            }
//        }
//        ifs.close();
//    }
//    double num_of_particle=10000;
//    vector<double> particle_x;
//    vector<double> particle_y;
//
//    std::srand( time(NULL) );
//    int count=0;
//    for(int i=0; i<num_of_particle; i++){
//        int tmp=rand()%300;
//        particle_x.push_back(double(double(tmp)/100.0));
//    }
//    for(int i=0; i<num_of_particle; i++){
//        int tmp=rand()%300;
//        particle_y.push_back(double(double(tmp)/100.0));
//    }
//
//    vector<double> x_t,y_t;
//    vector<vector<double>> u_t(time_step);
//    vector<vector<double>> v_t(time_step);
//
//    for(int i=0; i<y_point_num; i++){
//        for(int j=0; j<x_point_num; j++){
//            x_t.push_back(j*dx);
//            y_t.push_back(i*dy);
//        }
//    }
//
//    for(int i=0; i<time_step; i++){
//        for(int j=0; j<x_point_num; j++){
//            u_t[i].push_back(1.0);
//        }
//        for(int j=0; j<y_point_num-2; j++){
//            u_t[i].push_back(0.0);
//            for(int k=0; k<x_point_num-2; k++){
//                u_t[i].push_back(u[i][j*(y_point_num-2)+k]);
//            }
//            u_t[i].push_back(0.0);
//        }
//        for(int j=0; j<x_point_num; j++){
//            u_t[i].push_back(0.0);
//        }
//    }
//
//    for(int i=0; i<time_step; i++){
//        for(int j=0; j<x_point_num; j++){
//            v_t[i].push_back(0.0);
//        }
//        for(int j=0; j<y_point_num-2; j++){
//            v_t[i].push_back(0.0);
//            for(int k=0; k<x_point_num-2; k++){
//                v_t[i].push_back(v[i][j*(y_point_num-2)+k]);
//            }
//            v_t[i].push_back(0.0);
//        }
//        for(int j=0; j<x_point_num; j++){
//            v_t[i].push_back(0.0);
//        }
//    }
//
//    for(int i=0; i<time_step-1; i++){
//        for(int j=0;j<num_of_particle; j++){
//            int x_p=particle_x[j]/dx;
//            int y_p=particle_y[j]/dy;
//            int left_upper_node=y_p*x_point_num+x_p;
//            int right_upper_node=y_p*x_point_num+x_p+1;
//            int left_lower_node=(y_p+1)*x_point_num+x_p;
//            int right_lower_node=(y_p+1)*x_point_num+x_p+1;
//
//            double distance_bottom=y_t[left_lower_node]-particle_y[j];
//            double distance_left=particle_x[j]-x_t[left_lower_node];
//
//            double velocity_x=(1.0-distance_bottom)*(1.0-distance_left)*u_t[i][left_lower_node]
//            +(1.0-distance_bottom)*distance_left*u_t[i][right_lower_node]
//            +distance_bottom*(1.0-distance_left)*u_t[i][left_upper_node]
//            +distance_bottom*distance_left*u_t[i][right_upper_node];
//
//            double velocity_y=(1.0-distance_bottom)*(1.0-distance_left)*v_t[i][left_lower_node]
//            +(1.0-distance_bottom)*distance_left*v_t[i][right_lower_node]
//            +distance_bottom*(1.0-distance_left)*v_t[i][left_upper_node]
//            +distance_bottom*distance_left*v_t[i][right_upper_node];
//
//            //first next time step
//            double first_step_x=particle_x[j]+0.5*dt*velocity_x;
//            double first_step_y=particle_y[j]+0.5*dt*velocity_y;
//
//            x_p=first_step_x/dx;
//            y_p=first_step_y/dy;
//            left_upper_node=y_p*x_point_num+x_p;
//            right_upper_node=y_p*x_point_num+x_p+1;
//            left_lower_node=(y_p+1)*x_point_num+x_p;
//            right_lower_node=(y_p+1)*x_point_num+x_p+1;
//
//            distance_bottom=y_t[left_lower_node]-first_step_y;
//            distance_left=first_step_x-x_t[left_lower_node];
//            
//            double first_velocity_x=(1.0-distance_bottom)*(1.0-distance_left)*u_t[i][left_lower_node]
//            +(1.0-distance_bottom)*distance_left*u_t[i][right_lower_node]
//            +distance_bottom*(1.0-distance_left)*u_t[i][left_upper_node]
//            +distance_bottom*distance_left*u_t[i][right_upper_node];
//
//            double first_velocity_y=(1.0-distance_bottom)*(1.0-distance_left)*v_t[i][left_lower_node]
//            +(1.0-distance_bottom)*distance_left*v_t[i][right_lower_node]
//            +distance_bottom*(1.0-distance_left)*v_t[i][left_upper_node]
//            +distance_bottom*distance_left*v_t[i][right_upper_node];
//
//            //second next time step
//            double second_step_x=particle_x[j]+0.5*dt*first_velocity_x;
//            double second_step_y=particle_y[j]+0.5*dt*first_velocity_y;
//
//            x_p=second_step_x/dx;
//            y_p=second_step_y/dy;
//            left_upper_node=y_p*x_point_num+x_p;
//            right_upper_node=y_p*x_point_num+x_p+1;
//            left_lower_node=(y_p+1)*x_point_num+x_p;
//            right_lower_node=(y_p+1)*x_point_num+x_p+1;
//
//            distance_bottom=y_t[left_lower_node]-second_step_y;
//            distance_left=second_step_x-x_t[left_lower_node];
//            
//            double second_velocity_x=(1.0-distance_bottom)*(1.0-distance_left)*u_t[i][left_lower_node]
//            +(1.0-distance_bottom)*distance_left*u_t[i][right_lower_node]
//            +distance_bottom*(1.0-distance_left)*u_t[i][left_upper_node]
//            +distance_bottom*distance_left*u_t[i][right_upper_node];
//
//            double second_velocity_y=(1.0-distance_bottom)*(1.0-distance_left)*v_t[i][left_lower_node]
//            +(1.0-distance_bottom)*distance_left*v_t[i][right_lower_node]
//            +distance_bottom*(1.0-distance_left)*v_t[i][left_upper_node]
//            +distance_bottom*distance_left*v_t[i][right_upper_node];
//
//            //third time step
//            double third_step_x=particle_x[j]+0.5*dt*second_velocity_x;
//            double third_step_y=particle_y[j]+0.5*dt*second_velocity_y;
//
//            x_p=third_step_x/dx;
//            y_p=third_step_y/dy;
//            left_upper_node=y_p*x_point_num+x_p;
//            right_upper_node=y_p*x_point_num+x_p+1;
//            left_lower_node=(y_p+1)*x_point_num+x_p;
//            right_lower_node=(y_p+1)*x_point_num+x_p+1;
//
//            distance_bottom=y_t[left_lower_node]-third_step_y;
//            distance_left=third_step_x-x_t[left_lower_node];
//            
//            double third_velocity_x=(1.0-distance_bottom)*(1.0-distance_left)*u_t[i][left_lower_node]
//            +(1.0-distance_bottom)*distance_left*u_t[i][right_lower_node]
//            +distance_bottom*(1.0-distance_left)*u_t[i][left_upper_node]
//            +distance_bottom*distance_left*u_t[i][right_upper_node];
//
//            double third_velocity_y=(1.0-distance_bottom)*(1.0-distance_left)*v_t[i][left_lower_node]
//            +(1.0-distance_bottom)*distance_left*v_t[i][right_lower_node]
//            +distance_bottom*(1.0-distance_left)*v_t[i][left_upper_node]
//            +distance_bottom*distance_left*v_t[i][right_upper_node];
//
//            //update
//            particle_x[j]=particle_x[j]+(1.0/6.0)*dt*velocity_x+(1.0/3.0)*dt*first_velocity_x+
//            (1.0/3.0)*dt*second_velocity_x+(1.0/6.0)*dt*third_velocity_x;
//            
//            particle_y[j]=particle_y[j]+(1.0/6.0)*dt*velocity_y+(1.0/3.0)*dt*first_velocity_y+
//            (1.0/3.0)*dt*second_velocity_y+(1.0/6.0)*dt*third_velocity_y;
//
//            if(particle_x[j]<0.0) particle_x[j]=0.0;
//            if(particle_y[j]<0.0) particle_y[j]=0.0;
//            if(particle_x[j]>3.0) particle_x[j]=3.0;
//        }
//        string filename;
//        if(i%10==0){
//            filename="particle_file/particle"+to_string(i/10)+".vtk";
//            cout << filename << endl;
//            export_particle(filename, particle_x, particle_y);
//        }
//    }
//
//}