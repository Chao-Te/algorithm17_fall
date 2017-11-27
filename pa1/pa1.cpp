#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>

using namespace std;

int main(int argc, char** argv)
{
    // declare varaible
    string str;
    string in_file_name = argv[1];
    string out_file_name = argv[2];

    int n_control_pts = 0;
    int n_sample_pts = 0;
    int white_pos = 0;
    int upper_base = 0;
    int now_base = 0;

    float sample_space =0 ; // = 1/(n_sample_pts-1)
    float t0 = 0;           //t0 =t
    float t1 = 0;           // t1 = 1-t
    float upper_x = 0;
    float upper_y = 0;
    float upper_right_x = 0;
    float upper_right_y =0 ;

    float* dyn_table = NULL;
    float* out_pts = NULL;

    // read file
    ifstream infile;
    infile.open(in_file_name); //string to float stof string to int atoi
    getline(infile, str);// get first line
    n_control_pts = atoi(str.c_str());
    dyn_table = new float[n_control_pts*n_control_pts*2];
    for(int i = 0; i < n_control_pts; i++){//set hte first row of table
        getline(infile, str);
        white_pos = str.find(' ');
        dyn_table[i*2] = stof(str.substr(0,white_pos).c_str());
        dyn_table[i*2+1] = stof(str.substr(white_pos+1, str.size()).c_str());
    }
    getline(infile, str);
    n_sample_pts = atoi(str.c_str());
    infile.close();

    sample_space = 1.0/float(n_sample_pts-1);
    out_pts = new float[n_sample_pts*2];

    // dynamic programming 
    for(int k = 0; k < n_sample_pts ; k++){
        t0 = k*sample_space;
        t1 = 1-t0;
        for(int row = 1 ; row < n_control_pts ; row++) {
            upper_base = (row-1)*n_control_pts*2;
            now_base = row*n_control_pts*2;
            for(int col=0 ; col < (n_control_pts - row) ; col++){
                upper_x = dyn_table[upper_base + col*2];
                upper_y = dyn_table[upper_base + col*2 + 1];
                upper_right_x = dyn_table[upper_base + (col+1)*2];
                upper_right_y = dyn_table[upper_base + (col+1)*2 +1];
                dyn_table[now_base + col*2] = t1*upper_x + t0*upper_right_x;
                dyn_table[now_base + col*2 +1] = t1*upper_y + t0*upper_right_y;
            }
        }
        out_pts[k*2] = dyn_table[n_control_pts*(n_control_pts-1)*2];
        out_pts[k*2+1] =  dyn_table[n_control_pts*(n_control_pts-1)*2+1];
    }
    
    //write out results
    ofstream outfile;
    outfile.open(out_file_name);
    for(int i = 0; i < n_sample_pts; i++){
        outfile << setprecision(2) << fixed << out_pts[i*2] << '	'
                << setprecision(2) << fixed << out_pts[i*2+1]
                <<endl;
    }
    outfile.close();
    //release memory
    delete[] dyn_table;
    delete[] out_pts;

    return 0;
}