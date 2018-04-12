#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <bitset>
#include <cstdio>

using namespace std;

/*  
            N
        |   |   |
        |   |   |   
    -------------------
        | 1 | 2 |
W   -------------------  E
        | 4 | 8 |
    -------------------
        |   |   |
        |   |   |
            S
*/

// use DP and inversion deletion

#define MAXNUM 10000000
#define LOCAL_NUM_CARS 51


class Car
{
public:
    //member function
    Car(const char& to_direction, const char& occupy_blocks){
        m_to_direction = to_direction;
        m_occupy_blocks = occupy_blocks;
    }
    ~Car(){
    }

    //data member
    char m_to_direction;
    char m_occupy_blocks;
};


class Solver
{
public:
    //constructor
    Solver(){
        m_in_cars = new vector<Car*>[4];
        m_c_null = new Car('0',0);

        m_TraceBack = 0;
        m_Time = 0;
        m_num_cars = 0;
        

        // allocate memory for m_max_score
        m_score = new float[5]; // max_score,x,y,z,u

        //allocate memory for m_Corner
        m_Corner = new int[4];

        for(int i = 0; i<5; i++)
        	m_score[i] = 0;
        for(int i = 0; i<4; i++)
        	m_Corner[i] = 0;

        //memset(m_max_score, 0, sizeof(float)*5);
        //memset(m_Corner, 0, sizeof(int)*4);

        for(int i = 0; i<4; i++)
        	last_car_index[i] = 0;

        m_direction[0] = 1;	// N
        m_direction[1] = 2;	// E
        m_direction[2] = 4; // S
        m_direction[3] = 8; // W

        m_direction[4] = 1+2;
        m_direction[5] = 1+4;
        m_direction[6] = 1+8;
        m_direction[7] = 2+4;
        m_direction[8] = 2+8;
        m_direction[9] = 4+8;

        m_direction[10] = 1+2+4;
        m_direction[11] = 1+2+8;
        m_direction[12] = 1+4+8;
        m_direction[13] = 2+4+8;
        m_direction[14] = 1+2+4+8;

        x_dim = LOCAL_NUM_CARS;
        y_dim = LOCAL_NUM_CARS;
        z_dim = LOCAL_NUM_CARS;
        u_dim = LOCAL_NUM_CARS;
    }
    
    ~Solver(){
        for (int i = 0; i < 4; i++){
            for(int j = 0; j < m_in_cars[0].size(); j++){
                if(m_in_cars[i][j] != m_c_null)
                    delete m_in_cars[i][j];
            }
        }
        if (m_c_null != NULL)
            delete m_c_null;

        // release m_TraceBack
        if(m_TraceBack != 0)
        {
            for (int i = 0; i < x_dim; i++)
            {
                for(int j = 0; j < y_dim; j++)
                {
                    for(int k = 0; k < z_dim; k++)
                    {
                        delete[] m_TraceBack[i][j][k];
                        delete[] m_Time[i][j][k];
                        delete[] m_num_cars[i][j][k];
                    }
                    delete[] m_TraceBack[i][j];
                    delete[] m_Time[i][j];
                    delete[] m_num_cars[i][j];
                }
                delete[] m_TraceBack[i];
                delete[] m_Time[i];
                delete[] m_num_cars[i];
            }
 
            delete[] m_TraceBack;
            delete[] m_Time;
            delete[] m_num_cars;
        }
        //delete[] m_Score;
        delete[] m_score;
        delete[] m_Corner;
    }

    void readfile(const string& in_file_name)
    {
        
        // read input file
        string str[4];
        string str1;
        int num_car = 0;
        ifstream infile;
        infile.open(in_file_name.c_str());
        for (int i = 0; i < 4; i++)
            getline(infile, str[i]); // get N,E,S,W
        infile.close();

        // find number of car in a row, and reserve memory for vector for speed issue
        num_car = (str[0].size()-2)/3;
        for (int i = 0; i < 4; i++)
        {
            m_in_cars[i].reserve(num_car);
            cout<<str[i]<<endl;
        }

        Car* c;

        // for car from N
        str1 = str[0].substr(2, str[0].size()-2);
        //cout<<str1<<endl;
        for(int i = 2; i < str1.size(); i+=3){ 
            if (str1[i] == 'S')
                c= new Car('S', 5);
            else if (str1[i] == 'W')
                c = new Car('W', 1);
            else if (str1[i] == 'E')
                c = new Car('E', 13);
            else if (str1[i] == '0')
                c = m_c_null;
            m_in_cars[0].push_back(c);
            if(str1[i] != '0')
            	last_car_index[0] = m_in_cars[0].size() - 1;
        }

        // for car from E
        str1 = str[1].substr(2, str[1].size()-2);
        //cout<<str1<<endl;
        for (int i = 2; i < str1.size(); i+=3){
            if(str1[i] == 'N')
                c = new Car('N', 2);
            else if (str1[i] == 'S')
                c = new Car('S',7);
            else if (str1[i] == 'W')
                c = new Car('W', 3);
            else if (str1[i]== '0')
                c = m_c_null;
            m_in_cars[1].push_back(c);
            if(str1[i] != '0')
            	last_car_index[1] = m_in_cars[1].size() - 1;
        }

        //for car from S
        str1 = str[2].substr(2, str[2].size()-2);
        //cout<<str1<<endl;
        for (int i = 2; i < str1.size(); i+=3)
        {
            if (str1[i] == 'N')
                c = new Car('N',10);
            else if (str1[i] == 'E')
                c = new Car('E', 8);
            else if (str1[i] == 'W')
                c = new Car('W', 11);
            else if (str1[i] == '0')
                c = m_c_null;
            m_in_cars[2].push_back(c);
            if(str1[i] != '0')
            	last_car_index[2] = m_in_cars[2].size() - 1;
        }

        // car from W
        str1 = str[3].substr(2, str[3].size()-2);
        //cout<<str1<<endl;
        for (int i = 2; i < str1.size(); i+=3){
            if (str1[i] == 'N')
                c = new Car('N', 14);
            else if (str1[i] == 'S')
                c = new Car('S', 4);
            else if (str1[i] == 'E')
                c = new Car('E', 12);
            else if (str1[i] == '0')
                c = m_c_null;
            m_in_cars[3].push_back(c);

            if(str1[i] != '0')
            	last_car_index[3] = m_in_cars[3].size() - 1;
            
        }
    }  

    char choose_direction(Car* c_n, Car* c_e, Car* c_s, Car* c_w, int* dst, int& min_dst, short* index, int& max_index)
    {
    	// record whether the current score can come from i (0~15) direction
        bool record[15];
        char dir = 0;
        int max_num = -1;
        max_index = 0;
        //int  max_index = -1;
        //memset(record, 0, sizeof(bool)*15);
        
        for(int i = 0; i < 4; i++)
        { //single car
            record[i] = true;
        }


        record[4] = !(c_n->m_occupy_blocks & c_e->m_occupy_blocks);// n e
        //cout<<"                N E checked"<<endl;

        record[5] = !(c_n->m_occupy_blocks & c_s->m_occupy_blocks);// n s
        //cout<<"                N S checked"<<endl;

        record[6] = !(c_n->m_occupy_blocks & c_w->m_occupy_blocks);// n w
        //cout<<"                N W checked"<<endl;

        record[7] = !(c_e->m_occupy_blocks & c_s->m_occupy_blocks);// e s
        //cout<<"                E S checked"<<endl;

        record[8] = !(c_e->m_occupy_blocks & c_w->m_occupy_blocks);// e w
        //cout<<"                E W checked"<<endl;

        record[9] = !(c_s->m_occupy_blocks & c_w->m_occupy_blocks);// s w
        //cout<<"                S W checked"<<endl;
        
        record[10] = record[4] & record[5] & record[7]; // n e s
        record[11] = record[4] & record[6] & record[8]; // n e w
        record[12] = record[5] & record[6] & record[9]; // n s w
        record[13] = record[7] & record[8] & record[9]; // s e w
        record[14] = record[10] & record[11] & record[12] & record[13]; // n e s w

        for(int i = 0; i < 15; i++)
        {
            if(min_dst > dst[i] && record[i] && max_num <= index[i])
            {
                dir = m_direction[i];
                min_dst = dst[i];
                max_num = index[i];
                max_index = i;
            }
        }
        return dir;
    }


    void initDP(const int& x_max, const int& y_max, const int& z_max, const int& u_max)
    {

        x_dim = x_max;
        y_dim = y_max;
        z_dim = z_max;
        u_dim = u_max;

        m_TraceBack = new char***[x_dim];
        m_Time = new short***[x_dim];
        m_num_cars = new short***[x_dim];
        
        // allocate memory m_TraceBack
        for(int i = 0; i < x_dim; i++)
        {
            m_TraceBack[i] =    new char**[y_dim];
            m_Time[i] =         new short**[y_dim];
            m_num_cars[i] =     new short**[y_dim];
            for (int j = 0; j < y_dim; j++)
            {
                m_TraceBack[i][j] = new char*[z_dim];
                m_Time[i][j] =      new short*[z_dim];
                m_num_cars[i][j] =  new short*[z_dim];
                for (int k = 0; k < z_dim; k++)
                {
                    m_TraceBack[i][j][k] =  new char[u_dim];
                    m_Time[i][j][k] =       new short[u_dim];
                    m_num_cars[i][j][k] =   new short[u_dim];
                }
            }    
        }

		m_TraceBack[0][0][0][0] = 0;
        m_Time[0][0][0][0]      = 0;
        m_num_cars[0][0][0][0]  = 0;

        for(int x = 1; x < x_max; x++)
        {
        	m_TraceBack[x][0][0][0] = 1;
        	m_Time[x][0][0][0] = x;     
        }
        for(int y = 1; y < y_max; y++)
        {
        	m_TraceBack[0][y][0][0] = 2;
        	m_Time[0][y][0][0] = y;         
        }
        for(int z = 1; z < z_max; z++)
        {
        	m_TraceBack[0][0][z][0] = 4;
        	m_Time[0][0][z][0] = z;
        }
        for(int u = 1; u < u_max; u++)
        {
        	m_TraceBack[0][0][0][u] = 8;
        	m_Time[0][0][0][u] = u;
        }

        initDP_numCars(x_max, y_max, z_max, u_max);
    }

    void initDP_numCars(const int& x_max, const int& y_max, const int& z_max, const int& u_max)
    {
        for(int x = 1; x < x_max; x++)
        {
            if(m_in_cars[0][m_Corner[0] + x - 1]->m_to_direction !='0')
                m_num_cars[x][0][0][0] = m_num_cars[x-1][0][0][0] + 1;
            else 
                m_num_cars[x][0][0][0] = m_num_cars[x-1][0][0][0];
        }
        for(int y = 1; y < y_max; y++)
        {
            if(m_in_cars[1][m_Corner[1] + y - 1]->m_to_direction !='0')
                m_num_cars[0][y][0][0] = m_num_cars[0][y-1][0][0] + 1;
            else
                m_num_cars[0][y][0][0] = m_num_cars[0][y-1][0][0];
        }
        for(int z = 1; z < z_max; z++)
        {
            if(m_in_cars[2][m_Corner[2] + z - 1]->m_to_direction != '0')
                m_num_cars[0][0][z][0] = m_num_cars[0][0][z-1][0] + 1;
            else 
                m_num_cars[0][0][z][0] = m_num_cars[0][0][z-1][0];
        }
        for(int u = 1; u < u_max; u++)
        {
            if(m_in_cars[3][m_Corner[3] + u - 1]->m_to_direction != '0')
                m_num_cars[0][0][0][u] = m_num_cars[0][0][0][u-1] + 1;
            else
                m_num_cars[0][0][0][u] = m_num_cars[0][0][0][u-1];
        }
    }

    void resetDP(const int& x_max, const int& y_max, const int& z_max, const int& u_max)
    {
        if(m_Time != 0)
        {
    		for (int i = 0; i < y_dim; i++)
    		{
                for(int j = 0; j < z_dim; j++)
                {
                    for(int k = 0; k < u_dim; k++)
                    {
                        delete [] m_TraceBack[i][j][k];
                        delete [] m_Time[i][j][k];
                    }
                    delete [] m_TraceBack[i][j];
                    delete [] m_Time[i][j];
                }
                delete [] m_TraceBack[i];
                delete [] m_Time[i];
            }
            delete [] m_Time;
            delete [] m_TraceBack;
        }

        initDP(x_max, y_max, z_max, u_max);

    }

    void run() //bool writeFile, const string& filename)
    {
        // from direction -> m_traceBack value
        // N              -> 1
        // E              -> 2
        // S              -> 4
        // W              -> 8

    	vector<char> tb_path;					// record the path of traceback
    	vector<char> schedule;

    	schedule.reserve(m_in_cars[0].size()*4);
        if(check_rest_size())
    	   initDP(x_dim, y_dim, z_dim, u_dim); 	// boundary condition

        int step = 0;
        while (check_rest_size())
        {
            initDP_numCars(x_dim, y_dim, z_dim, u_dim);

            cout<<"Step: "<<step<<endl;
        	cout<<"Fixed size window"<<endl;
        	cout<<"    car index (N, E, S, W): ";
        	cout<<"("<<m_Corner[0]<<", "<<m_Corner[1]<<", "<<m_Corner[2]<<", "<<m_Corner[3]<<")"<<endl;

        	//calculate score for DP table (x_dim, y_dim, z_dim, u_dim)
        	// store the index needed in m_Score
        	windowDP_score(LOCAL_NUM_CARS, LOCAL_NUM_CARS, LOCAL_NUM_CARS, LOCAL_NUM_CARS);
        	
        	DP_traceback(m_score[1], m_score[2], m_score[3], m_score[4], tb_path);
        	//DP_traceback(x_dim -1, y_dim -1, z_dim -1, u_dim -1, tb_path);

        	
        	for(int i = 0; i < tb_path.size(); i++)
        		schedule.push_back(tb_path[tb_path.size()-1-i]); // reverse
        	
            int start = traj[0].size();
            updateTrajectory(tb_path);

            // optimize trajectory starting from the last position before update
            // use the information of m_Corner and x_dim, y_dim, z_dim, u_dim
            // reorder such that the flow of multiple cars is placed first
            optimizeTrajectory(start);

        	//printTB(schedule, false, " ");
            //printTrajectory();

            for(int i = 0; i<4; i++)
                m_Corner[i] = m_Corner[i] + m_score[i+1];

            

            step++;

        }

        if (!check_rest_size())
        {
        	cout<<" Dynamic sized window"<<endl;
			cout<<"    car index (N, E, S, W): ";
        	cout<<"("<<m_Corner[0]<<", "<<m_Corner[1]<<", "<<m_Corner[2]<<", "<<m_Corner[3]<<")"<<endl;
        	
        	// dimension required for DP (=car num + 1)
        	int x_max = last_car_index[0] - m_Corner[0] + 1 + 1;
        	int y_max = last_car_index[1] - m_Corner[1] + 1 + 1;
        	int z_max = last_car_index[2] - m_Corner[2] + 1 + 1;
        	int u_max = last_car_index[3] - m_Corner[3] + 1 + 1;

        	// reset DP table size and set boundary condition
        	cout<<"        Reset DP table"<<endl;
        	resetDP(x_max, y_max, z_max, u_max);

        	cout<<"            "<<x_dim<<" "<<y_dim<<" "<<z_dim<<" "<<u_dim<<endl;

        	cout<<"        Compute DP table"<<endl;

        	windowDP_score(x_dim, y_dim, z_dim, u_dim);

            m_score[1] = x_dim-1;
            m_score[2] = y_dim-1;
            m_score[3] = z_dim-1;
            m_score[4] = u_dim-1;

        	cout<<"        Traceback DP"<<endl;
        	DP_traceback(x_dim -1, y_dim -1, z_dim -1, u_dim -1, tb_path);

        	//for(int i = 0; i < tb_path.size(); i++)
        	//	schedule.push_back(tb_path[tb_path.size()-1-i]); // reverse

        	cout<<"        Print schedule"<<endl;
        	//printTB(schedule, true, filename);

            int start = traj[0].size();
            updateTrajectory(tb_path);
            //printTrajectory();

            optimizeTrajectory(start);
            //printTrajectory();

            for(int i = 0; i<4; i++)
                m_Corner[i] = m_Corner[i] + m_score[i+1];
        	
        }


        //run algorithm
        //Trace Back
    }


    void windowDP_score(const int& x_max, const int& y_max, const int& z_max, const int& u_max)
    {
		int* dst = new int[15];		// dst better use time
		short num_cars[15];
        short acc_cars[15];
        char from_direction = 0;
        int min_dst = MAXNUM;
        int max_index = 0;
        float heuristic_value = 0;
        /* 
        direction to distance
            N(x-1) : dst[0],  E(y-1) : dst[1],  S(z-1) : dst[2],  W(u-1) : dst[3]
            
            N_E    : dst[4],  N_S    : dst[5],  N_W    : dst[6],  E_S    : dst[7], E_W    : dst[8], S_W     : dst[9]
            
            N_E_S  : dst[10], N_E_W  : dst[11], N_S_W  : dst[12], E_S_W  : dst[13]
            
            N_E_S_W: dst[14]    
        */

        for(int x = 0; x < x_max; x++)
        {	//from north
            for(int y = 0; y < y_max; y++)
            {	//from east
            	for(int z = 0; z < z_max; z++)
            	{	//from south
            		for(int u = 0; u < u_max; u++)
            		{	//from west

            			//cout<<"    ("<<x<<", "<<y<<", "<<z<<", "<<u<<"): "<<endl;
			                
			            if(x == 0 && y == 0 && z == 0)
			            {
			            	//cout<<m_Time[x][y][z][u]<<endl;
			              	continue;
			            }
			            else if(x == 0 && y == 0 && u == 0)
			            {
			            	//cout<<m_Time[x][y][z][u]<<endl;
			               	continue;
			            }
			            else if(x == 0 && z == 0 && u == 0)
			            {
			            	//cout<<m_Time[x][y][z][u]<<endl;
			               	continue;
			            }
			            else if(y == 0 && z == 0 && u == 0)
			            {
			            	//cout<<m_Time[x][y][z][u]<<endl;
			               	continue;
			            }

						min_dst = MAXNUM;
			                // reset 15 neighboring values each iteration
			            for (int k = 0; k < 15; k ++)
			            {
			                dst[k] = -1;
			                num_cars[k] = 0;
			            }
			                
			                // not boundary condition, from 15 directions
			                // direction = 0 -> dst[that direction] = MAXNUM
			                // the car from that direction will be replaced with m_c_null
			            Car* c_n = m_c_null;
			            Car* c_e = m_c_null;
			            Car* c_s = m_c_null;
			            Car* c_w = m_c_null;

			            if(x != 0)
			              	c_n = m_in_cars[0][m_Corner[0] + x -1];
			            else
			            {
			             	dst[0]  = MAXNUM;
			               	dst[4]  = MAXNUM;
			               	dst[5]  = MAXNUM;
			               	dst[6]  = MAXNUM;
			               	dst[10] = MAXNUM;
			               	dst[11] = MAXNUM;
			               	dst[12] = MAXNUM;
			               	dst[14] = MAXNUM;
		                }
		                if(c_n->m_to_direction != '0')
		                {
		                	num_cars[0]  += 1;
			               	num_cars[4]  += 1;
			               	num_cars[5]  += 1;
			               	num_cars[6]  += 1;
			               	num_cars[10] += 1;
			               	num_cars[11] += 1;
			               	num_cars[12] += 1;
			               	num_cars[14] += 1;
		                }

			            if(y != 0)
			               	c_e = m_in_cars[1][m_Corner[1] + y -1];
			            else
			            {
			              	dst[1]  = MAXNUM;
			               	dst[4]  = MAXNUM;
			               	dst[7]  = MAXNUM;
			               	dst[8]  = MAXNUM;
			               	dst[10] = MAXNUM;
			               	dst[11] = MAXNUM;
			               	dst[13] = MAXNUM;
			               	dst[14] = MAXNUM;
			            }
			            if(c_e->m_to_direction != '0')
			            {
			              	num_cars[1]  += 1;
			                num_cars[4]  += 1;
			                num_cars[7]  += 1;
			                num_cars[8]  += 1;
			                num_cars[10] += 1;
			                num_cars[11] += 1;
			                num_cars[13] += 1;
			                num_cars[14] += 1;
			            }

			            if(z != 0)
			               	c_s = m_in_cars[2][m_Corner[2] + z -1];
			            else
			            {
		                	dst[2]  = MAXNUM;
		                	dst[5]  = MAXNUM;
		                	dst[7]  = MAXNUM;
		                	dst[9]  = MAXNUM;
		                	dst[10] = MAXNUM;
		                	dst[12] = MAXNUM;
		                	dst[13] = MAXNUM;
		                	dst[14] = MAXNUM;
		                }
		                if(c_s->m_to_direction != '0')
		                {
			             	num_cars[2]  += 1;			              	
							num_cars[5]  += 1;
			              	num_cars[7]  += 1;
			              	num_cars[9]  += 1;
			              	num_cars[10] += 1;
			              	num_cars[12] += 1;
			              	num_cars[13] += 1;
			              	num_cars[14] += 1;
		                }

			            if(u != 0)
			              	c_w = m_in_cars[3][m_Corner[3] + u -1];
			            else
		                {
		                	dst[3]  = MAXNUM;
		                	dst[6]  = MAXNUM;
		                	dst[8]  = MAXNUM;
		                	dst[9]  = MAXNUM;
		                	dst[11] = MAXNUM;
		                	dst[12] = MAXNUM;
		                	dst[13] = MAXNUM;
		                	dst[14] = MAXNUM;
		                }
		                if(c_w->m_to_direction != '0')
		                {
		                	num_cars[3]  += 1;
		                	num_cars[6]  += 1;
		                	num_cars[8]  += 1;
		                	num_cars[9]  += 1;
		                	num_cars[11] += 1;
		                	num_cars[12] += 1;
		                	num_cars[13] += 1;
		                	num_cars[14] += 1;
		                }

		                //cout<<"        find value"<<endl;
			            for(int i = 0; i<15; i++)
			            {
			            	//cout<<"            "<<i<<endl;
			              	if(dst[i] == -1)
			              	{

		               			switch(i)
		                		{
		                			case 0:
		                				dst[i] = m_Time[x-1][y][z][u] + 1; 
                                        acc_cars[i] = m_num_cars[x-1][y][z][u];
                                        break;
		                			case 1:
		                				dst[i] = m_Time[x][y-1][z][u] + 1; 
                                        acc_cars[i] = m_num_cars[x][y-1][z][u];
                                        break;
		                			case 2:
		                				dst[i] = m_Time[x][y][z-1][u] + 1; 
                                        acc_cars[i] = m_num_cars[x][y][z-1][u];
                                        break;
		                			case 3:
		                				dst[i] = m_Time[x][y][z][u-1] + 1; 
                                        acc_cars[i] = m_num_cars[x][y][z][u-1];
                                        break;
		                			case 4:
		                				dst[i] =      m_Time[x-1][y-1][z][u] + 1; 
                                        acc_cars[i] = m_num_cars[x-1][y-1][z][u];
                                        break;
		                			case 5:
		                				dst[i] =      m_Time[x-1][y][z-1][u] + 1; 
                                        acc_cars[i] = m_num_cars[x-1][y][z-1][u];
                                        break;
		                			case 6:
		                				dst[i] =      m_Time[x-1][y][z][u-1] + 1; 
                                        acc_cars[i] = m_num_cars[x-1][y][z][u-1];
                                        break;
		                			case 7:
		                				dst[i] = m_Time[x][y-1][z-1][u] + 1; 
                                        acc_cars[i] = m_num_cars[x][y-1][z-1][u];
                                        break;
		                			case 8:
		                				dst[i] = m_Time[x][y-1][z][u-1] + 1; 
                                        acc_cars[i] = m_num_cars[x][y-1][z][u-1];
                                        break;
		                			case 9:
		                				dst[i] = m_Time[x][y][z-1][u-1] + 1; 
                                        acc_cars[i] = m_num_cars[x][y][z-1][u-1];
                                        break;
		                			case 10:
		                				dst[i] = m_Time[x-1][y-1][z-1][u] + 1; 
                                        acc_cars[i] = m_num_cars[x-1][y-1][z-1][u];
                                        break;
		                			case 11:
		                				dst[i] = m_Time[x-1][y-1][z][u-1] + 1; 
                                        acc_cars[i] = m_num_cars[x-1][y-1][z][u-1];
                                        break;
		                			case 12:
		                				dst[i] = m_Time[x-1][y][z-1][u-1] + 1;
                                        acc_cars[i] = m_num_cars[x-1][y][z-1][u-1]; 
                                        break;
		                			case 13:
		                				dst[i] = m_Time[x][y-1][z-1][u-1] + 1; 
                                        acc_cars[i] = m_num_cars[x][y-1][z-1][u-1];
                                        break;
		                			case 14:
		                				dst[i] = m_Time[x-1][y-1][z-1][u-1] + 1;
                                        acc_cars[i] = m_num_cars[x-1][y-1][z-1][u-1];
		                			default:
		                				 ;
		                		}
			                }
			            }

			            from_direction = choose_direction(c_n, c_e, c_s, c_w, dst, min_dst, num_cars, max_index);
			            m_Time[x][y][z][u] = min_dst;
			            m_TraceBack[x][y][z][u] = from_direction;
                        m_num_cars[x][y][z][u] = acc_cars[max_index] + num_cars[max_index];

			            //cout<<endl;
                        /*
			            for(int i = 0; i<15; i++)
			            {
			            	if(from_direction == m_direction[i])
                            {
			            		//cout<<"        from direction: "<<i<<endl;
                            }
			            }
                        */
			            //cout<<endl<<"        "<<m_Time[x][y][z][u]<<endl;

			            //if (m_score[0] < float(x*x+y*y+z*z+u*u)/float(min_dst))
			            //{
			            //    m_score[0] = float(x*x+y*y+z*z+u*u)/float(min_dst);

                        if( (x >= x_dim/2) && (x==y) && (y==z) && (z==u) )
                        {
                            if(float(m_num_cars[x][y][z][u])/float(m_Time[x][y][z][u]) >= heuristic_value)
                            {
                                //cout<<"num of cars = "<<m_num_cars[x][y][z][u]<<", time = "<<m_Time[x][y][z][u]<<endl;
    			                m_score[1] = x;
    			                m_score[2] = y;
    			                m_score[3] = z;
    			                m_score[4] = u;
                                heuristic_value = float(m_num_cars[x][y][z][u])/float(m_Time[x][y][z][u]);
                                //cout<<"heuristic value update: "<<heuristic_value<<endl;
                            }
                        }
			            //}  
			        }
			    }
			}
        }
        delete [] dst;
    }


    void DP_traceback(const int& x_start, const int& y_start, const int& z_start, const int& u_start, vector<char>& tb_path)
    {
    	int temp = x_dim + y_dim + z_dim + u_dim;
    	int id_x = x_start;
    	int id_y = y_start;
    	int id_z = z_start;
    	int id_u = u_start;
    	char ptr; 

    	tb_path.clear();
    	tb_path.reserve(temp);

    	for(int i = 0; i < temp; i++)
    	{
    		ptr = m_TraceBack[id_x][id_y][id_z][id_u];
    		if(ptr != 0)
    		{
    			//cout<<"            "<<(ptr & 1)<<" "<<((ptr & 2)>>1)<<" ";
    			//cout<<((ptr & 4)>>2)<<" "<<((ptr & 8)>>3)<<endl;
    			tb_path.push_back(ptr);
    			id_x = id_x - (ptr & 1);
    			id_y = id_y - ((ptr & 2)>>1);
    			id_z = id_z - ((ptr & 4)>>2);
    			id_u = id_u - ((ptr & 8)>>3);
    		}
    		else
    		{
    			break;
    		}
    	}
    }


    // debugging function
    void printInCars(){
        /*
        print the read in cars
        */
        // print N
        for (int i = 0; i < 4; i++){
            if(i == 0)
                cout << "From North : (\t"<<last_car_index[0]<<")"<<endl;
            else if(i == 1)
                cout << "From East  : (\t"<<last_car_index[1]<<")"<<endl;
            else if(i == 2)
                cout << "From South : (\t"<<last_car_index[2]<<")"<<endl;
            else
                cout << "From West  : (\t"<<last_car_index[3]<<")"<<endl;
            for(int j = 0; j < m_in_cars[i].size(); j++){
                if (m_in_cars[i][j] == m_c_null)
                    cout << 0 << " ";
                else 
                    cout << m_in_cars[i][j]->m_to_direction << " ";

            }
            cout << endl;
        }
    }

    bool check_rest_size()
    {
        /*
        check if the reset of some direction < 5
        */
        for(int i = 0; i < 4; i++){
            if(m_Corner[i] + LOCAL_NUM_CARS-1 >= last_car_index[i] + 1)
                return false; 
                //should create a table with size
                // (m_in_cars[0]-m_corner[0])*(m_in_cars[1]-m_corner[1])
                //*(m_in_cars[2]-m_corner[2])*(m_in_cars[3]-m_corner[3])
        }
        return true;
    }

    void printTB(const vector<char>& path, bool writeFile, const string& filename)
    {
    	// convert char in path into N E S W form
    	// require the information of m_in_cars

    	// record the index of cars go through
    	int id_N = 0, id_E = 0, id_S = 0, id_W = 0;

    	// record the car destination for 4 lanes at different times
    	vector<char> trajectory[4];
    	for(int i = 0; i < path.size(); i++)
    	{
    		if((path[i] & 1) == 1)
    		{
    			trajectory[0].push_back(m_in_cars[0][id_N]->m_to_direction);
    			id_N += 1;
    		}
    		else
    			trajectory[0].push_back('0');
    		if((path[i] & 2) == 2)
    		{
    			trajectory[1].push_back(m_in_cars[1][id_E]->m_to_direction);
    			id_E += 1;
    		}
    		else
    			trajectory[1].push_back('0');
    		if((path[i] & 4) == 4)
    		{
    			trajectory[2].push_back(m_in_cars[2][id_S]->m_to_direction);
    			id_S += 1;
    		}
    		else
    			trajectory[2].push_back('0');
    		if((path[i] & 8) == 8)
    		{
    			trajectory[3].push_back(m_in_cars[3][id_W]->m_to_direction);
    			id_W += 1;
    		}
    		else
    			trajectory[3].push_back('0');
    	}
    	for(int i = 0; i < 4; i++)
    	{
    		if(i == 0)
    			cout<<"    N:"<<endl;
    		else if(i == 1)
    			cout<<"    E:"<<endl;
    		else if(i == 2)
    			cout<<"    S:"<<endl;
    		else
    			cout<<"    W:"<<endl;
    		cout<<"        ";
    		for(int j = 0; j < trajectory[0].size(); j++)
    		{
    			cout<<trajectory[i][j]<<" ";
    		}
    		cout<<endl;
    	}
    	cout<<"        Trajectory printed"<<endl;

    	if(writeFile)
    	{
    		writeSchedule(filename, trajectory);
    		cout<<"        Trajectory written to file"<<endl;
    	}
    }


    void updateTrajectory(const vector<char>& tb_path)
    {
        // convert char in path into N E S W form
        // require the information of m_in_cars

        // record the index of cars go through
        int id_N = m_Corner[0];
        int id_E = m_Corner[1];
        int id_S = m_Corner[2];
        int id_W = m_Corner[3];

        // record the car destination for 4 lanes at different times in traj
        for(int i = tb_path.size()-1; i >= 0; i--)
        {
            if((tb_path[i] & 1) == 1)
            {
                traj[0].push_back(m_in_cars[0][id_N]->m_to_direction);
                id_N += 1;
            }
            else
                traj[0].push_back('0');
            if((tb_path[i] & 2) == 2)
            {
                traj[1].push_back(m_in_cars[1][id_E]->m_to_direction);
                id_E += 1;
            }
            else
                traj[1].push_back('0');
            if((tb_path[i] & 4) == 4)
            {
                traj[2].push_back(m_in_cars[2][id_S]->m_to_direction);
                id_S += 1;
            }
            else
                traj[2].push_back('0');
            if((tb_path[i] & 8) == 8)
            {
                traj[3].push_back(m_in_cars[3][id_W]->m_to_direction);
                id_W += 1;
            }
            else
                traj[3].push_back('0');
        }
    }


    void optimizeTrajectory(const int& start)
    {
        // place flow of multiple cars in front of that of a single car

        // record the time when the latest cars enter the intersection
        
        int first_cars[4];

        bool current_cars[4];
        int current_car_flow;
        bool prev_cars[4];
        int prev_car_flow;

        // first_cars originally == time when the first car enter the intersection
        for(int i = 0; i<4; i++)
        {
            int temp = m_Corner[i] + m_score[i+1];
            first_cars[i] = temp;
            for(int j = m_Corner[i]; j < temp; j++)
            {
                if(m_in_cars[i][j]->m_to_direction != '0')
                {    
                    first_cars[i] = j;
                    break;
                }
            }
        }


        // optimize the trajectory from "start" poistion to end
        for(int i = start; i < traj[0].size(); i++)
        {


            //printTrajectory();
            // initialize trajectory at the current moment
            //cout<<"    i = "<<i<<endl;
            current_car_flow = 0;

            bool exchange = false;
            int  exchange_idx = 0;

            for(int j = 0; j<4; j++)
            {
                if(traj[j][i] != '0')
                {
                    current_car_flow += 1;
                    current_cars[j] = true;
                }
                else
                    current_cars[j] = false;
                prev_cars[j] = false;



            }

            
            // check trajectory of previous time
            for(int j = i-1; j>= start; j--)
            {

                bool check_yield = false;
                prev_car_flow = 0;
                for(int k = 0; k<4; k++)
                {
                    if(traj[k][j] != '0')
                    {
                        prev_car_flow += 1;
                        //if(prev_cars[k])
                        //    check_yield = true;
                        prev_cars[k] = true;
                    }
                    else
                        prev_cars[k] = false;
                    check_yield = ( check_yield | (prev_cars[k]&current_cars[k]) );
                }

               // cout<<"           j = "<<j<<", flow = "<<prev_car_flow<<endl;

                if(check_yield)
                {
                    //cout<<"        "<<i<<" yields "<<j<<endl;
                    break; 
                }     
                    // means the current car in some lane must yield some car in the same lane
                
                if(prev_car_flow < current_car_flow)
                {
                    for(int k = 0; k<4; k++)
                    {
                        if(current_cars[k]) 
                        {

                            if(j <first_cars[k])
                            {
                                //cout<<"               cannot exchange: j = "<<j;
                                //cout<<", first_car["<<k<<"] = "<<first_cars[k]<<endl;
                                break;
                            }
                        }  
                        if(k == 3)
                        {
                            exchange = true;
                            exchange_idx = j;
                        }  
                    }
                }
            }

            char char_temp;
            for(int k = 0; k<4; k++)
            {
                
                if(exchange)
                {
                    //cout<<"        exchange with "<<exchange_idx<<endl;
                    char_temp = traj[k][i];
                    for(int t = i; t > exchange_idx; t--)
                    {
                        traj[k][t] = traj[k][t-1];
                    }
                    traj[k][exchange_idx] = char_temp;
                    
                }

                // update the first car not passing 
                if(current_cars[k])
                {
                    //cout<<k<<endl;
                    int max_temp = m_Corner[k] + m_score[k+1];
                    //first_cars[k] = max_temp;
                    
                    for(int m = first_cars[k]+1; m < max_temp; m++)
                    {
                        
                        if(m_in_cars[k][m]->m_to_direction != '0')
                        {    

                            first_cars[k] = m;

                            //if(k == 3)
                            //    cout<<cout<<"                 = "
                            break;
                        }
                        if(m == max_temp - 1)
                            first_cars[k] = m;
                    }
                }

            }


        }
        
        

    }

    void writeSchedule(const string& filename, const vector<char>* trajectory)
    {
        FILE* fout;
        fout = fopen(filename.c_str(), "w");
        
        for(int i = 0; i<4; i++)
        {
        	if(i == 0)
    			fprintf(fout, "N:");
    		else if(i == 1)
    			fprintf(fout, "E:");
    		else if(i == 2)
    			fprintf(fout, "S:");
    		else
    			fprintf(fout, "W:");
    		
    		for(int j = 0; j < trajectory[0].size(); j++)
    		{
    			if(trajectory[i][j] != '0')
    				fprintf(fout, " 1%c", trajectory[i][j]);
    			else
    				fprintf(fout, " 00");
    		}
    		fprintf(fout, "\n");
        }
    }


    // use member variable traj
    void printTrajectory()
    {
        for(int i = 0; i < 4; i++)
        {
            if(i == 0)
                cout<<"    N:"<<endl;
            else if(i == 1)
                cout<<"    E:"<<endl;
            else if(i == 2)
                cout<<"    S:"<<endl;
            else
                cout<<"    W:"<<endl;
            cout<<"        ";
            for(int j = 0; j < traj[0].size(); j++)
            {
                cout<<traj[i][j]<<" ";
            }
            cout<<endl;
        }
        cout<<"        Trajectory printed"<<endl;
    }

    // use member variable traj
    void writeTrajectory(const string& filename)
    {
        FILE* fout;
        fout = fopen(filename.c_str(), "w");
        
        for(int i = 0; i<4; i++)
        {
            if(i == 0)
                fprintf(fout, "N:");
            else if(i == 1)
                fprintf(fout, "E:");
            else if(i == 2)
                fprintf(fout, "S:");
            else
                fprintf(fout, "W:");
            
            for(int j = 0; j < traj[0].size(); j++)
            {
                if(traj[i][j] != '0')
                    fprintf(fout, " 1%c", traj[i][j]);
                else
                    fprintf(fout, " 00");
            }
            fprintf(fout, "\n");
        }
    }


private:
    char m_direction[15];
    vector<Car*>* m_in_cars;

    Car*    	m_c_null;       // No car!
    char****    m_TraceBack;   // Tensor for trace back [6][6][6][6]
    short****   m_Time;        // Tensor to record Time [2][6][6][6]

    short****   m_num_cars;    // record the number of cars having gone through 

    float*      m_score;   //record maximum score and direction of local DP table [5]

    int*        m_Corner;      // the corner of the starting point [4], on the tensor n*n*n*n

    int  		x_dim;			// record the size of DP table (x_dim*y_dim*z_dim*u_dim)
    int 		y_dim;
    int 		z_dim;
    int 		u_dim;

    int 		last_car_index[4];

    vector<char>        traj[4];        // record the schedule of 4 lanes

};

int main(int argc, char** argv)
{
    string in_file_name = argv[1];
    string out_file_name = argv[2];

    vector<Car*>* in_cars = new vector<Car*>[4];
    vector<Car*>* out_cars = new vector<Car*>[4];
    
    //Car* c_null;
    Solver solver;
    solver.readfile(in_file_name);

    cout<<endl;
    solver.printInCars();

    solver.run();

    solver.writeTrajectory(out_file_name);
    solver.printTrajectory();

    cout<<"Solved"<<endl;
}
