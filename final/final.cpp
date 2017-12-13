#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <bitset>

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

class Car
{
public:
    //member function
    Car(const char& to_direction, const char& occupy_blocks){
        m_to_direction = to_direction;
        m_occupy_blocks = occupy_blocks;
    }
    ~Car();

    //data member
    char m_to_direction;
    char m_occupy_blocks;
};
Car* readfile(const string& in_file_name, vector<Car*>* in_cars){
    Car* c_null = new Car(0,0);
    // read input file
    string str[4];
    string str1;
    int num_car = 0;
    ifstream infile;
    infile.open(in_file_name.c_str());
    for (int i = 0; i < 4; i++)
        getline(infile, str[i]); // get N,E,S,W
    infile.close();
    
    // debugging message
    for (int i = 0; i < 4; i++)
        cout << str[i] <<endl;

    // find number of car in a row, and reserve memory for vector for speed issue
    num_car = (str[0].size()-2)/3;
    for (int i = 0; i < 4; i++)
        in_cars[i].reserve(num_car);

    Car* c;
    // for car from N
    str1 = str[0].substr(2, str[0].size());
    for(int i = 2; i < str1.size(); i+=3){ 
        if (str1[i] == 'S')
            c= new Car('S', 5);
        else if (str[i] == 'W')
            c = new Car('W', 1);
        else if (str[i] == 'E')
            c = new Car('E', 13);
        else if (str[i] == '0',0)
            c = c_null;
        in_cars[0].push_back(c);
    }
    // for car from E
    str1 = str[1].substr(2, str[1].size());
    for (int i = 0; i < str1.size(); i+=3){
        if(str1[i] == 'N')
            c = new Car('N', 2);
        else if (str1[i] == 'S')
            c = new Car('S',7);
        else if (str1[i] == 'W')
            c = new Car('W', 3);
        else if (str1[i]== '0')
            c = c_null;
        in_cars[1].push_back(c);
    }
    //for car from S
    str1 = str[2].substr(2, str[2].size());
    for (int i = 0; i < str1.size(); i+=3){
        if (str1[i] == 'N')
            c = new Car('N',10);
        else if (str[i] == 'E')
            c = new Car('E', 8);
        else if (str[i] == 'W')
            c = new Car('W', 11);
        else if (str[i] == '0')
            c = c_null;
    }
    // car from W
    str1 = str[3].substr(2, str[3].size());
    for (int i = 0; i < str1.size(); i+=3){
        if (str1[i] == 'N')
            c = new Car('N', 14);
        else if (str1[i] == 'S')
            c = new Car('S', 4);
        else if (str1[i] == 'E')
            c = new Car('E', 12);
        else if (str[i] == '0')
            c = c_null
        in_cars.push_back(c);
    }
    return c_null;
}
void outfile(const string& out_file_name, vector<Car*>* out_cars){
    ofstream outfile;
    outfile.open(out_file_name.c_str());
    for (int d = 0; d < 4; d++){
        if(d == 0)        outfile << "N:";
        else if (d == 1)  outfile << "E:";
        else if (d == 2)  outfile << "S:";
        else              outfile << "W:";

        for (int i = 0; i < out_cars[0].size() ; i++){
            if out_cars[0][i].m_to_direction == '0':
                outfile << " 00";
            else{
                coutfile << " 1" << out_cars[0][i].m_to_direction;
                delete out_cars[0][i];
            } 
        }
        outfile << endl;
    }
     
}
int bit2int(char blocks){
    int num = 0;
    for (int i = 0; i < 4; i++){
        if (blocks%2)
            num ++;
        blocks = blocks >> 1;
    }
    return num;
}
void check_2_cars(Car* cars,int* records){
    // i = 0, j = 1, N&E records[1]
    // i = 0, j = 2, N&S records[2]
    // i = 0, j = 3, N&W records[3]
    // i = 1, j = 2, E&S records[4]
    // i = 1, j = 3, E&W records[5]
    // i = 2, j = 3, S&W records[6]
    int r_idx = 1;
    for(int i = 0; i < 4; i++){
        for (int j = i+1; j < 4; j++){
            if(int(cars[i].m_occupy_blocks & cars[j].m_occupy_blocks) > 0)  
                records[r_idx] = -1;
            else    
                records[r_idx] = bit2int(cars[i].m_occupy_blocks | cars[j].m_occupy_blocks)
        }
    }
    // check N & E & S
    if (records[1] == -1 || records[2] == -1 || records[4] == -1) records[6] = -1;
    // check N & E & W||
    if (records[1] == -1 || records[3] == -1 || records[5] == -1) records[7] = -1;
    // check N & S & W||
    if (records[2] == -1 || records[3] == -1 || records[6] == -1) records[8] = -1;
    //check E & S & W||
    if (records[4] == -1 || records[5] == -1 || records[6] == -1) records[9] = -1;   
}
void check_3_cars(Cars* cars, int* records){
    // N & E & S
    if (records[6] != -1)
        records[6] = bit2int(cars[0]| cars[1] | cars[2]);
    if (records[7] != -1)
        records[7] = bit2int(cars[0]| cars[1] | cars[3]);
    if (records[8] != -1)
        records[8] = bit2int(cars[0]| cars[2] | cars[3]);
    if (records[9] != -1)
        records[9] = bit2int(cars[1]| cars[2] | cars[3]);
    for (int i = 0; i < 4 ; i++){
        if (records[i] == -1){
            records[10] = -1;
            break;
        }
    }
}
int choose_cars(int * records){
    // 1 : North
    // 2 : East
    // 4 : Sorth
    // 8 : West
    if (records[10] != -1)
        return (1+2+4+8);
    if (records[9] != -1)
        return (2+4+8);
    if (records[8] != -1)
        return (1+4+8);
    if (records[7] != -1)
        return (1+2+8);
    if (records[6] != -1)
        return (1+2+4);
    if (records[5] != -1)
        return (4+8);
    if (records[4] != -1)
        return (2+8);
    if (records[3] != -1)
        return (2+4);
    if (records[2] != -1)
        return (1+8);
    if (records[1] != -1)
        return (1+4);
    if (records[0] != -1)
        return (1+2);
    
}
void alg1(vector<Car*>* in_cars, vector<Car*>* out_cars, Car* c_null){

    int records[11] ;// C^4_2 = 6 + C^4_3 = 4 + C^4_4 = 1  = 11
    // record[10] = N & E & S & W, 
    // record[0] : N & E, record[1] : N & S, record[2] : N & W
    // reocrd[3] : E & S, record[4] : E & W, record[5] : S & W
    // record[6] : N & E & S, record[7] : N & E & W
    // record[8] : N & S & W, record[9] : E & S & W
    // -1 means those cars can not pass together at the same time
    // other number means those cars can pass at the same time 
    // and the number of blocks they occupied

    int num_cars = in_cars[0].size();
    int t_cars = 0;// cars to move forward
    Car cars[4]; // cars from north, east, south, west
    bool cars_empty[4];
    bool all_empty = false;
    for (int i = 0; i < 4; i++)
        cars_empty[i] = false;
    
    while(!all_empty){

        for (int i = 0; i < 4; i++){
            if (cars_empty[i])
                cars[i] = c_null;
            else
                cars[i] = in_cars[i][0]
        }
    
        check_2_cars(cars, records);
        // empty lane(no cars now) should be move forward immediately
        check_3_cars(cars, records);
        // check_4_cars
        if (records[10]!= -1){
            records[10] = bit2int(cars[0].m_occupy_blocks | cars[1].m_occupy_blocks | 
                                  cars[2].m_occupy_blocks | cars[3].m_occupy_blocks);
        }

        t_cars = choose_cars(records);

        for (int i = 0; i < 4; i++){
            if (cars_empty[i]){ // no more car in line
                out_cars[i].append(c_null);
            }
            else if (in_cars[i][0] == c_null || t_cars%2) { //// no car waiting there now
                in_cars[i].erase(0);
                out_cars[i].append(cars[i]);
            }
            else{
                cout << "something wrong when assigning cars" <<endl;
            }
            t_cars >> 1;
        }

        // clear record, and reset memory
        all_empty = true;
        for (int i = 0; i < 11; i++)
            records[i] = 0;
            c_empty[i] = in_cars[i].empty();
            if (!c_empty[i])
                all_empty = false;
        t_cars = 0;
    }   
}

int main(int argc, char** argv)
{
    string in_file_name = argv[1];
    string out_file_name = argv[2];

    vector<Car*>* in_cars = new vector<Car*>[4];
    vector<Car*>* out_cars = new vector<Car*>[4];
    
    Car* c_null;

    c_null = readfile(in_file_name, in_cars)
    
    // reserve memory for vector
    for (int i = 0; i < 4: i++)
        out_cars.reserve(in_cars[0].size()*4);

    // do algo
    alg1(in_cars, out_cars, c_null);

    // write result to output file
    outfile(out_file_name, out_cars);

    // release memory 
    delete c_null;
    for(int i = 0; i < 4; i++){
        in_cars[i].clear();
        out_cars[i].clear();
    }

}
