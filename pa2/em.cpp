#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>

// test1.txt : optimal = 48
using namespace std;

typedef multimap<const int, const int>   edgeMap;
typedef pair<const int, const int>  edgePair;

class Bin
{
public:
    Bin(){
        m_idx = 0;
        m_value = 0;
    }
    ~Bin(){}
    int m_idx;
    int m_value;

};
class Node
{
public:
    // member function
    Node(const int& x, const int& y, const int& flow){
        m_x = x;
        m_y = y;
        m_flow = flow;
        m_size = flow;
    }
    ~Node(){
    }
    void printNode(){
        cout << "(x, y, f, s) = ("<< m_x << ", " 
             << m_y << ", " << m_flow << ")" <<endl;// ", "<< m_size <<")" << endl;
    }
    //data member
    int m_x;
    int m_y;
    int m_flow;
    int m_size;

};

class Edge
{
public:
    // member function
    Edge(){
        m_capacity = 0;
        m_size = 0;
    }

    ~Edge(){
    }

    void setEdge(const Node& source, const Node& sink){
        m_capacity = min(source.m_flow, sink.m_flow);
        m_size = 0;
        m_dst = abs(source.m_x - sink.m_x) + abs(source.m_y - sink.m_y);
        m_forward = 0;
        m_backward = 0;
    }
    void printEdge(){
        cout << "(capacity, size, dist, forward, backward) = ("
             << m_capacity << ", " << m_size << ", " << m_dst << ", " 
             << m_forward << ", " << m_backward << ")" <<endl;
    }
    int getArea(){
        return m_size*m_dst;
    }
    void set1res(){
        m_backward = m_size;
        m_forward = m_capacity - m_size;
    }

    void updateflow(const int& min_flow, bool forward){
        if(forward)
            m_size += min_flow;
        else
            m_size-=min_flow;
        m_forward = m_capacity - m_size;
        m_backward = m_size;
    }
    // data member
    int m_capacity;
    int m_size;
    int m_dst;
    int m_forward;  // source -> sink
    int m_backward; // sink   -> source
};

class Graph
{
public:
    Graph(){
        m_n_sources = 0;
        m_n_sinks = 0;
        m_total_flow = 0;
        m_G = NULL;
        m_dyn_table = NULL;
    }

    ~Graph(){
        if(m_G != NULL)
            delete[] m_G;
        if(m_dyn_table != NULL)
            delete[] m_dyn_table;
    }

    void constructGraph(const string& in_file_name){
        string str;
        int n_nodes = 0;
        int x_cord = 0;
        int y_cord = 0;
        int flow = 0;
        int pos;
        int row = 0;
        edgeMap emap;

        ifstream infile;
        infile.open(in_file_name.c_str()); //string to float stof string to int atoi
        getline(infile, str);// get first line
        n_nodes = atoi(str.c_str());
        //cout << "number of nodes = " << n_nodes <<endl;
        // reserve memory for speed issue
        m_sources.reserve(n_nodes);
        m_sinks.reserve(n_nodes);
        for (int i = 0; i < n_nodes; i++){
            infile >> x_cord >> y_cord >> flow;
            //cout << "x = " << x_cord << ", y = " << y_cord << ", flow = " << flow <<endl;
            if (flow > 0 ){// read a source node
                m_sources.push_back(Node(x_cord, y_cord, flow));
                m_total_flow += flow;
            }
            else {
                m_sinks.push_back(Node(x_cord, y_cord, abs(flow)));
            }
        }
        

        m_n_sources = m_sources.size();
        m_n_sinks = m_sinks.size();
        m_G = new Edge [m_n_sources*m_n_sinks];
        for(int i = 0; i < m_n_sources; i++){
            for(int j = 0; j < m_n_sinks; j++){
                pos = row+j;
                m_G[pos].setEdge(m_sources[i], m_sinks[j]);
                emap.insert(edgePair(m_G[pos].m_dst, pos));
            }
            row += m_n_sinks;
        }

        // set up initial flow 
        for (edgeMap::iterator it=emap.begin(); it!=emap.end(); ++it){
            pos = it->second;
            y_cord = pos/m_n_sinks;
            x_cord = pos%m_n_sinks;
            flow = m_G[pos].m_capacity - m_G[pos].m_size;//the capacity of edge now
            //print Map for debugging
            //cout << "source : " << y_cord << "; sink : " << x_cord << "; dst = " << it->first <<endl; 

            if (m_sources[y_cord].m_size == 0 || m_sinks[x_cord].m_size == 0 || flow ==0) 
                continue;// source out of water, drain full of water, edge full
            else{// push water into edge
                flow = min(flow, m_sinks[x_cord].m_size);
                flow = min(flow, m_sources[y_cord].m_size);
                m_sinks[x_cord].m_size   -= flow;
                m_sources[y_cord].m_size -= flow;
                m_G[pos].m_size          += flow;
            }
        }

        // set initial residual graph
        for(row = 0; row < m_n_sources * m_n_sinks; row += m_n_sinks){
            for(int j = 0; j < m_n_sinks; j++){
                m_G[row+j].set1res();
            }
        }

        // construct dynamic table for finding negative cycle
        m_dyn_table = new Bin[(m_n_sources+m_n_sinks+1)*(m_n_sources+m_n_sinks+2)];
    }
    
    void run(){
        cout << "in run" <<endl;
        int n_row = 2 + m_n_sources + m_n_sinks;
        int n_col = 1 + m_n_sources + m_n_sinks;
        int row = 0;
        int pos = 0;
        int min_value;
        int nei_value;
        int idx;
        vector<int> cycle_node;
        bool neg_cycle = true;
        bool* hits = new bool[1+m_n_sources+m_n_sinks];
        
        //bell-forman, matrix.shape = (m_n_sources+m_n_sinks+1, (m_n_sources+m_n_sinks+2))
        while(neg_cycle){
            for (int i = 0; i < 1+m_n_sources+m_n_sinks; i++)
                hits[i] = false;
            
            neg_cycle = false;

            cycle_node.clear();
            // construct table
            for (int i = 2; i < n_row; i++){
                // check sources( col = [1, m_n_sources+1) 
                for(int j = 1; j < 1+m_n_sources; j++){
                    pos = (i-1)*n_col + j; // position for previous m_dyn_table
                    min_value = m_dyn_table[pos].m_value;
                    idx = j;
                    pos = (j-1)*m_n_sinks;// position for m_G
                    for(int k = 0; k < m_n_sinks; k++){
                        if (m_G[pos].m_forward != 0){
                            nei_value = m_dyn_table[(i-1)*n_col + (1 + m_n_sources + k)].m_value + m_G[pos].m_dst;
                            if (nei_value < min_value){
                                min_value = nei_value;
                                idx = 1 + m_n_sources + k;
                            }
                        }
                        pos ++;
                    }
                    pos = i*n_col + j;
                    m_dyn_table[pos].m_value = min_value;
                    m_dyn_table[pos].m_idx = idx;
                }

                //check sinks (row = [m_n_source+1, m_source+m_n_sinks+1)
                for(int j = 1+m_n_sources; j < n_col; j++){
                    pos = (i-1)*n_col + j;
                    min_value = m_dyn_table[pos].m_value;
                    idx = j;
                    for(int k = 0; k < m_n_sources; k++){

                        pos = k*m_n_sinks + (j-1-m_n_sources);
                        if (m_G[pos].m_backward != 0){
                            nei_value = m_dyn_table[(i-1)*n_col + (k+1)].m_value - m_G[pos].m_dst;
                            if (nei_value < min_value){
                                min_value = nei_value;
                                idx = 1+k;
                            }
                        }
                    }
                    pos = i*n_col + j;
                    m_dyn_table[pos].m_value = min_value;
                    m_dyn_table[pos].m_idx = idx;
                }

            }
            printTable();
            min_value = 10000000;
            // check negative cycle
            for(int j = 0; j < n_col; j++){
                if (m_dyn_table[(n_row-2)*n_col+j].m_value != m_dyn_table[(n_row-1)*n_col+j].m_value){
                    // negative cycle exits
                    neg_cycle = true;
                    cout  << "Found negative cycle" <<endl;
                    // trace back
                    pos = j;
                    
                    for (int i = m_n_sources+m_n_sinks+1; i>=1; i--){
                        // debugging
                        if(pos == 0){
                            cout << "something wierd" << endl;
                        }
                        cout << "pos = " << pos << endl;
                        hits[pos] = true;
                        row = i*n_col;
                        cycle_node.push_back(pos-1);
                        if (pos <= m_n_sources ){ // it is a sources, and forward
                            cout <<  "it is source" <<endl;
                            if (min_value > m_G[(pos-1)*m_n_sinks + m_dyn_table[row + pos].m_idx-1-m_n_sources].m_forward)
                                min_value = m_G[(pos-1)*m_n_sinks + m_dyn_table[row + pos].m_idx-1-m_n_sources].m_forward;
                        }
                        else{// it is a sink and backward
                            cout << "it is sink" <<endl; 
                            if (min_value > m_G[(m_dyn_table[row + pos].m_idx-1)*m_n_sinks + pos-1-m_n_sources].m_backward)
                                min_value = m_G[(m_dyn_table[row + pos].m_idx-1)*m_n_sinks + pos-1-m_n_sources].m_backward;
                        }
                        pos = m_dyn_table[row + pos].m_idx;
                        
                        if (hits[pos]){
                            cout << "time to eliminate negative cycle" <<endl;
                            cycle_node.push_back(pos-1);
                            updateflow(cycle_node, min_value);
                        }
                    }

                    break;
                }
                pos += (m_n_sources+m_n_sinks+2);
            }
            // debugging message
            break;
        }   
        delete[] hits;
        
        
    } 
    void writeGraph(const string& out_file_name){
        int row = 0;
        m_total_flow = 0;
        ofstream outfile;
        outfile.open(out_file_name.c_str());
        outfile << computeArea() <<endl;
        for(int i=0; i < m_n_sources; i++){
            for(int j = 0; j < m_n_sinks; j++){
                if (m_G[row+j].m_size != 0){
                    m_total_flow += m_G[row+j].m_size;
                    outfile << m_sources[i].m_x << ' ' << m_sources[i].m_y << ' '
                            << m_sinks[j].m_x   << ' ' << m_sinks[j].m_y   << ' '
                            << m_G[row+j].m_size << endl;
                }
            }
            row += m_n_sinks;
        }
    }

    void updateflow(const vector<int>& cycle_node, const int& min_flow){// update flow
        int head = 0;
        int tail = 0; 
        for (int i = 0; i < cycle_node.size(); i++){
            cout << cycle_node[i] << endl;
        }
        //cout << "min_flow = " << min_flow <<endl;
        for(int i = 0; i < cycle_node.size()-1; i++ ){
            if (cycle_node[i] < m_n_sources){//source and forward
                
                tail = cycle_node[i];// tail = source
                head = cycle_node[i+1] - m_n_sources;// head = sinks
                //cout << "i = " << i << " : source -> sink" << "(" << tail << "," << head << ")"<< endl;
                m_G[tail*m_n_sinks+head].updateflow(min_flow, true);
            }
            else{// sink and backward
                
                tail = cycle_node[i] - m_n_sources; //sink
                head = cycle_node[i+1]; //source
                //cout << "i = " << i << " : sink -> source" << "(" << tail << "," << head << ")"<<endl;
                m_G[head*m_n_sinks+tail].updateflow(min_flow, false);
            }
        }

    }

    void printGraph(){// DEBUGGING FUNCTION
        // print Node
        cout << "number of nodes = " << m_n_sources+m_n_sinks <<endl;
        cout << "*****************print sources Node*****************" <<endl;
        for(int i = 0; i < m_n_sources; i ++){
            m_sources[i].printNode();
        }
        cout << "*****************print sinks Node ******************" <<endl;
        for(int i = 0; i < m_n_sinks; i++)
            m_sinks[i].printNode();
        int row = 0;
        cout << "*****************print Graph matrix*****************" << endl;
        for(int i = 0; i < m_n_sources; i++){
            for(int j = 0; j< m_n_sinks; j++){
                m_G[row+j].printEdge();
            }
            row += m_n_sinks;
        }

    }
    
    int computeArea(){
        int row = 0;
        int area = 0;
        for (int i = 0; i < m_n_sources; i++){
            for(int j = 0; j < m_n_sinks; j++){
                area += m_G[row+j].getArea();
            }
            row+=m_n_sinks;
        }
        return area;
    }
    void check(){// DEBUGGING FUNCTION
        int* cal_source = new int[m_n_sources];
        int* cal_sink   = new int[m_n_sinks];
        int row = 0;
        bool res = true;
        for(int i = 0; i < m_n_sources; i++)
            cal_source[i] = 0;
        for(int i = 0; i < m_n_sinks; i++)
            cal_sink[i] = 0;
        for(int i = 0; i < m_n_sources; i++){
            for(int j = 0; j < m_n_sinks; j++){
                cal_source[i] += m_G[row+j].m_size;
                cal_sink[j] += m_G[row+j].m_size;
            }
            row += m_n_sinks;
        } 
        // check source
        for(int i = 0; i < m_n_sources; i++){
            if (cal_source[i] != m_sources[i].m_flow)
                cout << "error for flow out of source " << i << endl;
        }

        // check sinks
        for(int i = 0; i < m_n_sinks; i++){
            if (cal_sink[i] != m_sinks[i].m_flow)
                cout << "error for flow out of sink " << i << endl;
        }

        delete[] cal_source;
        delete[] cal_sink;
    }
    void printTable(){// DEBUGGING FUNCTION
        cout << "-------------- print distance ----------------" <<endl;
        for(int i = 0; i < m_n_sources+m_n_sinks+2; i++){
            for(int j = 0; j < m_n_sources+m_n_sinks+1; j++){
                cout.width(4);
                cout << right << m_dyn_table[i*(m_n_sources+m_n_sinks+1) + j].m_value;
            }
            cout << endl;
        }
        cout << "-------------- print idx      ----------------" << endl;
        for(int i = 0; i < m_n_sources+m_n_sinks+2; i++){
            for(int j = 0; j < m_n_sources+m_n_sinks+1; j++){
                cout.width(4);
                cout << right << m_dyn_table[i*(m_n_sources+m_n_sinks+1) + j].m_idx;
            }
            cout << endl;
        }


    }
    // public data member
    int m_total_flow;
private:
    vector<Node> m_sources;
    vector<Node> m_sinks;

    int m_n_sources;// number of row for graph
    int m_n_sinks;// number of column for graph 
    
    Bin* m_dyn_table; // for finding negative cycle
    // 0 : distance to NULL node
    // 1 : the row of it comes from(Node)
    // #row = m_n_sources + m_n_sinks + 1
    // #col = m_n_sources + m_n_sinks + 2

    Edge* m_G;// row is source, column is sink
    /*   t1  t2   t3 ...
    s1  e11  e12  e13
    s2  e21  e22  e23
    s3
    .
    .
    .
    */
    
};

int main(int argc, char** argv)
{
    // declare variable
    string in_file_name = argv[1];
    string out_file_name = argv[2];

    Graph G;

    G.constructGraph(in_file_name);
    cout << "total flow when reading Graph : " << G.m_total_flow <<endl;
    G.printGraph(); // correct
    G.check();
    G.run();
    G.writeGraph(out_file_name);
    cout << "total flow after ALGORITHM Graph : " << G.m_total_flow <<endl;
    return 0;
}
