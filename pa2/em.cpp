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
typedef pair<const int, const bool> cyclePair;
typedef pair<int, int> Bin;

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
    bool operator<(const Node& r_node)const{
        if (m_x < r_node.m_x)
            return true;
        else if(m_x > r_node.m_x)
            return false;
        else if(m_y < r_node.m_y)
            return true;
        else if(m_y > r_node.m_y)
            return false;
        else if(m_flow > r_node.m_flow)
            return true;
        else
            return false; 
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
        m_dst = 0;
        m_forward = 0;
        m_s_x = 0;
        m_s_y = 0;
        m_t_x = 0;
        m_t_y = 0; 
    }
    //Edge(const Node& source, const Node& sink){
    //    m_capacity = 0;
    //    m_size = 0;
    //    m_dst = 0;
    //    m_forward = 0;
    //    m_s_x = 0;
    //    m_s_y = 0;
    //    m_t_x = 0;
    //    m_t_y = 0; 
    //}

    ~Edge(){
    }
    void setEdge(const Node& source, const Node& sink){
        m_capacity = min(source.m_flow, sink.m_flow);
        m_size = 0;
        m_dst = abs(source.m_x - sink.m_x) + abs(source.m_y - sink.m_y);
        m_forward = 0;
        m_s_x = source.m_x;
        m_s_y = source.m_y;
        m_t_x = sink.m_x;
        m_t_y = sink.m_y;
    }
    void printEdge(){
        cout << "(capacity, size, dist, forward, backward, m_s_id, m_t_id) = ("
             << m_capacity << ", " << m_size << ", " << m_dst << ", " 
             << m_forward << ", " << m_size <<  ")" << endl;
    }
    int getArea(){
        return m_size*m_dst;
    }
    void set1res(){
        //m_backward = m_size;
        m_forward = m_capacity - m_size;
    }

    void updateflow(const int& min_flow, bool forward){
        if(forward)
            m_size += min_flow;
        else
            m_size-=min_flow;
        m_forward = m_capacity - m_size;
        //m_backward = m_size;
    }
    bool operator<(const Edge& r_edge) const{
        if(m_s_x < r_edge.m_s_x)
            return true;
        else if(m_s_x > r_edge.m_s_x)
            return false;
        else if(m_s_y < r_edge.m_s_y)
            return true;
        else if(m_s_y > r_edge.m_s_y)
            return false;
        else if(m_t_x < r_edge.m_t_x)
            return true;
        else if(m_t_x > r_edge.m_t_x)
            return false;
        else if(m_t_y < r_edge.m_t_y)
            return true;
        else if(m_t_y > r_edge.m_t_y)
            return false;
        else if(m_size < r_edge.m_size)
            return true;
        else 
            return false;
    }
    // data member
    int m_capacity;
    int m_size;
    int m_forward;  // source -> sink
    int m_dst;
    int m_s_x;
    int m_s_y;
    int m_t_x;
    int m_t_y;
    //int m_backward; // sink   -> source
};

class Graph
{
public:
    Graph(){
        m_n_sources = 0;
        m_n_sinks = 0;
        m_total_flow = 0;
        m_dyn_table = NULL;
    }

    ~Graph(){
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
        
        // reserve memory for speed issue
        m_sources.reserve(n_nodes);
        m_sinks.reserve(n_nodes);
        for (int i = 0; i < n_nodes; i++){
            infile >> x_cord >> y_cord >> flow;
            if (flow > 0 ){// read a source node
                m_sources.push_back(Node(x_cord, y_cord, flow));
            }
            else {
                m_sinks.push_back(Node(x_cord, y_cord, abs(flow)));
            }
        }
        sort(m_sources.begin(), m_sources.end());
        sort(m_sinks.begin(), m_sinks.end());
        // setup edge graph
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
        
        // set up initial flow, greedy
        for (edgeMap::iterator it=emap.begin(); it!=emap.end(); ++it){
            pos = it->second;
            y_cord = pos/m_n_sinks;
            x_cord = pos%m_n_sinks;
            flow = m_G[pos].m_capacity - m_G[pos].m_size;//the capacity of edge now

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
        for(int i = 0; i < m_n_sources * m_n_sinks; i++){    
            m_G[i].set1res();
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
        int pos2 = 0;
        int min_value;
        int nei_value;
        bool neg_cycle = true;
        bool* hits = new bool[1+m_n_sources+m_n_sinks];
        vector<int> cycle_node;
        
        //bell-forman, matrix.shape = (m_n_sources+m_n_sinks+1, (m_n_sources+m_n_sinks+2))
        while(neg_cycle){
            for (int i = 0; i < 1+m_n_sources+m_n_sinks; i++)
                hits[i] = false;
            neg_cycle = false;
            cycle_node.clear();
        
            for (int i = 2; i < n_row; i++){
                row = (i-1)*n_col;
                pos = row +n_col;
                for (int j = 1; j < n_col; j++){
                    m_dyn_table[pos+j].first = m_dyn_table[row+j].first;
                    m_dyn_table[pos+j].second = j;
                }
                pos = 0;
                for (int r = 0; r < m_n_sources; r++){
                    for(int c = 0; c < m_n_sinks; c++){
                        //source - > sink
                        pos=r*m_n_sinks + c;
                        if (m_G[pos].m_forward != 0){
                            //cout << "In forward, r = " << r << "; c = " << c <<endl;
                            pos2 = row + n_col + (1 + r);
                            nei_value = m_dyn_table[row + (1 + m_n_sources + c)].first + m_G[pos].m_dst;
                            if (m_dyn_table[pos2].first > nei_value){
                                m_dyn_table[pos2].first = nei_value;
                                m_dyn_table[pos2].second = 1 + m_n_sources + c;
                            }
                        }
                        //sink->source
                        if(m_G[pos].m_size != 0){
                            //cout << "In backward, r = " << r << "; c = " << c <<endl;
                            pos2 = row + n_col + (1 + m_n_sources + c);
                            nei_value = m_dyn_table[row + (1 + r)].first - m_G[pos].m_dst;
                            if (m_dyn_table[pos2].first > nei_value){
                                m_dyn_table[pos2].first = nei_value;
                                m_dyn_table[pos2].second = 1 + r;
                            }
                        }
                    }
                    
                }
            }
            
            // check negative cycle
            for(int j = 1; j < n_col; j++){
                if (m_dyn_table[(n_row-2)*n_col+j].first != m_dyn_table[(n_row-1)*n_col+j].first){
                    // negative cycle exits
                    neg_cycle = true;
                    pos = j;
                    for (int i = m_n_sources+m_n_sinks+1; i>=1; i--){
                        hits[pos] = true;
                        row = i*n_col;
                        cycle_node.push_back(pos-1);
                        
                        pos = m_dyn_table[row + pos].second;
                        
                        if (hits[pos]){
                            cycle_node.push_back(pos-1);
                            while(cycle_node[cycle_node.size()-1] != cycle_node[0]){
                                cycle_node.erase(cycle_node.begin());
                            }
                            updateflow(cycle_node);
                            break;
                        }
                    }	
                    if(neg_cycle)
                        break;
                }
                pos += (m_n_sources+m_n_sinks+2);
            }
            //cout << "Area = " << computeArea() << endl;
        }   
        delete[] hits;
        
        
    }

    void writeGraph(const string& out_file_name){
        int row = 0;
        //m_total_flow = 0;
        ofstream outfile;
        outfile.open(out_file_name.c_str());
        outfile << computeArea() <<endl;
        // sort edge
        vector<Edge> sortG(m_G, m_G + m_n_sources*m_n_sinks); 
        sort (sortG.begin(), sortG.end());

        for(int i = 0; i < sortG.size(); i++){
            if (sortG[i].m_size != 0){
                //m_total_flow += sortG[i].m_size;
                outfile << sortG[i].m_s_x << ' ' << sortG[i].m_s_y << ' '
                        << sortG[i].m_t_x   << ' ' << sortG[i].m_t_y   << ' '
                        << sortG[i].m_size << endl;
            }
        }
    }

    void updateflow(const vector<int>& cycle_node){
        //find min value
        vector<cyclePair> cycle_edge;
        //vector<bool> fflow;
        int min_flow = 1000000000;
        int pos = 0;
        cycle_edge.reserve(cycle_node.size()-1);
        for(int i = 0; i < cycle_node.size()-1; i++){
            if (cycle_node[i] < m_n_sources){//source and forward
                pos = cycle_node[i]*m_n_sinks+(cycle_node[i+1] - m_n_sources);
                cycle_edge.push_back(edgePair(pos, true));
                if (m_G[pos].m_forward < min_flow)
                    min_flow = m_G[pos].m_forward;
            }
            else {
                pos = cycle_node[i+1]*m_n_sinks+(cycle_node[i] - m_n_sources);
                cycle_edge.push_back(edgePair(pos, false));
                if (m_G[pos].m_size < min_flow)
                    min_flow = m_G[pos].m_size;
            }
        }
        for(int i = 0; i < cycle_edge.size(); i++){
            if(cycle_edge[i].second)
                m_G[cycle_edge[i].first].updateflow(min_flow, true);
            else    
                m_G[cycle_edge[i].first].updateflow(min_flow, false);
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
        int area = 0;
        for(int i = 0; i < m_n_sources*m_n_sinks; i++){
            area += m_G[i].getArea();
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
                if(m_G[row+j].m_size<0)
                    cout << "flow < 0 !!!" <<endl;
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
    /*
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


    }*/
    
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
    cout << "total area after reading Graph : " << G.computeArea() <<endl;
    //G.printGraph(); // correct
    G.run();
    //G.printGraph();
    G.check();
    G.writeGraph(out_file_name);
    cout << "total flow after algorithm Graph : " << G.m_total_flow <<endl;
    return 0;
}
