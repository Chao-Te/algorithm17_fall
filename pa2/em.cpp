#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>

using namespace std;

typedef multimap<const int, int>   edgeMap;
typedef pair<const int, int>  edgePair;

class Node
{
public:
    // member function
    Node(const int& x, const int& y, const int& flow){
        m_x = x;
        m_y = y;
        m_flow = flow;
        //m_size = 0;
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
    //int m_size;

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

            if (m_sources[y_cord].m_flow == 0 || m_sinks[x_cord].m_flow == 0 || flow ==0) 
                continue;// source out of water, drain full of water, edge full
            else{// push water into edge
                flow = min(flow, m_sinks[x_cord].m_flow);
                flow = min(flow, m_sources[y_cord].m_flow);
                m_sinks[x_cord].m_flow   -= flow;
                m_sources[y_cord].m_flow -= flow;
                m_G[pos].m_size          += flow;
            }
        }
        // set initial residual graph

        // construct dynamic table for finding negative cycle
        m_dyn_table = new int[(m_n_sources+m_n_sinks+1)*(m_n_sources+m_n_sinks+2)*2];
    }

    void writeGraph(const string& out_file_name){
        int row = 0;
        ofstream outfile;
        outfile.open(out_file_name.c_str());
        outfile << computeArea() <<endl;
        for(int i=0; i < m_n_sources; i++){
            for(int j = 0; j < m_n_sinks; j++){
                if (m_G[row+j].m_size != 0){
                    outfile << m_sources[i].m_x << ' ' << m_sources[i].m_y << ' '
                            << m_sinks[j].m_x   << ' ' << m_sinks[j].m_y   << ' '
                            << m_G[row+j].m_size << endl;
                }
            }
            row += m_n_sinks;
        }
    }


    void updateflow(){// update size
    }

    void updateres(){// update forward and backward

    }

    void printGraph(){// debugging function
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

private:
    vector<Node> m_sources;
    vector<Node> m_sinks;

    int m_n_sources;// number of row for graph
    int m_n_sinks;// number of column for graph 

    int* m_dyn_table; // for finding negative cycle
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
    //G.printGraph(); // correct
    G.writeGraph(out_file_name);
    return 0;
}