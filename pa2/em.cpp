#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <math>
#include <algorithm>

using namespace std;

class Node
{
public:
    // member function
    Node(const int& x, const int& y, const int& flow){
        m_x = x;
        m_y = y;
        m_flow = flow;
    }
    ~Node(){
    }
    void printNode(){
        cout << "(x, y, f) = ("<< m_x << ", " 
             << m_y << ", " << f << ")" << endl;
    }
    //data member
    int m_x;
    int m_y;
    int m_flow;

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
        m_capacity = max(source.m_flow, sink.m_flow);
        m_size = 0;
        m_dst = abs(source.m_x - sink.m_x) + abs(source.m_y - sink.m_y);
        m_forward = 0;
        m_backward = 0;
    }
    void printEdge(){
        cout << "(capacity, size, dist, forward, backward) = ("
             << m_capacity << ", " << m_size << ", " << m_dist << ", " 
             << m_forward << ", " << m_backward << ")" <<endl;
    }
    // data member
    int m_capacity;
    int m_size;
    int m_dist;
    int m_forward;
    int m_backward;
};
class Graph
{
public:
    Graph(){
        m_n_sources = 0;
        m_n_sinks = 0;
        m_G = NULL;
    }
    ~Graph(){
        if(m_G != NULL)
            delete[] m_G;
    }
    void constructGraph(const string& in_file_name){
        string str;
        int n_nodes = 0;
        int x_cord = 0;
        int y_cord = 0;
        int flow = 0;
        int row = 0;
        ifstream infile;
        infile.open(in_file_name.c_str()); //string to float stof string to int atoi
        getline(infile, str);// get first line
        n_nodes = atoi(str.c_str());
        cout << "number of nodes = " << n_nodes <<endl;
        // reserve memory for speed issue
        m_sources.reserve(n_nodes);
        m_sinks.reserve(n_nodes);
        for (int i = 0; i < n_nodes; i++){
            infile >> x_cord >> y_cord >> flow;
            cout << "x = " << x_cord << ", y = " << y_cord << ", flow = " << flow <<endl;
            if (flow > 0 ){// read a source node
                m_sources.append(Node(x_cord, y_cord, flow);
            }
            else {
                m_sinks.append(Node(x_cord, y_cord, abs(flow)));
            }
        }
        m_n_sources = m_sources.size();
        m_n_sinks = m_sinks.size();
        m_G = new Edge [m_n_sources*m_n_sinks];
        for(int i = 0; i < m_n_sources; i++){
            for(int j = 0; j < m_n_sinks; j++){
                m_G[row+j].setEdge(m_sources[i], m_sinks[j]);
            }
            row += m_n_sources;
        }
    }
    void printGraph(){
        // print Node
        cout << "*****************print sources Node*****************" <<endl;
        for(int i = 0; i < m_n_sources; i ++){
            m_sources[i].printNode();
        }
        cout << "*****************print sinks Node ******************" <<endl;
        for(int i = 0; i < m_n_sinks; i++)
            m_sinks[i].printNode();
        int row = 0;
        for(int i = 0; i < m_n_sources; i++){
            for(int j = 0; j< m_n_sinks; j++){
                m_G[i][j].printEdge();
            }
            row += m_n_sources;
        }

    }

private:
    vector<Node> m_sources;
    vector<Node> m_sinks;

    int m_n_sources;// number of row for graph
    int m_n_sinks;// number of column for graph 

    Edge* m_G;// row is source, column is sink
    /*

    */
};

int main(int argc, char** argv)
{
    // declare variable
    string in_file_name = argv[1];
    string out_file_name = argv[2];

    Graph G;

    G.constructGraph(in_file_name);
    G.printGraph();
    
    return 0;
    
}