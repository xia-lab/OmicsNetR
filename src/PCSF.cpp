
// ************************  Prize-collecting Steiner Forest (PCSF) *********************//
// This is a C++ implementation of our heuristic algorithm for the PCSF (Akhmedov et al 2017).
// The heuristic is developed for functional analyses of large biological networks in a reasonable time.


// Loading the required libraries
#include <Rcpp.h>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/random.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <time.h>
#include <stdio.h>
#include <boost/graph/graphviz.hpp>
#include <fstream>
#include <vector>
#include <utility>
#include <string>
#include <math.h>
#include <iomanip>
#include <boost/limits.hpp>
#include <queue>

using namespace boost;
using namespace Rcpp;
using namespace std;


// Node class
class Node
{  public:
  vector <int> children;
  int father;
  int size;
  double price;
  Node(){};
  Node(const Node & other){
    father=other.father;
    children= other.children;}

  Node & operator= (const Node & other){
    father=other.father;
    children= other.children;
    return *this;
  }

};

// Properties of Vertices in the network
struct VertexProperties  {
  VertexProperties() :c(0) {}
  VertexProperties(string const & name) : name(name),c(0){}
  string name;
  double c;
};

// Properties of the graph that is used within the BOOST
typedef adjacency_list<vecS, vecS, undirectedS, VertexProperties, property < edge_weight_t, double > > GraphBase;
typedef graph_traits<GraphBase>::vertex_iterator vertex_iterator;
typedef graph_traits<GraphBase>::out_edge_iterator edge_iterator;
typedef graph_traits<GraphBase>::edge_iterator graph_edge_iterator;
typedef graph_traits<GraphBase>::edge_descriptor Edge;
typedef graph_traits<GraphBase>::vertex_descriptor Vertex;


GraphBase g, g_adjusted; GraphBase G; GraphBase G_pruned;
property_map<GraphBase, edge_weight_t>::type weight_g;
property_map<GraphBase, edge_weight_t>::type weight_g_adjusted;
property_map<GraphBase, edge_weight_t>::type weight_G;
property_map<GraphBase, edge_weight_t>::type weight_G_pruned;

int Root = -1;
static map <string, int> g_map;
static map <string, int> G_map;
static map <string, int> G_pruned_map;

void clear_variables(){
  g.clear();
  g_adjusted.clear();
  G.clear();
  G_pruned.clear();
  g_map.clear();
  G_map.clear();
  G_pruned_map.clear();
}

// Map
int idx_g(String const & id)
{
  map<string, int>::iterator mit = g_map.find(id);
  if (mit == g_map.end())
    return g_map[id] = add_vertex(VertexProperties(id), g);
  return mit->second;
}
int idx_G(string const & id)
{
  map<string, int>::iterator mit = G_map.find(id);
  if (mit == G_map.end())
    return G_map[id] = add_vertex(VertexProperties(id), G);
  return mit->second;
}
int idx_G_pruned(string const & id)
{
  map<string, int>::iterator mit = G_pruned_map.find(id);
  if (mit == G_pruned_map.end())
    return G_pruned_map[id] = add_vertex(VertexProperties(id), G_pruned);
  return mit->second;
}


// Reading the input network
void read_input_graph(CharacterVector from, CharacterVector to,  NumericVector cost, CharacterVector prize, NumericVector prize_v)
{
  for(int i=0; i < from.size(); i++){
    add_edge(vertex(idx_g(from[i]), g), vertex(idx_g(to[i]), g), cost[i], g);
  }
  for(int i=0; i<prize.size(); i++){
    g[idx_g(prize[i])].c = prize_v[i];
  }

  // cerr << num_edges(g) << " edges, " << num_vertices(g) << " vertices" << endl;

  g_adjusted = g;
  weight_g_adjusted = get(edge_weight, g_adjusted);
  graph_edge_iterator ei, ei_end; double penalty; Vertex sour, tar;
  for(tie(ei, ei_end) = edges(g_adjusted); ei != ei_end; ++ei){
    sour = source(*ei,g_adjusted); tar = target(*ei,g_adjusted); penalty = 0;
    if(g_adjusted[sour].c < 0 && g_adjusted[tar].c < 0){
      penalty = g_adjusted[sour].c + g_adjusted[tar].c;
    } else if( g_adjusted[sour].c < 0 ){
      penalty = g_adjusted[sour].c;
    } else if( g_adjusted[tar].c < 0 ){
      penalty = g_adjusted[tar].c;
    }
    weight_g_adjusted[*ei] = weight_g_adjusted[*ei] + abs(penalty);
  }
}


// A function to dynamically remove the leaf node if its prize smaller than connection cost.
// It is used called within the process_leafs() function.
void clear(vector <Node> & predecessor, int & current_node){
  int node=current_node;
  for(unsigned int j=0; j<predecessor[node].children.size(); j++){
    clear(predecessor, predecessor[node].children[j]);
  }
  predecessor[node].father=-1;
}

// A function to dynamically remove the leaf node if its prize smaller than connection cost.
void process_leafs(vector <Node> & predecessor, int & current_node, Edge &e, bool &found){
  for(unsigned int i=0; i<predecessor[current_node].children.size(); i++){
    process_leafs(predecessor, predecessor[current_node].children[i], e, found);
  }
  int node = current_node;
  if(node != predecessor[node].father){
    boost::tuples::tie(e,found) = edge( vertex(predecessor[node].father, G_pruned),vertex(node, G_pruned) , G_pruned);
    if(predecessor[node].price - weight_G_pruned[e] <= 0) {
      clear(predecessor, node);
    }
  }

}

// A function to dynamically sum up the prizes of vertices.
void price_collect(vector <Node> & predecessor, int & current_node, Edge &e, bool &found){
  for(unsigned int i=0; i<predecessor[current_node].children.size(); i++){
    price_collect(predecessor, predecessor[current_node].children[i], e, found);
  }
  int node = current_node;
  if(node != predecessor[node].father){
    boost::tuples::tie(e,found) = edge(vertex(predecessor[node].father, G_pruned),vertex(node, G_pruned) , G_pruned);
    if(predecessor[node].price - weight_G_pruned[e] > 0)
      predecessor[predecessor[node].father].price = predecessor[predecessor[node].father].price + predecessor[node].price - weight_G_pruned[e];
  }
}




// After reading the input network information from the input file, the algorithm constructs a
vector< Vertex > constructG(vector<int> & terminals, int &Root){

  // Distance: all-pairs-shortest-path distance matrix
  // perPath: List of arcs in all-pairs-shortest-path distance matrix
  vector <vector <vector<int> > > perPath;
  vector <vector<double> > Distance;
  perPath.resize (terminals.size());
  Distance.resize (terminals.size());
  for (unsigned int i = 0; i < terminals.size(); ++i) {
    perPath [i].resize(terminals.size());
    Distance [i].resize(terminals.size());
  }

  // Computing all-pairs-shortest-path distance matrix
  Vertex from; int current, pred, outer = 0, inner;
  std::vector<Vertex> p = vector<Vertex> (num_vertices(g));
  std::vector<double> d = vector<double> (num_vertices(g));

  for (std::vector<int>::iterator first=terminals.begin(); first!=terminals.end(); ++first){

    from = vertex(*first, g_adjusted);
    dijkstra_shortest_paths(g_adjusted, from,
                            predecessor_map(boost::make_iterator_property_map(p.begin(), get(boost::vertex_index, g_adjusted))).
                              distance_map(boost::make_iterator_property_map(d.begin(), get(boost::vertex_index, g_adjusted))));

    inner = outer + 1;
    for (std::vector<int>::iterator second = first+1; second!=terminals.end(); ++second) {
      current=*second; pred=p[current];
      Distance[outer][inner] = d[current];
      Distance[inner][outer] = d[current];
      while(pred!=current){perPath[outer][inner].push_back(current); current=pred; pred=p[current]; } perPath[outer][inner].push_back(current);
      inner++;
    }
    outer++;
  }




  // Heuristic Clustering, given large input network, the algorithm clusters input network into
  // smaller clusters, and solves the MST afterwards
  set<int> V;
  set<int> D;
  unsigned int root_index = -1;
  vector<int> node_labels;
  node_labels.resize(terminals.size());
  for(unsigned int i=0; i< node_labels.size(); i++){
    node_labels[i] = 0;}
  for(unsigned int i=0; i<terminals.size(); i++){
    if(terminals[i] == Root){
      root_index = i;}
    else{
      V.insert(i);}
  }
  node_labels[root_index] = INT_MAX;
  int clusterID=0; int targ;

  while(!V.empty()){
    clusterID++;  current = *V.begin();
    node_labels[current] =clusterID;
    V.erase(current);

    D.clear();
    for (unsigned int i=0; i < terminals.size(); i++) {
      targ = i;
      if(node_labels[targ] == 0 && current != targ && i != root_index){
        if(g[terminals[current]].c >= Distance[current][targ] && g[terminals[targ]].c >= Distance[current][targ]){
          if(g[terminals[targ]].c > 0){
            D.insert(targ); V.erase(targ);
          }
          node_labels[targ]=clusterID;
        }
      }
    }

    while(!D.empty()){
      current = *D.begin();  D.erase(current);
      for (unsigned int i=0; i < terminals.size(); i++) {
        targ = i;
        if(node_labels[targ] == 0 && current != targ && i != root_index){
          if(g[terminals[current]].c >= Distance[current][targ] && g[terminals[targ]].c >= Distance[current][targ]){
            if(g[terminals[targ]].c > 0){ D.insert(targ); V.erase(targ);}
            node_labels[targ]=clusterID;
          }
        }
      }
    }

  }

  // Identfying the vertex membership with respect to clusters
  vector<vector <int> > clusters(clusterID+1);
  for(unsigned int i=0; i< node_labels.size(); i++){
    for(int j=0; j<=clusterID; j++){
      if(node_labels[i] == j){
        clusters[j].push_back(i);
      }
    }
  }

  std::vector<int>::iterator it, itt; int num_clusters=0;
  for(unsigned int i=1; i<clusters.size(); i++){
    if(clusters[i].size() > 1) num_clusters++;
  }

  if(num_clusters == 0){
    //cout<<"There is no tree in construct G ()"<<endl;
    //return 0;
  }


  // Regrouping the singletone and dobletone clusters after clustering
  unsigned int threshold_num = 2; int min_index; double min_distance;
  for(unsigned int i=1; i<clusters.size(); i++){
    if( clusters[i].size() <= threshold_num){
      for (it=clusters[i].begin(); it!=clusters[i].end(); ++it){
        min_index = -1; min_distance = DBL_MAX;
        for(unsigned int j=1; j<clusters.size(); j++){
          if( clusters[j].size() > threshold_num){
            for (itt=clusters[j].begin(); itt!=clusters[j].end(); ++itt){
              if(min_distance > Distance[*it][*itt] -g[terminals[*it]].c - g[terminals[*itt]].c ){
                min_distance = Distance[*it][*itt] -g[terminals[*it]].c - g[terminals[*itt]].c; min_index = j;
              }
            }
          }
        }


        if (min_index != -1){
          clusters[min_index].push_back(*it);
          *it = -1;
        }
      }
    }

  }



  // Construct an artificial graph G, which is composed of all clusters determined
  // from Heuristic Clustering phase

  string str; unsigned int index1=-1, index2=-1;

  for (unsigned int l = 0; l<terminals.size(); l++){
    str=to_string(l);
    index1=idx_G(str);
    G[index1].c = g[terminals[l]].c;
    G[index1].name = g[terminals[l]].name;
  }


  for (unsigned int l = 0; l<terminals.size(); l++){
    if(l != root_index){
      str=to_string(l);
      index1=idx_G(str);
      add_edge(root_index, index1, Distance[root_index][index1], G);
    }
  }

  for(unsigned int i = 1; i < clusters.size(); i++){
    for (it=clusters[i].begin(); it!=clusters[i].end(); ++it){
      if(*it != -1){
        str=to_string(*it); index1=idx_G(str);
        for (itt=it+1; itt!=clusters[i].end(); ++itt){
          if(*itt != -1){
            str=to_string(*itt); index2=idx_G(str);
            add_edge(index1, index2, Distance[*it][*itt], G);
          }
        }
      }
    }
  }



  weight_G = get(edge_weight, G);
  vector < Vertex > spanning_tree_G(num_vertices(G));
  prim_minimum_spanning_tree(G, & spanning_tree_G[0]);


  Edge beg; Vertex sour, tar; double cost;
  Edge e; bool found;

  weight_g = get(edge_weight, g);

  vector<int> path; index1=0; index2=0;
  edge_iterator out_i, out_end; int add=0;

  // Solving the Minimum Spanning Tree on G
  for(unsigned int i = 0; i < spanning_tree_G.size(); ++i ){

    if(spanning_tree_G[i]!=i ){

      if(i> spanning_tree_G[i]) path=perPath[spanning_tree_G[i]][i];
      else path=perPath[i][spanning_tree_G[i]];

      for(unsigned int j=0; j<path.size()-1; j++){

        sour=vertex(path[j], g); tar= vertex(path[j+1], g);
        boost::tuples::tie(beg, found) = edge(sour, tar,g);
        cost=get(weight_g, beg);

        index1=idx_G_pruned(to_string(path[j])); index2=idx_G_pruned(to_string(path[j+1]));
        add=0;
        for (boost::tuples::tie(out_i, out_end) = out_edges(vertex(index1,G_pruned), G_pruned); out_i != out_end; ++out_i) {
          if(target(*out_i, G_pruned)==index2) add++;
        }
        if(!add) add_edge(index1, index2, cost, G_pruned);
      }
  }
  }

  weight_G_pruned = get(edge_weight, G_pruned);

  vector< Vertex >spanning_tree_G_pruned(num_vertices(G_pruned));
  prim_minimum_spanning_tree(G_pruned, &spanning_tree_G_pruned[0]);


  double total1=0;
  for(unsigned int i = 0; i < spanning_tree_G_pruned.size(); ++i ){
    if(spanning_tree_G_pruned[i] != i){
      sour= vertex(i,G_pruned); tar=vertex(spanning_tree_G_pruned[i], G_pruned);
      boost::tuples::tie(beg, found) = edge(sour, tar,G_pruned);
      total1+=get(weight_G_pruned, beg);
    }
  }


  return spanning_tree_G_pruned;

  }


// After obtaining MST tree, the algorithm prunes the leaf nodes
// which have prizes smaller than connection cost
double cut(int &Root, vector< Vertex > &spanning_tree_G_pruned,  vector< string > &tree_from,  vector< string > &tree_to,  vector< double > &tree_cost, map < string, double > &tree_terminals){

  weight_G_pruned = get(edge_weight, G_pruned);

  Edge e; bool found;

  Edge beg; Vertex sour, tar; int ancestor=-1; bool select=false;


  int root= -1;
  if(Root == -1){
    double max=0;
    vertex_iterator ei, ef;
    for(tie(ei, ef)= vertices(G_pruned); ei!=ef; ei++){
      if(g[boost::lexical_cast<int>(G_pruned[*ei].name)].c > max){
        root = *ei;
        max = g[boost::lexical_cast<int>(G_pruned[*ei].name)].c;
      }
    }
  } else {root = idx_G_pruned(to_string(Root));}



  select = true;

  bool ancestor_changed=true; unsigned int father, temp; ancestor = root;
  if (select){
    father=spanning_tree_G_pruned[ancestor];
    spanning_tree_G_pruned[ancestor]=ancestor;
    while(ancestor_changed){
      if(spanning_tree_G_pruned[father]==father){
        ancestor_changed=false;
        spanning_tree_G_pruned[father]=ancestor;
      }else{
        temp=spanning_tree_G_pruned[father];
        spanning_tree_G_pruned[father]=ancestor;
        ancestor=father; father=temp;
      }

    }

  }else{
    //cout <<"There is no tree"<<endl;
    return 0;
  }



  vector<Node> predecessor(num_vertices(G_pruned));
  if(select){
    for(unsigned int i = 0; i < spanning_tree_G_pruned.size(); ++i ){
      if(spanning_tree_G_pruned[i]!=i){
        predecessor[i].father=spanning_tree_G_pruned[i];
        predecessor[spanning_tree_G_pruned[i]].children.push_back(i);
      }else{predecessor[i].father=i;}
    }
  }


  for(unsigned int i = 0; i < predecessor.size(); ++i ){
    predecessor[i].size=predecessor[i].children.size();
    predecessor[i].price=g[boost::lexical_cast<int>(G_pruned[i].name)].c;

  }


  price_collect(predecessor, root, e, found);

  process_leafs(predecessor, root, e, found);

  weight_g = get(edge_weight, g);

  // Tree
  for(unsigned int i = 0; i < predecessor.size(); ++i ){
    if(predecessor[i].father != -1 &&  predecessor[i].father != (int) i ){
      sour= vertex(boost::lexical_cast<int>(G_pruned[i].name),g); tar = vertex(boost::lexical_cast<int>(G_pruned[predecessor[i].father].name),g);
      boost::tuples::tie(beg, found) = edge(sour, tar,g);
      tree_from.push_back(g[boost::lexical_cast<int>(G_pruned[i].name)].name);
      tree_to.push_back(g[boost::lexical_cast<int>(G_pruned[predecessor[i].father].name)].name);
      tree_cost.push_back(weight_g[beg]);
      tree_terminals[g[boost::lexical_cast<int>(G_pruned[i].name)].name] = g[boost::lexical_cast<int>(G_pruned[i].name)].c;
      tree_terminals[g[boost::lexical_cast<int>(G_pruned[predecessor[i].father].name)].name] = g[boost::lexical_cast<int>(G_pruned[predecessor[i].father].name)].c;
    }
  }



  double total = 0, lostPrice =0;
  int uncovered_nodes = 0;
  for(unsigned int i = 0; i < predecessor.size(); ++i ){
    if(predecessor[i].father != -1){
      if(predecessor[i].father != (int) i){
        sour= vertex(i,G_pruned); tar = vertex(predecessor[i].father,G_pruned);
        boost::tuples::tie(beg, found) = edge(sour, tar,G_pruned);
        total+=get(weight_G_pruned, beg);
      }
    }
  }


  // Lsit of nodes that are outside of final tree
  vector<int> calculatecost(num_vertices(g));
  for(unsigned int i = 0; i < predecessor.size(); ++i ){
    if(predecessor[i].father != -1 && predecessor[i].father != (int) i){
      sour= vertex(boost::lexical_cast<int>(G_pruned[i].name),g); tar = vertex(boost::lexical_cast<int>(G_pruned[predecessor[i].father].name),g);
      calculatecost[sour]=1; calculatecost[tar]=1;
    }
  }

  // Uncovered nodes
  for(unsigned int i = 0; i < num_vertices(g); ++i ){
    if(calculatecost[i] == 0 && (int) i != root ){
      lostPrice += g[i].c;
      uncovered_nodes++;
    }
  }

  // The list of Nodes in the final Tree
  for(unsigned int i = 0; i < num_vertices(g); ++i ){
    if(calculatecost[i] == 1){
    }
  }

  // Objective value
  return total + lostPrice;

}


//' Internal function \code{call_sr}
//'
//' This function is internally used to solve the PCST.
//'
//' @keywords internal
//'
//' @param from  A \code{CharacterVector} that corresponds to \code{head} nodes of the edges.
//' @param to A \code{CharacterVector} that corresponds the \code{tail} nodes of the edges.
//' @param cost A \code{NumericVector} which represents the edge weights.
//' @param node_names A \code{CharacterVector} demonstrates the names of the nodes.
//' @param node_prizes A \code{NumericVector} which corresponds to the node prizes.
//' @author Murodzhon Akhmedov
//'
// [[Rcpp::export]]
List call_sr(CharacterVector from, CharacterVector to,  NumericVector cost, CharacterVector node_names, NumericVector node_prizes)
{
  clear_variables();

  vector <int> terminals;

  read_input_graph(from, to, cost, node_names, node_prizes);

  Root = idx_g("DUMMY");

  double max_price=0; int max_price_index = -1;
  for(unsigned int i=0; i<num_vertices(g); i++){
    if (g[i].c > max_price){
      max_price = g[i].c;
      max_price_index = i;
    }
  }

  if(Root != -1){
    for(unsigned int i=0; i<num_vertices(g); i++){
      if( (int) i != Root && g[i].c >0){
        terminals.push_back(i);
      }
    }
  }else{
    Root = max_price_index;
    for(unsigned int i=0; i<num_vertices(g); i++){
      if( (int) i != Root && g[i].c >0){
        terminals.push_back(i);
      }
    }
  }

  terminals.push_back(Root);


  if(terminals.size() <=1){
    // There is no tree
    //return 0;
  }


  vector< Vertex > spanning_tree;
  spanning_tree = constructG(terminals, Root);

  vector< string > tree_from;
  vector< string > tree_to;
  vector< double > tree_cost;
  map < string, double > tree_terminals;
  double obj = cut(Root, spanning_tree, tree_from, tree_to, tree_cost, tree_terminals); if (obj == 0) return 0;

  CharacterVector tree_f(tree_from.size());
  CharacterVector tree_t(tree_to.size());
  NumericVector tree_c(tree_cost.size());
  CharacterVector tree_ter(tree_terminals.size());
  NumericVector tree_ter_p(tree_terminals.size());

  for(unsigned int i=0; i<tree_from.size(); i++){
    tree_f[i]=tree_from[i];
    tree_t[i]=tree_to[i];
    tree_c[i]=tree_cost[i];
  }

  int counter = 0;
  for (std::map<string, double>::iterator it=tree_terminals.begin(); it!=tree_terminals.end(); ++it){
    tree_ter[counter] = it->first;
    tree_ter_p[counter] = it->second;
    counter++;
  }

  List tree = List::create(tree_from, tree_to, tree_cost, tree_ter, tree_ter_p);


  return tree;
}






