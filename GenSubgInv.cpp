#include "GraphMatrix.h"
#include "Common.h"
#include "Timer.h"
#include <cmath>
#include <sys/types.h>
#include <unistd.h>
#include "GraphUtils.h"
#include "DynamicGraph.h"
#include "Fase.h"
#include "Random.h"
#include <vector>

using namespace std;

#define SWITCH_FAC 10
#define ITER_FACTOR 1
#define T_INIT 0.001
#define T_FACTOR 0.9
#define TEMP_LIMIT -1000.0
#define E_THRES 0.0001 // 0.1% dif count maximum
#define S1 0
#define T1 1
#define S2 2
#define T2 3

bool directed = 0;
vector<pair<int,string>> real_count;
vector<pair<int,string>> random_count;
int subg_size = 3;
double e_threshold = E_THRES;
double t_init = T_INIT;
double time_count_acc = 0.0;
int switch_factor = SWITCH_FAC;
bool use_er = false;        // use classical ER to generate a random graph
bool keep_deg_seq = true;   // use find_3 instead of find_4
bool switching = true;      // random graph based on markov-chain switching a real network
Fase *f_random;
Fase *f_real;
int t_factor = 0;
double cool_factor = T_FACTOR;
double acc_time_count = 0.0;
int tries_max = 15;
long n_iters = 0;

FILE *ofile;

/***
Print graph2:

Print in the same format as input graph (1 to N)

***/

void print_graph2(Graph *g) {
  for (int i = 0; i<g->numNodes(); i++) {
    vector<int> neighs = g->outEdges(i)[0];
    for (int j = 0; j<neighs.size(); j++) {
      if (directed || neighs.at(j) > i)
	fprintf("%d %d\n",i+1,neighs.at(j)+1);
    }
  }
}

/***
Print graph:

For each node v in graph g, takes list l of the out edges of v and prints:
v connected to: l

***/

void print_graph(Graph *g) {
  for (int i = 0; i<g->numNodes(); i++) {
    vector<int> neighs = g->outEdges(i)[0];
    cout << i << " connected to: ";
    for (int j = 0; j<neighs.size(); j++) {
      cout << neighs.at(j) << " ";
    }
    cout << endl;
  }
}


/***
Print count:

Input: An array of pairs<string,int>, representing subgraphs and their frequency

Prints the array given as input.
***/

void print_subgraphs(vector<pair<int,string>> counts) {

  for (int i = 0; i<counts.size(); i++)
    fprintf(f_freq, "%s %d\n",counts[i].second.c_str(), counts[i].first);
  
}

/***
Copy graph:

For each edge (u,v) in graph source, creates edge (u,v) in graph target

***/

void copy_graph(Graph *source, Graph *target) {

  for (int i = 0; i<source->numNodes(); i++) {
    vector<int> neighs = source->outEdges(i)[0];
    for (int j = 0; j<neighs.size(); j++) {
      target->addEdge(i,neighs[j]);
    }
  }
  
}


/***
Find 3:

Find 3 nodes (s1, t1, s2) such that:
- s1 != t1, t1 != s2, s1 != s2
- s1 has edge to t1
- s2 does not have dge to t1

Changes connection s1->t1 to s2->t1
***/

int find_3(Graph *g, int *nodes) {
  
  int s1,s2,t1;
  int n = g->numNodes();

  s1 = rand()%n;
  while (g->outEdges(s1)[0].size() == 0)
    s1 = rand()%n;

  t1= g->outEdges(s1)[0].at(rand()%g->outEdges(s1)[0].size());

  s2 = rand()%n;
  while (s2 == s1 || s2 == t1 || g->hasEdge(s2,t1))
    s2 = rand()%n;

  nodes[S1] = s1;
  nodes[T1] = t1;
  nodes[S2] = s2;
  return 1;
}

/***
Find 4:

Input: Graph g, array of int nodes
Output: changed array of int nodes, 0 if unsuccessful, 1 if successful

Finds 4 nodes (s1, t1, s2, t2) in g, such that:
- s1 has edge to t1
- s2 has edge to t2
- s1, t1, s2, t2 all different
- s1 does not have edge to t2
- s2 does not have edge to t1

Changes connection s1->t1 and s2->t2 to s1->t2 and s2->t1

Nodes are returned in the array of int nodes:
- s1 in position 0 (macro'ed as S1) 
- t1 in position 1 (macro'ed as T1)
- s2 in position 2 (macro'ed as S2)
- t2 in position 3 (macro'ed as T2)

***/

int find_4(Graph *g, int* nodes) {

  int s1,s2,t1,t2;
  int numNodes = g->numNodes();
  int s1_out_edges, s2_out_edges, t1_out_edges;
  
  /* Find first source (s1):
     - Must have at least 1 neighbour to generate the first target (t1)
     - Cannot be linked to all other nodes, to generate the second target (t2)
   */
  s1 = rand()%numNodes;
  s1_out_edges=g->outEdges(s1)[0].size();
  while (s1_out_edges == 0 ||
	 s1_out_edges >= numNodes-1) {
    s1 = rand()%numNodes;
    s1_out_edges=g->outEdges(s1)[0].size();
  }
  
  /*
    Find first target (t1):
    - Cannot be linked with all other nodes, to generate the second source (s2)
    - If s1 does not have a neighbour that respects the first condition, we fail and retry 
   */
  t1= g->outEdges(s1)[0].at(rand()%s1_out_edges);
  int while_breaker = 0;
  t1_out_edges=g->outEdges(t1)[0].size();
  while (t1_out_edges >= numNodes - 1 && while_breaker < 2*s1_out_edges) {
    t1= g->outEdges(s1)[0].at(rand()%s1_out_edges);
    t1_out_edges=g->outEdges(t1)[0].size();
    while_breaker++;
  }

  if (while_breaker == 2*s1_out_edges) {
    return 0;
  }
  
  /*
    Find second source (s2):
    - Must be different from s1
    - Must be different from t1 and cannot be linked with it
    - Must have at least 1 neighbour to generate second target (t2)
    - If it has only 1 neighbour, this neighbour must be different from s1
 */  
  s2 = rand() % numNodes;
  s2_out_edges = g->outEdges(s2)[0].size();
  while_breaker=0;
  while ((s2 == s1 ||
	 s2 == t1 ||
	 g->hasEdge(s2,t1) ||
	 s2_out_edges == 0 ||
	(s2_out_edges == 1 && g->outEdges(s2)[0].at(0) == s1)) && while_breaker < 4*numNodes) {
    s2 = rand() % numNodes;
    s2_out_edges = g->outEdges(s2)[0].size();
    while_breaker++;
  }

  if (while_breaker == 4*numNodes)
    return 0;
  
  /*
    Find second target (t2):
    - Must be different from s1 and cannot be linked with it
    - If s2 does not have a neighbour that respects the first condition, we fail and retry 
   */
  while_breaker = 0;
  t2 = g->outEdges(s2)[0].at(rand()%s2_out_edges);
  while ((t2 == s1 || g->hasEdge(s1,t2)) &&
	 while_breaker < 2*s2_out_edges) {
    t2 = g->outEdges(s2)[0].at(rand()%s2_out_edges);
    while_breaker++;
  }

  if (while_breaker == 2*s2_out_edges) {
    return 0;
  }

  nodes[S1] = s1;
  nodes[T1] = t1;
  nodes[S2] = s2;
  nodes[T2] = t2;

  return 1;
}

/***
Randomize graph:

Implements a Markov-Chain Algorithm of edge swapping in Graph g.
Performs exactly number_of_switches swaps in G, using edges found through the Find 4 function.
Maintains degree sequence.
***/

void randomize_graph_deg_seq(Graph *g, int number_of_switches) {
  int k=0;
  int *nodes = new int[4];
  while (k<number_of_switches) {
    if (find_4(g, nodes)) {
      g->rmEdge(nodes[S1],nodes[T1]);
      g->rmEdge(nodes[S2],nodes[T2]);
      g->addEdge(nodes[S1],nodes[T2]);
      g->addEdge(nodes[S2],nodes[T1]);
      if (!directed) {
	g->rmEdge(nodes[T1],nodes[S1]);
	g->rmEdge(nodes[T2],nodes[S2]);
	g->addEdge(nodes[T2],nodes[S1]);
	g->addEdge(nodes[T1],nodes[S2]);
      }
      k++;
    }
  }
}

/***
Randomize graph:

Implements a Markov-Chain Algorithm of edge swapping in Graph g.
Performs exactly number_of_switches swaps in G, using edges found through the Find 3 function.
Does not maintain degree sequence.
***/

void randomize_graph_3(Graph *g, int number_of_switches) {
  int k=0;
  int *nodes = new int[3];
  while (k<number_of_switches) {
    if (find_3(g, nodes)) {
      g->rmEdge(nodes[S1],nodes[T1]);
      g->addEdge(nodes[S2],nodes[T1]);
      if (!directed) {
	g->rmEdge(nodes[T1],nodes[S1]);
	g->addEdge(nodes[T1],nodes[S2]);
      }
      k++;
    }
  }
}

/***

Energy calculation based on mfinder (modified 2):

(Sum(k in K) [ abs(V_real,k - V_random,k)/(V_real,k + V_random,k) ]) / |K|

Calculates the average of the deviation of each type of subgraph.

***/

double energy_m2(bool p, vector<pair<int,string>> new_count, vector<pair<int,string>> r_count) {

  double totenergy = 0.0;
  double subE=0.0, addE,div;
  int i=0,j=0;
  int diff_subgs = 0;
  while (i < new_count.size() && j < r_count.size()) {
    if (new_count[i].second.compare(r_count[j].second) == 0) {
      addE = (double)(new_count[i].first + r_count[j].first);
      if (addE != 0.0) {
	subE = (double)(abs(r_count[j].first - new_count[i].first));
	div = subE/addE;
	totenergy += div;
	diff_subgs++;
      }
      i++;
      j++;
    }
    else {
      if (new_count[i].second.compare(r_count[j].second) < 0) {
	if (new_count[i].second.length() > 0 && new_count[i].first > 0) {
	  totenergy += 1.0;
	  diff_subgs++;
	}
	i++;
      }
      else {
	if (r_count[j].second.length() > 0 && r_count[j].first > 0) {
	  totenergy += 1.0;
	  diff_subgs++;
	}
	j++;
      }
    }
  }
  while (i < new_count.size() && new_count[i].first > 0) {
    i++;
    totenergy += 1.0;
    diff_subgs++;
  }
  while (j < r_count.size() && r_count[j].first > 0) {
    j++;
    totenergy += 1.0;
    diff_subgs++;
  }
  return totenergy/diff_subgs;
}


/***

Metropolis acceptance probability.
Oracle:
- Returns "Accept" if new state is closer to the target.
- Returns "Reject" if temperature is too low.
- Calculates bolzmann factor, and takes a random value. Returns "Accept" if random value is lesser than the factor, "Reject" otherwise.

***/

int metropolis(double deltaE, double t) {

  double bolz;
  
  bolz = exp(-deltaE/t);
  if (deltaE < 0.0)
    return 1;
  if (-(deltaE/t) < TEMP_LIMIT) {
    return 0;
  }
  else {
    double randvar = (double)rand()/RAND_MAX;
    return randvar < bolz;
  }
}


/***
Add edge and remove edge with fase on the fly counting
***/

void addEdge(Graph *g, int a, int b) {
  
  f_random->countAddEdge(a,b);
  g->addEdge(a,b);
  
}

void rmEdge(Graph *g, int a, int b) {
  
  f_random->countRemoveEdge(a,b);
  g->rmEdge(a,b);
  
}


void change_state(Graph *g, int *nodes) {
  
  if (keep_deg_seq) {
    while (!find_4(g, nodes));
    rmEdge(g,nodes[S1],nodes[T1]);
    rmEdge(g,nodes[S2],nodes[T2]);
    addEdge(g,nodes[S1],nodes[T2]);
    addEdge(g,nodes[S2],nodes[T1]);;
    
    if (!directed) {
      rmEdge(g,nodes[T1],nodes[S1]);
      rmEdge(g,nodes[T2],nodes[S2]);
      addEdge(g,nodes[T2],nodes[S1]);
      addEdge(g,nodes[T1],nodes[S2]);
    }
  }
  else {
    find_3(g, nodes);
    
    rmEdge(g,nodes[S1],nodes[T1]);
    addEdge(g,nodes[S2],nodes[T1]);
    
    if (!directed) {
      rmEdge(g,nodes[T1],nodes[S1]);
      addEdge(g,nodes[T1],nodes[S2]);
    }
  } 
}

void change_state_back(Graph *g, int *nodes) {

  if (keep_deg_seq) {  
    addEdge(g,nodes[S1],nodes[T1]);
    addEdge(g,nodes[S2],nodes[T2]);
    rmEdge(g,nodes[S1],nodes[T2]);
    rmEdge(g,nodes[S2],nodes[T1]);
    if (!directed) {
      addEdge(g,nodes[T1],nodes[S1]);
      addEdge(g,nodes[T2],nodes[S2]);
      rmEdge(g,nodes[T2],nodes[S1]);
      rmEdge(g,nodes[T1],nodes[S2]);
    }
  }
  else {
    addEdge(g,nodes[S1],nodes[T1]);
    rmEdge(g,nodes[S2],nodes[T1]);
    if (!directed) {
      addEdge(g,nodes[T1],nodes[S1]);
      rmEdge(g,nodes[T1],nodes[S2]);
    }
  }
}
  
/***

Main function to apply random switches to the graph.

Based on a Metropolis Monte Carlo Markov chain.

***/

bool metropolis_switches(Graph *random_g, int ms, int mss) {

  double t = t_init*t_factor;
  double e1, e2, deltaE;
  int succ_count = 0, swi_count = 0;
  int *nodes = new int[4];
  int change;
  int n = 0;
  double last_e;
  e1 = energy_m2(false, random_count, real_count);
  
  while (e1 > e_threshold) {

    if (succ_count > mss || swi_count > ms) {
      succ_count = 0;
      swi_count = 0;
      t = t*cool_factor;
      if (last_e == e1)
	n++;
      else {
	last_e = e1;
	n=0;
      }
      if (n >= 8) {
        return true;
      }
    }

    swi_count++;

    change_state(random_g, nodes);

    random_count = f_random->subgraphCount();

    e2 = energy_m2(false,random_count,real_count);
    deltaE = e2-e1;
    change=metropolis(deltaE,t);
    if (change) {
      e1 = e2;
      succ_count++;
    }
    else
      change_state_back(random_g, nodes);
  }
  return false;
}



/***

Main function.

To distribute:
- Create a .h, only this function set as public
- Set invariants (deg seq, subgraphs of size k)
- Set methodology (random graph to originate from switches to the real network or using well known null models (er, kronecker), switching factor, enery threshold, initial temperature, temperature decreasing factor)

***/


void generate_random_graph(Graph *real_g, Graph **random_g, int dir, int sbg) {

  Graph *new_g;
  int number_of_switches;
  int max_success_changes, max_changes;
  int i,j;
  directed = dir;
  subg_size = sbg;
  bool c;
  /* Initialization of real_count array:
     - Takes real_g
     - Uses FaSE to count subgraph occurences
  */
  
  f_real = new Fase(real_g,directed);
  f_real->runCensus(subg_size);
  real_count = f_real->subgraphCount();
  
  /* Initialization of the new graph: */  

  do {
    t_factor++;
    
    new_g = new GraphMatrix();
    new_g->createGraph(real_g->numNodes(), real_g->type());
    
    number_of_switches=switch_factor*real_g->numEdges();
    
    if (switching) {
      copy_graph(real_g, new_g);
      
      /* Randomizing graph, using a Markov Chain Algorithm:
	 - Number of switches = factor * nEdges + E(X)
	 - X rand var, uniform dist between 0 and factor * nEdges
	 - factor default = 10
	 - Number of switches <= 2 * factor * nEdges
	 - factor is tunable: [1,inf[
      */
  
      if (keep_deg_seq)
	randomize_graph_deg_seq(new_g,number_of_switches);
      else
	randomize_graph_3(new_g,number_of_switches);
    }
    else {
      if (use_er) { 
	while (new_g->numEdges() < real_g->numEdges()) {
	  i = rand() % new_g->numNodes();
	  j = rand() % new_g->numNodes();
	  if (i!=j && !(new_g->hasEdge(i,j))) {
	    new_g->addEdge(i,j);
	    if (!directed)
	      new_g->addEdge(j,i);
	  }
	}
      }
    }
 
    /* 
       Initialization of random_count array
    */
    
    f_random = new Fase(new_g,directed);
    
    f_random->runCensus(subg_size);
    
    random_count = f_random->subgraphCount();
    
    /* Metropolis Monte-Carlo Algorithm. 
       - Temperature change breakpoints depend on number of switches*/
    
    max_changes=number_of_switches;
    max_success_changes=max_changes/10;
    c=metropolis_switches(new_g,max_changes,max_success_changes);
  } while (c && t_factor < tries_max);
  *random_g = new_g;
}

void print_help() {
  cout << "----------- GENK USAGE ----------" << endl;
  cout << "Main settings: ./genk -s <Subgraph Size> -g <input graph file> -o <output graph file> [arguments]\n\n\n";
  cout << "All commands:\n";
  cout << "-h: Displays this menu, ignores all other commands.\n";
  cout << "-s: Sets the subgraph size whose frequencies to preserve.\n";
  cout << "-g: input graph file with subgraph counts to replicate. Format should be a list of edges, one edge per line, with origin and destination nodes separated by a single space. The nodes should be numbered from 1 to the total number of nodes.\n";
  cout << "-o: outfile graph file with the random graph generated by the method.\n";
  cout << "-d: Directed graph (Default is undirected).\n";
  cout << "-w: Weighted graph (Default is unweighted). Weighted graphs are simply treated as unweighted by our method.\n";
  cout << "-t: Initial temperature for the simulated annealing procedure (Default is 0.001).\n";
  cout << "-e: Energy threshold to stop the simulated annealing procedure (Default is 0.0001).\n";
  cout << "-switches: Number of switches per iteration. This is a factor, f, multiplied by the number of edges in the graph, so the number of switches is f * |E|. Integer number is expected (Default is 1).\n";
  cout << "-cool: Cooling factor for the simulated annealing procedure (Default is 0.9).\n";
  cout << "-maxtries: Number of tries to approximate the energy under the specified threshold (Default is 15).\n";
  cout << "-degseq: Generates a random graph with a different degree sequence that the original one.\n";
  cout << "-er: Uses a Erdos-Renyi as the original random graph.\n";
  cout << "-seed: Initializes a random seed (Default is getpid()).\n";
  cout << "---------------------------------" << endl;
}
  
  

int main(int argc, char **argv) {

  int dir=0;
  int sgs;
  string graph_file;
  string outfile("Results.txt");
  int i;
  bool weigh = false;
  double r_seed = getpid();
  for (i = 1; i<argc; i++) {
    if (!strcmp("-s",argv[i]))
      sgs = atoi(argv[++i]);
    if (!strcmp("-d",argv[i]))
      dir = 1;
    if (!strcmp("-g",argv[i]))
      graph_file = argv[++i];
    if (!strcmp("-t",argv[i]))
      t_init = atof(argv[++i]);
    if (!strcmp("-e",argv[i]))
      e_threshold = atof(argv[++i]);
    if (!strcmp("-switches",argv[i])) {
      switch_factor = atoi(argv[++i]);
      switch_factor = switch_factor < 1 ? 1 : switch_factor;
    }
    if (!strcmp("-degseq",argv[i]))
      keep_deg_seq = false;
    if (!strcmp("-er",argv[i])) {
      use_er = true;
      switching = false;
    }
    if (!strcmp("-w",argv[i]))
      weigh=true;
    if (!strcmp("-cool",argv[i]))
      cool_factor=atof(argv[++i]);
   
    if (!strcmp("-maxtries",argv[i]))
      tries_max = atoi(argv[++i]);
    if (!strcmp("-o",argv[i]))
      out_file=argv[++i];
    if (!strcmp("-h",argv[i])) {
      print_help();
      return 0;
    }
    if (!strcmp("-seed",argv[i]))
      r_seed = atof(argv[++i]);
  }

  srand(r_seed);
  Random::init(r_seed);

  ofile = fopen(outfile.c_str(),"w");
  
  sgs = sgs < 3? 3:sgs;
 
  Graph *real_g = new GraphMatrix();
  Graph *random_g;

  GraphUtils::readFileTxt(real_g,graph_file.c_str(),dir,weigh,1);

  generate_random_graph(real_g, &random_g, dir, sgs);

  print_graph2(random_g);
  
  return 0;
}
  
