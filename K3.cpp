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

int main(int argc, char **argv) {

  string graph_file;
  graph_file = argv[1];
  
  FILE *f_freq = fopen("3dir.csv", "a");
 
  Graph *g = new GraphMatrix();

  GraphUtils::readFileTxt(g,graph_file.c_str(),1,0,1);
  vector<pair<int,string>> counts;
  Fase *final_fase = new Fase(g,1);
  final_fase->runCensus(3);
  counts = final_fase->subgraphCount();

  for (int i = 0; i<counts.size(); i++)
    fprintf(f_freq, "%s,%s,%d\n",argv[2],counts[i].second.c_str(), counts[i].first);
    
  return 0;
}
  
