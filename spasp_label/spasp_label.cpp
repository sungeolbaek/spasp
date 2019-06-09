#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include "graph.h"
using namespace std;

static void usage_desc() 
{
  printf("\nUsage:\n");
  printf("	spasp_label filename vertex_number [-c] [-k num_level] [-m memory]\n");
  printf("Description:\n");
  printf("	filename	The graph file.\n");
  printf("	vertex_number	The number of vertices of the graph.\n");
  printf("	-c	Complete all the levels to make G_k to be empty. \n");
  printf("	-k	Set num_level to be the number of levels k for k-level hierarchy.\n");
  printf("	-m	Set memory to be the size of memory (in MB) to be used. The default value is 4096. \n");
}

int main(int argc, char* argv[])
{
  if (argc == 1 || argc == 2) 
  {
	usage_desc();
	return 1;
  }
	
  int i = 3;
  bool flag_c=false;
  bool flag_k=false;
  int mem = 4096;
  int round_number; 

  while (i < argc) 
  {
	if (strcmp("-c", argv[i]) == 0) 
    {
      i++;
      flag_c=true;
    }
    else if (strcmp("-k", argv[i]) == 0)
    {
      i++;
      flag_k = true;
      round_number = atoi(argv[i++])-1;
    }
    else if (strcmp("-m", argv[i]) == 0)
    {
      i++;
      mem = atoi(argv[i++]);
    } 
  }
	
  Graph g;
  g.set_mem(mem);
  g.set_vertex_num(atoi(argv[2]));
  g.run(argv[1], flag_c, flag_k, round_number);

  return 0;
}
