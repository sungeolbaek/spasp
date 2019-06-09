#ifndef GRAPH_H_HHWU
#define GRAPH_H_HHWU

#include <stdio.h>
#include <vector>
#include <string.h>
#include <algorithm>
#include <set>
#include <queue>
#include <iostream>
#include "Timer.h"
#include "io.h"

using namespace std;
typedef int vertex_id;

const int BLK_SZ = 1024*1024;
const int FILE_NUM = 1024;
const int max_level=1023;
const int SZ_PTR=sizeof(void*);
const int LABEL_INNER_NODE_NUMBER= 1024*1024;
const int LABEL_INNER_EDGE_NUMBER= LABEL_INNER_NODE_NUMBER *6;

struct edgepair 
{
  vertex_id node1;
  vertex_id node2;
  int weight;
};

struct offsetInfo 
{
  long buff_offset;
  char level;
};


const int SZ_VERTEX = sizeof(int);
const int VERTEX_PER_BLK = BLK_SZ/SZ_VERTEX;
const int SZ_OFFSET = sizeof(offsetInfo);
const int OFFSET_PER_BLK = BLK_SZ/SZ_OFFSET;
const int SZ_EDGE = sizeof(edgepair);
const int EDGE_PER_BLK = BLK_SZ/SZ_EDGE;

struct noderange
{
  vertex_id nodeid;
  vertex_id degree;
  int begin;
};

bool compare(const noderange a, const noderange b)
{
  if(a.degree == b.degree)
    return a.nodeid < b.nodeid;

  return a.degree < b.degree;
}

bool compare_edge(const edgepair& a, const edgepair& b)
{
  if(a.node1==b.node1)
  {
	if(a.node2 == b.node2)
	  return a.weight < b.weight;
	else
	  return a.node2 < b.node2;

  }
  return a.node1 < b.node1;
}

struct minvertex
{
  vertex_id nodeid;
  int degree;
  int fileindex;
  friend bool operator < (minvertex a, minvertex b)
  {
    if(a.degree == b.degree)
    {
      return a.nodeid > b.nodeid;
	}
	return (a.degree > b.degree); 
  }
};

struct minedge
{
  vertex_id node1;
  vertex_id node2;
  int weight;
  int fileindex;
  friend bool operator < (const minedge &a, const minedge &b)
  {
    if(a.node1 == b.node1)
      return a.node2 > b.node2;
    return a.node1 > b.node1;
  }
};

struct weightedge
{	
  vertex_id nodeid;
  int weight;
};

bool compare_weightedge(const weightedge& a, const weightedge& b)
{
  return a.nodeid < b.nodeid;
}

struct minvid 
{
  vertex_id nodeid;
  int levelindex;
};

bool compare_vid(const minvid& a, const minvid& b)
{
  return a.nodeid < b.nodeid;
}

class Graph
{
public:
  Graph(); 
  void set_vertex_num(int v_num); 
  void set_mem(int mem_size);
  void initial(); 
  void write_graph(int &ptr_node_index, vector <noderange> & indexNode, vector <weightedge> & e, WriteBuffer & write_io); // write the nodes and edges by write_io
  void write_edge(int & ptr_edge_index, vector <edgepair> & edge, WriteBuffer & write_io); // write the edges into file by write_io
  void sortFile(const char* filePath1, const char* filePath2); //sort the graph according to the vertex degree
  void computeIS(const char* filePath1, const char* filePath2, const char* filePath3); //computing independent set, getting the new graph
  
  void labelInit(int &u, int &u_deg, vector <weightedge> & neighbor, WriteBuffer & write_label); // the initial labels when doing independent set
  bool label_loop(int & ptr_node_index, vector <noderange> & indexNode, vector <weightedge> & e, vector < vector <weightedge> > & e2, vector <int> & curr, vector <int> & pos_node_inner, int &ptr_node_index_inner, int &ptr_e_inner, vector <noderange> & indexNode_inner, vector <weightedge> & e_inner, long & label_num); // when the memory for the inner loop is full, updating the labels
  int labelUpdate_inner_loop(int & ptr_node_index, int & ptr_e, vector <noderange> & indexNode, vector <weightedge> & e, vector < vector <weightedge> > & e2, vector <int> & curr, FILE *sortfile, int lev, WriteBuffer &write_level,  WriteBuffer &write_label, long & curr_offset, vector <long> & offset); // when the memory for outer loop is full, call the function to update the labels
  void labelUpdate(int lev, WriteBuffer &write_level, WriteBuffer &write_label, long & curr_offset, vector <long> & offset); // update the labels of level lev
  void offsetUpdate(const char* filePath1, const char* filePath2); //inputpath, offsetfile
  void adjLabel(const char* adjfilePath, const char* adjoffsetPath);
  void topdownLabel(const char* adjfilePath, const char* adjoffsetPath, const char* filePath1, const char* filePath2); //doing top down labeling process
  void labelinttobyte(const char* filePath1, const char* filePath2);
  void print_statistics(const char* filePath); //write out the statistics of the remaining graph
  
  void run(char* filePath, bool flag_c, bool flag_k, int round_number); 

public:
  long V, E, mem;
  long max_node_m, max_edge_m;
  long  max_node_label_m, max_edge_label_m;
  FILE ** level;
  int round;
  long v_num_old, e_num_old;
  long v_num_new, e_num_new;
  long label_number_total;
  double ratio;
  long limit_label_m;
  bool stopis;
  int max_file_num;
  vector <char> levelInfo;
};


Graph::Graph()
{
  mem = 4096;
}


void Graph::set_vertex_num(int v_num)
{
  V = v_num;
}


void Graph::set_mem(int mem_size)
{
  mem = mem_size;
}


void Graph::initial()
{
  stopis= false;
  level = (FILE**)malloc((max_level+1)*SZ_PTR);
  if(level == NULL)
  {
    printf("[%d] malloc fail!\n", __LINE__);
    exit(1);
  }
  
  max_node_m = mem*1024*4;  
  max_edge_m = max_node_m*6;
  
  max_node_label_m= mem*256;
  max_edge_label_m= max_node_label_m*12;
  
  limit_label_m = max_edge_label_m*32;
  
  levelInfo.resize(V+1);
	
  max_file_num=256;
}

void filecopy(FILE *dest, FILE *src)
{
  const int size = 16384; // 16*1024
  char buffer[size];

  while (!feof(src))
  {
    int n = fread(buffer, 1, size, src);
    fwrite(buffer, 1, n, dest);
  }

  fflush(dest);
}

void inputFileCopy(char * dest, char * src)
{
  FILE * infile  = fopen(src,  "rb");
  FILE * outfile = fopen(dest, "wb");

  filecopy(outfile, infile);

  fclose(infile);
  fclose(outfile);
}

// the start point of running
void Graph::run(char* filePath, bool flag_c, bool flag_k, int round_number)
{
  Timer tt;
  tt.start();
  
  initial();
  	
  char SORT_INPUT_NAME[100], REMAIN_GRAPH_NAME[100], HDS_GRAPH_NAME[100];
  char LABEL_FILE_NAME[100], FINAL_LABEL_FILE_NAME[100], OFFSET_FILE_NAME[100];
  char GRAPH_GK_info[100], HDS_OFFSET_FILE_NAME[100];
  char ADJ_LABEL_FILE_NAME[100], ADJ_OFFSET_FILE_NAME[100];
  
  strcpy(SORT_INPUT_NAME,filePath);
  strcat(SORT_INPUT_NAME,".sort");
  strcpy(LABEL_FILE_NAME, filePath);
  strcat(LABEL_FILE_NAME,".label_temp");
  strcpy(FINAL_LABEL_FILE_NAME, filePath);
  strcat(FINAL_LABEL_FILE_NAME,".label");
  strcpy(OFFSET_FILE_NAME, filePath);
  strcat(OFFSET_FILE_NAME,".offset");
  strcpy(GRAPH_GK_info, filePath);
  strcat(GRAPH_GK_info,".info");
  strcpy(ADJ_LABEL_FILE_NAME, filePath);
  strcat(ADJ_LABEL_FILE_NAME,".adjlabel");
  strcpy(ADJ_OFFSET_FILE_NAME, filePath);
  strcat(ADJ_OFFSET_FILE_NAME,".adjoffset");
  
  ratio=0.0;
  round = 0;

  sprintf(HDS_GRAPH_NAME, "%s.%d.hds", filePath, round);
  inputFileCopy(HDS_GRAPH_NAME, filePath);
  sprintf(HDS_OFFSET_FILE_NAME, "%s.%d.hds.offset", filePath, round);	
  offsetUpdate(filePath, HDS_OFFSET_FILE_NAME);
  
  printf("round: %d\n", round+1);
  sortFile(filePath, SORT_INPUT_NAME);
  sprintf(REMAIN_GRAPH_NAME, "%s.%d.hds", filePath, round+1);  
  computeIS(filePath, SORT_INPUT_NAME, REMAIN_GRAPH_NAME);	
  
  sprintf(HDS_OFFSET_FILE_NAME, "%s.%d.hds.offset", filePath, round+1);  
  offsetUpdate(REMAIN_GRAPH_NAME, HDS_OFFSET_FILE_NAME);
  
  if(flag_c)
  { // complete the round
  	for (round=1; ; round++)
    {
      printf("round: %d\n", round+1);
      sortFile(REMAIN_GRAPH_NAME, SORT_INPUT_NAME);
      sprintf(HDS_GRAPH_NAME, "%s.%d.hds", filePath, round+1); 

      computeIS(REMAIN_GRAPH_NAME, SORT_INPUT_NAME, HDS_GRAPH_NAME); 
      if(stopis)
        break;
      
      sprintf(HDS_OFFSET_FILE_NAME, "%s.%d.hds.offset", filePath, round+1);  
      offsetUpdate(HDS_GRAPH_NAME, HDS_OFFSET_FILE_NAME);
      strcpy(REMAIN_GRAPH_NAME, HDS_GRAPH_NAME);
    }
  }
  else if(flag_k)
  { // set the number of levels
    printf("[%d] total round_number: %d\n", __LINE__, round_number);
  	for (round=1; round<round_number; round++)
    {
      printf("round: %d\n", round+1);
      sortFile(REMAIN_GRAPH_NAME, SORT_INPUT_NAME);
      sprintf(HDS_GRAPH_NAME, "%s.%d.hds", filePath, round+1); 

      computeIS(REMAIN_GRAPH_NAME, SORT_INPUT_NAME, HDS_GRAPH_NAME);  
      if(stopis)
        break;

      sprintf(HDS_OFFSET_FILE_NAME, "%s.%d.hds.offset", filePath, round+1);  
      offsetUpdate(HDS_GRAPH_NAME, HDS_OFFSET_FILE_NAME);
      strcpy(REMAIN_GRAPH_NAME, HDS_GRAPH_NAME);
    }
  }
  else 
  { //automatically
  	for (round=1; ratio<=0.95; round++)
    {
      printf("round: %d\n", round+1);
      sortFile(REMAIN_GRAPH_NAME, SORT_INPUT_NAME);
      sprintf(HDS_GRAPH_NAME, "%s.%d.hds", filePath, round+1); 

      computeIS(REMAIN_GRAPH_NAME, SORT_INPUT_NAME, HDS_GRAPH_NAME);  
      if(stopis)
        break;

      sprintf(HDS_OFFSET_FILE_NAME, "%s.%d.hds.offset", filePath, round+1);  
      offsetUpdate(HDS_GRAPH_NAME, HDS_OFFSET_FILE_NAME);
      strcpy(REMAIN_GRAPH_NAME, HDS_GRAPH_NAME);
    }
  }

  printf("<creat topdown label file> ");
  
  topdownLabel(ADJ_LABEL_FILE_NAME, ADJ_OFFSET_FILE_NAME, LABEL_FILE_NAME, OFFSET_FILE_NAME);
  labelinttobyte(LABEL_FILE_NAME, FINAL_LABEL_FILE_NAME);
  print_statistics(GRAPH_GK_info);
  
  remove(SORT_INPUT_NAME);
  remove(LABEL_FILE_NAME);

  free(level);
  tt.stop();
  printf("The total time: %lf\n", tt.GetRuntime());
} 

// write the nodes in indexNode and edges in e by write_io
void Graph::write_graph(int &ptr_node_index, vector <noderange> & indexNode, vector <weightedge> & e, WriteBuffer & write_io)
{
  for(int i = 0; i < ptr_node_index; i++)
  {
    write_io.write(&indexNode[i].nodeid, 1);
    write_io.write(&indexNode[i].degree, 1);

    //printf("u:%d, deg:%d\n", indexNode[i].nodeid, indexNode[i].degree);
    for(int j=0; j<indexNode[i].degree; j++)
    {
      write_io.write(&e[j+indexNode[i].begin].nodeid, 1);
      write_io.write(&e[j+indexNode[i].begin].weight, 1);
      //printf("v:%d, weight:%d\n", 
      //	e[j+indexNode[i].begin].nodeid, e[j+indexNode[i].begin].weight);
	  }
  }
  write_io.flush();
  ptr_node_index = 0;
}

// write the edges in edge into file by write_io
void Graph::write_edge(int & ptr_edge_index, vector <edgepair> & edge, WriteBuffer & write_io)
{	
  edgepair previous;
  previous.node1 = previous.node2 = -1;

  for(int i=0; i<ptr_edge_index; i++)
  {
    if(edge[i].node1 != previous.node1 || edge[i].node2 != previous.node2)
    {
      write_io.write(&edge[i].node1, 1);
      write_io.write(&edge[i].node2, 1);
      write_io.write(&edge[i].weight, 1);
      
      previous.node1=edge[i].node1;
      previous.node2=edge[i].node2;
    }
  }
  write_io.flush();
  ptr_edge_index = 0;
}

//sort the graph according to the vertex degree
void Graph::sortFile(const char* filePath1, const char* filePath2)
{
  Timer tt;
  tt.start();

  FILE* file = fopen(filePath1,"rb");
  if(file == NULL)
  {
    printf("Cannot open source file %d\n", __LINE__);
    exit(1); 
  }
  FILE* sortfile = fopen(filePath2,"wb");
  if(sortfile == NULL)
  {
    printf("Cannot open source file %d\n", __LINE__);
    exit(1); 
  }
	
  int ptr_file_number = 0;
  /* seed[i] is ordered as degree size */
  FILE **seed = (FILE **)malloc(FILE_NUM*SZ_PTR);
  if(seed == NULL)
  {
    printf("[%d] malloc fail!\n", __LINE__);
    exit(1);
  }
	
  vector <noderange> indexNode; 
  vector <weightedge> e;

  indexNode.resize(max_node_m+1);
  e.resize(max_edge_m+1);
	
  int ptr_node_index = 0, ptr_e = 0;
  long node_number_total =0, edge_number_total = 0;

  ReadBuffer read_io(file);	
  int u, u_deg;
  int v, weight;

  /* file format : Vid, degree, (edge(vid1), w, .. ), (edge(vid2), w, .. ) ... */
  while(!read_io.isend)
  {   
    read_io.read(&u);
    read_io.read(&u_deg);
    if(ptr_node_index >= max_node_m || (ptr_e + u_deg) >= max_edge_m)
    {
      //sort according to the node degree 
	    std::sort(indexNode.begin(), indexNode.begin()+ptr_node_index, compare);
      /* Creates and opens a temporary file with unique auto-generated filename.  */
      seed[ptr_file_number] = tmpfile();
      if(seed[ptr_file_number]==NULL)
      {
        printf("ERROR: cannot create temp file\n");
        exit(1);
      }
			
      WriteBuffer write_io(seed[ptr_file_number]);
      
      write_graph(ptr_node_index, indexNode, e, write_io);
      write_io.flush();
      ptr_e = 0;
      ptr_file_number++;
    }
    
    indexNode[ptr_node_index].nodeid = u;
    indexNode[ptr_node_index].degree = u_deg;
    indexNode[ptr_node_index].begin = ptr_e;
    ptr_node_index++;
    
    node_number_total++;
    edge_number_total += u_deg;
    
    //vertex store in indexNode, edge store in e
    for(int i = 0; i < u_deg; i++)
    {
      read_io.read(&v);
      read_io.read(&weight);
      
      e[ptr_e].nodeid=v;
      e[ptr_e].weight=weight;
      ptr_e++;
    }
  }

  /* for the last part of the file. */
	
  if(ptr_node_index > 0)
  { //for the last part of the file
    std::sort(indexNode.begin(), indexNode.begin()+ptr_node_index, compare);
    seed[ptr_file_number] = tmpfile();
    if(seed[ptr_file_number]==NULL)
    {
      printf("ERROR: cannot create temp file\n");
      exit(1);
    }
    WriteBuffer write_io(seed[ptr_file_number]);
    write_graph(ptr_node_index, indexNode, e, write_io);
    write_io.flush();
    ptr_e = 0;
    ptr_file_number++;
  }

  fclose(file);

  /* creating seed[] is done. */
	
  v_num_old = node_number_total;
  e_num_old = edge_number_total;
	
  //merge Sort
  priority_queue<minvertex> a;
  minvertex new_vertex;
  ReadBuffer seed_io[ptr_file_number+1];
  for(int i = 0; i < ptr_file_number; i++)
  {
    seed_io[i].open(seed[i]);
    seed_io[i].read(&new_vertex.nodeid);
    seed_io[i].read(&new_vertex.degree);
    new_vertex.fileindex = i;
    a.push(new_vertex);
  }
	
  WriteBuffer write_io_1(sortfile);
  node_number_total = edge_number_total = 0;

  minvertex minone;
  while(!a.empty())
  {
    minone = a.top();
    a.pop();
    
    write_io_1.write(&minone.nodeid, 1);
    write_io_1.write(&minone.degree, 1);
    
    node_number_total++;
    edge_number_total += minone.degree;

    for(int i = 0; i < minone.degree; i++)
    {
      seed_io[minone.fileindex].read(&v);
      seed_io[minone.fileindex].read(&weight);
    		
      write_io_1.write(&v, 1);
      write_io_1.write(&weight, 1);
    }
		
    if(!seed_io[minone.fileindex].isend)
    {
      seed_io[minone.fileindex].read(&new_vertex.nodeid);
      seed_io[minone.fileindex].read(&new_vertex.degree);
      
      new_vertex.fileindex = minone.fileindex;
      a.push(new_vertex);
    }
  }
  write_io_1.flush();
	
  for(int i = 0; i < ptr_file_number; i++)
  {
    fclose(seed[i]);
  }
	
  free(seed);
  fclose(sortfile);
  tt.stop();
  printf("The time used for sorting: %lf\n");
}

//computing independent set, getting the new graph
void Graph::computeIS(const char* filePath1, const char* filePath2, const char* filePath3)
{
  Timer tt;
  tt.start();
  
  printf("Start computing independent set:");
  FILE* file = fopen(filePath1,"rb");
  if(file == NULL)
  {
    printf("Cannot open source file %d\n", __LINE__);
    exit(1); 
  }
  
  FILE* sortfile = fopen(filePath2,"rb");
  if(sortfile == NULL)
  {
    printf("Cannot open source file %d\n", __LINE__);
    exit(1); 
  }
  
  FILE* remainfile = fopen(filePath3,"wb");
  if(remainfile == NULL)
  {
    printf("Cannot open source file %d\n", __LINE__);
    exit(1); 
  }
  
  vector <bool> mark, deleteable; 
  vector <edgepair> edge;
  vector <weightedge> neighbor;
  
  mark.resize(V+1, false);
  deleteable.resize(V+1, true);
  neighbor.resize(V+1);

  ReadBuffer read_io(sortfile);
  int ptr_edge_index=0, mark_node_count=0;

  FILE ** seed = (FILE**)malloc(FILE_NUM*SZ_PTR);
  if(seed == NULL)
  {
    printf("[%d] malloc fail!\n", __LINE__);
    exit(1);
  }
  
  int ptr_file_number = 0;
  
  int u, u_deg;
  int v, weight;
  edgepair new_edgepair;

  while(!read_io.isend) //according to the sortfile, get the independent set
  {
    read_io.read(&u);
    read_io.read(&u_deg);
    
    for(int i = 0; i < u_deg; i++)
    {
      read_io.read(&neighbor[i].nodeid);
      read_io.read(&neighbor[i].weight);
    }
    
    if(deleteable[u] && !mark[u] )
    {
      if((u_deg-1)*(u_deg-1) + ptr_edge_index >= max_edge_m)
      {
        std::sort(edge.begin(), edge.begin()+ptr_edge_index, compare_edge);
        seed[ptr_file_number] = tmpfile();
        WriteBuffer write_io(seed[ptr_file_number]);
        write_edge(ptr_edge_index, edge, write_io);
        edge.clear();
        ptr_file_number ++;
        
        if(ptr_file_number>=max_file_num)
        {
          fprintf(stderr,"ERROR: Too Many Temp File, out of disk size\n");
          exit(1);
        }
      }
      
      mark[u]=true;
      mark_node_count++;//delete node
      
      for(int i=0; i<u_deg; i++)
      {
        deleteable[neighbor[i].nodeid]=false;
      }
	  
      for(int j=0; j<u_deg; j++)
      {
        for(int k=0; k<u_deg; k++)
        {
          if(neighbor[k].nodeid!=neighbor[j].nodeid)
          {
            new_edgepair.node1=neighbor[j].nodeid;
            new_edgepair.node2=neighbor[k].nodeid;
            new_edgepair.weight=neighbor[j].weight+neighbor[k].weight;
            edge.push_back(new_edgepair);
            
            ptr_edge_index++;
          }
        }
      }
    }
  }
  
  if(ptr_edge_index>0)
  {
    std::sort(edge.begin(), edge.begin()+ptr_edge_index, compare_edge);
    seed[ptr_file_number] = tmpfile();
    if(seed[ptr_file_number]==NULL)
    {
      printf("[%d] cannot create temp file!\n", __LINE__);
      exit(1);
    }
    WriteBuffer write_io(seed[ptr_file_number]);
    write_edge(ptr_edge_index, edge, write_io);
    edge.clear();
    ptr_file_number ++;
    
    if(ptr_file_number>=max_file_num)
    {
      fprintf(stderr,"ERROR: Too Many Temp File, out of disk size\n");
      exit(1);
    }
  }  
  vector<bool>().swap(deleteable);
  fclose(sortfile);
	
  printf("Getting augmenting edges:\n");
  // merge sort of the augmenting edges
  priority_queue<minedge> a;
  minedge new_edge;
  ReadBuffer seed_io[ptr_file_number];
  
  FILE * tempfile = tmpfile();
  if(tempfile == NULL)
  {
    printf("[%d] cannot create temp file!\n", __LINE__);
    exit(1);
  }
  WriteBuffer write_io(tempfile);
  
  for(int i=0;i<ptr_file_number;i++)
  {
  	seed_io[i].open(seed[i]);
  	seed_io[i].read(&new_edge.node1);
  	seed_io[i].read(&new_edge.node2);
  	seed_io[i].read(&new_edge.weight);
  	new_edge.fileindex = i;
  	a.push(new_edge);
  }
  
  ptr_edge_index = 0;
  int min_file_index=0;
  
  minedge minone_edge;
  while (!a.empty())
  {	
  	minone_edge = a.top();
  	a.pop();
  	
  	write_io.write(&minone_edge.node1, 1);
  	write_io.write(&minone_edge.node2, 1);
  	write_io.write(&minone_edge.weight, 1);
  	
  	min_file_index = minone_edge.fileindex;
  
  	if(!seed_io[minone_edge.fileindex].isend)
  	{
  	  seed_io[minone_edge.fileindex].read(&new_edge.node1);
  	  seed_io[minone_edge.fileindex].read(&new_edge.node2);
  	  seed_io[minone_edge.fileindex].read(&new_edge.weight);
  		
  	  new_edge.fileindex = minone_edge.fileindex;
  	  a.push(new_edge);
    }
  }
  write_io.flush();
  
  vector<edgepair>().swap(edge);
  
  for(int i = 0; i < ptr_file_number; i++)
  {
  	fclose(seed[i]);
  }

  free(seed);
  printf("Generating the new graph:\n");
  //scan original graph to get new graph		
  rewind(tempfile);
  WriteBuffer write_io_1(remainfile);
  
  vector <int> neighborlist;
  vector <int> pos;
  neighborlist.resize(V+1, -1);
  pos.resize(V+1, -1);
  
  int ptr_node_index=0, ptr_e=0;
  edgepair * edgebuff = (edgepair *)malloc(BLK_SZ);
  if(edgebuff == NULL)
  {
    printf("[%d] malloc fail!\n", __LINE__);
    exit(1);
  }

  // tempfile has augment edges set.
  int edge_read = fread(edgebuff,SZ_EDGE,EDGE_PER_BLK,tempfile);
  
  int ptr_edge_buff=0, remarknumber=0;
  long node_number_total=0, edge_number_total=0;
  
  ReadBuffer read_io_1(file); // original file or next graph file.
  
  level[round] = tmpfile(); //tmpfile()to store different level
  if(level[round] == NULL)
  {
    printf("[%d] cannot create temp file!\n", __LINE__);
    exit(1);
  }
  
  WriteBuffer write_label(level[round]);
  
  vector <noderange> indexNode;
  indexNode.resize(max_node_m+1);
  
  vector < weightedge > e;
  e.resize(max_edge_m+1);
  
  while(!read_io_1.isend)
  {
  	read_io_1.read(&u);
  	read_io_1.read(&u_deg);

  	if(!mark[u]) //u is not in independent set
  	{
      /* independent set is removed in next graph */
      levelInfo[u] = round+1;
  	  if(ptr_node_index >= max_node_m || (ptr_e+16*u_deg)>= max_edge_m)//need to flush out the edges
      {
        node_number_total+=ptr_node_index;
        for(int i=0; i<ptr_node_index; i++)
        {
          std::sort(e.begin()+indexNode[i].begin, e.begin()+indexNode[i].begin+indexNode[i].degree, compare_weightedge);
          edge_number_total += indexNode[i].degree;
        }
        //printf("[%d]write remain graph for %d\n", __LINE__, u);
        
        write_graph(ptr_node_index, indexNode, e, write_io_1);
        write_io_1.flush();
        ptr_e = 0;
      }
  
      indexNode[ptr_node_index].nodeid = u;
  	  indexNode[ptr_node_index].degree = u_deg;
  	  indexNode[ptr_node_index].begin = ptr_e;
  	  ptr_node_index++;
  
      for(int i = 0; i < u_deg; i++)
  	  {
  	    read_io_1.read(&v);
        read_io_1.read(&weight);
        if(!mark[v])
        {
          neighborlist[v] = u;
          pos[v] = ptr_e;
          e[ptr_e].nodeid = v;
          e[ptr_e].weight = weight;
          ptr_e++;
        }
        else 
          indexNode[ptr_node_index-1].degree--;
  	  }
  	}
  	else
  	{
  	  for(int i=0;i<u_deg;i++)
      {
        read_io_1.read(&neighbor[i].nodeid);
        read_io_1.read(&neighbor[i].weight);
  	  }
  	  remarknumber++;
  	  labelInit(u, u_deg, neighbor, write_label); //the inital labels in level L[round]
  	}

  	if(!mark[u])
  	{
  	  while (!feof(tempfile) || ptr_edge_buff!=edge_read) 
      {
        if(edgebuff[ptr_edge_buff].node1>u) 
          break;
        if(edgebuff[ptr_edge_buff].node1 == u)
        {
          if(neighborlist[edgebuff[ptr_edge_buff].node2] != u)
          {
            neighborlist[edgebuff[ptr_edge_buff].node2] = u;
            indexNode[ptr_node_index-1].degree++;
            e[ptr_e].nodeid=edgebuff[ptr_edge_buff].node2;
            e[ptr_e].weight=edgebuff[ptr_edge_buff].weight;
            pos[edgebuff[ptr_edge_buff].node2]=ptr_e;
            ptr_e++;
            
            if(ptr_e>=max_edge_m)
            {  //need to flush out the edges
              ptr_node_index--;
              int last_ptr_node = ptr_node_index; 
              
              node_number_total+=ptr_node_index;
              for(int i=0; i<ptr_node_index; i++)
              {
                std::sort(e.begin()+indexNode[i].begin, e.begin()+indexNode[i].begin+indexNode[i].degree, compare_weightedge);
                edge_number_total += indexNode[i].degree;
              }
              //printf("[%d]write remain graph for %d\n", __LINE__, u);
              
              write_graph(ptr_node_index, indexNode, e, write_io_1);
              write_io_1.flush();
              ptr_e = 0;
              
              if(indexNode[last_ptr_node].degree >= max_edge_m-1)
              {
                printf("ERROR: Too High Degree Vertex \n");
                exit(1);
              }
              
              for(int i=indexNode[last_ptr_node].begin; i<indexNode[last_ptr_node].begin+indexNode[last_ptr_node].degree; i++)
              {
                e[ptr_e].nodeid=e[i].nodeid;
                e[ptr_e].weight=e[i].weight;
                pos[e[i].nodeid]=ptr_e;
                ptr_e++;
              }
              
              ptr_node_index=1;
            }
          }
          else if(edgebuff[ptr_edge_buff].weight < e[pos[edgebuff[ptr_edge_buff].node2]].weight)
            e[pos[edgebuff[ptr_edge_buff].node2]].weight=edgebuff[ptr_edge_buff].weight;
        }
        ptr_edge_buff++;
        if(ptr_edge_buff==edge_read)
        {
          if(!feof(tempfile)) 
          {
            edge_read = fread(edgebuff,SZ_EDGE,EDGE_PER_BLK,tempfile);
            ptr_edge_buff=0;
          }
        }
      }
  	}
  }

  free(edgebuff);
  
  //for the last part
  if(ptr_node_index > 0)
  {
  	node_number_total+=ptr_node_index;
  	for(int i=0; i<ptr_node_index; i++)
    {
  	  std::sort(e.begin()+indexNode[i].begin, e.begin()+indexNode[i].begin+indexNode[i].degree, compare_weightedge);
  	  edge_number_total += indexNode[i].degree;
  	}
  	//printf("[%d]write remain graph for last part\n", __LINE__);
  	write_graph(ptr_node_index, indexNode, e, write_io_1);
  	write_io_1.flush();
  	ptr_e = 0;
  }
  
  fclose(tempfile);
  fclose(file);
  fclose(remainfile);
  
  if(node_number_total==0)
  	stopis=true;
  else 
  	stopis=false;
  
  printf("The number of deleted nodes:%d\n", remarknumber);
  printf("Edge number in total:%ld\n", edge_number_total);
  printf("Node number in total:%ld\n", node_number_total);
  
  tt.stop();
  printf("The time used for getting independent set and adding augmented edges:%lf\n", tt.GetRuntime());
}

//the inital labels, it is called when computing independent set 
void Graph::labelInit(int &u, int &u_deg, vector <weightedge> & neighbor, WriteBuffer & write_label)
{
  int u_deg_new = u_deg+1; // distance 가 0 인 node u 의 label 추가 위함.
  write_label.write(&u, 1);
  write_label.write(&u_deg_new, 1);

  bool finish = false;
  int dis = 0;
  int x = u;
  for(int i = 0; i < u_deg; i++)
  {
    if(!finish)
	  {
      if(neighbor[i].nodeid > u)
      {
        finish = true;
        write_label.write(&x, 1);
        write_label.write(&dis, 1);
      }
    }
  	write_label.write(&neighbor[i].nodeid, 1);
  	write_label.write(&neighbor[i].weight, 1);
  }
  if(!finish)
  {
    finish = true;
  	write_label.write(&x, 1);
  	write_label.write(&dis, 1);
  }
  write_label.flush();
}

// when the memory for the inner loop is full, updating the labels
bool Graph::label_loop(int & ptr_node_index, vector <noderange> & indexNode, vector <weightedge> & e, 
  vector < vector <weightedge> > & e2, vector <int> & curr, vector <int> & pos_node_inner, 
  int &ptr_node_index_inner, int &ptr_e_inner, vector <noderange> & indexNode_inner, 
  vector <weightedge> & e_inner, long & label_num)
{
  int lastid=indexNode_inner[ptr_node_index_inner-1].nodeid;
  for(int i=0; i<ptr_node_index; i++)
  {      
    while(curr[i]<indexNode[i].begin+indexNode[i].degree) 
    {
      if(e[curr[i]].nodeid > lastid) 
        break;
      if(pos_node_inner[e[curr[i]].nodeid] != -1) 
      { //e[curr[i]].nodeid is found
        int k=indexNode_inner[pos_node_inner[e[curr[i]].nodeid]].begin;
        int end=k+indexNode_inner[pos_node_inner[e[curr[i]].nodeid]].degree;
        int size=(int)e2[i].size();
        for(int j=0; j<size; j++)
        {
          while(k<end)
          {
            if(e_inner[k].nodeid>e2[i][j].nodeid) 
              break;
            if(e_inner[k].nodeid==e2[i][j].nodeid)
            {
              if(e2[i][j].weight>=e_inner[k].weight+e[curr[i]].weight)
                e2[i][j].weight=e_inner[k].weight+e[curr[i]].weight;
            }
            if(e_inner[k].nodeid<e2[i][j].nodeid)
            {  // add new label
              weightedge temp;
              temp.nodeid=e_inner[k].nodeid;
              temp.weight=e_inner[k].weight+e[curr[i]].weight;
              e2[i].push_back(temp);
              
              label_num++;
              if(label_num >= limit_label_m) // memory is full
                return true;
            }
            k++;
          }
          if(k==end) break;
        }
        while(k<end)
        { // add new label
          weightedge temp;
          temp.nodeid=e_inner[k].nodeid;
          temp.weight=e_inner[k].weight+e[curr[i]].weight;
          e2[i].push_back(temp);
          k++;
          
          label_num++;
          if(label_num >= limit_label_m) // memory is full
            return true;
        }
        sort(e2[i].begin(), e2[i].end(), compare_weightedge);
      }
      curr[i]++;
    }
  }

  ptr_node_index_inner=0;
  ptr_e_inner=0;
	
  return false;
}


int Graph::labelUpdate_inner_loop(int & ptr_node_index, int & ptr_e, vector <noderange> & indexNode,
  vector <weightedge> & e, vector < vector <weightedge> > & e2, vector <int> & curr, FILE *sortfile,
  int lev, WriteBuffer &write_level,  WriteBuffer &write_label, long & curr_offset, vector <long> & offset)
{	
  vector <int> pos_node_inner;
  pos_node_inner.resize(V+1, -1);
	
  vector <noderange> indexNode_inner;
  indexNode_inner.resize(LABEL_INNER_NODE_NUMBER+1);	
  vector <weightedge> e_inner;
  e_inner.resize(LABEL_INNER_EDGE_NUMBER+1);
	
  ReadBuffer read_inner(sortfile);
  int ptr_node_index_inner=0, ptr_e_inner=0;
  int u_inner, u_deg_inner;
	
  long label_num=ptr_e;
  bool label_exceed_mem = false;
  while(!read_inner.isend)
  {
    read_inner.read(&u_inner);
    read_inner.read(&u_deg_inner);
    if(ptr_node_index_inner >= LABEL_INNER_NODE_NUMBER || (ptr_e_inner+u_deg_inner)>= LABEL_INNER_EDGE_NUMBER) 
    { //the memory for the inner loop is full
      label_exceed_mem=label_loop(ptr_node_index, indexNode, e, e2, curr, 
        pos_node_inner, ptr_node_index_inner, ptr_e_inner, indexNode_inner, e_inner, label_num);
      if(label_exceed_mem) 
        return 0;	
    }
    
    indexNode_inner[ptr_node_index_inner].nodeid = u_inner;
    indexNode_inner[ptr_node_index_inner].degree = u_deg_inner;
    indexNode_inner[ptr_node_index_inner].begin = ptr_e_inner;
    pos_node_inner[u_inner] = ptr_node_index_inner;
    ptr_node_index_inner++;		
    
    for(int i=0;i<u_deg_inner;i++)
    {
      read_inner.read(&e_inner[ptr_e_inner].nodeid);
      read_inner.read(&e_inner[ptr_e_inner].weight);
      ptr_e_inner++;
    }
  }
  
  if(ptr_node_index_inner>0)
  { //handle the last part
    label_exceed_mem=label_loop(ptr_node_index, indexNode, e, e2, curr, 
	    pos_node_inner, ptr_node_index_inner, ptr_e_inner, indexNode_inner, e_inner, label_num);
    if(label_exceed_mem)
      return 0;
  }

  //write the labels into labelfile and a new file for level[lev]
  int size=0;
  for(int i=0; i<ptr_node_index; i++)
  {  
    size = (int) e2[i].size();
    write_label.write(&indexNode[i].nodeid, 1);
    write_label.write(&size, 1);
    write_level.write(&indexNode[i].nodeid, 1);
    write_level.write(&size, 1);
    
    offset[indexNode[i].nodeid]=curr_offset;
    curr_offset += 8+5*size;  // u:4, u_deg:4, v: 4, weight: 1

    for(int j=0; j<size; j++)
    {
      write_label.write(&e2[i][j].nodeid, 1);
      write_label.write(&e2[i][j].weight, 1);
      write_level.write(&e2[i][j].nodeid, 1);
      write_level.write(&e2[i][j].weight, 1);
    }
  }
  write_label.flush();
  write_level.flush();
  for(int i=0; i<ptr_node_index; i++)
  {
    vector<weightedge>().swap(e2[i]);
  }
	
  int num = ptr_node_index;
  ptr_node_index = 0;
  ptr_e = 0;
	
  return num;
}

// update the labels for level lev
void Graph::labelUpdate(int lev, WriteBuffer &write_level, WriteBuffer &write_label, long & curr_offset, vector <long> & offset)
{
  FILE * sortfile; 
  vector <minvid> vidarray;
  vidarray.resize(round+1);
	
  int u, u_deg;
  int v, weight;
	
  if(lev+1!=round-1) //merge and sort label L[lev+1...round-1] node into vidarrary
  {
    sortfile = tmpfile();
    WriteBuffer write_io(sortfile);
    ReadBuffer read_seed[round+1];

    for(int i=lev+1; i<round; i++) 
    {
      read_seed[i].open(level[i]);
      read_seed[i].read(&vidarray[i].nodeid);
      vidarray[i].levelindex=i;
    }
	
  	sort(vidarray.begin()+lev+1, vidarray.begin()+round, compare_vid);
  
  	int level_size = round-lev-1, min_level_index=0;
  	while(level_size!=0)
  	{
  	  min_level_index=vidarray[lev+1].levelindex;
  	  u = vidarray[lev+1].nodeid;
  	  read_seed[min_level_index].read(&u_deg);
  			
  	  write_io.write(&u, 1);
  	  write_io.write(&u_deg, 1);
  			
  	  for(int i=0; i<u_deg; i++)
      {
        read_seed[min_level_index].read(&v);
        read_seed[min_level_index].read(&weight);
        write_io.write(&v, 1);
        write_io.write(&weight, 1);
  	  }
  
  	  if(!read_seed[min_level_index].isend)
      { //sort
        read_seed[min_level_index].read(&vidarray[lev+1].nodeid);
        vidarray[lev+1].levelindex=min_level_index;			
        sort(vidarray.begin()+lev+1, vidarray.begin()+lev+1+level_size, compare_vid);
  	  }
  	  else 
      {
        level_size--;
        vidarray.erase(vidarray.begin()+lev+1);
  	  }
    }
    write_io.flush();
  }
  else
    sortfile = level[lev+1];

  vector <noderange> indexNode;
  indexNode.resize(max_node_label_m+1);	
  vector <weightedge> e;
  e.resize(max_edge_label_m+1);
  vector < vector <weightedge> > e2;
  e2.resize(max_node_label_m+1);
  int ptr_node_index=0, ptr_e=0;
	
  vector <int> curr;
  curr.resize(max_node_label_m+1);

  ReadBuffer read_lev(level[lev]);
	
  vector <noderange> indexNode_pack;
  indexNode_pack.resize(max_node_label_m+1);	
  vector <weightedge> e_pack;
  e_pack.resize(max_edge_label_m+1); 
	
  int ptr_node_index_pack=0;
	
  int ptr_total=0, ptr_start=0;
  int processed_num=0;
  while(!read_lev.isend)
  {
  	read_lev.read(&u);
  	read_lev.read(&u_deg);
  	if(ptr_node_index >= max_node_label_m || (ptr_e+u_deg)>= max_edge_label_m) 
  	{ // when the memory for outer loop is full
  	  indexNode_pack = indexNode;
  	  e_pack = e;
  	  ptr_node_index_pack = ptr_node_index;
  	  ptr_total=0;
  	  processed_num=0;
  			
  	  while(ptr_total < ptr_node_index_pack)
  	  {
  	    processed_num=labelUpdate_inner_loop(ptr_node_index, ptr_e, indexNode, e, e2, curr, 
  		  sortfile, lev, write_level, write_label, curr_offset, offset);
  				
        if(processed_num==0)
        {  //out of memory, retry
          for(int i=0; i<max_node_label_m; i++)
          {
            vector<weightedge>().swap(e2[i]);
          }
          max_node_label_m = max_node_label_m/2;
    					
          ptr_node_index=0;
          ptr_e=0;
    					
          for(int i=ptr_start; i<ptr_node_index_pack; i++)
          {
            if(ptr_node_index >= max_node_label_m || (ptr_e+indexNode_pack[i].degree)>= max_edge_label_m)
            {
              if(ptr_node_index<=1)
              {
                fprintf(stderr,"ERROR: Too High Degree");
                exit(1);
              } 
              break;	
            }
            indexNode[ptr_node_index].nodeid = indexNode_pack[i].nodeid; 
            indexNode[ptr_node_index].degree = indexNode_pack[i].degree; 
            indexNode[ptr_node_index].begin = ptr_e;
            curr[ptr_node_index] = ptr_e;
            ptr_node_index ++;
      				
            for(int j=indexNode_pack[i].begin; j<indexNode_pack[i].begin+indexNode_pack[i].degree; j++)
            {
              e[ptr_e]=e_pack[j];
              e2[ptr_node_index-1].push_back(e[ptr_e]);
              ptr_e++;
            }
          }
        }
        else 
        {
          ptr_total+=processed_num;
          ptr_start+=processed_num;
    					
          ptr_node_index=0;
          ptr_e=0;
          for(int i=ptr_start; i<ptr_node_index_pack; i++)
          {
            if(ptr_node_index >= max_node_label_m || (ptr_e+indexNode_pack[i].degree)>= max_edge_label_m)
            {
              if(ptr_node_index<=1)
              {
                fprintf(stderr,"ERROR: Too High Degree");
                exit(1);
              } 
              break;	
            }
            indexNode[ptr_node_index].nodeid = indexNode_pack[i].nodeid; 
            indexNode[ptr_node_index].degree = indexNode_pack[i].degree; 
            indexNode[ptr_node_index].begin = ptr_e;
            curr[ptr_node_index] = ptr_e;
            ptr_node_index ++;
    				
            for(int j=indexNode_pack[i].begin; j<indexNode_pack[i].begin+indexNode_pack[i].degree; j++)
            {
              e[ptr_e]=e_pack[j];
              e2[ptr_node_index-1].push_back(e[ptr_e]);
              ptr_e++;
            }
          }
        }
  	  }	
    }
  	
    indexNode[ptr_node_index].nodeid = u;
    indexNode[ptr_node_index].degree = u_deg;
    indexNode[ptr_node_index].begin = ptr_e;
    curr[ptr_node_index] = ptr_e;
    ptr_node_index++;
    for(int i=0;i<u_deg;i++)
    {
      read_lev.read(&e[ptr_e].nodeid);
      read_lev.read(&e[ptr_e].weight);
      e2[ptr_node_index-1].push_back(e[ptr_e]);
      ptr_e++;
    }
  }
  
  if(ptr_node_index>0)
  { //handle the last part
    indexNode_pack = indexNode;
    e_pack = e;
    ptr_node_index_pack = ptr_node_index;
    ptr_total=0;
    processed_num=0;
		
    while(ptr_total < ptr_node_index_pack)
    {
      processed_num=labelUpdate_inner_loop(ptr_node_index, ptr_e, indexNode, e, e2, curr, 
      sortfile, lev, write_level, write_label, curr_offset, offset);
      if(processed_num==0)
      {  //out of memory, retry
        for(int i=0; i<max_node_label_m; i++)
        {
          vector<weightedge>().swap(e2[i]);
        }
        max_node_label_m = max_node_label_m/2;
        
        ptr_node_index=0;
        ptr_e=0;
        
        for(int i=ptr_start; i<ptr_node_index_pack; i++)
        {
          if(ptr_node_index >= max_node_label_m || (ptr_e+indexNode_pack[i].degree)>= max_edge_label_m)
          {
            if(ptr_node_index<=1)
            {
              fprintf(stderr,"ERROR: Too High Degree");
              exit(1);
            } 
            break;	
          }
          indexNode[ptr_node_index].nodeid = indexNode_pack[i].nodeid; 
          indexNode[ptr_node_index].degree = indexNode_pack[i].degree; 
          indexNode[ptr_node_index].begin = ptr_e;
          curr[ptr_node_index] = ptr_e;
          ptr_node_index ++;
          
          for(int j=indexNode_pack[i].begin; j<indexNode_pack[i].begin+indexNode_pack[i].degree; j++)
          {
            e[ptr_e]=e_pack[j];
            e2[ptr_node_index-1].push_back(e[ptr_e]);
            ptr_e++;
          }
        }
      }
      else 
      {
        ptr_total+=processed_num;
        ptr_start+=processed_num;
        
        ptr_node_index=0;
        ptr_e=0;
        for(int i=ptr_start; i<ptr_node_index_pack; i++)
        {
          if(ptr_node_index >= max_node_label_m || (ptr_e+indexNode_pack[i].degree)>= max_edge_label_m)
          {
            if(ptr_node_index<=1)
            {
              fprintf(stderr,"ERROR: Too High Degree");
              exit(1);
            }
            break;	
          }
          indexNode[ptr_node_index].nodeid = indexNode_pack[i].nodeid; 
          indexNode[ptr_node_index].degree = indexNode_pack[i].degree; 
          indexNode[ptr_node_index].begin = ptr_e;
          curr[ptr_node_index] = ptr_e;
          ptr_node_index ++;
          
          for(int j=indexNode_pack[i].begin; j<indexNode_pack[i].begin+indexNode_pack[i].degree; j++)
          {
            e[ptr_e]=e_pack[j];
            e2[ptr_node_index-1].push_back(e[ptr_e]);
            ptr_e++;
          }
        }
      }
    }	
  }
	
  if(lev+1!=round-1)
    fclose(sortfile);
}

/* creating offset file */
void Graph::offsetUpdate(const char* filePath1, const char* filePath2) //inputpath, offsetfile 
{
  Timer tt;
  tt.start();
  
  printf("Start computing offsetUpdate:\n");
  vector <long> offset;
  offset.resize(V+1, -1);
  
  long curr_offset=0;
  
  FILE* levelGraphFile = fopen(filePath1,"rb");
  if(levelGraphFile == NULL)
  {
    printf("Cannot open source file %d\n", __LINE__);
    exit(1); 
  }
  
  ReadBuffer read_io(levelGraphFile);
  int u, u_deg;
  int v, weight;
  
  while(!read_io.isend)
  {
    read_io.read(&u);
    read_io.read(&u_deg);
    
    offset[u]=curr_offset;
    //curr_offset += 8+5*u_deg; /* (u:4, u_deg:4) + (v: 4, weight: 1) * u_deg */
    curr_offset += 2+2*u_deg; /* (u:1, u_deg:1) + (v: 1, weight: 1) * u_deg */
    
    //printf("[%d] u:%d, u_deg:%d\n", __LINE__, u, u_deg); 
    
    for(int i=0; i<u_deg; i++)
    {
      read_io.read(&v);
      read_io.read(&weight);
    }
  }
  
  fclose(levelGraphFile);
  
  //record the label offset and level info in offset file 
  FILE * offsetFile = fopen(filePath2,"wb");
  if(offsetFile == NULL)
  {
    printf("Cannot open source file %d\n", __LINE__);
    exit(1); 
  }
  
  offsetInfo * write_buff_offsetInfo = (offsetInfo *) calloc(SZ_OFFSET, BLK_SZ);
  if(write_buff_offsetInfo == NULL)
  {
    printf("[%d] malloc fail!\n", __LINE__);
    exit(1);
  }
  
  int ptr_write_buff_offset = 0;
  for(int i=0; i<V; i++)
  {
    write_buff_offsetInfo[ptr_write_buff_offset].buff_offset = offset[i];
    write_buff_offsetInfo[ptr_write_buff_offset].level = levelInfo[i];
    
    //printf("[%d] vId:%d, Lv:%d, offset:%d\n", __LINE__, 
    //  i, write_buff_offsetInfo[ptr_write_buff_offset].level,
    //  write_buff_offsetInfo[ptr_write_buff_offset].buff_offset);
    
    ptr_write_buff_offset++;
    
    if(ptr_write_buff_offset==OFFSET_PER_BLK)
    {
      fwrite(write_buff_offsetInfo,SZ_OFFSET,ptr_write_buff_offset,offsetFile);
      ptr_write_buff_offset = 0;
    }
  }
  
  if(ptr_write_buff_offset>0)
  {
    fwrite(write_buff_offsetInfo,SZ_OFFSET,ptr_write_buff_offset,offsetFile);
    ptr_write_buff_offset = 0;
  }		
  
  fclose(offsetFile);
  free(write_buff_offsetInfo);
  
  tt.stop();
  printf("The time used for labeling:%lf\n", tt.GetRuntime());
}

void Graph::adjLabel(const char* adjfilePath, const char* adjoffsetPath)
{
  vector <long> adjOffset;
  adjOffset.resize(V+1, -1); 

  long curr_adjOffset=0;
  int u, u_deg;
  int v, weight;

  FILE *adjLabelFile = fopen(adjfilePath,"wb");
  if(adjLabelFile == NULL)
  {
    printf("Cannot open source file %d\n", __LINE__);
    exit(1); 
  }

  FILE *adjOffsetFile = fopen(adjoffsetPath,"wb");
  if(adjOffsetFile == NULL)
  {
    printf("Cannot open source file %d\n", __LINE__);
    exit(1); 
  }

  WriteBuffer write_label(adjLabelFile);
  WriteBuffer write_offset(adjOffsetFile);

  for(int i=round-1; i>=0; i--)
  {
    ReadBuffer read_io(level[i]);
    
  	while(!read_io.isend)
    {
      read_io.read(&u);
      read_io.read(&u_deg);
      
      adjOffset[u]=curr_adjOffset;
      curr_adjOffset += 8 + 8*u_deg; /* (u:4, u_deg:4) + (v: 4, weight: 4) * u_deg */
      
      write_label.write(&u, 1);
      write_label.write(&u_deg, 1);
      
      //printf("[%d] u:%d, u_deg:%d\n", __LINE__, u, u_deg);
      
      for(int i=0; i<u_deg; i++)
      {
        read_io.read(&v);
        read_io.read(&weight);
        write_label.write(&v, 1);
        write_label.write(&weight, 1);
        
        //printf(" (u:%d, w:%d)", v, weight);
      }
      //printf("\n");
    }
  	write_label.flush();
  }

  fclose(adjLabelFile);

  offsetInfo * buff_adjOffsetInfo = (offsetInfo *) calloc(SZ_OFFSET, V);
  if(buff_adjOffsetInfo == NULL)
  {
    printf("[%d] malloc fail!\n", __LINE__);
    exit(1);
  }

  //printf("[%d] print ADJ offset.\n", __LINE__);
  
  int pBuff_adjOffset = 0;
  for(int i=0; i<V; i++)
  {
    buff_adjOffsetInfo[pBuff_adjOffset].buff_offset = adjOffset[i];
    buff_adjOffsetInfo[pBuff_adjOffset].level = levelInfo[i];
    
    #if 0
    printf("[%d] vId:%d, Lv:%d, offset:%d\n", __LINE__, 
    	i, buff_adjOffsetInfo[pBuff_adjOffset].level,
    	buff_adjOffsetInfo[pBuff_adjOffset].buff_offset);
    #endif
	
    pBuff_adjOffset++;    
  }

  if(pBuff_adjOffset>0)
  {
    fwrite(buff_adjOffsetInfo,SZ_OFFSET, pBuff_adjOffset, adjOffsetFile);
    pBuff_adjOffset = 0;
  }		
  
  fclose(adjOffsetFile);
}

//doing top down labeling 
void Graph::topdownLabel(const char* adjfilePath, const char* adjoffsetPath, 
  const char* filePath1, const char* filePath2) //adjlabelfile, oadjffsetfile, labelfile, offsetfile
{
  Timer tt;
  tt.start();
  
  /* create adj label for EIED */
  adjLabel(adjfilePath, adjoffsetPath);
  
  vector <long> offset;
  offset.resize(V+1, -1);
  
  long curr_offset=0;
  
  ReadBuffer read_io(level[round-1]); //level[h]->labelFile
  FILE * labelFile = fopen(filePath1,"wb");
  if(labelFile == NULL)
  {
    printf("Cannot open source file %d\n", __LINE__);
    exit(1); 
  }
  
  WriteBuffer write_label(labelFile);
  
  int u, u_deg;
  int v, weight;

  while(!read_io.isend)
  {
    read_io.read(&u);
    read_io.read(&u_deg);
    
    offset[u]=curr_offset;
    curr_offset += 8+5*u_deg; /* (u:4, u_deg:4) + (v: 4, weight: 1) * u_deg */
    
    write_label.write(&u, 1);
    write_label.write(&u_deg, 1);
    
    //printf("[%d] u:%d, u_deg:%d\n", __LINE__, u, u_deg);
    
    for(int i=0; i<u_deg; i++)
    {
      read_io.read(&v);
      read_io.read(&weight);
      write_label.write(&v, 1);
      write_label.write(&weight, 1);
    }
  }
  write_label.flush();
  
  FILE * new_levelfile;
  if(round >= 2)
  {	
    for(int i=round-2; i>=0; i--)
    { // from h-1 to 1
      new_levelfile = tmpfile();
      WriteBuffer write_level(new_levelfile);
      printf("update label i:%d\n");
      /* update the labels for level i */
      labelUpdate(i, write_level, write_label, curr_offset, offset); 
      fclose(level[i]);
      level[i] = new_levelfile;
    }
    for(int i=round-1; i>=0; i--)
    {
      fclose(level[i]);
    }
  }
	
  write_label.flush();
  fclose(labelFile);

  //record the label offset and level info in offset file 
  FILE * offsetFile = fopen(filePath2,"wb");
  if(offsetFile == NULL)
  {
    printf("Cannot open source file %d\n", __LINE__);
    exit(1); 
  }

  offsetInfo * write_buff_offsetInfo = (offsetInfo *) calloc(SZ_OFFSET, BLK_SZ);
  if(write_buff_offsetInfo == NULL)
  {
    printf("[%d] malloc fail!\n", __LINE__);
    exit(1);
  }

  //printf("[%d] print full label offset.\n", __LINE__);
	
  int ptr_write_buff_offset = 0;
  for(int i=0; i<V; i++)
  {
    write_buff_offsetInfo[ptr_write_buff_offset].buff_offset = offset[i];
    write_buff_offsetInfo[ptr_write_buff_offset].level = levelInfo[i];
    
    //printf("[%d] vId:%d, Lv:%d, offset:%d\n", __LINE__, 
    //	i, write_buff_offsetInfo[ptr_write_buff_offset].level,
    //	write_buff_offsetInfo[ptr_write_buff_offset].buff_offset);
    
    ptr_write_buff_offset++;
    
    if(ptr_write_buff_offset==OFFSET_PER_BLK)
    {
      fwrite(write_buff_offsetInfo,SZ_OFFSET,ptr_write_buff_offset,offsetFile);
      ptr_write_buff_offset = 0;
    }
  }
	
  if(ptr_write_buff_offset>0)
  {
    fwrite(write_buff_offsetInfo,SZ_OFFSET,ptr_write_buff_offset,offsetFile);
    ptr_write_buff_offset = 0;
  }		
  
  fclose(offsetFile);
  free(write_buff_offsetInfo);
  
  tt.stop();
  printf("The time used for labeling:%lf\n", tt.GetRuntime());
}

//change the label file from int to byte
void Graph::labelinttobyte(const char* filePath1, const char* filePath2) //labelfile, finallabelfile
{
  FILE * intfile = fopen(filePath1,"rb");
  if(intfile == NULL)
  {
    printf("Cannot open source file %d\n", __LINE__);
    exit(1); 
  }
  
  FILE * bytefile = fopen(filePath2,"wb");
  if(bytefile == NULL)
  {
    printf("Cannot open source file %d\n", __LINE__);
    exit(1); 
  }
  
  label_number_total=0;
  
  char * write_buff = (char *)malloc(BLK_SZ); 
  if(write_buff == NULL)
  {
    printf("[%d] malloc fail!\n", __LINE__);
    exit(1);
  }
  
  char * read_buff_1 = (char *)malloc(BLK_SZ);
  if(read_buff_1 == NULL)
  {
    printf("[%d] malloc fail!\n", __LINE__);
    exit(1);
  }
  
  vertex_id * buff = (vertex_id *) read_buff_1;
  
  int num_read = fread(buff,SZ_VERTEX,VERTEX_PER_BLK,intfile);
  int ptr_buff = 0;
  
  int num_write=0;
  while(!feof(intfile) || ptr_buff!=num_read)
  {
    int u = buff[ptr_buff];
    ptr_buff++;
    if(ptr_buff==num_read)
    {
      num_read = fread(buff,SZ_VERTEX,VERTEX_PER_BLK,intfile);
      ptr_buff = 0;
    }
    
    //printf("[%d] vId=%d, ", __LINE__, u);
    memcpy(write_buff+num_write, &u, sizeof(u));
    num_write += sizeof(u);
    if(num_write > BLK_SZ-sizeof(int))
    {
      fwrite(write_buff,1,num_write,bytefile);
      num_write = 0;
    }
    
    int u_deg = buff[ptr_buff];
    ptr_buff++;
    if(ptr_buff==num_read)
    {
      num_read = fread(buff,SZ_VERTEX,VERTEX_PER_BLK,intfile);
      ptr_buff = 0;
    }
    
    //printf("deg=%d : ", u_deg);
    memcpy(write_buff+num_write, &u_deg, sizeof(u_deg));
    num_write += sizeof(u_deg);
    if(num_write > BLK_SZ-sizeof(int))
    {
      fwrite(write_buff,1,num_write,bytefile);
      num_write = 0;
    }
    
    label_number_total += u_deg;
    
    for(int j=0; j<u_deg; j++)
    {
      int v = buff[ptr_buff];
      ptr_buff++;
      if(ptr_buff==num_read)
      {
        num_read = fread(buff,SZ_VERTEX,VERTEX_PER_BLK,intfile);
        ptr_buff = 0;
      }
      
      //printf("(%d, ", v);
      memcpy(write_buff+num_write, &v, sizeof(v));
      num_write += sizeof(v);
      if(num_write > BLK_SZ-sizeof(char))
      {
        fwrite(write_buff,1,num_write,bytefile);
        num_write = 0;
      }
      
      char distance = (char) buff[ptr_buff];
      ptr_buff++;
      if(ptr_buff==num_read)
      {
        num_read = fread(buff,SZ_VERTEX,VERTEX_PER_BLK,intfile);
        ptr_buff = 0; 
      }
      
      //printf("%d) ", distance);
      memcpy(write_buff+num_write, &distance, sizeof(distance));
      num_write += sizeof(distance);
      if(num_write > BLK_SZ-sizeof(int))
      {
        fwrite(write_buff,1,num_write,bytefile);
        num_write = 0;
      }
    }	
    //printf("\n");
  }
		
  if(num_write>0)
  {
    fwrite(write_buff,1,num_write,bytefile);
    num_write = 0;
  }
  
  fclose(bytefile);
  fclose(intfile);
  free(write_buff);
  free(read_buff_1);
}

//print out the statistics for remaing graph
void Graph::print_statistics(const char* filePath)
{
  FILE * graph_gk = fopen(filePath, "w");
  if(graph_gk == NULL)
  {
    printf("Cannot open source file %d\n", __LINE__);
    exit(1); 
  }
  
  fprintf(graph_gk, "%ld %ld %ld %d\n", V, v_num_new, e_num_new, round);
  fclose(graph_gk);
  
  printf("[%d] V:%ld, v_num_new:%ld, e_num_new:%ld, round:%d\n", __LINE__, 
    V, v_num_new, e_num_new, round);
}

#endif

