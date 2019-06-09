#include <iostream>
#include <ctime>
#include <vector>
#include <algorithm>
#include <queue>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "io.h"
#include "timecounter.h"

#define DEBUG 0
#define DEBUG_DIJKSTR 0

using namespace std;

typedef int vertex_id;
const long INF=0xFFFF; // 0xFFFF to denote infinity

struct IVCInfo //intermediate vertex candidates info
{
  int nodeId;
  char level;
};

struct compare_IVC
{
  bool operator() ( const IVCInfo& a, const  IVCInfo& b) const
  { return a.level < b.level; }
};

struct offsetInfo 
{
  long buff_offset;
  char level;
};

struct candiNodeInfo 
{
  int s_u_dist; /* distance between start & node. */
  int u_t_dist; /* distance between node & target. */
};

struct seedNodeInfo 
{
  int nodeId;
  int l_s; /* label_s index */
  int l_t; /* label_t index */
};

const int BLK_SZ = 1024*64;  //1024*64 sizeof blk for I/O
const int SZ_VERTEX = sizeof(int);
const int VERTEX_PER_BLK = BLK_SZ/SZ_VERTEX;
const int SZ_OFFSET = sizeof(offsetInfo);
const int OFFSET_PER_BLK = BLK_SZ/SZ_OFFSET;

FILE * labelFile; // store the label content
FILE * offsetFile; //store the position of label of each vertex
FILE * remainFile; // store the remaining graph G_k

/* each level info of vids  */
char * levelInfo;
/* candidate node set */
candiNodeInfo* candiNodeSet;
/* seed node set */
seedNodeInfo* seedNodeSet = NULL;
int seedNodeSet_ptr = 0;

/* maximum level among leveled graphs */
char max_level, tmp_max_level;  

long * offset;
int nGraph; /* the size of G */
int min_sum, old_min_sum, min_sum_id, s_u_dist, u_t_dist;
int max_node_m; // node number of G_i 
int max_edge_m; // edge number of G_i 
bool flag_label; //true: only use label; false: label + graphX

offsetInfo* buff_offsetInfo;
char* read_buff_2;

struct edgepair 
{
  vertex_id node1;
  vertex_id node2;
  int weight;
};

struct weightedge 
{ 
  vertex_id nodeid;
  char weight;
};

struct noderange
{
  vertex_id nodeid;
  vertex_id degree;
  int begin;
};

struct compare_weightedge
{
  bool operator() ( const weightedge& a, const  weightedge& b) const
  { 
    return a.weight > b.weight;
  }
};

int s_deg=0, t_deg=0;
weightedge *label_s, *label_t;

/* graphX */
vector< vector<pair<int, char> > > graphX; 

int inter_node = -1;

/* final result edge set */
vector <edgepair> v_resultSet;
/* seed edge set */
vector <edgepair> v_seedEdge;

char SVA_GRAPH_LABEL[100], SVA_GRAPH_OFFSET[100];
FILE *SVA_LABEL_FP, *SVA_OFFSET_FP;
offsetInfo* tmp_buff_offsetInfo;

long long checkSize(char * fileName)
{
  FILE * pFile = fopen (fileName,"rb");
  fseek (pFile, 0, SEEK_END); 
  long long ans =ftell(pFile);
  fclose(pFile);
  return ans;
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

void read_vertex(int & u, int & u_deg, int & num_read, int &ptr_buff, vertex_id * buff, FILE * file)
{
  u = buff[ptr_buff];
  ptr_buff++;
  if(ptr_buff==num_read)
  {
    if(!feof(file))
    {
      num_read = fread(buff,SZ_VERTEX,VERTEX_PER_BLK,file);
      ptr_buff = 0;
    }
  }
  
  u_deg = buff[ptr_buff];
  ptr_buff++;
  if(ptr_buff==num_read)
  {
    if(!feof(file))
    {
      num_read = fread(buff,SZ_VERTEX,VERTEX_PER_BLK,file);
      ptr_buff = 0;
    }
  }
}

void read_Gi()
{
  for(int i=0; i<nGraph; i++)
  {
    offset[i]=-1;
  }
  
  rewind(offsetFile);
  int num_read = fread(buff_offsetInfo,SZ_OFFSET,OFFSET_PER_BLK,offsetFile);
  int ptr_buff = 0;
  int ptr_offset = 0;
  
  while(!feof(offsetFile) || ptr_buff != num_read)
  {
    offset[ptr_offset]=buff_offsetInfo[ptr_buff].buff_offset;

    //printf("[%d] vId:%d, offset:%ld\n", __LINE__, ptr_offset, offset[ptr_offset]);
	  levelInfo[ptr_offset++]=buff_offsetInfo[ptr_buff++].level;
	
    if(ptr_buff==num_read)
    {
      if(!feof(offsetFile))
      {
        num_read = fread(buff_offsetInfo,SZ_OFFSET,OFFSET_PER_BLK,offsetFile);
        ptr_buff = 0;
      }
    }
  }
  fclose(offsetFile);
  
  if(!flag_label)
  {
    rewind(remainFile);
 
    vertex_id* buff_remain = (vertex_id *) buff_offsetInfo;
    int ptr_buff_remain = 0; 
    int num_read_remain = fread(buff_remain,SZ_VERTEX,VERTEX_PER_BLK,remainFile);
    
    while(!feof(remainFile) || ptr_buff_remain != num_read_remain)
    {
      int u;
      int u_deg;
      read_vertex(u, u_deg, num_read_remain, ptr_buff_remain, buff_remain, remainFile);
      
      for(int i=0;i<u_deg;i++)
	  {
        int v;
        int weight;
        read_vertex(v, weight, num_read_remain, ptr_buff_remain, buff_remain, remainFile);
       
        /* construct graph i edges for graphX that is built later */
        graphX[u].push_back(make_pair(v, weight));
      }
    }

    fclose(remainFile);
  }
}

double total_query_time;
Runtimecounter rt;	
double query_t;
double largest_time, cur_t_time;

void read_label(int v, int &deg, weightedge *label)
{ 
  char id[10];
  int ptr=0;
  /* covert char buffer from offsetInfo */
  char *buff = (char *) buff_offsetInfo;

  fseeko64(labelFile, offset[v], SEEK_SET);
  int num_read = fread(buff,1,BLK_SZ,labelFile);

  ptr+=sizeof(int);

  /* read node v's degree */
  memcpy(&deg, buff+ptr, sizeof(int));
  ptr+=sizeof(int);
  
  for(int i=0; i<deg; i++)
  {
    memcpy(id, buff+ptr, sizeof(char));
    ptr+=sizeof(char);
    
    if(ptr==num_read)
    {
      num_read = fread(buff,1,BLK_SZ,labelFile);
      ptr = 0;
    }
  
    memcpy(id+sizeof(char), buff+ptr, sizeof(char));
    ptr+=sizeof(char);
    if(ptr==num_read)
    {
      num_read = fread(buff,1,BLK_SZ,labelFile);
      ptr = 0;
    }
  
    memcpy(id+2*sizeof(char), buff+ptr, sizeof(char));
    ptr+=sizeof(char);
    if(ptr==num_read)
    {
      num_read = fread(buff,1,BLK_SZ,labelFile);
      ptr = 0;
    }
    
    memcpy(id+3*sizeof(char), buff+ptr, sizeof(char));
    ptr+=sizeof(char);
    if(ptr==num_read)
    {
      num_read = fread(buff,1,BLK_SZ,labelFile);
      ptr = 0;
    }
    
    memcpy(&label[i].nodeid, id, sizeof(int));
    
    memcpy(&label[i].weight, buff+ptr, sizeof(char));

    //printf("[%d](%d, %d, %d)\n", __LINE__, v, label[i].nodeid, label[i].weight);
    
    ptr+=sizeof(char);
   
    if(ptr==num_read)
    {
      num_read = fread(buff,1,BLK_SZ,labelFile);
      ptr = 0;
    }
  } 
}

void findAllpath (int t, vector< vector<weightedge> > &ASP_ary, vector<edgepair> &v_stPathset)
{
  edgepair tempPair;
  
  if(t == -1)
  { //we reached the initial pId -1 of start node
    return;
  }

  for(int i = 0; i < ASP_ary[t].size(); i++)
  {
    //set v_stPathset using ASP_ary[t][i]
    if(ASP_ary[t][i].nodeid != -1)
    {
      tempPair.node1 = ASP_ary[t][i].nodeid;
      tempPair.node2 = t;
      tempPair.weight = ASP_ary[t][i].weight;

      v_stPathset.push_back(tempPair);
    }
    
    findAllpath(ASP_ary[t][i].nodeid, ASP_ary, v_stPathset);
  }
}

/* remove duplicate pair */
void removeDuplicate(vector<edgepair> &v_set)
{ 
  edgepair tmpPair;
  std::vector <edgepair>::iterator iter;
  
  //printf("\n<<[%d] Start remove Duplicate>>\n", __LINE__);
  
  std::sort(v_set.begin(), v_set.end(), compare_edge);
  
  tmpPair = v_set[0];

  for(iter = v_set.begin()+1; iter != v_set.end(); )
  {
    if((tmpPair.node1 == (*iter).node1) && (tmpPair.node2 == (*iter).node2)) 
    {
      iter = v_set.erase(iter);
    }
    else
    {
      tmpPair = *iter;
      iter++;
    }
  }
}

void calc_label(int s, int t)
{
  int ptr_s=0, ptr_t=0;
  char id_s[10], id_t[10];

  old_min_sum = min_sum=INF; //0xFFFF denote they are not reached
  seedNodeSet_ptr = 0;
  s_u_dist = u_t_dist = 0;
  inter_node = -1;
  s_deg = t_deg = 0;
  
  if(offset[s] >= 0) // s is not in Gx
  { 
    /* change offsetInfo into char buff */
    char * buff_s = (char *) buff_offsetInfo;
    long offset_s=offset[s]; 

    fseeko64(labelFile,offset_s,SEEK_SET);
    /* read block */
    int num_read_s=fread(buff_s,1,BLK_SZ,labelFile);
    ptr_s=0;

    ptr_s+=sizeof(int);

    /* read deg of node s */
    memcpy(&s_deg, buff_s+ptr_s, sizeof(int));
    ptr_s+=sizeof(int);
    
    for(int i=0; i<s_deg; i++)
    {
      memcpy(id_s, buff_s+ptr_s, sizeof(char));
      ptr_s+=sizeof(char);
      
      if(ptr_s==num_read_s)
      {
        num_read_s=fread(buff_s,1,BLK_SZ,labelFile);
        ptr_s = 0;
      }
    
      memcpy(id_s+sizeof(char), buff_s+ptr_s, sizeof(char));
      ptr_s+=sizeof(char);
      if(ptr_s==num_read_s)
      {
        num_read_s=fread(buff_s,1,BLK_SZ,labelFile);
        ptr_s = 0;
      }
    
      memcpy(id_s+2*sizeof(char), buff_s+ptr_s, sizeof(char));
      ptr_s+=sizeof(char);
      if(ptr_s==num_read_s)
      {
        num_read_s=fread(buff_s,1,BLK_SZ,labelFile);
        ptr_s = 0;
      }
      
      memcpy(id_s+3*sizeof(char), buff_s+ptr_s, sizeof(char));
      ptr_s+=sizeof(char);
      if(ptr_s==num_read_s)
      {
        num_read_s=fread(buff_s,1,BLK_SZ,labelFile);
        ptr_s = 0;
      }
      
      memcpy(&label_s[i].nodeid, id_s, sizeof(int));
      
      memcpy(&label_s[i].weight, buff_s+ptr_s, sizeof(char));

      //printf("[%d](%d, %d, %d)\n", __LINE__, s, label_s[i].nodeid, label_s[i].weight);
      
      /* in case of the node wiht the highest level */
      if(label_s[i].nodeid == t)
	  {
	    if(min_sum > label_s[i].weight)
	    {
	      old_min_sum = min_sum = label_s[i].weight;
          min_sum_id = label_s[i].nodeid;
          s_u_dist = min_sum;
        }
      }  
      ptr_s+=sizeof(char);
      
      if(ptr_s==num_read_s)
      {
        num_read_s=fread(buff_s,1,BLK_SZ,labelFile);
        ptr_s = 0;
      }
    }
  }
  
  if(offset[t] >= 0) // t is not in Gx
  { 
    char * buff_t = (char *) read_buff_2;
    long offset_t=offset[t]; 
    
    fseeko64(labelFile,offset_t,SEEK_SET);
    int num_read_t=fread(buff_t,1,BLK_SZ,labelFile);
    ptr_t=0;
    
    ptr_t+=sizeof(int);
    
    memcpy(&t_deg, buff_t+ptr_t, sizeof(int));
    ptr_t+=sizeof(int);
    
    for(int i=0; i<t_deg; i++)
    {
      memcpy(id_t, buff_t+ptr_t, sizeof(char));
      ptr_t+=sizeof(char);
      
      if(ptr_t==num_read_t)
      {
        num_read_t=fread(buff_t,1,BLK_SZ,labelFile);
        ptr_t = 0;
      }
      
      memcpy(id_t+sizeof(char), buff_t+ptr_t, sizeof(char));
      ptr_t+=sizeof(char);
      if(ptr_t==num_read_t)
      {
        num_read_t=fread(buff_t,1,BLK_SZ,labelFile);
        ptr_t = 0;
      }
      
      memcpy(id_t+2*sizeof(char), buff_t+ptr_t, sizeof(char));
      ptr_t+=sizeof(char);
      if(ptr_t==num_read_t)
      {
        num_read_t=fread(buff_t,1,BLK_SZ,labelFile);
        ptr_t = 0;
      }
      
      memcpy(id_t+3*sizeof(char), buff_t+ptr_t, sizeof(char));
      ptr_t+=sizeof(char);
      if(ptr_t==num_read_t)
      {
        num_read_t=fread(buff_t,1,BLK_SZ,labelFile);
        ptr_t = 0;
      }
      
      memcpy(&label_t[i].nodeid, id_t, sizeof(int));
      
      memcpy(&label_t[i].weight, buff_t+ptr_t, sizeof(char));

      //printf("[%d](%d, %d, %d)\n", __LINE__, t, label_t[i].nodeid, label_t[i].weight);
      
      if(label_t[i].nodeid == s)
	  {
	    if(min_sum > label_t[i].weight)
	    {
	      old_min_sum = min_sum = label_t[i].weight;
          min_sum_id = label_t[i].nodeid;
          u_t_dist = min_sum;
        }
      }  
      ptr_t+=sizeof(char);
      
      if(ptr_t==num_read_t)
      {
        num_read_t=fread(buff_t,1,BLK_SZ,labelFile);
        ptr_t = 0;
      }
    }
  }

  ptr_s=0;
  ptr_t=0;
  
  /* entry node ids were orderd */
  for(ptr_s=0; ptr_s<s_deg; ptr_s++)
  {
    while(ptr_t<t_deg)
    {
      if(label_t[ptr_t].nodeid > label_s[ptr_s].nodeid)
        break;
      if(label_t[ptr_t].nodeid == label_s[ptr_s].nodeid)
      {
        /* update minimum distance */
        if(label_t[ptr_t].weight + label_s[ptr_s].weight <= min_sum)
        {
          min_sum = label_t[ptr_t].weight + label_s[ptr_s].weight;
  
          if(old_min_sum > min_sum)
          {
            /* reset count for seed nodes */
            seedNodeSet_ptr = 0;
            old_min_sum = min_sum;
          }
          
          seedNodeSet[seedNodeSet_ptr].nodeId = label_t[ptr_t].nodeid;
          seedNodeSet[seedNodeSet_ptr].l_s = ptr_s;
          seedNodeSet[seedNodeSet_ptr++].l_t = ptr_t;
        }
      }
      ptr_t++;
    }
    
    if(ptr_t==t_deg)
    {
      break;
    }
  }
  
  if(min_sum < INF)
  {
    if((0 == s_deg) || (0 == t_deg))
    {
      seedNodeSet[seedNodeSet_ptr++].nodeId = min_sum_id;
      candiNodeSet[min_sum_id].s_u_dist = s_u_dist;
      candiNodeSet[min_sum_id].u_t_dist = u_t_dist;
    }
    else 
    {
      for(int i=0; i<seedNodeSet_ptr; i++)
      {
        /* add candidiate node */
        candiNodeSet[seedNodeSet[i].nodeId].s_u_dist = label_s[seedNodeSet[i].l_s].weight;
        candiNodeSet[seedNodeSet[i].nodeId].u_t_dist = label_t[seedNodeSet[i].l_t].weight;
      }
    }
  }
}

void searchGraphX(int s, int t, vector<edgepair> &v_stPathset)
{
  #if DEBUG == 2
  printf("[%d] searchGraphX .. \n", __LINE__);
  #endif
  int tmpId=0;
  char tmpWeight;
  
  //read label
  if(offset[s] >= 0)
  {
    read_label(s, s_deg, label_s);
    
    for(int i=0; i<s_deg; i++)
    {
      //the vertex is in the remaining graph
      if(offset[label_s[i].nodeid] < 0)
      {
        tmpId = label_s[i].nodeid;
        tmpWeight = label_s[i].weight;
        
        #if DEBUG_DIJKSTR
        printf("[%d] pq_s.push (%d, %d, %d)\n", __LINE__, s, label_s[i].nodeid, label_s[i].weight);
        #endif
       
        graphX[s].push_back(make_pair(tmpId, tmpWeight));
        /* update graphX */
        graphX[tmpId].push_back(make_pair(s, tmpWeight));
      }
    }
  }

  if(offset[t] >= 0)
  {
    read_label(t, t_deg, label_t);

    for(int i=0; i<t_deg; i++)
    {
      //the vertex is in the remaining graph
      if(offset[label_t[i].nodeid] < 0)
      {
        tmpId = label_t[i].nodeid;
        tmpWeight = label_t[i].weight;
        
        #if DEBUG_DIJKSTR
        printf("[%d] pq_s.push (%d, %d, %d)\n", __LINE__, s, label_t[i].nodeid, label_t[i].weight);
        #endif
       
        graphX[t].push_back(make_pair(tmpId, tmpWeight));
        /* update graphX */
        graphX[tmpId].push_back(make_pair(t, tmpWeight));
      }
    }
  }

#if 0
  for (int i = 0; i < nGraph; i++)
  {
    printf("[%d] graphX [%d] : ", __LINE__, i);
    for (int j = 0; j < graphX[i].size(); j++)
    {
      printf("(%d, %d) ", graphX[i][j].first, graphX[i][j].second);
    }
    printf("\n");
  }
#endif
 
  vector<int> dist(nGraph, 0xFFFF);  
  priority_queue< pair<int, int> > pq;

  vector< vector<weightedge> > ASP_ary; /* store ASP edges */  
  weightedge tmpEdge;
 
  ASP_ary.resize(nGraph+1);
  tmpEdge.nodeid = -1;
  tmpEdge.weight = -1;
  
  ASP_ary[s].push_back(tmpEdge);
  dist[s] = 0; /* distance table */

  pq.push(make_pair(0, s)); 
  
  while (!pq.empty())
  {
    int cost = -pq.top().first; 
    int node = pq.top().second;

    pq.pop();

    //printf("[%d] node: %d\n", __LINE__, node);

    /* if target node is found */
    if(node == t)
    {
      //printf("[%d] min_sum : %d --> updated! --> min_sum : %d\n", __LINE__, min_sum, cost);
      min_sum = cost;
      break;
    }

    /* check adj nodes */
    for (int i = 0; i < graphX[node].size(); i++)
    {
      int v = graphX[node][i].first;
      int nextDist = cost + graphX[node][i].second;

      //printf("[%d] v:%d, weight:%d\n", __LINE__, v, weight);
      
      if (dist[v] > nextDist)
      {
        dist[v] = nextDist;

        tmpEdge.nodeid = node;
        tmpEdge.weight = graphX[node][i].second;
          
        ASP_ary[v].clear();
        ASP_ary[v].push_back(tmpEdge);
        pq.push(make_pair(-nextDist, v));
          
      }
      else if (dist[v] == nextDist)
      {
        tmpEdge.nodeid = node;
        tmpEdge.weight = graphX[node][i].second;
        
        ASP_ary[v].push_back(tmpEdge);     
      }
    }
  }

  findAllpath(t, ASP_ary, v_stPathset);
  /* where we remove duplicate edges */
  removeDuplicate(v_stPathset);

  /* re-init graphX */
  {
    if(offset[s] >= 0) 
    {
      for(int i=0; i<s_deg; i++)
      {
        //the vertex is in the remaining graph
        if(offset[label_s[i].nodeid] < 0)
        {
          tmpId = label_s[i].nodeid;
  
          graphX[s].erase(graphX[s].end());
          graphX[tmpId].erase(graphX[tmpId].end());
        }
      }
    }
  
    if(offset[t] >= 0)
    {
      for(int i=0; i<t_deg; i++)
      {
        //the vertex is in the remaining graph
        if(offset[label_t[i].nodeid] < 0)
        {
          tmpId = label_t[i].nodeid;
          
          graphX[t].erase(graphX[t].end());
          graphX[tmpId].erase(graphX[tmpId].end());
        }
      }
    }
  }
  
#if 0 //verify re-init result
  printf("[%d] re-init grapnX\n", __LINE__);
  for (int i = 0; i < nGraph; i++)
  {
    printf("[%d] graphX [%d] : ", __LINE__, i);
    for (int j = 0; j < graphX[i].size(); j++)
    {
      printf("(%d, %d) ", graphX[i][j].first, graphX[i][j].second);
    }
    printf("\n");
  }
#endif

#if 0 //print ASP on graphX
  for (int i = 0; i < v_stPathset.size(); i++)
  {
    printf("[%d] ASP (%d, %d, %d)\n", __LINE__, 
      v_stPathset[i].node1, v_stPathset[i].node2, v_stPathset[i].weight);
  }
  printf("\n");
#endif
  
}

void query(int s, int t, char fseed)
{
  #if DEBUG == 2
  printf("[%d] query s offset[%d]: %ld, query t offset[%d]: %ld\n", __LINE__, 
   s, offset[s], t, offset[s]);
  #endif
  
  min_sum=INF; //0xFFFF denote they are not reached
  old_min_sum = min_sum;
  seedNodeSet_ptr = 0;
  s_u_dist = u_t_dist = 0;
  inter_node = -1;
  s_deg = t_deg = 0;
  
  int ptr_s=0, ptr_t=0, num_edge=0;
  int node_s, node_t;
  char id_s[10];
  char id_t[10];
  
  if(offset[s] >= 0)
  { 
    char * buff_s = (char *) buff_offsetInfo;
    long offset_s=offset[s]; 

    fseeko64(labelFile,offset_s,SEEK_SET);
	  
    int num_read_s=fread(buff_s,1,BLK_SZ,labelFile);
    ptr_s=0;

    ptr_s+=sizeof(int);

    memcpy(&s_deg, buff_s+ptr_s, sizeof(int));
    ptr_s+=sizeof(int);
    
    for(int i=0; i<s_deg; i++)
    {
      memcpy(id_s, buff_s+ptr_s, sizeof(char));
      ptr_s+=sizeof(char);
      
      if(ptr_s==num_read_s)
      {
        num_read_s=fread(buff_s,1,BLK_SZ,labelFile);
        ptr_s = 0;
      }
    
      memcpy(id_s+sizeof(char), buff_s+ptr_s, sizeof(char));
      ptr_s+=sizeof(char);
      if(ptr_s==num_read_s)
      {
        num_read_s=fread(buff_s,1,BLK_SZ,labelFile);
        ptr_s = 0;
      }
    
      memcpy(id_s+2*sizeof(char), buff_s+ptr_s, sizeof(char));
      ptr_s+=sizeof(char);
      if(ptr_s==num_read_s)
      {
        num_read_s=fread(buff_s,1,BLK_SZ,labelFile);
        ptr_s = 0;
      }
      
      memcpy(id_s+3*sizeof(char), buff_s+ptr_s, sizeof(char));
      ptr_s+=sizeof(char);
      if(ptr_s==num_read_s)
      {
        num_read_s=fread(buff_s,1,BLK_SZ,labelFile);
        ptr_s = 0;
      }
      
      memcpy(&label_s[i].nodeid, id_s, sizeof(int));
      
      memcpy(&label_s[i].weight, buff_s+ptr_s, sizeof(char));

      //printf("[%d](%d, %d, %d)\n", __LINE__, s, label_s[i].nodeid, label_s[i].weight);
      
      if(label_s[i].nodeid == t)
	  {
	    if(min_sum > label_s[i].weight)
	    {
	      old_min_sum = min_sum = label_s[i].weight;
          min_sum_id = label_s[i].nodeid;
          s_u_dist = min_sum;
        }
      }  
      ptr_s+=sizeof(char);
      
      if(ptr_s==num_read_s)
      {
        num_read_s=fread(buff_s,1,BLK_SZ,labelFile);
        ptr_s = 0;
      }
    }
    
  }
  
  if(offset[t] >= 0)
  { // t is not in remaining graph
    char * buff_t = (char *) read_buff_2;
    long offset_t=offset[t]; 
    
    fseeko64(labelFile,offset_t,SEEK_SET);
    int num_read_t=fread(buff_t,1,BLK_SZ,labelFile);
    ptr_t=0;
    
    ptr_t+=sizeof(int);
    
    memcpy(&t_deg, buff_t+ptr_t, sizeof(int));
    ptr_t+=sizeof(int);
    
    for(int i=0; i<t_deg; i++)
    {
      memcpy(id_t, buff_t+ptr_t, sizeof(char));
      ptr_t+=sizeof(char);
      
      if(ptr_t==num_read_t)
      {
        num_read_t=fread(buff_t,1,BLK_SZ,labelFile);
        ptr_t = 0;
      }
      
      memcpy(id_t+sizeof(char), buff_t+ptr_t, sizeof(char));
      ptr_t+=sizeof(char);
      if(ptr_t==num_read_t)
      {
        num_read_t=fread(buff_t,1,BLK_SZ,labelFile);
        ptr_t = 0;
      }
      
      memcpy(id_t+2*sizeof(char), buff_t+ptr_t, sizeof(char));
      ptr_t+=sizeof(char);
      if(ptr_t==num_read_t)
      {
        num_read_t=fread(buff_t,1,BLK_SZ,labelFile);
        ptr_t = 0;
      }
      
      memcpy(id_t+3*sizeof(char), buff_t+ptr_t, sizeof(char));
      ptr_t+=sizeof(char);
      if(ptr_t==num_read_t)
      {
        num_read_t=fread(buff_t,1,BLK_SZ,labelFile);
        ptr_t = 0;
      }
      
      memcpy(&label_t[i].nodeid, id_t, sizeof(int));
      
      memcpy(&label_t[i].weight, buff_t+ptr_t, sizeof(char));

      //printf("[%d](%d, %d, %d)\n", __LINE__, t, label_t[i].nodeid, label_t[i].weight);
      
      if(label_t[i].nodeid == s)
	  {
	    if(min_sum > label_t[i].weight)
	    {
	      old_min_sum = min_sum = label_t[i].weight;
          min_sum_id = label_t[i].nodeid;
          u_t_dist = min_sum;
        }
      }  
      ptr_t+=sizeof(char);
      
      if(ptr_t==num_read_t)
      {
        num_read_t=fread(buff_t,1,BLK_SZ,labelFile);
        ptr_t = 0;
      }
    }
  }

  ptr_s=0;
  ptr_t=0;
  
  for(ptr_s=0; ptr_s<s_deg; ptr_s++)
  {
    while(ptr_t<t_deg)
    {
      if(label_t[ptr_t].nodeid > label_s[ptr_s].nodeid)
        break;
      if(label_t[ptr_t].nodeid == label_s[ptr_s].nodeid)
      {
        if(label_t[ptr_t].weight + label_s[ptr_s].weight <= min_sum)
        {
          min_sum = label_t[ptr_t].weight + label_s[ptr_s].weight;

          if(true == fseed)
          {    
            if(old_min_sum > min_sum)
            {
              seedNodeSet_ptr = 0;
              old_min_sum = min_sum;
            }
            
            seedNodeSet[seedNodeSet_ptr].nodeId = label_t[ptr_t].nodeid;
            seedNodeSet[seedNodeSet_ptr].l_s = ptr_s;
            seedNodeSet[seedNodeSet_ptr++].l_t = ptr_t;
          }
        }

      }
      ptr_t++;
    }
    
    if(ptr_t==t_deg)
    {
      break;
    }
  }
  
  if((true == fseed) && (min_sum < INF))
  {
    if((0 == s_deg) || (0 == t_deg))
    {
      seedNodeSet[seedNodeSet_ptr++].nodeId = min_sum_id;
      candiNodeSet[min_sum_id].s_u_dist = s_u_dist;
      candiNodeSet[min_sum_id].u_t_dist = u_t_dist;
    }
    else 
    {
      for(int i=0; i<seedNodeSet_ptr; i++)
      {
        candiNodeSet[seedNodeSet[i].nodeId].s_u_dist = label_s[seedNodeSet[i].l_s].weight;
        candiNodeSet[seedNodeSet[i].nodeId].u_t_dist = label_t[seedNodeSet[i].l_t].weight;
      }
    }
  }

  if(!flag_label)
  {
    int tmpId=0;
    char tmpWeight;
    
    if(offset[s] >= 0)
    {      
      for(int i=0; i<s_deg; i++)
      {
        //the vertex is in the remaining graph
        if(offset[label_s[i].nodeid] < 0)
        {
          tmpId = label_s[i].nodeid;
          tmpWeight = label_s[i].weight;
            
          graphX[s].push_back(make_pair(tmpId, tmpWeight));
          graphX[tmpId].push_back(make_pair(s, tmpWeight));
        }
      }
    }
  
    if(offset[t] >= 0)
    {
      for(int i=0; i<t_deg; i++)
      {
        //the vertex is in the remaining graph
        if(offset[label_t[i].nodeid] < 0)
        {
          tmpId = label_t[i].nodeid;
          tmpWeight = label_t[i].weight;
          
          graphX[t].push_back(make_pair(tmpId, tmpWeight));
          graphX[tmpId].push_back(make_pair(t, tmpWeight));
        }
      }
    }

    vector<int> dist(nGraph, 0xFFFF);  
    priority_queue< pair<int, int> > pq;
   
    dist[s] = 0;
  
    pq.push(make_pair(0, s));
    
    while (!pq.empty())
    {
      int cost = -pq.top().first; 
      int node = pq.top().second;
  
      pq.pop();
 
      if(node == t)
      {
        //printf("[%d] min_sum : %d, current_min_sum : %d\n", __LINE__, min_sum, cost);
        if(cost < min_sum)
          min_sum = cost;
        
        break;
      }
  
      for (int i = 0; i < graphX[node].size(); i++)
      {
        int v = graphX[node][i].first;
        int nextDist = cost + graphX[node][i].second;
  
        //printf("[%d] v:%d, weight:%d\n", __LINE__, v, weight);
        
        if (dist[v] > nextDist)
        {
          dist[v] = nextDist;            
          pq.push(make_pair(-nextDist, v));     
        }
      }
    }
  
    /* re-init graphX */
    {
      if(offset[s] >= 0) 
      {
        for(int i=0; i<s_deg; i++)
        {
          //the vertex is in the remaining graph
          if(offset[label_s[i].nodeid] < 0)
          {
            tmpId = label_s[i].nodeid;
    
            graphX[s].erase(graphX[s].end());
            graphX[tmpId].erase(graphX[tmpId].end());
          }
        }
      }
    
      if(offset[t] >= 0)
      {
        for(int i=0; i<t_deg; i++)
        {
          //the vertex is in the remaining graph
          if(offset[label_t[i].nodeid] < 0)
          {
            tmpId = label_t[i].nodeid;
            
            graphX[t].erase(graphX[t].end());
            graphX[tmpId].erase(graphX[tmpId].end());
          }
        }
      }
    }
  }

  #if DEBUG == 1
  printf("[%d] s_deg:%d, t_deg:%d, min_sum_id:%d, min_sum:%d, seedNodeSet_ptr:%d\n", __LINE__, 
    s_deg, t_deg, min_sum_id, min_sum, seedNodeSet_ptr);
  #endif
}

int vertex_s, vertex_t, dist_st;
char *gAnsVid;

char getAdjNodes(ReadBuffer &read_io, edgepair *label_node, 
	  	 int &label_node_ptr, int AnsVid)
{
  int u, u_deg, v, weight;
  
  read_io.read(&u);
  read_io.read(&u_deg);

  //printf("[%d] u:%d, u_deg:%d\n", __LINE__, u, u_deg);
  for(int i=0; i<u_deg; i++)
  {
    read_io.read(&v);

    if(v == AnsVid)
    {
      return false;
    }
	
    read_io.read(&weight);

    //printf("[%d] v:%d, weight:%d\n", __LINE__, v, weight);
    
    label_node[label_node_ptr].node1 = u;
    label_node[label_node_ptr].node2 = v;
    label_node[label_node_ptr++].weight = weight;
  }

  return true;
}


void start_1NIS(void)
{
  int tmp_vSize;
  edgepair curEdgePair;
  long tmp_offset;
  
  edgepair *label_node1 = (edgepair *) malloc(sizeof(edgepair) * (nGraph+1));
  if(label_node1 == NULL)
  {
    printf("[%d] calloc fail!\n", __LINE__);
    exit(1);
  }
  
  edgepair *label_node2 = (edgepair *) malloc(sizeof(edgepair) * (nGraph+1));
  if(label_node2 == NULL)
  {
    printf("[%d] calloc fail!\n", __LINE__);
    exit(1);
  }

  ReadBuffer read_io;

  tmp_vSize = v_seedEdge.size();

  while(tmp_vSize > 0)
  { 
    char find_additional_node = false, chk_no_pair_node = false;
    int label_node1_ptr = 0, label_node2_ptr = 0;
  
    curEdgePair = v_seedEdge.front();
    
    v_seedEdge.erase(v_seedEdge.begin());

    #if DEBUG == 2
    printf("<<[%d] Start searching for q_answer pair (%d, %d, %d)>>\n", __LINE__, 
      curEdgePair.node1, curEdgePair.node2, curEdgePair.weight);
    #endif
    
    //printf("[%d] internode processing for (%d, %d, %d)\n", __LINE__, 
    //  curEdgePair.node1, curEdgePair.node2, curEdgePair.weight);

    if((0 <= tmp_buff_offsetInfo[curEdgePair.node1].buff_offset) &&
       (0 <= tmp_buff_offsetInfo[curEdgePair.node2].buff_offset))
    {
      /* get the node1's offset */
      tmp_offset = tmp_buff_offsetInfo[curEdgePair.node1].buff_offset;
  	  read_io.open(SVA_LABEL_FP, tmp_offset);
  
      /* get the adj set of node1 */
  	  chk_no_pair_node = getAdjNodes(read_io, label_node1, label_node1_ptr, curEdgePair.node2);
  
      //printf("[%d] chk_no_pair_node:%d, label_node1_ptr:%d)\n", __LINE__, 
      //  chk_no_pair_node, label_node1_ptr);
      
  	  if(true == chk_no_pair_node)
      {
        tmp_offset = tmp_buff_offsetInfo[curEdgePair.node2].buff_offset;
        read_io.open(SVA_LABEL_FP, tmp_offset);
  
        chk_no_pair_node = getAdjNodes(read_io, label_node2, label_node2_ptr, curEdgePair.node1);
  
        //printf("[%d] chk_no_pair_node:%d, label_node2_ptr:%d)\n", __LINE__, 
        //  chk_no_pair_node, label_node2_ptr);
  
        for(int i=0; i<label_node1_ptr; i++)
        {
          if(true == gAnsVid[label_node1[i].node2]) 
            continue;
          
          for(int j=0; j<label_node2_ptr; j++)
          {
            if(label_node2[j].node2 > label_node1[i].node2)
              break;
            
            if(true == gAnsVid[label_node2[j].node2])
              continue;
            
            if((label_node2[j].node2 == label_node1[i].node2) && 
               (curEdgePair.weight == (label_node2[j].weight + label_node1[i].weight)))
            {
              edgepair tmpPair1, tmpPair2;
  
              //printf("[%d] we find inter node %d\n", __LINE__, label_node2[j].node2);
              
              tmpPair1 = label_node1[i];
              tmpPair2 = label_node2[j];

              /* classify ASP path results */
              if (tmpPair1.weight > 1)
              { 
                #if DEBUG
                printf("[%d] go to answer pair (%d, %d, %d)\n", __LINE__, 
                  tmpPair1.node1, tmpPair1.node2, tmpPair1.weight);
                #endif
                v_seedEdge.push_back(tmpPair1);
              }
              else if (tmpPair1.weight == 1)
              {
                #if DEBUG
                printf("[%d] go to q_reseltSet pair (%d, %d, %d)\n", __LINE__, 
                  tmpPair1.node1, tmpPair1.node2, tmpPair1.weight);
                #endif
                v_resultSet.push_back(tmpPair1);
              }
  
              if (tmpPair2.weight > 1)
              {
                #if DEBUG
                printf("\n[%d] go to answer pair (%d, %d, %d)\n", __LINE__, 
                  tmpPair2.node1, tmpPair2.node2, tmpPair2.weight);
                #endif
                v_seedEdge.push_back(tmpPair2);
              }
              else if (tmpPair2.weight == 1)
              {
                #if DEBUG
                printf("\n[%d] go to q_reseltSet pair (%d, %d, %d)\n", __LINE__, 
                  tmpPair2.node1, tmpPair2.node2, tmpPair2.weight);
                #endif
                v_resultSet.push_back(tmpPair2);
              }

              find_additional_node = true;
            }
          }
        }
      }
    }
    
    /* if the paths are not found, move the next level Gi. */
    if(false == find_additional_node)
    {
      //printf("[%d] go to answer pair again (%d, %d, %d)\n", __LINE__, 
      //  curEdgePair.node1, curEdgePair.node2, curEdgePair.weight);
      v_seedEdge.push_back(curEdgePair);
    }

    tmp_vSize--;
  }

  if(v_seedEdge.size() > 0)
    removeDuplicate(v_seedEdge);

  free(label_node1);
  free(label_node2);

  #if DEBUG == 1
  printf("\n<<[%d] End calcPath.>>\n", __LINE__);
  #endif
}

void generateSeedEdges(char tmpLevel, priority_queue <IVCInfo, vector<IVCInfo>, compare_IVC> &pq_IVC)
{
  int u, u_deg, v, weight; 
  long tmp_offset;
  ReadBuffer read_io;
  edgepair curEdgePair;

  /* get the graph gi */
  while(tmpLevel == pq_IVC.top().level && (!pq_IVC.empty()))
  {
    int s_seed_dist, t_seed_dist, dist_s = -1, dist_t = -1;
    int IVCNode = pq_IVC.top().nodeId;
    pq_IVC.pop();

    //printf("[%d] IVCNode : %d\n", __LINE__, IVCNode);
    
    tmp_offset = tmp_buff_offsetInfo[IVCNode].buff_offset;

    read_io.open(SVA_LABEL_FP, tmp_offset);
    read_io.read(&u);
    read_io.read(&u_deg);

    //printf("[%d] u : %d, u_deg : %d\n", __LINE__, u, u_deg);
        
    for(int i=0; i<u_deg; i++)
    {
      read_io.read(&v);

      if(u == v)
        continue;
      
      read_io.read(&weight);
      
      /* check triangle equality for seedNodeSet */
      if(true == gAnsVid[v])
      {
        //printf("[%d] v : %d, weight : %d\n", __LINE__, v, weight);
        
        if(-1 == dist_s && -1 == dist_t)
        {
          query(vertex_s, u, false);
          dist_s = min_sum;

          if(dist_st < dist_s)
          {
            continue;
          }

          query(u, vertex_t, false);
          dist_t = min_sum;
        }
        
        s_seed_dist = candiNodeSet[v].s_u_dist;
        t_seed_dist = candiNodeSet[v].u_t_dist;
          
        #if DEBUG == 2
        printf("[%d] #%d, weight:%d, dist_st:%d, dist_s:%d, dist_t:%d\n", __LINE__, 
          i, weight, dist_st, dist_s, dist_t);
        #endif
        
        /* onSP check */
        if(dist_st == (dist_s + dist_t))
        {
          char isOnSp = false;
          
          if(dist_s > s_seed_dist)
          {
            if((dist_s == s_seed_dist + weight) && (dist_t == t_seed_dist - weight))
              isOnSp = true;
          }
          else if(dist_s < s_seed_dist)
          {
            if((dist_s == s_seed_dist - weight) && (dist_t == t_seed_dist + weight))
              isOnSp = true;
          }

          if(true == isOnSp)
          {
            curEdgePair.node1 = u;
            curEdgePair.node2 = v;
            curEdgePair.weight = weight;

            if (curEdgePair.weight > 1)
            {
              #if DEBUG
              printf("\n[%d] go to answer pair (%d, %d, %d)\n", __LINE__, u, v, weight);
              #endif
              v_seedEdge.push_back(curEdgePair);
            }
            else if (curEdgePair.weight == 1)
            {
              #if DEBUG
              printf("\n[%d] weight is 1!. go to result pair: (%d, %d, %d)\n", __LINE__, u, v, weight);
              #endif
              v_resultSet.push_back(curEdgePair);
            }

            if(false == gAnsVid[u])
            {
              candiNodeSet[u].s_u_dist = dist_s;
              candiNodeSet[u].u_t_dist = dist_t;
            
              gAnsVid[u] = true;
            } 
          } 
        }
      }
    }
  }
}

void extendSeedNode(char* fileName, int u_s_deg, int u_t_deg, 
  weightedge *u_lable_s, weightedge *u_label_t, char u_max_level)
{
  char tmpLevel, *tmpCheck;
  priority_queue <IVCInfo, vector<IVCInfo>, compare_IVC> pq_IVC;
  IVCInfo temp;  

  tmp_max_level = u_max_level;
  
  for(int i=0; i<seedNodeSet_ptr; i++)
  { 
    //printf("[%d] seed node %d..\n", __LINE__, seedNodeSet[i].nodeId);
    gAnsVid[seedNodeSet[i].nodeId] = true;
    tmpLevel = levelInfo[seedNodeSet[i].nodeId];
    
    if(tmp_max_level < tmpLevel)
      tmp_max_level = tmpLevel;
  }
  
  tmpLevel = tmp_max_level;

  if((u_s_deg > 0) || (u_t_deg > 0))
  {
    tmpCheck = (char *) calloc(sizeof(char), nGraph+1);
    if(NULL == tmpCheck)
    {
      printf("Calloc fail %d\n", __LINE__);
      exit(1); 
    }

    for(int i=0; i<u_s_deg; i++)
    {
      if((false == tmpCheck[u_lable_s[i].nodeid]) && 
         (levelInfo[u_lable_s[i].nodeid] <= tmpLevel))
      {
        temp.nodeId = u_lable_s[i].nodeid;
        temp.level = levelInfo[temp.nodeId];
  
        pq_IVC.push(temp);
        tmpCheck[u_lable_s[i].nodeid] = true;
      }
    }
    
    for(int i=0; i<u_t_deg; i++)
    {
      if((false == tmpCheck[u_label_t[i].nodeid]) && 
         (levelInfo[u_label_t[i].nodeid] <= tmpLevel))
      {
        temp.nodeId = u_label_t[i].nodeid;
        temp.level = levelInfo[temp.nodeId];
    
        pq_IVC.push(temp);
        tmpCheck[u_label_t[i].nodeid] = true;
      }
    }                

    free(tmpCheck);
  
    if(pq_IVC.size() > 0)
    {
      while((tmpLevel <= pq_IVC.top().level) && (pq_IVC.size() > 0))
      {
        //printf("[%d] remove pq_IVC.top().nodeId : %d, pq_IVC.top().level : %d\n", __LINE__, 
          //pq_IVC.top().nodeId, pq_IVC.top().level);
        pq_IVC.pop();
      }
    }
  }
  
  //printf("[%d] pq_IVC size : %d\n", __LINE__, pq_IVC.size());
  
  tmpLevel--;
  
  while(tmpLevel >= 0)
  {
    /* access target graph file */
    sprintf(SVA_GRAPH_LABEL, "%s.%d.hds", fileName, tmpLevel);
    sprintf(SVA_GRAPH_OFFSET, "%s.%d.hds.offset", fileName, tmpLevel);
    
    SVA_LABEL_FP = fopen(SVA_GRAPH_LABEL,"rb"); 
    if(NULL == SVA_LABEL_FP)
    {
      printf("Cannot open source file %d\n", __LINE__);
      exit(1); 
    }
    
    SVA_OFFSET_FP = fopen(SVA_GRAPH_OFFSET,"rb");  
    if(NULL == SVA_OFFSET_FP)
    {
      printf("Cannot open source file %d\n", __LINE__);
      exit(1); 
    }
    
    #if DEBUG
    printf("\n[%d]open %s\n", __LINE__, SVA_GRAPH_OFFSET);
    #endif
    fread(tmp_buff_offsetInfo, sizeof(offsetInfo), nGraph, SVA_OFFSET_FP);

    /* =============== 1-hop neighbor intersection search ================ */
    start_1NIS();

    /* ================== seed edges generation ==================== */
    if(pq_IVC.size() > 0)
      generateSeedEdges(tmpLevel, pq_IVC);

    fclose(SVA_LABEL_FP);
    fclose(SVA_OFFSET_FP);
    
    tmpLevel--;
  }
}

void query_path(char* fileName)
{  
  if(!flag_label)
  {
    weightedge *u_label_s = NULL, *u_label_t = NULL;
    int u_s_deg = 0, u_t_deg = 0;

    vector<edgepair> v_stPathset;
      
    searchGraphX(vertex_s, vertex_t, v_stPathset);
    
    u_s_deg = s_deg; u_t_deg = t_deg;
      
    if(u_s_deg > 0)
    {
      u_label_s = (weightedge *) calloc(sizeof(weightedge), (u_s_deg+1));
      if(NULL == u_label_s)
      {
        printf("Calloc fail %d\n", __LINE__);
        exit(1); 
      }

      for(int i=0; i<u_s_deg; i++)
      {
        u_label_s[i] = label_s[i];
      }
    }

    if(u_t_deg > 0)
    {
      u_label_t = (weightedge *) calloc(sizeof(weightedge), (u_t_deg+1));
      if(NULL == u_label_t)
      {
        printf("Calloc fail %d\n", __LINE__);
        exit(1); 
      }

      for(int i=0; i<u_t_deg; i++)
      {
        u_label_t[i] = label_t[i];
      }
    }
  
    #if DEBUG == 1
    printf("[%d] first query result : inter_node:%d, min_sum:%d, seedNodeSet_ptr:%d\n", __LINE__, 
      inter_node, min_sum, seedNodeSet_ptr);
    #endif

    if(0 == min_sum || INF == min_sum)
    {
      printf("[%d] path is not found.\n", __LINE__);
      return;
    }
    
    dist_st = min_sum;

    {
      edgepair tmpPair;

      for(int i=0; i<v_stPathset.size(); i++)
      {
        tmpPair = v_stPathset[i];
  
        if((offset[tmpPair.node1] >= 0) || (offset[tmpPair.node2] >= 0))
        { 
          /* the node in Gx is added to seedNodeSet */
          if(offset[tmpPair.node1] < 0)
          {
            seedNodeSet[seedNodeSet_ptr].nodeId = tmpPair.node1;
            seedNodeSet_ptr++;
          }
          else if(offset[tmpPair.node2] < 0)
          {
            seedNodeSet[seedNodeSet_ptr].nodeId = tmpPair.node2;
            seedNodeSet_ptr++;
          }
  
          if(vertex_s == tmpPair.node1)
          {
            candiNodeSet[tmpPair.node2].s_u_dist = tmpPair.weight;
            candiNodeSet[tmpPair.node2].u_t_dist = dist_st - tmpPair.weight;
          }
          else if(vertex_s == tmpPair.node2)
          {
            candiNodeSet[tmpPair.node1].s_u_dist = tmpPair.weight;
            candiNodeSet[tmpPair.node1].u_t_dist = dist_st - tmpPair.weight;
          }
  
          if(vertex_t == tmpPair.node1)
          {
            candiNodeSet[tmpPair.node2].s_u_dist = dist_st - tmpPair.weight;
            candiNodeSet[tmpPair.node2].u_t_dist = tmpPair.weight; 
          }
          else if(vertex_t == tmpPair.node2)
          {
            candiNodeSet[tmpPair.node1].s_u_dist = dist_st - tmpPair.weight;
            candiNodeSet[tmpPair.node1].u_t_dist = tmpPair.weight; 
          }
  
        }

        if((offset[tmpPair.node1] < 0) && (offset[tmpPair.node2] < 0))
        {
          if (tmpPair.weight > 1)
          { 
            #if DEBUG
            printf("[%d] go to answer pair (%d, %d, %d)\n", __LINE__, 
              tmpPair.node1, tmpPair.node2, tmpPair.weight);
            #endif
            v_seedEdge.push_back(tmpPair);
          }
          else if (tmpPair.weight == 1)
          {
            #if DEBUG
            printf("[%d] go to q_reseltSet pair (%d, %d, %d)\n", __LINE__, 
              tmpPair.node1, tmpPair.node2, tmpPair.weight);
            #endif
            v_resultSet.push_back(tmpPair);
          }
        }
      }
    }
     
    extendSeedNode(fileName, u_s_deg, u_t_deg, u_label_s, u_label_t, max_level);

    if(NULL != u_label_s)
      free(u_label_s);

    if(NULL != u_label_t)
      free(u_label_t);
  }
  else // use only 2 hop labels index
  {  
    query(vertex_s, vertex_t, true);

    #if DEBUG
    printf("[%d]vertex_s:%d, vertex_t:%d, min_sum:%d\n", __LINE__, vertex_s, vertex_t, min_sum);
    #endif
    
    dist_st = min_sum;

    #if DEBUG
    printf("\n<<[%d] Start finding via node pair based on num %d of Nodes.>>\n\n", __LINE__, seedNodeSet_ptr);
    #endif

    extendSeedNode(fileName, s_deg, t_deg, label_s, label_t, 0);
  }
  
  if(0 == v_resultSet.size())
  {
    #if DEBUG
    printf("[%d] we does not found additional nodes\n", __LINE__);
    #endif
    return;
  }
}

void reInit()
{
  for (int i=0; i<nGraph+1; i++)
  {
    candiNodeSet[i].s_u_dist = 0;
    candiNodeSet[i].u_t_dist = 0;
    gAnsVid[i] = 0;
    tmp_buff_offsetInfo[i].buff_offset = 0;
  }
  dist_st = 0;
}

void usage_desc() 
{
  printf("\nUsage description:\n");
  printf("	filename [-m s_node t_node]\n");
  printf("	or\n");
  printf("	filename [-n the num of query]\n");
  printf("	or\n");
  printf("	filename [-l s_node t_node]\n");
  printf("	(-m	s_node t_node : ASP between s_node and t_node using label + graphX.)\n");
  printf("	(-n	num_query : random query using the number of num_query.)\n");
  printf("	(-l	s_node t_node : ASP between s_node and t_node using label only.)\n");
}

int main(int argc, char* argv[])
{
  if (argc == 1) 
  {
    usage_desc();
    return 1;
  }
  
  flag_label=false;
  int num_query=10;
  int i = 2;

  while (i < argc)
  {
    if (strcmp("-n", argv[i]) == 0) 
	{
      i++;
      num_query = atoi(argv[i++]);
    }
    else if (strcmp("-l", argv[i]) == 0) 
    {
      i++;
      flag_label=true;
      vertex_s = atoi(argv[i++]);
      vertex_t = atoi(argv[i++]);
      num_query = 1;
    }
    else if (strcmp("-m", argv[i]) == 0) 
	{
      i++;
      vertex_s = atoi(argv[i++]);
      vertex_t = atoi(argv[i++]);
      num_query = 1;
    }
  }

  //printf("[%d] vertex_s:%d, vertex_t:%d, flag_label:%d\n", __LINE__, vertex_s, vertex_t, flag_label);
  
  char GRAPH_GK_NAME[100], FINAL_LABEL_FILE_NAME[100], OFFSET_FILE_NAME[100];
  char GRAPH_GK_info[100];

  strcpy(FINAL_LABEL_FILE_NAME, argv[1]);
  strcpy(OFFSET_FILE_NAME, argv[1]);
  strcpy(GRAPH_GK_info, argv[1]);

  strcat(FINAL_LABEL_FILE_NAME,".label");
  strcat(OFFSET_FILE_NAME,".offset");
  strcat(GRAPH_GK_info,".info");

  labelFile = fopen(FINAL_LABEL_FILE_NAME,"rb");
  if(NULL == labelFile)
  {
    printf("[%d] open fail\n", __LINE__);
    exit(1);
  }
    
  offsetFile = fopen(OFFSET_FILE_NAME,"rb");
  if(NULL == offsetFile)
  {
    printf("[%d] open fail\n", __LINE__);
    exit(1);
  }

  FILE * graph_gk=fopen(GRAPH_GK_info, "r");
  if(NULL == graph_gk)
  {
    printf("[%d] open fail\n", __LINE__);
    exit(1);
  }
  
  int x=fscanf(graph_gk, "%d%d%d%d", &nGraph, &max_node_m, &max_edge_m, &max_level);
  if(x!=4)
  {
    printf("[%d] error! graph G_k info file error\n", __LINE__);
  }
  fclose(graph_gk);

  printf("[%d]<< Basic info. nGraph:%d, max_node_m:%d, max_edge_m:%d, max_level:%d>>\n", __LINE__, 
    nGraph, max_node_m, max_edge_m, max_level);
  
  //initialize several buffers
  buff_offsetInfo = (offsetInfo *)calloc(SZ_OFFSET, BLK_SZ);
  if(buff_offsetInfo == NULL)
  {
    printf("[%d] calloc fail!\n", __LINE__);
    exit(1);
  }
  
  read_buff_2 = (char *)malloc(BLK_SZ);
  
  // for query
  offset =  (long *) malloc(sizeof(long) * (nGraph+1));
  if(offset == NULL)
  {
    printf("[%d] calloc fail!\n", __LINE__);
    exit(1);
  }
  
  levelInfo = (char *) malloc(sizeof(char) * (nGraph+1));
  if(levelInfo == NULL)
  {
    printf("[%d] calloc fail!\n", __LINE__);
    exit(1);
  }
  
  label_s= (weightedge *) malloc(sizeof(weightedge) * (nGraph+1));
  if(label_s == NULL)
  {
    printf("[%d] calloc fail!\n", __LINE__);
    exit(1);
  }
  
  label_t= (weightedge *) malloc(sizeof(weightedge) * (nGraph+1));
  if(label_t == NULL)
  {
    printf("[%d] calloc fail!\n", __LINE__);
    exit(1);
  }
  
  //for grapX
  if(!flag_label)
  {
    /* set graphX size */
    graphX.resize(nGraph+1);
  
    sprintf(GRAPH_GK_NAME, "%s.%d.hds", argv[1], max_level);
    
    remainFile = fopen(GRAPH_GK_NAME,"rb");
  }

  read_Gi();
  
  int index_s=0;
  int index_t=0;
  
  total_query_time=0.0;
  largest_time=0;	
  
  /* create candidate node set */
  candiNodeSet = (candiNodeInfo*)calloc(sizeof(candiNodeInfo), (nGraph+1));
  if(candiNodeSet == NULL)
  {
    printf("[%d] calloc fail!\n", __LINE__);
    exit(1);
  }

  gAnsVid = (char*)calloc(sizeof(char), (nGraph+1));
  if(gAnsVid == NULL)
  {
    printf("[%d] calloc fail!\n", __LINE__);
    exit(1);
  }

  tmp_buff_offsetInfo = (offsetInfo *)calloc(sizeof(offsetInfo), nGraph+1);
  if(tmp_buff_offsetInfo == NULL)
  {
    printf("[%d] calloc fail!\n", __LINE__);
    exit(1);
  }

  /* create seed node set */
  seedNodeSet = (seedNodeInfo*)calloc(sizeof(seedNodeInfo), (nGraph+1));
  if(seedNodeSet == NULL)
  {
    printf("[%d] calloc fail!\n", __LINE__);
    exit(1);
  }

  for (int j=0; j<num_query; j++)
  {
    reInit();
    
    if(num_query > 1)
    {
      vertex_s=rand()%nGraph;
      vertex_t=rand()%nGraph;
    }

    //printf("[%d]<< #%d, Start doing query %d <--> %d>>\n", __LINE__, j, vertex_s, vertex_t);
    
    rt.start();

    if(vertex_s != vertex_t)
      query_path(argv[1]);
    else
      dist_st = 0;
    
    rt.stop();

    query_t=rt.GetRuntime(); 
    total_query_time+=query_t;
    cur_t_time+=query_t;
    
    if(cur_t_time > largest_time)
      largest_time=cur_t_time;

    //printf("[%4d] The query time: %lf\n", j, cur_t_time);
    
    cur_t_time = 0;

#if 1
    //print result

    if(num_query == 1)
    {
      if(dist_st == INF)
        printf("[%d]<< manual query %d <--> %d : path is not found>>\n", __LINE__, vertex_s, vertex_t);
      else
        printf("[%d]<< manual query %d <--> %d : %d>>\n", __LINE__, vertex_s, vertex_t, dist_st);
    }
    
    if(v_resultSet.size() > 0)
      removeDuplicate(v_resultSet);
    
    edgepair tmpPair;
    for(int i=0; i<v_resultSet.size(); i++)
    {
      tmpPair = v_resultSet[i];
    
      printf("%2d, edge (%d - %d)\n", i+1, tmpPair.node1, tmpPair.node2);
    }
#endif

    /* clear result */
    v_resultSet.clear();
    v_seedEdge.clear();
  }
 
  printf("\nThe average query time: %lf\n",total_query_time/num_query);
  printf("The largest query time: %lf\n", largest_time);

  free(seedNodeSet);
  free(candiNodeSet);
  free(tmp_buff_offsetInfo);
  free(gAnsVid);
  free(label_s);
  free(label_t);  
  free(offset);
  free(levelInfo);
  free(buff_offsetInfo);
  free(read_buff_2);
  
  fclose(labelFile);
  
  return 0;
}

