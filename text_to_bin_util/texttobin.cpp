#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
using namespace std;

typedef int vertex_id;

const long int BLK_SZ = 4194304*8;  //1024*1024*4 sizeof blk for I/O
const int SZ_VERTEX = sizeof(int);
const int VERTEX_PER_BLK = BLK_SZ/SZ_VERTEX;

long no_fread = 0;
long no_fwrite = 0;

int fread_c(void * ptr, size_t size, size_t count,FILE * stream)
{
  int returnV = fread (ptr,size,count,stream);
  no_fread++;
  
  return returnV;
}

int fwrite_c ( const void * ptr, size_t size, size_t count, FILE * stream )
{
  int returnV = fwrite(ptr,size,count,stream);
  no_fwrite++;
  
  return returnV;
}

static void usage_desc() 
{
  printf("\nUsage:\n");
  printf("	text_to_bin filename1 filename2\n");
  printf("Description:\n");
  printf("	filename1 : text file (input file name, eg. input.txt).\n");
  printf("	filename2 : bin file (output file name. eg. output.bin\n");
}

//argv[1]:text file
//argv[2]:binary file
int main(int argc, char* argv[])
{
  if (argc == 1 || argc == 2) 
  {
    usage_desc();
    return 1;
  }
  vertex_id * write_buff = (vertex_id *)malloc(SZ_VERTEX*(VERTEX_PER_BLK+2));

  FILE * textFile = fopen(argv[1],"r");
  FILE * bfile = fopen(argv[2], "wb");
 
  if(textFile == NULL)
  {
    fprintf(stderr,"ERROR: cannot open file %s \n",argv[1]);
    exit(-1);
  }

  vertex_id node, degree; 
  vertex_id weight, node2; 

  long node_number_total=0;
  long edge_number_total=0;
	
  int x;

  int read_number=0;
  while(!feof(textFile))
  {
    x=fscanf(textFile, "%d", &node);
    write_buff[read_number++]=node;
    if(read_number == VERTEX_PER_BLK)
    {
      fwrite_c(write_buff,SZ_VERTEX,read_number,bfile);
      read_number = 0;
    }
		
    x=fscanf(textFile, "%d", &degree);
    write_buff[read_number++]=degree;
    if(read_number == VERTEX_PER_BLK)
    {
      fwrite_c(write_buff,SZ_VERTEX,read_number,bfile);
      read_number = 0;
    }

    node_number_total++;
    edge_number_total += degree;

    for(int j=0; j<degree; j++)
    {
      x=fscanf(textFile, "%d", &node2);
      write_buff[read_number++]=node2;
      if(read_number == VERTEX_PER_BLK)
      {
        fwrite_c(write_buff,SZ_VERTEX,read_number,bfile);
        read_number = 0;
      }
			
      x=fscanf(textFile, "%d", &weight);
      write_buff[read_number++]=weight;
      if(read_number == VERTEX_PER_BLK)
      {
        fwrite_c(write_buff,SZ_VERTEX,read_number,bfile);
        read_number = 0;
      }
    }

    x=fscanf(textFile, "\n");
  }

		
  if(read_number>0)
  {
    fwrite_c(write_buff,SZ_VERTEX,read_number,bfile);
    read_number = 0;
  }
	
  printf("node number total:%ld\n", node_number_total);
  printf("edge number total:%ld\n", edge_number_total);

  fclose(bfile);
  fclose(textFile);
  free(write_buff);
}


