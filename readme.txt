# An implementation of SPASP Algorithm (CPP Version)
===========================
The code is an implementation of the SPASP algorithm in the following Information Science 2019 paper:
Efficient Single-Pair All-Shortest-Path Query Processing for Massive Dynamic Networks. 

(1) Input file format:
The input file format of SPASP algorithm is a binary file format.  
An user can use "text_to_bin_util" to convert text format into binary format. 
The format of the input file consists of the following: 
(vertex u) (degree of u) (v_1 1) (v_2 1) ... (v_n 1), where v_1 to v_n are a set of adjacent vertices of u.
You can refer to "graph11.txt" in spasp_sample folder.

(2) Indexing: 
An user can build 2-hop labels as the following:
#spasp_label filename vertex_number [-c] [-k num_level]
Description:
  filename        The graph file.
  vertex_number   The number of vertices of the graph.
  (-c      Complete all the levels to make G_k to be empty.)
  (-k      Set num_level to be the number of levels k for k-level hierarchy.)
  (-m      Set memory to be the size of memory (in MB) to be used. The default value is 4096.)
Some inforamtion about the indexing process can be found in the screen, e.g. indexing time, # of iterations. 
We modified the indexing method presented by "Is-label: an independent-set based labeling scheme for point-to-point distance querying" paper.

(3) Implementation: 
An user can implement SPASP query processing as the following:
#spasp_query filename [-m s_node t_node]
 or
 spasp_query filename [-n the num of query]
 or
 spasp_query filename [-l s_node t_node]
 (-m     s_node t_node : ASP between s_node and t_node using label + graphX.)
 (-n     num_query : random query using the number of num_query.)
 (-l     s_node t_node : ASP between s_node and t_node using label only.)
Some inforamtion about the SPASP process can be found in the screen, e.g. processing time, ASP inforamtion. 

