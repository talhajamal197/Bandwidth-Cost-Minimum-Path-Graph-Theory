#include <fstream>
#include <conio.h>
#include <iostream>
#include<bits/stdc++.h> 
#include <cmath>
#include <graphics.h> 
using namespace std;

# define INF 0x3f3f3f3f 
////////////////////////prims
#define primsV 100  
  int N;
// A utility function to find the vertex with  
// minimum key value, from the set of vertices  
// not yet included in MST  
int minKey(int key[], bool mstSet[],int n)  
{  
    // Initialize min value  
    int min = INT_MAX, min_index;  
  
    for (int v = 0; v < n; v++)  
        if (mstSet[v] == false && key[v] < min)  
            min = key[v], min_index = v;  
  
    return min_index;  
}  
  
// A utility function to print the  
// constructed MST stored in parent[]  
//void printMST(int parent[], int graph[primsV][primsV],int n)  
//{  
//    cout<<"Edge \tWeight\n";  
//    for (int i = 1; i < n; i++)  
//        cout<<parent[i]<<" - "<<i<<" \t"<<graph[i][parent[i]]<<" \n";  
//}  
  
// Function to construct and print MST for  
// a graph represented using adjacency  
// matrix representation  
void primMST(int *mat,int n,int *x,int *y)  
{  
    int graph[n][n];
    for(int i=0;i<n;i++)
	 {
	 	for(int j=0;j<n;j++)
	 	{
	 		graph[i][j]=*((mat+i*n) + j);
		 }
		 
	 }
    // Array to store constructed MST  
    int parent[n];  
      
    // Key values used to pick minimum weight edge in cut  
    int key[n];  
      
    // To represent set of vertices not yet included in MST  
    bool mstSet[n];  
  
    // Initialize all keys as INFINITE  
    for (int i = 0; i < n; i++)  
        key[i] = INT_MAX, mstSet[i] = false;  
  
    // Always include first 1st vertex in MST.  
    // Make key 0 so that this vertex is picked as first vertex.  
    key[0] = 0;  
    parent[0] = -1; // First node is always root of MST  
  
    // The MST will have V vertices  
    for (int count = 0; count < n - 1; count++) 
    {  
        // Pick the minimum key vertex from the  
        // set of vertices not yet included in MST  
        int u = minKey(key, mstSet,n);  
  
        // Add the picked vertex to the MST Set  
        mstSet[u] = true;  
  
        // Update key value and parent index of  
        // the adjacent vertices of the picked vertex.  
        // Consider only those vertices which are not  
        // yet included in MST  
        for (int v = 0; v < n; v++)  
  
            // graph[u][v] is non zero only for adjacent vertices of m  
            // mstSet[v] is false for vertices not yet included in MST  
            // Update the key only if graph[u][v] is smaller than key[v]  
            if (graph[u][v] && mstSet[v] == false && graph[u][v] < key[v])  
                parent[v] = u, key[v] = graph[u][v];  
    }  
  
    // print the constructed MST  
    int gd = DETECT, gm; 
    
    
    
      char snum[5];
      char d[5];
     initwindow(1000, 1500,"VCR"); 
     setbkcolor(DARKGRAY);
    for(int i=0; i<n;i++)
    {
    	
	
       circle(x[i], y[i], 7);
       itoa(i, snum, 10);
       outtextxy(x[i], y[i], snum);
    
	}
    
    // initgraph initializes the graphics system 
    // by loading a graphics driver from disk 
    getch(); 
     cout<<"Edge \tWeight\n";  
    for (int i = 1; i < n; i++)
	{  setcolor(YELLOW);
       cout<<parent[i]<<" - "<<i<<" \t"<<graph[i][parent[i]]<<" \n";   
	   delay(80); 
       line(x[parent[i]],y[parent[i]], x[i], y[i]);
        
    	circle(x[i], y[i], 9);

	 }  
	  setcolor(WHITE);
        setbkcolor(RED);
        settextstyle(6, 0, 4);
         outtextxy(200,600, "PRIM'S ALGORITHM : ");
         setbkcolor(DARKGRAY);
	 getch();
	 closegraph();
    
} 
//////////
typedef  pair<int, int> iPair; 
  
// Structure to represent a GraphKruskal 
struct GraphKruskal 
{ 
    int V, E; 
    vector< pair<int, iPair> > edges; 
  
    // Constructor 
    GraphKruskal(int V, int E) 
    { 
        this->V = V; 
        this->E = E; 
    } 
  
    // Utility function to add an edge 
    void addEdge(int u, int v, float w) 
    { 
        edges.push_back({w, {u, v}}); 
    } 
  
    // Function to find MST using Kruskal's 
    // MST algorithm 
    float kruskalMST(int *x,int *y,int n); 
}; 
  
// To represent Disjoint Sets 
struct DisjointSets 
{ 
    int *parent, *rnk; 
    int n; 
  
    // Constructor. 
    DisjointSets(int n) 
    { 
        // Allocate memory 
        this->n = n; 
        parent = new int[n+1]; 
        rnk = new int[n+1]; 
  
        // Initially, all vertices are in 
        // different sets and have rank 0. 
        for (int i = 0; i <= n; i++) 
        { 
            rnk[i] = 0; 
  
            //every element is parent of itself 
            parent[i] = i; 
        } 
    } 
  
    // Find the parent of a node 'u' 
    // Path Compression 
    int find(int u) 
    { 
        /* Make the parent of the nodes in the path 
           from u--> parent[u] point to parent[u] */
        if (u != parent[u]) 
            parent[u] = find(parent[u]); 
        return parent[u]; 
    } 
  
    // Union by rank 
    void merge(int x, int y) 
    { 
        x = find(x), y = find(y); 
  
        /* Make tree with smaller height 
           a subtree of the other tree  */
        if (rnk[x] > rnk[y]) 
            parent[y] = x; 
        else // If rnk[x] <= rnk[y] 
            parent[x] = y; 
  
        if (rnk[x] == rnk[y]) 
            rnk[y]++; 
    } 
}; 
  
 /* Functions returns weight of the MST*/ 
  
float GraphKruskal::kruskalMST(int *x,int *y,int n) 
{ 
    float mst_wt = 0; // Initialize result 
  
    // Sort edges in increasing order on basis of cost 
    sort(edges.begin(), edges.end()); 
  
    // Create disjoint sets 
    DisjointSets ds(V); 
  
    // Iterate through all sorted edges 
    vector< pair<int, iPair> >::iterator it; 
      int gd = DETECT, gm; 
    
    
    
      char snum[5];
   
    initwindow(1000, 1500,"VCR"); 
     
    for(int i=0; i<n;i++)
    {
	   setbkcolor(DARKGRAY);
       circle(x[i], y[i], 7);
       itoa(i, snum, 10);
       
       settextstyle(8, 0, 1); 
       outtextxy(x[i], y[i], snum);
    
	}
    
    // initgraph initializes the graphics system 
    // by loading a graphics driver from disk 
    getch(); 
    
    
    
    
    int i=1;
     
    for (it=edges.begin(); it!=edges.end(); it++) 
    { 
        int u = it->second.first; 
        int v = it->second.second; 
  
        int set_u = ds.find(u); 
        int set_v = ds.find(v); 
  
        // Check if the selected edge is creating 
        // a cycle or not (Cycle is created if u 
        // and v belong to same set) 
        if (set_u != set_v) 
        {   setcolor(YELLOW); 
            i++;

             
            // Current edge will be in the MST 
            // so print it 
            cout << u << " - " << v << endl; 
            delay(80);
            line(x[u],y[u], x[v], y[v]);;
            // Update MST weight 
            mst_wt += it->first; 
  
            // Merge two sets 
            ds.merge(set_u, set_v); 
        } 
    } 
     setcolor(WHITE);
        setbkcolor(RED);
        settextstyle(6, 0, 4);
         outtextxy(200,600, " KRUSKAL'S ALGORITHM ");
    getch(); 
  
    // closegraph function closes the graphics 
    // mode and deallocates all memory allocated 
    // by graphics system . 
    closegraph();
  
    return mst_wt; 
} //////////////////////KRuskal End
  
class Graph 
{ 
    int V;    // No. of vertices 
  
    // In a weighted graph, we need to store vertex  
    // and weight pair for every edge 
    list< pair<int, int> > *adj; 
  
public: 
    Graph(int V);  // Constructor 
  
    // function to add an edge to graph 
    void addEdge(int u, int v, int w); 
  
    // prints shortest path from s 
    void shortestPath(int s,int *x,int *y,int n); 
}; 
  
// Allocates memory for adjacency list 
Graph::Graph(int V) 
{ 
    this->V = V; 
    adj = new list< pair<int, int> >[V]; 
} 
  
void Graph::addEdge(int u, int v, int w) 
{ 
    adj[u].push_back(make_pair(v, w)); 
    adj[v].push_back(make_pair(u, w)); 
} 
  
// Prints shortest paths from src to all other vertices 
void Graph::shortestPath(int src,int *x,int *y,int n) 
{ 

    set< pair<int, int> > setds; 

    vector<int> dist(V, INF); 
    setds.insert(make_pair(0, src)); 
    dist[src] = 0; 

    while (!setds.empty()) 
    { 

        pair<int, int> tmp = *(setds.begin()); 
        setds.erase(setds.begin()); 
        int u = tmp.second; 

        list< pair<int, int> >::iterator i; 
        for (i = adj[u].begin(); i != adj[u].end(); ++i) 
        { 

            int v = (*i).first; 
            int weight = (*i).second; 

            if (dist[v] > dist[u] + weight) 
            { 
   
                if (dist[v] != INF) 
                    setds.erase(setds.find(make_pair(dist[v], v))); 
  
                // Updating distance of v 
                dist[v] = dist[u] + weight; 
                setds.insert(make_pair(dist[v], v)); 
            } 
        } 
    } 
    int gd = DETECT, gm; 
    
    
    
      char snum[5];
   
     initwindow(1000, 1500,"VCR"); 
    for(int i=0; i<n;i++)
    {
	
       circle(x[i], y[i], 7);
       itoa(i, snum, 10);
       outtextxy(x[i], y[i], snum);
    
	}
    
    // initgraph initializes the graphics system 
    // by loading a graphics driver from disk 
    getch(); 
   
  
    // closegraph function closes the graphics 
    // mode and deallocates all memory allocated 
    // by graphics system . 
    
  
    // Print shortest distances stored in dist[] 
    printf("Vertex   Distance from Source\n"); 
    int dij[2][V];
   int sum =0;
    for (int i = 0; i < V; ++i) 
    {
    	int d=dist[i];
    	dij[0][i]=i;
		dij[1][i]=dist[i];
    	
    	printf("%d \t\t %d\n", i, dist[i]);
    //	line(x[src],y[src], x[i], y[i]);
    	circle(x[i], y[i], 9);
    	sum=sum+dist[i];
		 
    	
	}
	int k[V],t1,t2;
	for (int i = 0; i < V; ++i) 
    {
		 for(int j=0;j<V-1;j++)
		 {
		 	if(dij[1][j]>dij[1][j+1])
		 	{
		 		t1=dij[0][j];
		 		t2=dij[1][j];
		 		dij[0][j]=dij[0][j+1];
		 		dij[1][j]=dij[1][j+1];
		 		dij[0][j+1]=t1;
		 		dij[1][j+1]=t2;
			 }
		 }
    	
	}
	for (int i = 0; i < V; ++i) 
    {
		k[i]=dij[0][i];
    	
	}
	int prv=1;
		for (int i = 0; i < V; ++i) 
    {    setcolor(YELLOW);
	    delay(50);
    	line(x[prv],y[prv], x[k[i]], y[k[i]]);
    	prv=k[i];
	}


	setbkcolor(RED);
    settextstyle(6, 0, 4);
   outtextxy(200,650, "DIJSTRA ALGORITHM ");
	cout<<"\nShortest Path through Dijstra Is :"<<sum;
	getch();
	
	closegraph();

	
        
        
} 








//////////////////////////////////////////////belmen///////////////////////////////////////////////


struct Edge { 
    int src, dest, weight; 
}; 
  
// a structure to represent a connected, directed and 
// weighted graph 
struct bGraph { 
    // V-> Number of vertices, E-> Number of edges 
    int V, E; 
  
    // graph is represented as an array of edges. 
    struct Edge* edge; 
}; 
  
// Creates a graph with V vertices and E edges 
struct bGraph* createGraph(int V, int E) 
{ 
    struct bGraph* graph = new bGraph; 
    graph->V = V; 
    graph->E = E; 
    graph->edge = new Edge[E]; 
    return graph; 
} 
  
// A utility function used to print the solution 
void printArr(int dist[], int n,int *x,int *y,int N) 
{ 
     int gd = DETECT, gm; 
    
    
    
      char snum[5];
   
     initwindow(1000, 1500,"VCR"); 
    for(int i=0; i<N;i++)
    {
	
       circle(x[i], y[i], 7);
       itoa(i, snum, 10);
       settextstyle(8, 0, 2); 
       outtextxy(x[i], y[i], snum);
    
	}
    
    // initgraph initializes the graphics system 
    // by loading a graphics driver from disk 
    getch(); 
    printf("Vertex   Distance from Source\n"); 
    int dij[2][n];
    for (int i = 0; i < n; ++i) 
    {
    	setcolor(YELLOW);
    	dij[0][i]=i;
		dij[1][i]=dist[i];
        printf("%d \t\t %d\n", i, dist[i]); 
      //  line(x[1],y[1], x[i], y[i]);
    	circle(x[i], y[i], 9);
        
    }
    	int k[n],t1,t2;
	for (int i = 0; i < n; ++i) 
    {
		 for(int j=0;j<n-1;j++)
		 {
		 	if(dij[1][j]>dij[1][j+1])
		 	{
		 		t1=dij[0][j];
		 		t2=dij[1][j];
		 		dij[0][j]=dij[0][j+1];
		 		dij[1][j]=dij[1][j+1];
		 		dij[0][j+1]=t1;
		 		dij[1][j+1]=t2;
			 }
		 }
}
	for (int i = 0; i < n; ++i) 
    {
		k[i]=dij[0][i];
    	
	}
	int prv=1;
		for (int i = 0; i < n; ++i) 
    {
	    delay(70);
    	line(x[prv],y[prv], x[k[i]], y[k[i]]);
    	prv=k[i];
	}
     	setcolor(WHITE);
        setbkcolor(RED);
        settextstyle(8, 0, 2);
         outtextxy(100,600, "BELMENFORD  ALGORITHM ");
     getch();
	closegraph();

} 
//int dij[2][V];
//
//    for (int i = 0; i < V; ++i) 
//    {
//    	int d=dist[i];
//    	dij[0][i]=i;
//		dij[1][i]=dist[i];
//    	
//    	printf("%d \t\t %d\n", i, dist[i]);
//    //	line(x[src],y[src], x[i], y[i]);
//    	circle(x[i], y[i], 9);
//		 
//    	
//	}
//	int k[V],t1,t2;
//	for (int i = 0; i < V; ++i) 
//    {
//		 for(int j=0;j<V-1;j++)
//		 {
//		 	if(dij[1][j]>dij[1][j+1])
//		 	{
//		 		t1=dij[0][j];
//		 		t2=dij[1][j];
//		 		dij[0][j]=dij[0][j+1];
//		 		dij[1][j]=dij[1][j+1];
//		 		dij[0][j+1]=t1;
//		 		dij[1][j+1]=t2;
//			 }
//		 }
//}
//	for (int i = 0; i < V; ++i) 
//    {
//		k[i]=dij[0][i];
//    	
//	}
//	int prv=1;
//		for (int i = 0; i < V; ++i) 
//    {
//	    
//    	line(x[prv],y[prv], x[k[i]], y[k[i]]);
//    	prv=k[i];
//	}
//     getch();
//	closegraph();




  
// The main function that finds shortest distances from src to 
// all other vertices using Bellman-Ford algorithm.  The function 
// also detects negative weight cycle 
void BellmanFord(struct bGraph* graph, int src,int *x,int *y,int n) 
{ 
    int V = graph->V; 
    int E = graph->E; 
    int dist[V]; 
  
    // Step 1: Initialize distances from src to all other vertices 
    // as INFINITE 
    for (int i = 0; i < V; i++) 
        dist[i] = INT_MAX; 
    dist[src] = 0; 
  
    // Step 2: Relax all edges |V| - 1 times. A simple shortest 
    // path from src to any other vertex can have at-most |V| - 1 
    // edges 
    for (int i = 1; i <= V - 1; i++) { 
        for (int j = 0; j < E; j++) { 
            int u = graph->edge[j].src; 
            int v = graph->edge[j].dest; 
            int weight = graph->edge[j].weight; 
            if (dist[u] != INT_MAX && dist[u] + weight < dist[v]) 
                dist[v] = dist[u] + weight; 
        } 
    } 
  
    // Step 3: check for negative-weight cycles.  The above step 
    // guarantees shortest distances if graph doesn't contain 
    // negative weight cycle.  If we get a shorter path, then there 
    // is a cycle. 
    for (int i = 0; i < E; i++) { 
        int u = graph->edge[i].src; 
        int v = graph->edge[i].dest; 
        int weight = graph->edge[i].weight; 
        if (dist[u] != INT_MAX && dist[u] + weight < dist[v]) { 
            printf("Graph contains negative weight cycle"); 
            return; // If negative cycle is detected, simply return 
        } 
    } 
  
    printArr(dist, V,x,y,n); 
  
    return; 
} 



















//////////////////Floyd///









void floyds(int *mat,int n,int src,int *x,int *y)
{
	int b[n][n];
	for(int i=0;i<n;i++)
	 {
	 	for(int j=0;j<n;j++)
	 	{
	 		b[i][j]=*((mat+i*n) + j);
		 }
		 
	 }
	 
    int i, j, k;
    for (k = 0; k < n; k++)
    {
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                if ((b[i][k] * b[k][j] != 0) && (i != j))
                {
                    if ((b[i][k] + b[k][j] < b[i][j]) || (b[i][j] == 0))
                    {
                        b[i][j] = b[i][k] + b[k][j];
                    }
                }
            }
        }
    }
  //  for (i = 0; i < n; i++)
    //{'
        int gd = DETECT, gm; 
    
    
    
      char snum[5];
   
       initwindow(1000, 1500,"VCR"); 
       for(int i=0; i<n;i++)
      {
	
         circle(x[i], y[i], 7);
         itoa(i, snum, 10);
         outtextxy(x[i], y[i], snum);
    
       }
      
    
    // initgraph initializes the graphics system 
    // by loading a graphics driver from disk 
         getch(); 
        int floyd[n][2];
        cout<<"\nMinimum Cost With Respect to Node:"<<src<<endl;
        for (j = 0; j < n; j++)
        {
        	setcolor(WHITE);
        
            cout<<b[1][j]<<"\t";
            floyd[j][0]=j;
            floyd[j][1]=b[1][j];  
           // line(x[src],y[src], x[b[1][j]], y[b[1][j]]);
    	    circle(x[i], y[i], 9);
        }
        int t1,t2;
        	cout<<"\n";
     		cout<<"\n";
//      for (int i = 0; i < n; ++i) 
//       {
//       	  cout<<floyd[i][0]<<" "<<floyd[i][1]<<"\n";
//       }
        
        for (int i = 0; i < n; ++i) 
       {
		 for(int j=0;j<n-1;j++)
		 {
		 	if(floyd[j][1]>floyd[j+1][1])
		 	{
		 		t1=floyd[j][0];
		 		t2=floyd[j][1];
		 		floyd[j][0]=floyd[j+1][0];
		 		floyd[j][1]=floyd[j+1][1];
		 		floyd[j+1][0]=t1;
		 		floyd[j+1][1]=t2;
			 }
		 }
    	
     	}
     	cout<<"\n";
     	int prev=1;
      for (int i = 0; i < n; ++i) 
       {
       	  if(floyd[i][1]==0)
       	  {
       	  	i++;
			 }
			 setcolor(YELLOW);
			 	delay(70);
       	  line(x[prev],y[prev],x[floyd[i][0]],y[floyd[i][0]]);
       	  prev=floyd[i][0];
       }
        

        setcolor(WHITE);
        setbkcolor(RED);
        settextstyle(6, 0, 2);
         outtextxy(100,600, "Floydwarshal Algorthim ");
       getch();
       closegraph();
    //}
}





///////////////////////floyd end




////////////
int floydwarshalAlgorithm(int *mat,int n,int src, int *x,int *y)
{
	
	floyds((int *)mat,n,src,x,y);
}









  

int DijstraAlgorithm(int *mat,int n,int *x,int *y)
{ 
     int V = n; 
    Graph g(V); 
  
	int k[n][n];
	 for(int i=0;i<n;i++)
	 {
	 	for(int j=0;j<n;j++)
	 	{
	 		k[i][j]=*((mat+i*n) + j);
		 }
		 
	 } for(int i=0;i<n;i++)
	 {
	 	for(int j=0;j<n;j++)
	 	{
	 		if(k[i][j]!=0)
	 		 g.addEdge(i, j, k[i][j]);
	 		
		 }
		 
	 }
	g.shortestPath(1,x,y,n); 
  
}
int KruskalAlgorithm(int *mat,int n,int *x,int *y,int Edges)
{
	int V = n; 
    GraphKruskal g(V,Edges); 
  
	int k[n][n];
	 for(int i=0;i<n;i++)
	 {
	 	for(int j=0;j<n;j++)
	 	{
	 		k[i][j]=*((mat+i*n) + j);
		 }
		 
	 }
	for(int i=0;i<n;i++)
   {
   	  for(int j=0;j<n;j++)
   	  {
   	  	 if(k[i][j]!=0)
   	  	 g.addEdge(i, j, k[i][j]); ; 
	  }
	  
   }
   int mst_wt = g.kruskalMST(x,y,n);
   cout<<"Total Cost by KRUSKAL is:"<<mst_wt;
}








/////////////////////////////belmen ford
int BellmanfordAlgorithm(int *mat,int n,int *x,int *y,int Edges)
{ 
     int V = n; 
     
  struct bGraph* graph = createGraph(V,Edges); 
	int k[n][n];
	 for(int i=0;i<n;i++)
	 {
	 	for(int j=0;j<n;j++)
	 	{
	 		k[i][j]=*((mat+i*n) + j);
		 }
		 
	 } 
	 int e=0;
	 for(int i=0;i<n;i++)
     {
   	  for(int j=0;j<n;j++)
   	   {
   	  	 if(k[i][j]!=0)
   	  	 {
   	  	 	cout<<"E:"<<e<<" SRC: "<<i<<" DEST: "<<j<<"w:"<<k[i][j]<<"\n";
   	  		graph->edge[e].src = i; 
            graph->edge[e].dest = j; 
            graph->edge[e].weight = k[i][j];
			e++; 
			
		 }
	   }
	  }
   BellmanFord(graph,1,x,y,n ); 
  
}
 
















int main () 
{
	
     
	 cout<<"\n\n\n\t\t\t\tPlease wait while loading\n\n";
	 char a=177, b=219;
	 cout<<"\t\t\t\t";
	 for (int i=0;i<=30;i++)
	 cout<<a;
	 cout<<"\r";
	 cout<<"\t\t\t\t";
	 for (int i=0;i<=30;i++)
	 {
	  cout<<b;
	  for (int j=0;j<=1e7;j++); //You can also use sleep function instead of for loop
	 }
	
	 system("cls");

	 int gd = DETECT, gm; 
    
    

      char snum[5];
   
     initwindow(1500, 1000,"VCR"); 
     setcolor(WHITE);
        setbkcolor(BLUE);
        settextstyle(6, 0, 4);
         outtextxy(100,100, "Welcome To MY Data Analysis And Algorithm Project : ");
        setbkcolor(RED);
        delay(700);
        settextstyle(6, 0, 4);
         outtextxy(410,300, "Author :Talha Jamal (17K-3877)");
        setbkcolor(RED);
         delay(700);
        settextstyle(6, 0, 4);
         outtextxy(480,450, "Topic  : Graph Thoery ");
         delay(500);
           outtextxy(480,550, "Press Enter To Continue: ");
         setbkcolor(DARKGRAY);
         
       getch();
       closegraph();
     
      cout<<"\n\n          	***********THE FOLDER INCLUDES TOTAL 10 INPUT TEXT FILE *******************\n";
      	 cout<<"             ________________________________________________________________________________\n";
		   cout<<"\n\n\n i.e               input10.txt  input20.txt  input30.txt  input40.txt  input50.txt\n\n                   input60.txt  input70.txt  input80.txt  input90.txt  input100.txt\n\n";
          cout<<"           ________________________________________________________________________________\n";
   int choice;
   
   cout<<"\n\n\tENTER OF HOW MUCH NODES FILE YOU WANT TO USE AS INPUT :";
   cin>>choice;
   float fileIn;
   ifstream infile; 
   if(choice ==10)
   infile.open("input10.txt");
   else if(choice==20)
    infile.open("input20.txt");
    else if(choice ==30)
   infile.open("input30.txt");
   else if(choice==40)
    infile.open("input40.txt");
    else if(choice ==50)
   infile.open("input50.txt");
   else if(choice==60)
    infile.open("input60.txt");
    else if(choice ==70)
   infile.open("input70.txt");
   else if(choice==80)
    infile.open("input80.txt");
    else if(choice ==90)
   infile.open("input90.txt");
   else if(choice==100)
    infile.open("input100.txt");
    else {
    	cout<<"\n NO SUCH FILE EXIST !!!";
    	return 0;
	}
    cout << "Reading from the file" << endl; 
    string gar;
    infile>>gar;
    int n;
    cout<<gar<<"\n";
   infile >> n; 
   int v[n],x[n],y[n];
   for(int i=0;i<n;i++)
   {
   	  infile >> fileIn;
   	  v[i]=fileIn;
   	  infile >> fileIn;
	  x[i]=fileIn*900;
	  infile >> fileIn;
   	  y[i]=fileIn*900;
   }
//    for(int i=0;i<n;i++)
//   {
//   	 
//   	  cout<<v[i]<<" "<<x[i]<<" "<<y[i];
//   	  cout<<"\n";
//   }
  
   float connect[n][22];
   int from[n];
   float garbage;
    for(int i=0;i<n;i++)
   {
   	  
   	  
   	  for(int j=0;j<22;j++)
   	  {
   	     	
   	     	connect[i][j]=0;
   	     	
	  }
	  
   }
  for(int i=0;i<n;i++)
   {
   	  infile>>fileIn;
   	  from[i]=fileIn;
   	  
   	  for(int j=0;j<2*from[i];j++)
   	  {
   	     	infile>>fileIn;
   	     	connect[i][j]=fileIn;
   	     	infile>>fileIn;
   	     	garbage=fileIn;
   	     	if(connect[i][j]>=10&&j%2!=0)
   	     	{
   	     		while(connect[i][j]>=10)
   	     		connect[i][j]=connect[i][j]/10;
   	     		
   	     		
			}
   	     	
   	     	
	  }
	  
}
  


infile.close();
  
   
   
   
 int sumOfEdges=0;
 for(int i=0;i<n;i++)
 {
 	sumOfEdges=sumOfEdges+from[i];
 }
 cout<<sumOfEdges;
 
 
   float adjMat[n][n];
    for(int i=0;i<n;i++)
   {
   	  for(int j=0;j<n;j++)
   	  {
   	  	adjMat[i][j]=0;
	  }
   }
   for(int i=0;i<n;i++)
   {
   	  for(int j=0;j<22;j++)
   	  {
   	  	 adjMat[i][(int)(connect[i][j])]=connect[i][j+1];
   	  	 j++;
	  }
   }
   
  
   int mat[n][n];
   for(int i=0;i<n;i++)
   {
   	  for(int j=0;j<n;j++)
   	  {
   	  	mat[i][j]=round(adjMat[i][j]);
	  }

   }
 
    int an;
    cout<<"\n\n\tENTER WHICH ALGORITHM YOU WANT TO APPLY :";
     cout<<"\n\n\n i.e  1)DIJIKRA    2)PRIMS  3)KRUSKAL  4)BELMEN  5)FLOYD WARSHAL  6)CLUSTERING COEFFICIENT\n\n\t ENTER ALGORITHM YOU WANT TO APPLY :";
     cin>>an;
     
     if(an==1)
     {
     	   DijstraAlgorithm((int *)mat,n,&x[0],&y[0]);
	 }
	 else if(an==2)
	 {
	 	primMST((int *)mat,n,&x[0],&y[0]);
	 }else if(an==3)
	 {
	 	int mst_wt=KruskalAlgorithm((int *)mat,n,&x[0],&y[0],sumOfEdges);
//	 	cout << "\nWeight of MST is " << mst_wt;
	 }else if(an==4)
	 {
	 	BellmanfordAlgorithm((int *)mat,n,&x[0],&y[0],sumOfEdges);
	 }else if(an==5)
	 {
	 	floydwarshalAlgorithm((int *)mat,n,1,&x[0],&y[0]);
	 }else if(an==6)
	 {
	 	float deg;
	 	float C;
	    float sumC=0;
	 	float clusteringCoefficiect[n] ;
	 	for(int i=0 ;i<n;i++)
	 	{
	 		deg=from[i];
	 	   C=(2*(deg))/(deg*(deg-1));
	 	   clusteringCoefficiect[i]=C;
			sumC=sumC+clusteringCoefficiect[i];
		 }
		 cout<<"Local Clustering Coefficient the Graph with input:"<<n<<".txt"<<" file is :"<<sumC/n;
		 
	 	
	 }
	 system("cls");
        
        
	     gd = DETECT, gm; 
         initwindow(1500, 1000,"VCR"); 
	     setcolor(WHITE);
         setbkcolor(BLACK);
         settextstyle(6, 0, 7);
         outtextxy(300,100, "THANK YOU !! ");
         getch();
       closegraph();
   return 0;
}

