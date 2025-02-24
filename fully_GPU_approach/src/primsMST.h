/* Image-driven minimum-spanning-tree calcuation using Prim's algorithm.
 Average run-time should be of n*log(n) complexity.
 Uses heap data structure to speed-up finding the next lowest edge-weight.
 Edge-weights 
*/

class Edge
{
public:
	double weight;
	int vert1;
	int vert2;
	Edge(double w=0,int v1=0,int v2=0);
	bool operator<(const Edge & b) const;
	void print();
};

Edge::Edge(double w,int v1,int v2){
	weight=w;
	vert1=v1;
	vert2=v2;
}

bool Edge::operator<(const Edge & b) const{
	return (this->weight>b.weight);
}

int newEdge(Edge edge1,Edge& edgeout,bool* vertices){
	bool new1=vertices[edge1.vert1]; 
	bool new2=vertices[edge1.vert2];
	int out1;
	if(new1^new2){
		if(new1){
			out1=edge1.vert2;
			edgeout=Edge(edge1.weight,edge1.vert1,edge1.vert2);
		}
		else {
			out1=edge1.vert1;
			edgeout=Edge(edge1.weight,edge1.vert2,edge1.vert1);
		}
	}
	else{
		out1=-1;
	}
	return out1;
}

float edgecost2weight(float val,float meanim){
    return exp(-val/meanim);
}

void primsGraph(float* im1,int* ordered,int* parents,float* edgemst,int step1,int m2,int n2,int o2){

	
	int m=m2/step1;
	int n=n2/step1;
	int o=o2/step1;
	
	int num_vertices=m*n*o; int sz=num_vertices;
	int len=m*n*o;
	timeval time1,time2;
	int num_neighbours=6;
	float* edgecost=new float[num_vertices*num_neighbours]; 
	int* index_neighbours=new int[num_vertices*num_neighbours];
	for(int i=0;i<num_vertices*num_neighbours;i++){
		edgecost[i]=0.0;
		index_neighbours[i]=-1;
	}
	
	int dx[6]={-1,1,0,0,0,0};
	int dy[6]={0,0,-1,1,0,0};
	int dz[6]={0,0,0,0,-1,1};
	int xx,yy,zz,xx2,yy2,zz2;
	//calculate edge-weights based on SAD of groups of voxels (for each control-point)
	for(int k=0;k<o;k++){
		for(int j=0;j<n;j++){
			for(int i=0;i<m;i++){
				for(int nb=0;nb<num_neighbours;nb++){
					if((i+dy[nb])>=0&(i+dy[nb])<m&(j+dx[nb])>=0&(j+dx[nb])<n&(k+dz[nb])>=0&(k+dz[nb])<o){
						index_neighbours[i+j*m+k*m*n+nb*num_vertices]=i+dy[nb]+(j+dx[nb])*m+(k+dz[nb])*m*n;
						//float randv=((float)rand()/float(RAND_MAX));
						//edgecost[i+j*m+k*m*n+nb*num_vertices]=randv;
						for(int k1=0;k1<step1;k1++){
							for(int j1=0;j1<step1;j1++){
								for(int i1=0;i1<step1;i1++){
									xx=j*step1+j1;
									yy=i*step1+i1;
									zz=k*step1+k1;
									xx2=(j+dx[nb])*step1+j1;
									yy2=(i+dy[nb])*step1+i1;
									zz2=(k+dz[nb])*step1+k1;
									edgecost[i+j*m+k*m*n+nb*num_vertices]+=fabs(im1[yy+xx*m2+zz*m2*n2]-im1[yy2+xx2*m2+zz2*m2*n2]);
								}
							}
						}
					}
				}
			}
		}
	}
    float meanim=0.0;
    for(int i=0;i<m2*n2*o2;i++){
        meanim+=im1[i];
    }
    meanim/=(float)(m2*n2*o2);
    float stdim=0.0;
    for(int i=0;i<m2*n2*o2;i++){
        stdim+=pow(im1[i]-meanim,2);
    }
    stdim=sqrt(stdim/(float)(m2*n2*o2));
    
    for(int i=0;i<sz*6;i++){
        edgecost[i]/=(float)pow(step1,3);
    }
    for(int i=0;i<sz*6;i++){
        edgecost[i]=-edgecost2weight(edgecost[i],2.0f*stdim);
    }
	
	float centrex=n/2;
	float centrey=m/2;
	float centrez=o/2;
	
	int root=m/2+n/2*m+o/2*m*n;
	
	vector<Edge> priority;
	bool* vertices=new bool[num_vertices];
	int* level=new int[num_vertices];
	for(int i=0;i<num_vertices;i++){
		vertices[i]=false;
		parents[i]=-1;
	}
	//int root=0;
	level[root]=0;
	int last=root;
	vertices[root]=true;
	Edge edgeout=Edge(0.0,-1,-1);
	Edge minedge=Edge(0.0,-1,-1);
	float cost=0.0;
	gettimeofday(&time1, NULL);
	
	for(int i=0;i<num_vertices-1;i++){ //run n-1 times to have all vertices added
		//add edges of new vertex to priority queue
		for(int j=0;j<num_neighbours;j++){
			int n=index_neighbours[last+j*num_vertices];
			if(n>=0){
				priority.push_back(Edge(edgecost[last+j*num_vertices],last,n));
				push_heap(priority.begin(),priority.end());
			}
		}
		last=-1;
		//find valid edge with lowest weight (main step of Prim's algorithm)
		while(last==-1){
			minedge=priority.front();
			pop_heap(priority.begin(),priority.end());
			priority.pop_back();
			bool new1=vertices[minedge.vert1]; //is either vertex already part of MST?
			bool new2=vertices[minedge.vert2];
			last=newEdge(minedge,edgeout,vertices); //return next valid vertex or -1 if edge exists already
		}
		cost+=edgeout.weight;
		vertices[last]=true;
		level[edgeout.vert2]=level[edgeout.vert1]+1;
		parents[edgeout.vert2]=edgeout.vert1;
	}
	
	//find correct ordering in constant time
	int maxlevel=0;
	for(int i=0;i<num_vertices;i++){
		if(level[i]>maxlevel)
			maxlevel=level[i];
	}
	maxlevel++;
	int* leveloffset=new int[maxlevel];
	int* levelcount=new int[maxlevel];
	for(int i=0;i<maxlevel;i++){
		leveloffset[i]=0;
		levelcount[i]=0;
	}
	for(int i=0;i<num_vertices;i++){
		if(level[i]<maxlevel-1)
			leveloffset[level[i]+1]++; //counting number of vertices in each level
	}
	for(int i=1;i<maxlevel;i++){
		leveloffset[i]+=leveloffset[i-1]; //cumulative sum
	}
	for(int i=0;i<num_vertices;i++){
		int num=leveloffset[level[i]]+levelcount[level[i]];
		levelcount[level[i]]++;
		ordered[num]=i;
	}

    
	gettimeofday(&time2, NULL);
	double timeAll=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
    nth_element(levelcount,levelcount+maxlevel/2,levelcount+maxlevel);
	//printf("Prims algorithm with %d levels finished in %f secs.\nMaximum %d, minimum %d, mean %d, and median %d width of tree.\n",
          // maxlevel,timeAll,*max_element(levelcount,levelcount+maxlevel),*min_element(levelcount,levelcount+maxlevel),(int)(num_vertices/maxlevel),levelcount[maxlevel/2]);
	for(int i=0;i<sz;i++){
        edgemst[i]=0.0f;
    }
    for(int i=1;i<sz;i++){
		int ochild=ordered[i];
		int oparent=parents[ordered[i]];
        for(int nb=0;nb<num_neighbours;nb++){
            int z=ochild/(m*n); int x=(ochild-z*m*n)/m; int y=ochild-z*m*n-x*m;
            int index=y+dy[nb]+(x+dx[nb])*m+(z+dz[nb])*m*n;
            if(index==oparent){
                edgemst[ochild]=-edgecost[ochild+nb*sz];
            }
        }
    }
	priority.clear();
	
	delete edgecost;
	delete index_neighbours;
	delete levelcount;
	delete leveloffset;
	delete vertices;
	delete level;
	
	
	
}



/*
__global__ void minKey_gpu()
{
    // Initialize min value
    int min = INT_MAX;

 
    for (int v = 0; v < array_sz; v++){
        if (min>key_g[v]){
            min = key_g[v]; 
            min_index_g = v;
        }
    }
 
    
}
*/

 #define SHMEM_SIZE 512

__global__ void primsGraph_1(float* edgecost,float* im1,int* index_neighbours,int num_vertices,int step1,int m,int n,int o,int m2,int n2,int o2){

    int num_neighbours=6;
    int dx[6]={-1,1,0,0,0,0};
    int dy[6]={0,0,-1,1,0,0};
    int dz[6]={0,0,0,0,-1,1};
    int xx,yy,zz,xx2,yy2,zz2;

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    if(i<m && j<n && k<o){
                for(int nb=0;nb<num_neighbours;nb++){
                    if((i+dy[nb])>=0&(i+dy[nb])<m&(j+dx[nb])>=0&(j+dx[nb])<n&(k+dz[nb])>=0&(k+dz[nb])<o){
                        index_neighbours[i+j*m+k*m*n+nb*num_vertices]=i+dy[nb]+(j+dx[nb])*m+(k+dz[nb])*m*n;
                        //float randv=((float)rand()/float(RAND_MAX));
                        //edgecost[i+j*m+k*m*n+nb*num_vertices]=randv;
                        for(int k1=0;k1<step1;k1++){
                            for(int j1=0;j1<step1;j1++){
                                for(int i1=0;i1<step1;i1++){
                                    xx=j*step1+j1;
                                    yy=i*step1+i1;
                                    zz=k*step1+k1;
                                    xx2=(j+dx[nb])*step1+j1;
                                    yy2=(i+dy[nb])*step1+i1;
                                    zz2=(k+dz[nb])*step1+k1;
                                    edgecost[i+j*m+k*m*n+nb*num_vertices]+=fabsf(im1[yy+xx*m2+zz*m2*n2]-im1[yy2+xx2*m2+zz2*m2*n2]);
                                }
                            }
                        }
                    }
                }
    
    }
}


 __global__ void minReduction11(int* key_g,int* key_r,int* indx_r,int sz) {
	
	// Allocate shared memory
	__shared__ int partial_min[SHMEM_SIZE];
	__shared__ int partial_ind[SHMEM_SIZE];

	// Calculate thread ID
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if(tid<sz){
	// Load elements into shared memory
	partial_min[threadIdx.x] = key_g[tid];
	partial_ind[threadIdx.x] = tid;
	__syncthreads();

	// Iterate of log base 2 the block dimension
	for (int s = 1; s < blockDim.x; s *= 2) {
		// Reduce the threads performing work by half previous the previous
		// iteration each cycle
		if (threadIdx.x % (2 * s) == 0) {
			if(partial_min[threadIdx.x] >partial_min[threadIdx.x + s]){
				partial_min[threadIdx.x] =partial_min[threadIdx.x + s];
				partial_ind[threadIdx.x] =partial_ind[threadIdx.x + s];
			}
		}
		__syncthreads();
	}

	// Let the thread 0 for this block write it's result to main memory
	// Result is inexed by this block
	if (threadIdx.x == 0) {
		key_r[blockIdx.x] = partial_min[0];
		indx_r[blockIdx.x] = partial_ind[0];
	}
	}
}

__global__ void minReduction22(int* key_r,int* indx_r,int sz) {
	// Allocate shared memory
	__shared__ int partial_min[SHMEM_SIZE];
	__shared__ int partial_ind[SHMEM_SIZE];

	// Calculate thread ID
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if(tid<sz){
	// Load elements into shared memory
	partial_min[threadIdx.x] = key_r[tid];
	partial_ind[threadIdx.x] = indx_r[tid];
	__syncthreads();

	// Iterate of log base 2 the block dimension
	for (int s = 1; s < blockDim.x; s *= 2) {
		// Reduce the threads performing work by half previous the previous
		// iteration each cycle
		if (threadIdx.x % (2 * s) == 0) {
			if(partial_min[threadIdx.x] >partial_min[threadIdx.x + s]){
				partial_min[threadIdx.x] =partial_min[threadIdx.x + s];
				partial_ind[threadIdx.x] =partial_ind[threadIdx.x + s];
			}
		}
		__syncthreads();
	}

	// Let the thread 0 for this block write it's result to main memory
	// Result is inexed by this block
	if (threadIdx.x == 0) {
		key_r[blockIdx.x] = partial_min[0];
		indx_r[blockIdx.x] = partial_ind[0];
	}
	}
}
/*
 __global__ void min1Reduction() {
     // Allocate shared memory
     __shared__ int partial_key[SHMEM_SIZE];
     __shared__ int partial_index[SHMEM_SIZE];
 
     // Calculate thread ID
     int tid = blockIdx.x * blockDim.x + threadIdx.x;
 
     // Load elements into shared memory
     partial_key[threadIdx.x] = key_g[tid];
     partial_index[threadIdx.x] = indx[tid];
     __syncthreads();
 
     // Iterate of log base 2 the block dimension
     for (int s = 1; s < blockDim.x; s *= 2) {
         // Reduce the threads performing work by half previous the previous
         // iteration each cycle
         if (threadIdx.x % (2 * s) == 0){

          if(partial_key[threadIdx.x] > partial_key[threadIdx.x + s]) {
             partial_key[threadIdx.x] = partial_key[threadIdx.x + s];
             partial_index[threadIdx.x] = partial_index[threadIdx.x + s];
         }
        }
         __syncthreads();
     }
 
     // Let the thread 0 for this block write it's result to main memory
     // Result is inexed by this block
     if (threadIdx.x == 0 ) {
         key_r[blockIdx.x] = partial_key[0];
         indx_r[blockIdx.x] = partial_index[0];
     }
 }

 __global__ void minReduction2() {
    // Allocate shared memory
    __shared__ int partial_key[SHMEM_SIZE];
    __shared__ int partial_index[SHMEM_SIZE];

    // Calculate thread ID
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    // Load elements into shared memory
    partial_key[threadIdx.x] = key_r[tid];
    partial_index[threadIdx.x] = indx_r[tid];
    __syncthreads();

    // Iterate of log base 2 the block dimension
    for (int s = 1; s < blockDim.x; s *= 2) {
        // Reduce the threads performing work by half previous the previous
        // iteration each cycle
        if (threadIdx.x % (2 * s) == 0){

         if(partial_key[threadIdx.x] > partial_key[threadIdx.x + s]) {
            partial_key[threadIdx.x] = partial_key[threadIdx.x + s];
            partial_index[threadIdx.x] = partial_index[threadIdx.x + s];
        }
    }
        __syncthreads();
    }

    // Let the thread 0 for this block write it's result to main memory
    // Result is inexed by this block
    if (threadIdx.x == 0 ) {
        key_r[blockIdx.x] = partial_key[0];
        indx_r[blockIdx.x] = partial_index[0];
    }
}
 */
__global__ void primsGraph_2_1_gpu(float* d_edgecost,int* d_index_neighbours,bool* d_vertices,int* d_level,int* d_parents ,int num_vertices,int root,int* d_key_g,int* d_key_r,bool* d_mstSet_g,bool* d_mstSet_r,int* d_indx_g,int* d_indx_r){


    int u = d_indx_r[0];
    d_mstSet_g[u] = true;
    d_key_g[u] = INT_MAX;
        
    int num_neighbours = 6;

    int v = threadIdx.x + blockIdx.x*blockDim.x;

    //for (int v = 0; v < num_neighbours; v++){
 
            // graph[u][v] is non zero only for adjacent
            // vertices of m mstSet[v] is false for vertices
            // not yet included in MST Update the key only
            // if graph[u][v] is smaller than key[v]
            int inx_neig = d_index_neighbours[u+v*num_vertices];
            if (d_edgecost[u+v*num_vertices] && d_mstSet_g[inx_neig] == false && d_edgecost[u+v*num_vertices] < d_key_g[inx_neig] && inx_neig>=0){
                    d_parents[inx_neig] = u, d_key_g[inx_neig] = d_edgecost[u+v*num_vertices];
                    d_level[inx_neig]=d_level[u]+1;
                    //cost+=edgecost[u+v*V];
                    d_vertices[inx_neig] = true;
            }
        //}
} 

__global__ void initialize_aray(int* d_key_g,bool* d_mstSet_g){

    int i = threadIdx.x + blockIdx.x*blockDim.x;


    //for (int i = 0; i < array_sz; i++)
        d_key_g[i] = INT_MAX;
        d_mstSet_g[i] = false;
    

}

__global__ void initialize_key(int* d_key_g){
    d_key_g[0] = 0;

}
void primsGraph_2_gpu(float* d_edgecost,int* d_index_neighbours,bool* d_vertices,int* d_level,int* d_parents ,int num_vertices,int root,int* d_key_g,int* d_key_r,bool* d_mstSet_g,bool* d_mstSet_r,int* d_indx_g,int* d_indx_r)
{

    
    cudaMemset(&d_level[root],0,sizeof(int));
    
    cudaMemset(&d_vertices[root],true,sizeof(bool));
    
    int num_neighbours = 6;

    initialize_aray<<<num_vertices,1>>>(d_key_g,d_mstSet_g);

    // Always include first 1st vertex in MST.
    // Make key 0 so that this vertex is picked as first
    // vertex.

    // TB Size
	const int TB_SIZE = 256;

	// Grid Size (No padding)
	int GRID_SIZE = floor_div(num_vertices , TB_SIZE);


    initialize_key<<<1,1>>>(d_key_g);


    cudaMemset(&d_parents[0],-1,sizeof(int)); // First node is always root of MST

    for (int count = 0; count < num_vertices - 1; count++) {
        // Pick the minimum key vertex from the
        // set of vertices not yet included in MST
        
        //minKey_gpu<<<1,1>>>();
 
        /*                
        int min = INT_MAX;
     
        for (int v = 0; v < array_sz ; v++){
            if (key_g[v] < min){
                min = key_g[v]; 
                u = v;
                }
            }
        */

        
        //find_minimum_kernel<<<256,256>>>(num_vertices);
        minReduction11<<<GRID_SIZE,TB_SIZE>>>(d_key_g,d_key_r,d_indx_r,num_vertices);
        cudaDeviceSynchronize();
        
        minReduction22<<<1,GRID_SIZE>>>(d_key_r,d_indx_r,GRID_SIZE);
        cudaDeviceSynchronize();
        //u = thrust::min_element(thrust::device, key_g, key_g + num_vertices) - key_g;
        //u=min_index;
        // Add the picked vertex to the MST Set
        //printf("%i ",u);
        //u = indx_r[0];
        //mstSet_g[u] = true;
        //key_g[u] = INT_MAX;
        // Update key value and parent index of
        // the adjacent vertices of the picked vertex.
        // Consider only those vertices which are not
        // yet included in MST
        /*
        for (int v = 0; v < num_neighbours; v++){
 
            // graph[u][v] is non zero only for adjacent
            // vertices of m mstSet[v] is false for vertices
            // not yet included in MST Update the key only
            // if graph[u][v] is smaller than key[v]
            int inx_neig = d_index_neighbours[u+v*num_vertices];
            if (d_edgecost[u+v*num_vertices] && mstSet_g[inx_neig] == false && d_edgecost[u+v*num_vertices] < key_g[inx_neig] && inx_neig>=0){
                    d_parents[inx_neig] = u, key_g[inx_neig] = d_edgecost[u+v*num_vertices];
                    d_level[inx_neig]=d_level[u]+1;
                    //cost+=edgecost[u+v*V];
                    d_vertices[inx_neig] = true;
            }
        }
        */
        primsGraph_2_1_gpu<<<6,1>>>(d_edgecost,d_index_neighbours,d_vertices,d_level,d_parents,num_vertices,root,d_key_g,d_key_r,d_mstSet_g,d_mstSet_r,d_indx_g,d_indx_r);
        cudaDeviceSynchronize();

    }
 
    // print the constructed MST
    //printMST(parent, graph);
}



int minKey(int *key, bool *mstSet,int num_vertices)
{
    // Initialize min value
    int min = INT_MAX;
    int min_index;

 
    for (int v = 0; v < num_vertices; v++){
        if (mstSet[v] == false && key[v] < min){
            min = key[v]; 
            min_index = v;
        }
    }
    return min_index;
 
    
}
                                          

/*
void primsGraph_2(float* d_edgecost,int* d_index_neighbours,bool* d_vertices,int* d_level,int* d_parents ,int num_vertices,int root)
{

    int* mutex;
    cudaMalloc((void**)&mutex,sizeof(int));
    mutex = 0;
    // Array to store constructed MST
    //int parent[V];
    int num_neighbours = 6;
    // Key values used to pick minimum weight edge in cut
    int key[num_vertices];
    
    // To represent set of vertices included in MST
    bool mstSet[num_vertices];
 
    // Initialize all keys as INFINITE
    for (int i = 0; i < num_vertices; i++)
        key[i] = INT_MAX, mstSet[i] = false;
 
    // Always include first 1st vertex in MST.
    // Make key 0 so that this vertex is picked as first
    // vertex.
    key[root] = 0;
    d_parents[root] = -1; // First node is always root of MST
    int u;
    // The MST will have V vertices
    for (int count = 0; count < num_vertices - 1; count++) {
        // Pick the minimum key vertex from the
        // set of vertices not yet included in MST
        
        //u = minKey(key, mstSet,num_vertices);
        //u=min_index;
        // Add the picked vertex to the MST Set
        //auto u = thrust::min_element(thrust::host, key, key + num_vertices) - key;
        int min = INT_MAX;
        for (int v = 0; v < num_vertices; v++){
            if (key[v] < min){
                min = key[v]; 
                u = v;
            }
        }
        
        mstSet[u] = true;
        key[u] = INT_MAX;
        // Update key value and parent index of
        // the adjacent vertices of the picked vertex.
        // Consider only those vertices which are not
        // yet included in MST
        for (int v = 0; v < num_neighbours; v++){
        
            // graph[u][v] is non zero only for adjacent
            // vertices of m mstSet[v] is false for vertices
            // not yet included in MST Update the key only
            // if graph[u][v] is smaller than key[v]
            int inx_neig = d_index_neighbours[u+v*num_vertices];
            if (d_edgecost[u+v*num_vertices] && mstSet[inx_neig] == false && d_edgecost[u+v*num_vertices] < key[inx_neig] && inx_neig>=0){
                    d_parents[inx_neig] = u, key[inx_neig] = d_edgecost[u+v*num_vertices];
                    d_level[inx_neig]=d_level[u]+1;
                    //cost+=edgecost[u+v*V];
                    d_vertices[inx_neig] = true;
            }
        }
    }
 
    // print the constructed MST
    //printMST(parent, graph);
}
*/
__global__ void primsGraph_3(float* edgecost,int step1,int sz,float stdim){

    int i = threadIdx.x + blockIdx.x*blockDim.x;

    //for(int i=0;i<sz*6;i++){
        edgecost[i]/=(float)pow(step1,3);
        edgecost[i]=-exp(-edgecost[i]/(2.0f*stdim));
    //}
    
}

__global__ void primsGraph_4(float* edgecost,float stdim,int sz){

    int i = threadIdx.x + blockIdx.x*blockDim.x;

    //for(int i=0;i<sz*6;i++){
        edgecost[i]=-exp(-edgecost[i]/(2.0f*stdim));


    //}

}

__global__ void primsGraph_5(float* im1,float* meanim,int sz){

    int i = threadIdx.x+blockIdx.x*blockDim.x;
    //for(int i=0;i<m2*n2*o2;i++){
        atomicAdd(meanim,im1[i]/sz);
    //}
}

__global__ void primsGraph_6(float* im1,float* stdim,float* meanim){

    int i = threadIdx.x+blockIdx.x*blockDim.x;
    //for(int i=0;i<m2*n2*o2;i++){
        atomicAdd(stdim,pow(im1[i]-*meanim,2));
    //}
}

__global__ void primsGraph_7(int* level,int* maxlevel,int num_vertices){

    for(int i=0;i<num_vertices;i++){
        if(level[i]>*maxlevel)
            *maxlevel=level[i];
    }
    *maxlevel++;

}

__global__ void primsGraph_8(int* level,int* leveloffset,int maxlevel,int num_vertices){

    //int i = threadIdx.x+ blockIdx.x*blockDim.x;
    for(int i=0;i<num_vertices;i++){
        if(level[i]<maxlevel-1)
            leveloffset[level[i]+1]++; //counting number of vertices in each level
    }

}
__global__ void primsGraph_9(int* leveloffset,int maxlevel){

    for(int i=1;i<maxlevel;i++){
        leveloffset[i]+=leveloffset[i-1]; //cumulative sum
    }
}

__global__ void primsGraph_10(int* leveloffset,int* levelcount,int* ordered,int* level,int num_vertices){

    for(int i=0;i<num_vertices;i++){
        int num=leveloffset[level[i]]+levelcount[level[i]];
        levelcount[level[i]]++;
        ordered[num]=i;
    }

}
__global__ void primsGraph_11(float* edgemst,float* edgecost,int* ordered,int* parents,int sz,int num_neighbours,int m,int n){

    int dx[6]={-1,1,0,0,0,0};
    int dy[6]={0,0,-1,1,0,0};
    int dz[6]={0,0,0,0,-1,1};
    
    int i = threadIdx.x + blockIdx.x*blockDim.x+1;
    //for(int i=1;i<sz;i++){
        if(i<sz){
        int ochild=ordered[i];
        int oparent=parents[ordered[i]];
        for(int nb=0;nb<num_neighbours;nb++){
            int z=ochild/(m*n); int x=(ochild-z*m*n)/m; int y=ochild-z*m*n-x*m;
            int index=y+dy[nb]+(x+dx[nb])*m+(z+dz[nb])*m*n;
            if(index==oparent){
                edgemst[ochild]=-edgecost[ochild+nb*sz];
            }
        }
    }

}
void primsGraph_GPU(float* d_im1,int* d_ordered,int* d_parents,float* d_edgemst,int step1,int m2,int n2,int o2,int sz1){

    
    int m=m2/step1;
    int n=n2/step1;
    int o=o2/step1;
    
    int num_vertices=m*n*o; int sz=num_vertices;
    int len=m*n*o;
    timeval time1,time2;
    int num_neighbours=6;
    
    float* edgecost =new float[num_vertices*num_neighbours]; 
    int* index_neighbours =new int[num_vertices*num_neighbours];
    
    float* d_edgecost;//=new float[num_vertices*num_neighbours]; 
    int* d_index_neighbours;//=new int[num_vertices*num_neighbours];
    cudaMalloc((void**)&d_edgecost,num_vertices*num_neighbours*sizeof(float));
    cudaMalloc((void**)&d_index_neighbours,num_vertices*num_neighbours*sizeof(int));

    int* d_key_g;int* d_key_r;
    bool* d_mstSet_g;bool* d_mstSet_r;
    int* d_indx_g;int* d_indx_r;
    cudaMalloc((void**)&d_key_g,num_vertices*sizeof(int));
    cudaMalloc((void**)&d_key_r,num_vertices*sizeof(int));
    cudaMalloc((void**)&d_mstSet_g,num_vertices*sizeof(bool));
    cudaMalloc((void**)&d_mstSet_r,num_vertices*sizeof(bool));
 	cudaMalloc((void**)&d_indx_g,num_vertices*sizeof(int));
 	cudaMalloc((void**)&d_indx_r,num_vertices*sizeof(int));
    /*
    for(int i=0;i<num_vertices*num_neighbours;i++){
        edgecost[i]=0.0;
        index_neighbours[i]=-1;
    }
    */

    cudaMemset(d_edgecost,0.0,num_vertices*num_neighbours*sizeof(float));
    cudaMemset(d_index_neighbours,-1,num_vertices*num_neighbours*sizeof(int));
    dim3 grid(floor_div(m,8),floor_div(n,8),floor_div(o,4));
    primsGraph_1<<<grid,dim3(8,8,4)>>>(d_edgecost,d_im1,d_index_neighbours,num_vertices,step1,m,n,o,m2,n2,o2);
    cudaDeviceSynchronize();
    
    //cudaMemcpy(edgecost,d_edgecost,num_neighbours*num_vertices*sizeof(float),cudaMemcpyDeviceToHost);
    //cudaMemcpy(index_neighbours,d_index_neighbours,num_neighbours*num_vertices*sizeof(int),cudaMemcpyDeviceToHost);
    
    int dx[6]={-1,1,0,0,0,0};
    int dy[6]={0,0,-1,1,0,0};
    int dz[6]={0,0,0,0,-1,1};
    int xx,yy,zz,xx2,yy2,zz2;
    //calculate edge-weights based on SAD of groups of voxels (for each control-point)
    /*
    for(int k=0;k<o;k++){
        for(int j=0;j<n;j++){
            for(int i=0;i<m;i++){
                for(int nb=0;nb<num_neighbours;nb++){
                    if((i+dy[nb])>=0&(i+dy[nb])<m&(j+dx[nb])>=0&(j+dx[nb])<n&(k+dz[nb])>=0&(k+dz[nb])<o){
                        index_neighbours[i+j*m+k*m*n+nb*num_vertices]=i+dy[nb]+(j+dx[nb])*m+(k+dz[nb])*m*n;
                        //float randv=((float)rand()/float(RAND_MAX));
                        //edgecost[i+j*m+k*m*n+nb*num_vertices]=randv;
                        for(int k1=0;k1<step1;k1++){
                            for(int j1=0;j1<step1;j1++){
                                for(int i1=0;i1<step1;i1++){
                                    xx=j*step1+j1;
                                    yy=i*step1+i1;
                                    zz=k*step1+k1;
                                    xx2=(j+dx[nb])*step1+j1;
                                    yy2=(i+dy[nb])*step1+i1;
                                    zz2=(k+dz[nb])*step1+k1;
                                    edgecost[i+j*m+k*m*n+nb*num_vertices]+=fabs(im1[yy+xx*m2+zz*m2*n2]-im1[yy2+xx2*m2+zz2*m2*n2]);
                                }
                            }
                        }

                    }
                }
            }
        }
    }
    */

    float* d_meanim;
    cudaMalloc((void**)&d_meanim,sizeof(float));
    cudaMemset(d_meanim,0.0,sizeof(float));

    float meanim=0.0;
    /*
    for(int i=0;i<m2*n2*o2;i++){
        meanim+=im1[i];
    }
    */
    primsGraph_5<<<m2*n2*o2,1>>>(d_im1,d_meanim,m2*n2*o2);


    //meanim/=(float)(m2*n2*o2);
    cudaMemcpy(&meanim,d_meanim,sizeof(float),cudaMemcpyDeviceToHost);

    float* d_stdim;
    cudaMalloc((void**)&d_stdim,sizeof(float));
    cudaMemset(d_stdim,0.0,sizeof(float));

    float stdim=0.0;
    /*
    for(int i=0;i<m2*n2*o2;i++){
        stdim+=pow(im1[i]-meanim,2);
    }
    */
    primsGraph_6<<<m2*n2*o2,1>>>(d_im1,d_stdim,d_meanim);
    cudaMemcpy(&stdim,d_stdim,sizeof(float),cudaMemcpyDeviceToHost);

    stdim=sqrt(stdim/(float)(m2*n2*o2));
    
    primsGraph_3<<<sz*6,1>>>(d_edgecost,step1,sz,stdim);
    cudaDeviceSynchronize();
    /*
    for(int i=0;i<sz*6;i++){
        edgecost[i]/=(float)pow(step1,3);
    }
    */
    
    //primsGraph_4<<<1,1>>>(d_edgecost,stdim,sz);
    //cudaDeviceSynchronize();

    /*
    for(int i=0;i<sz*6;i++){
        edgecost[i]=-edgecost2weight(edgecost[i],2.0f*stdim);
    }
    */
    
    float centrex=n/2;
    float centrey=m/2;
    float centrez=o/2;
    
    int root=m/2+n/2*m+o/2*m*n;
    
    vector<Edge> priority;
    bool* vertices=new bool[num_vertices];
    int* level=new int[num_vertices];
    
    bool* d_vertices;
    int* d_level;
    cudaMalloc((void**)&d_vertices,num_vertices*sizeof(bool));
    cudaMalloc((void**)&d_level,num_vertices*sizeof(int));

    cudaMemset(d_vertices,false,num_vertices*sizeof(bool));
    cudaMemset(d_parents,-1,num_vertices*sizeof(int));
    
    /*
    for(int i=0;i<num_vertices;i++){
        vertices[i]=false;
        parents[i]=-1;
    }
    */
    
    //int root=0;
    level[root]=0;
    int last=root;
    vertices[root]=true;
    //cudaMemcpy(&d_level[root],level[root],sizeof(int),cudaMemcpyHostToDevice);
    //cudaMemcpy(&d_vertices[root],vertices[root],sizeof(bool),cudaMemcpyHostToDevice);


    Edge edgeout=Edge(0.0,-1,-1);
    Edge minedge=Edge(0.0,-1,-1);
    float cost=0.0;
    gettimeofday(&time1, NULL);
    /*
    for(int i=0;i<num_vertices-1;i++){ //run n-1 times to have all vertices added
        //add edges of new vertex to priority queue
        for(int j=0;j<num_neighbours;j++){
            int n=index_neighbours[i+j*num_vertices];
                priority.push_back(Edge(edgecost[last+j*num_vertices],i,n));
        }
    }
    */
    //primsGraph_2(edgecost,index_neighbours,vertices,level,parents ,num_vertices,root);
    
    primsGraph_2_gpu(d_edgecost,d_index_neighbours,d_vertices,d_level,d_parents ,num_vertices,root,d_key_g,d_key_r,d_mstSet_g,d_mstSet_r,d_indx_g,d_indx_r);
    //cudaDeviceSynchronize();
    //cudaMemcpy(vertices,d_vertices,num_vertices*sizeof(bool),cudaMemcpyDeviceToHost);
    //cudaMemcpy(level,d_level,num_vertices*sizeof(int),cudaMemcpyDeviceToHost);
    //cudaMemcpy(parents,d_parents,num_vertices*sizeof(int),cudaMemcpyDeviceToHost);
    //cudaMemcpy(edgecost,d_edgecost,num_neighbours*num_vertices*sizeof(float),cudaMemcpyDeviceToHost);
    
    /*
    for(int i=0;i<num_vertices-1;i++){ //run n-1 times to have all vertices added
        //add edges of new vertex to priority queue
        for(int j=0;j<num_neighbours;j++){
            int n=index_neighbours[last+j*num_vertices];
            if(n>=0){
                priority.push_back(Edge(edgecost[last+j*num_vertices],last,n));
                push_heap(priority.begin(),priority.end());
            }
        }
        last=-1;
        //find valid edge with lowest weight (main step of Prim's algorithm)
        while(last==-1){
            minedge=priority.front();
            pop_heap(priority.begin(),priority.end());
            priority.pop_back();
            bool new1=vertices[minedge.vert1]; //is either vertex already part of MST?
            bool new2=vertices[minedge.vert2];
            last=newEdge(minedge,edgeout,vertices); //return next valid vertex or -1 if edge exists already
        }
        cost+=edgeout.weight;
        vertices[last]=true;
        level[edgeout.vert2]=level[edgeout.vert1]+1;
        parents[edgeout.vert2]=edgeout.vert1;
    }
    */
    
    //find correct ordering in constant time
    int maxlevel=0;
    int* d_maxlevel;
    cudaMalloc((void**)&d_maxlevel,sizeof(int));
    primsGraph_7<<<1,1>>>(d_level,d_maxlevel,num_vertices);
    cudaDeviceSynchronize();
    /*
    for(int i=0;i<num_vertices;i++){
        if(level[i]>maxlevel)
            maxlevel=level[i];
    }
    */
    cudaMemcpy(&maxlevel,d_maxlevel,sizeof(int),cudaMemcpyDeviceToHost);
    cout << "maxlevel : " << maxlevel << endl;
    maxlevel++;
    int* leveloffset=new int[maxlevel];
    int* levelcount=new int[maxlevel];

    int* d_levelcount;
    int* d_leveloffset;
    cudaMalloc((void**)&d_levelcount,maxlevel*sizeof(int));
    cudaMalloc((void**)&d_leveloffset,maxlevel*sizeof(int));

    cudaMemset(d_leveloffset,0,maxlevel*sizeof(int));
    cudaMemset(d_levelcount,0,maxlevel*sizeof(int));

    /*
    for(int i=0;i<maxlevel;i++){
        leveloffset[i]=0;
        levelcount[i]=0;
    }
    */
    primsGraph_8<<<1,1>>>(d_level,d_leveloffset,maxlevel,num_vertices);
    cudaDeviceSynchronize();
    /*
    for(int i=0;i<num_vertices;i++){
        if(level[i]<maxlevel-1)
            leveloffset[level[i]+1]++; //counting number of vertices in each level
    }
    */
    
    primsGraph_9<<<1,1>>>(d_leveloffset,maxlevel);
    cudaDeviceSynchronize();
    /*
    for(int i=1;i<maxlevel;i++){
        leveloffset[i]+=leveloffset[i-1]; //cumulative sum
    }
    */
    //cudaMemcpy(levelcount,d_levelcount,maxlevel*sizeof(int),cudaMemcpyDeviceToHost);
    //cudaMemcpy(leveloffset,d_leveloffset,maxlevel*sizeof(int),cudaMemcpyDeviceToHost);
    
    
    primsGraph_10<<<1,1>>>(d_leveloffset,d_levelcount,d_ordered,d_level,num_vertices);
    cudaDeviceSynchronize();
    /*
    for(int i=0;i<num_vertices;i++){
        int num=leveloffset[level[i]]+levelcount[level[i]];
        levelcount[level[i]]++;
        ordered[num]=i;
    }
    */

    gettimeofday(&time2, NULL);
    double timeAll=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
    //nth_element(levelcount,levelcount+maxlevel/2,levelcount+maxlevel);
    //printf("Prims algorithm with %d levels finished in %f secs.\nMaximum %d, minimum %d, mean %d, and median %d width of tree.\n",
          // maxlevel,timeAll,*max_element(levelcount,levelcount+maxlevel),*min_element(levelcount,levelcount+maxlevel),(int)(num_vertices/maxlevel),levelcount[maxlevel/2]);
    /*
    for(int i=0;i<sz;i++){
        edgemst[i]=0.0f;
    }
    */
    cudaMemset(d_edgemst,0.0f,sz*sizeof(float));
    /*
    for(int i=1;i<sz;i++){
        int ochild=ordered[i];
        int oparent=parents[ordered[i]];
        for(int nb=0;nb<num_neighbours;nb++){
            int z=ochild/(m*n); int x=(ochild-z*m*n)/m; int y=ochild-z*m*n-x*m;
            int index=y+dy[nb]+(x+dx[nb])*m+(z+dz[nb])*m*n;
            if(index==oparent){
                edgemst[ochild]=-edgecost[ochild+nb*sz];
            }
        }
    }
    */
    primsGraph_11<<<sz,1>>>(d_edgemst,d_edgecost,d_ordered,d_parents,sz,num_neighbours,m,n);
    cudaDeviceSynchronize();

    priority.clear();
    
    delete edgecost;
    delete index_neighbours;
    delete levelcount;
    delete leveloffset;
    delete vertices;
    delete level;
    cudaFree(d_mstSet_g);cudaFree(d_key_g);cudaFree(d_mstSet_r);cudaFree(d_key_r);
    cudaFree(d_indx_g);cudaFree(d_indx_g);
    
    
}
