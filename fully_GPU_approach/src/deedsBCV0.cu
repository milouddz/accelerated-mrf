#include <iostream>
#include <fstream>
#include <math.h>
#include <sys/time.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <map>
#include <functional>
#include <string.h>
#include <sstream>
#include <x86intrin.h>
#include <pthread.h>
#include <thread>
#include <cstddef>   
#include "zlib.h"
#include <sys/stat.h>

using namespace std;

//some global variables
int RAND_SAMPLES; //will all be set later (if needed)
int image_m; int image_n; int image_o; int image_d=1;
float SSD0=0.0; float SSD1=0.0; float SSD2=0.0; float distfx_global; float beta=1;


//float SIGMA=8.0;
int qc=1;

//struct for multi-threading of mind-calculation
struct mind_data{
	float* im1;
    float* d1;
    uint64_t* mindq;
    int qs;
    int ind_d1;
};

struct parameters{
    float alpha; int levels; bool segment,affine,rigid;
    vector<int> grid_spacing; vector<int> search_radius;
    vector<int> quantisation;
    string fixed_file,moving_file,output_stem,moving_seg_file,affine_file,deformed_file;
};

#include "imageIOgzType.h"
#include "transformations.h"
#include "primsMST.h"
#include "regularisation.h"
#include "MINDSSCbox.h"
#include "dataCostD.h"
#include "parseArguments.h"


int main (int argc, char * const argv[]) {
    
    //PARSE INPUT ARGUMENTS
    
    if(argc<4||argv[1][1]=='h'){
        cout<<"=============================================================\n";
        cout<<"Usage (required input arguments):\n";
        cout<<"./deedsBCV -F fixed.nii.gz -M moving.nii.gz -O output\n";
        cout<<"optional parameters:\n";
        cout<<" -a <regularisation parameter alpha> (default 1.6)\n";
        cout<<" -l <number of levels> (default 5)\n";
        cout<<" -G <grid spacing for each level> (default 8x7x6x5x4)\n";
        cout<<" -L <maximum search radius - each level> (default 8x7x6x5x4)\n";
        cout<<" -Q <quantisation of search step size> (default 5x4x3x2x1)\n";
        cout<<" -S <moving_segmentation.nii> (short int)\n";
        cout<<" -A <affine_matrix.txt> \n";
        cout<<"=============================================================\n";
        return 1;
    }
    parameters args;
    //defaults
    args.grid_spacing={8,7,6,5,4};
    args.search_radius={8,7,6,5,4};
    args.quantisation={5,4,3,2,1};
    args.levels=5;
    parseCommandLine(args, argc, argv);

    size_t split_fixed=args.fixed_file.find_last_of("/\\");
    if(split_fixed==string::npos){
        split_fixed=-1;
    }
    size_t split_moving=args.moving_file.find_last_of("/\\");
    if(split_moving==string::npos){
        split_moving=-1;
    }
    
    
    if(args.fixed_file.substr(args.fixed_file.length()-2)!="gz"){
        cout<<"images must have nii.gz format\n";
        return -1;
    }
    if(args.moving_file.substr(args.moving_file.length()-2)!="gz"){
        cout<<"images must have nii.gz format\n";
        return -1;
    }

    cout<<"Starting registration of "<<args.fixed_file.substr(split_fixed+1)<<" and "<<args.moving_file.substr(split_moving+1)<<"\n";
    cout<<"=============================================================\n";

    string outputflow;
    outputflow.append(args.output_stem);
    outputflow.append("_displacements.dat");
    string outputfile;
    outputfile.append(args.output_stem);
    outputfile.append("_deformed.nii.gz");
    
                      
    //READ IMAGES and INITIALISE ARRAYS
    
    timeval time1,time2,time1a,time2a;

	RAND_SAMPLES=1; //fixed/efficient random sampling strategy
	
	float* im1; float* im1b;
    int M,N,O,P; //image dimensions
    
    //==ALWAYS ALLOCATE MEMORY FOR HEADER ===/
	char* header=new char[352];
    
	readNifti(args.fixed_file,im1b,header,M,N,O,P);
    image_m=M; image_n=N; image_o=O;

	readNifti(args.moving_file,im1,header,M,N,O,P);
	

    if(M!=image_m|N!=image_n|O!=image_o){
        cout<<"Inconsistent image sizes (must have same dimensions)\n";
        return -1;
    }
	
	int m=image_m; int n=image_n; int o=image_o; int sz=m*n*o;
    
    float* d_im1; float* d_im1b;
    cudaMalloc((void**)&d_im1,m*n*o*sizeof(float));
    cudaMalloc((void**)&d_im1b,m*n*o*sizeof(float));
    
    //assume we are working with CT scans (add 1024 HU)
    float thresholdF=-1024; float thresholdM=-1024;

    for(int i=0;i<sz;i++){
        im1b[i]-=thresholdF;
        im1[i]-=thresholdM;
    }
    cudaMemcpy(d_im1,im1,sz*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(d_im1b,im1b,sz*sizeof(float),cudaMemcpyHostToDevice);

	float *warped1=new float[m*n*o];
    float* d_warped1;
    cudaMalloc((void**)&d_warped1,m*n*o*sizeof(float));

    //READ AFFINE MATRIX from linearBCV if provided (else start from identity)
    
    float* X=new float[16];
    float* d_X;
    cudaMalloc((void**)&d_X,16*sizeof(float));

    if(args.affine){
        size_t split_affine=args.affine_file.find_last_of("/\\");
        if(split_affine==string::npos){
            split_affine=-1;
        }

        cout<<"Reading affine matrix file: "<<args.affine_file.substr(split_affine+1)<<"\n";
        ifstream matfile;
        matfile.open(args.affine_file);
        for(int i=0;i<4;i++){
            string line;
            getline(matfile,line);
            sscanf(line.c_str(),"%f  %f  %f  %f",&X[i],&X[i+4],&X[i+8],&X[i+12]);
        }
        matfile.close();


    }
    else{
        cout<<"Starting with identity transform.\n";
        fill(X,X+16,0.0f);
        X[0]=1.0f; X[1+4]=1.0f; X[2+8]=1.0f; X[3+12]=1.0f;
    }
    
    for(int i=0;i<4;i++){
        printf("%+4.3f | %+4.3f | %+4.3f | %+4.3f \n",X[i],X[i+4],X[i+8],X[i+12]);//X[i],X[i+4],X[i+8],X[i+12]);
        
    }
    cudaMemcpy(d_X,X,16*sizeof(float),cudaMemcpyHostToDevice);
    //PATCH-RADIUS FOR MIND/SSC DESCRIPTORS

    vector<int> mind_step;
    for(int i=0;i<args.quantisation.size();i++){
        mind_step.push_back(floor(0.5f*(float)args.quantisation[i]+1.0f));
    }
    printf("MIND STEPS, %d, %d, %d, %d, %d \n",mind_step[0],mind_step[1],mind_step[2],mind_step[3],mind_step[4]);

	
    
	int step1; int hw1; float quant1;
	
	//set initial flow-fields to 0; i indicates backward (inverse) transform
	//u is in x-direction (2nd dimension), v in y-direction (1st dim) and w in z-direction (3rd dim)
	float* ux=new float[sz]; float* vx=new float[sz]; float* wx=new float[sz];
    float*d_ux;float*d_vx;float*d_wx;
    cudaMalloc((void**)&d_ux,sz*sizeof(float));
    cudaMalloc((void**)&d_vx,sz*sizeof(float));
    cudaMalloc((void**)&d_wx,sz*sizeof(float));

	for(int i=0;i<sz;i++){
		ux[i]=0.0; vx[i]=0.0; wx[i]=0.0;
	}

    cudaMemcpy(d_ux,ux,sz*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(d_vx,vx,sz*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(d_wx,wx,sz*sizeof(float),cudaMemcpyHostToDevice);

	int m2,n2,o2,sz2;
	int m1,n1,o1,sz1;
	m2=m/args.grid_spacing[0]; n2=n/args.grid_spacing[0]; o2=o/args.grid_spacing[0]; sz2=m2*n2*o2;
	float* u1=new float[sz2]; float* v1=new float[sz2]; float* w1=new float[sz2];
	float* u1i=new float[sz2]; float* v1i=new float[sz2]; float* w1i=new float[sz2];

    float* d_u1;float* d_v1;float* d_w1;
    float* d_u1i;float* d_v1i;float* d_w1i; 
    cudaMalloc((void**)&d_u1,sz2*sizeof(float));
    cudaMalloc((void**)&d_v1,sz2*sizeof(float));
    cudaMalloc((void**)&d_w1,sz2*sizeof(float));
    cudaMalloc((void**)&d_u1i,sz2*sizeof(float));
    cudaMalloc((void**)&d_v1i,sz2*sizeof(float));
    cudaMalloc((void**)&d_w1i,sz2*sizeof(float));

	for(int i=0;i<sz2;i++){		
		u1[i]=0.0; v1[i]=0.0; w1[i]=0.0;
		u1i[i]=0.0; v1i[i]=0.0; w1i[i]=0.0;
	}
	cudaMemcpy(d_u1,u1,sz2*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(d_v1,v1,sz2*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(d_w1,w1,sz2*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(d_u1i,u1i,sz2*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(d_v1i,v1i,sz2*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(d_w1i,w1i,sz2*sizeof(float),cudaMemcpyHostToDevice);
    int u1_u1i_size = sz2;
    float* warped0=new float[m*n*o];
    float* d_warped0;
    cudaMalloc((void**)&d_warped0,m*n*o*sizeof(float));

    float* d_SSD0;float* d_SSD1;
    cudaMalloc((void**)&d_SSD0,sizeof(float));
    cudaMalloc((void**)&d_SSD1,sizeof(float));
    cudaMemcpy(d_SSD0,&SSD0,sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(d_SSD1,&SSD1,sizeof(float),cudaMemcpyHostToDevice);

    //warpAffine(warped0,im1,im1b,X,ux,vx,wx);
    warpAffine_GPU<<<dim3(floor_div(m,8),floor_div(n,8),floor_div(o,4)),dim3(8,8,4)>>>(d_warped0,d_im1,d_im1b,d_X,d_ux,d_vx,d_wx,m,n,o,d_SSD0,d_SSD1);
    cudaDeviceSynchronize();
    cudaMemcpy(warped0,d_warped0,m*n*o*sizeof(float),cudaMemcpyDeviceToHost);
    cudaMemcpy(&SSD0,d_SSD0,sizeof(float),cudaMemcpyDeviceToHost);
    cudaMemcpy(&SSD1,d_SSD1,sizeof(float),cudaMemcpyDeviceToHost);
    uint64_t* im1_mind=new uint64_t[m*n*o];
    uint64_t* im1b_mind=new uint64_t[m*n*o];
    uint64_t* warped_mind=new uint64_t[m*n*o];
	
    uint64_t* d_im1_mind;uint64_t* d_im1b_mind;uint64_t* d_warped_mind;
    cudaMalloc((void**)&d_im1_mind,m*n*o*sizeof(uint64_t));
    cudaMalloc((void**)&d_im1b_mind,m*n*o*sizeof(uint64_t));
    cudaMalloc((void**)&d_warped_mind,m*n*o*sizeof(uint64_t));
    

	gettimeofday(&time1a, NULL);
    float timeDataSmooth=0;
	//==========================================================================================
	//==========================================================================================

    int* d_ordered;int* d_parents; float* d_edgemst;
    float* d_costall; 
    float* d_u0;float* d_v0;float* d_w0;

    for(int level=0;level<args.levels;level++){
        quant1=args.quantisation[level];
		
        float prev=mind_step[max(level-1,0)];//max(min(label_quant[max(level-1,0)],2.0f),1.0f);
        float curr=mind_step[level];//max(min(label_quant[level],2.0f),1.0f);
        
        float timeMIND=0; float timeSmooth=0; float timeData=0; float timeTrans=0;
		
        if(level==0|prev!=curr){
            gettimeofday(&time1, NULL);
            //descriptor(im1_mind,warped0,m,n,o,mind_step[level]);//im1 affine
            //descriptor(im1b_mind,im1b,m,n,o,mind_step[level]);
            cudaMemcpy(d_warped0,warped0,sz*sizeof(float),cudaMemcpyHostToDevice);
            descriptor_GPU(d_im1_mind,d_warped0,m,n,o,mind_step[level]);//im1 affine
            descriptor_GPU(d_im1b_mind,d_im1b,m,n,o,mind_step[level]);
            //cudaMemcpy(im1_mind,d_im1_mind,m*n*o*sizeof(uint64_t),cudaMemcpyDeviceToHost);
            //cudaMemcpy(im1b_mind,d_im1b_mind,m*n*o*sizeof(uint64_t),cudaMemcpyDeviceToHost);

            gettimeofday(&time2, NULL);
            timeMIND+=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
		}
		
		step1=args.grid_spacing[level];
		hw1=args.search_radius[level];
		
		int len3=pow(hw1*2+1,3);
		m1=m/step1; n1=n/step1; o1=o/step1; sz1=m1*n1*o1;
        
        float* costall=new float[sz1*len3]; 
        float* u0=new float[sz1]; float* v0=new float[sz1]; float* w0=new float[sz1];
		int* ordered=new int[sz1]; int* parents=new int[sz1]; float* edgemst=new float[sz1];
		
        
        cudaMalloc((void**)&d_u0,sz1*sizeof(float));
        cudaMalloc((void**)&d_v0,sz1*sizeof(float));
        cudaMalloc((void**)&d_w0,sz1*sizeof(float));

        cudaMalloc((void**)&d_ordered,sz1*sizeof(int));
        cudaMalloc((void**)&d_parents,sz1*sizeof(int));
        cudaMalloc((void**)&d_edgemst,sz1*sizeof(float));
        cudaMalloc((void**)&d_costall,sz1*len3*sizeof(float));
        
        cout<<"==========================================================\n";
		cout<<"Level "<<level<<" grid="<<step1<<" with sizes: "<<m1<<"x"<<n1<<"x"<<o1<<" hw="<<hw1<<" quant="<<quant1<<"\n";
		cout<<"==========================================================\n";
		
        //FULL-REGISTRATION FORWARDS
        gettimeofday(&time1, NULL);
		
        cudaMemcpy(d_u1,u1,u1_u1i_size*sizeof(float),cudaMemcpyHostToDevice);
        cudaMemcpy(d_v1,v1,u1_u1i_size*sizeof(float),cudaMemcpyHostToDevice);
        cudaMemcpy(d_w1,w1,u1_u1i_size*sizeof(float),cudaMemcpyHostToDevice);
        upsampleDeformationsCL_GPU(d_u0,d_v0,d_w0,d_u1,d_v1,d_w1,m1,n1,o1,m2,n2,o2);
        upsampleDeformationsCL_GPU(d_ux,d_vx,d_wx,d_u0,d_v0,d_w0,m,n,o,m1,n1,o1);
        /*
        cudaMemcpy(ux,d_ux,sz*sizeof(float),cudaMemcpyDeviceToHost);
        cudaMemcpy(vx,d_vx,sz*sizeof(float),cudaMemcpyDeviceToHost);
        cudaMemcpy(wx,d_wx,sz*sizeof(float),cudaMemcpyDeviceToHost);
        
        cudaMemcpy(u0,d_u0,sz1*sizeof(float),cudaMemcpyDeviceToHost);
        cudaMemcpy(v0,d_v0,sz1*sizeof(float),cudaMemcpyDeviceToHost);
        cudaMemcpy(w0,d_w0,sz1*sizeof(float),cudaMemcpyDeviceToHost);
        */
        //upsampleDeformationsCL(u0,v0,w0,u1,v1,w1,m1,n1,o1,m2,n2,o2);
        //upsampleDeformationsCL(ux,vx,wx,u0,v0,w0,m,n,o,m1,n1,o1);
        
        //float dist=landmarkDistance(ux,vx,wx,m,n,o,distsmm,casenum);
        cudaMemset(d_SSD0,0.0,sizeof(float));
        cudaMemset(d_SSD1,0.0,sizeof(float));
		warpAffine_GPU<<<dim3(floor_div(m,8),floor_div(n,8),floor_div(o,4)),dim3(8,8,4)>>>(d_warped1,d_im1,d_im1b,d_X,d_ux,d_vx,d_wx,m,n,o,d_SSD0,d_SSD1);
        cudaDeviceSynchronize();
        cudaMemcpy(warped1,d_warped1,m*n*o*sizeof(float),cudaMemcpyDeviceToHost);
        cudaMemcpy(&SSD0,d_SSD0,sizeof(float),cudaMemcpyDeviceToHost);
        cudaMemcpy(&SSD1,d_SSD1,sizeof(float),cudaMemcpyDeviceToHost);
    
		u1=new float[sz1]; v1=new float[sz1]; w1=new float[sz1];
        
        
        cudaFree(d_u1);cudaFree(d_v1);cudaFree(d_w1);
        cudaMalloc((void**)&d_u1,sz1*sizeof(float));
        cudaMalloc((void**)&d_v1,sz1*sizeof(float));
        cudaMalloc((void**)&d_w1,sz1*sizeof(float));

        gettimeofday(&time2, NULL);
		timeTrans+=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
        cout<<"T"<<flush;
        gettimeofday(&time1, NULL);
		//descriptor(warped_mind,warped1,m,n,o,mind_step[level]);
        
        cudaMemcpy(d_warped1,warped1,sz*sizeof(float),cudaMemcpyHostToDevice);
        cudaMemcpy(d_warped_mind,warped_mind,sz*sizeof(uint64_t),cudaMemcpyHostToDevice);
        descriptor_GPU(d_warped_mind,d_warped1,m,n,o,mind_step[level]);
        
        cudaMemcpy(warped_mind,d_warped_mind,sz*sizeof(uint64_t),cudaMemcpyDeviceToHost);

        gettimeofday(&time2, NULL);
		timeMIND+=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
        cout<<"M"<<flush;
        gettimeofday(&time1, NULL);
        
        //cudaMemcpy(d_im1b_mind,im1b_mind,sz*sizeof(uint64_t),cudaMemcpyHostToDevice);
        //cudaMemcpy(d_warped_mind,warped_mind,sz*sizeof(uint64_t),cudaMemcpyHostToDevice);
        dataCostCL_GPU(d_im1b_mind,d_warped_mind,d_costall,m,n,o,len3,step1,hw1,quant1,args.alpha,RAND_SAMPLES);
        cudaMemcpy(costall,d_costall,sz1*len3*sizeof(float),cudaMemcpyDeviceToHost);
        
        //dataCostCL(im1b_mind,warped_mind,costall,m,n,o,len3,step1,hw1,quant1,args.alpha,RAND_SAMPLES);
        
        gettimeofday(&time2, NULL);

		timeData+=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
        cout<<"D"<<flush;
        gettimeofday(&time1, NULL);
        cudaMemcpy(u0,d_u0,sz1*sizeof(float),cudaMemcpyDeviceToHost);
        cudaMemcpy(v0,d_v0,sz1*sizeof(float),cudaMemcpyDeviceToHost);
        cudaMemcpy(w0,d_w0,sz1*sizeof(float),cudaMemcpyDeviceToHost);

        primsGraph_GPU(d_im1b,d_ordered,d_parents,d_edgemst,step1,m,n,o,sz);

        cudaMemcpy(ordered,d_ordered,sz1*sizeof(int),cudaMemcpyDeviceToHost);
        cudaMemcpy(parents,d_parents,sz1*sizeof(int),cudaMemcpyDeviceToHost);
        cudaMemcpy(edgemst,d_edgemst,sz1*sizeof(float),cudaMemcpyDeviceToHost);

        regularisationCL(costall,u0,v0,w0,u1,v1,w1,hw1,step1,quant1,ordered,parents,edgemst);

        cudaMemcpy(u1,d_u1,u1_u1i_size*sizeof(float),cudaMemcpyHostToDevice);
        cudaMemcpy(v1,d_v1,u1_u1i_size*sizeof(float),cudaMemcpyHostToDevice);
        cudaMemcpy(w1,d_w1,u1_u1i_size*sizeof(float),cudaMemcpyHostToDevice);
        cudaMemcpy(d_costall,costall,sz1*len3*sizeof(float),cudaMemcpyHostToDevice);
        
        gettimeofday(&time2, NULL);
		timeSmooth+=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
        cout<<"S"<<flush;
		
        //FULL-REGISTRATION BACKWARDS
        gettimeofday(&time1, NULL);
		
        
        //cudaMemcpy(d_u1i,u1i,u1_u1i_size*sizeof(float),cudaMemcpyHostToDevice);
        //cudaMemcpy(d_v1i,v1i,u1_u1i_size*sizeof(float),cudaMemcpyHostToDevice);
        //cudaMemcpy(d_w1i,w1i,u1_u1i_size*sizeof(float),cudaMemcpyHostToDevice);
        
        upsampleDeformationsCL_GPU(d_u0,d_v0,d_w0,d_u1i,d_v1i,d_w1i,m1,n1,o1,m2,n2,o2);
        upsampleDeformationsCL_GPU(d_ux,d_vx,d_wx,d_u0,d_v0,d_w0,m,n,o,m1,n1,o1);
		
        //cudaMemcpy(u0,d_u0,sz1*sizeof(float),cudaMemcpyDeviceToHost);
        //cudaMemcpy(v0,d_v0,sz1*sizeof(float),cudaMemcpyDeviceToHost);
        //cudaMemcpy(w0,d_w0,sz1*sizeof(float),cudaMemcpyDeviceToHost);
        
        //cudaMemcpy(ux,d_ux,sz*sizeof(float),cudaMemcpyDeviceToHost);
        //cudaMemcpy(vx,d_vx,sz*sizeof(float),cudaMemcpyDeviceToHost);
        //cudaMemcpy(wx,d_wx,sz*sizeof(float),cudaMemcpyDeviceToHost);
        
        //upsampleDeformationsCL(u0,v0,w0,u1i,v1i,w1i,m1,n1,o1,m2,n2,o2);
        //upsampleDeformationsCL(ux,vx,wx,u0,v0,w0,m,n,o,m1,n1,o1);
        cudaMemset(d_SSD0,0.0,sizeof(float));
        cudaMemset(d_SSD1,0.0,sizeof(float));
        warpImageCL_GPU(d_warped1,d_im1b,d_warped0,d_ux,d_vx,d_wx,m,n,o,d_SSD0,d_SSD1);
		cudaMemcpy(warped1,d_warped1,sz*sizeof(float),cudaMemcpyDeviceToHost);
        cudaMemcpy(&SSD0,d_SSD0,sizeof(float),cudaMemcpyDeviceToHost);
        cudaMemcpy(&SSD1,d_SSD1,sizeof(float),cudaMemcpyDeviceToHost);

        u1i=new float[sz1]; v1i=new float[sz1]; w1i=new float[sz1];
        u1_u1i_size = sz1;
        cudaFree(d_u1i);cudaFree(d_v1i);cudaFree(d_w1i);
        cudaMalloc((void**)&d_u1i,sz1*sizeof(float));
        cudaMalloc((void**)&d_v1i,sz1*sizeof(float));
        cudaMalloc((void**)&d_w1i,sz1*sizeof(float));

        gettimeofday(&time2, NULL);
		timeTrans+=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
        cout<<"T"<<flush;
        gettimeofday(&time1, NULL);
		
        cudaMemcpy(d_warped1,warped1,sz*sizeof(float),cudaMemcpyHostToDevice);
        descriptor_GPU(d_warped_mind,d_warped1,m,n,o,mind_step[level]);
        cudaMemcpy(warped_mind,d_warped_mind,sz*sizeof(uint64_t),cudaMemcpyDeviceToHost);

        gettimeofday(&time2, NULL);
		timeMIND+=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
        cout<<"M"<<flush;
        gettimeofday(&time1, NULL);
        
                
        //cudaMemcpy(d_im1_mind,im1_mind,sz*sizeof(uint64_t),cudaMemcpyHostToDevice);
        //cudaMemcpy(d_warped_mind,warped_mind,sz*sizeof(uint64_t),cudaMemcpyHostToDevice);
        //cudaMemcpy(d_costall,costall,sz1*len3*sizeof(float),cudaMemcpyHostToDevice);
        dataCostCL_GPU(d_im1_mind,d_warped_mind,d_costall,m,n,o,len3,step1,hw1,quant1,args.alpha,RAND_SAMPLES);
        cudaMemcpy(costall,d_costall,sz1*len3*sizeof(float),cudaMemcpyDeviceToHost);
        

        //dataCostCL(im1_mind,warped_mind,costall,m,n,o,len3,step1,hw1,quant1,args.alpha,RAND_SAMPLES);
        gettimeofday(&time2, NULL);
		timeData+=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
        cout<<"D"<<flush;
        gettimeofday(&time1, NULL);
        cudaMemcpy(u0,d_u0,sz1*sizeof(float),cudaMemcpyDeviceToHost);
        cudaMemcpy(v0,d_v0,sz1*sizeof(float),cudaMemcpyDeviceToHost);
        cudaMemcpy(w0,d_w0,sz1*sizeof(float),cudaMemcpyDeviceToHost);
        
        //primsGraph(warped0,ordered,parents,edgemst,step1,m,n,o);
        primsGraph_GPU(d_warped0,d_ordered,d_parents,d_edgemst,step1,m,n,o,sz1);
        
        cudaMemcpy(ordered,d_ordered,sz1*sizeof(int),cudaMemcpyDeviceToHost);
        cudaMemcpy(parents,d_parents,sz1*sizeof(int),cudaMemcpyDeviceToHost);
        cudaMemcpy(edgemst,d_edgemst,sz1*sizeof(float),cudaMemcpyDeviceToHost);

        regularisationCL(costall,u0,v0,w0,u1i,v1i,w1i,hw1,step1,quant1,ordered,parents,edgemst);
        
        gettimeofday(&time2, NULL);
		timeSmooth+=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
        cout<<"S"<<flush;
		
        cout<<"\nTime: MIND="<<timeMIND<<", data="<<timeData<<", MST-reg="<<timeSmooth<<", transf.="<<timeTrans<<"\n speed="<<2.0*(float)sz1*(float)len3/(timeData+timeSmooth)<<" dof/s\n";
        
        gettimeofday(&time1, NULL);
        

        cudaMemcpy(d_u1i,u1i,sz1*sizeof(float),cudaMemcpyHostToDevice);
        cudaMemcpy(d_v1i,v1i,sz1*sizeof(float),cudaMemcpyHostToDevice);
        cudaMemcpy(d_w1i,w1i,sz1*sizeof(float),cudaMemcpyHostToDevice);
        cudaMemcpy(d_u1,u1,sz1*sizeof(float),cudaMemcpyHostToDevice);
        cudaMemcpy(d_v1,v1,sz1*sizeof(float),cudaMemcpyHostToDevice);
        cudaMemcpy(d_w1,w1,sz1*sizeof(float),cudaMemcpyHostToDevice);
        
        //consistentMappingCL_(u1,v1,w1,u1i,v1i,w1i,m1,n1,o1,step1);
        consistentMappingCL_GPU(d_u1,d_v1,d_w1,d_u1i,d_v1i,d_w1i,m1,n1,o1,step1);
        
        
        cudaMemcpy(d_costall,costall,sz1*len3*sizeof(float),cudaMemcpyHostToDevice);
        
        gettimeofday(&time2, NULL);
        float timeMapping=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
		
        //cout<<"Time consistentMapping: "<<timeMapping<<"  \n";
		
		//upsample deformations from grid-resolution to high-resolution (trilinear=1st-order spline)
		cudaMemcpy(u1i,d_u1i,sz1*sizeof(float),cudaMemcpyDeviceToHost);
        cudaMemcpy(v1i,d_v1i,sz1*sizeof(float),cudaMemcpyDeviceToHost);
        cudaMemcpy(w1i,d_w1i,sz1*sizeof(float),cudaMemcpyDeviceToHost);
        cudaMemcpy(u1,d_u1,sz1*sizeof(float),cudaMemcpyDeviceToHost);
        cudaMemcpy(v1,d_v1,sz1*sizeof(float),cudaMemcpyDeviceToHost);
        cudaMemcpy(w1,d_w1,sz1*sizeof(float),cudaMemcpyDeviceToHost);
        
        float jac = jacobian(u1,v1,w1,m1,n1,o1,step1);
		
        
        cout<<"SSD before registration: "<<SSD0<<" and after "<<SSD1<<"\n";
		m2=m1; n2=n1; o2=o1;
		cout<<"\n";
       
		delete u0; delete v0; delete w0;
		delete costall;
        
		delete parents; delete ordered;
        cudaFree(d_costall);
        cudaFree(d_u0);cudaFree(d_v0);cudaFree(d_w0);		
        cudaFree(d_parents);cudaFree(d_ordered);cudaFree(d_edgemst);
	}
    delete im1_mind;
    delete im1b_mind;
	//==========================================================================================
	//==========================================================================================
	
    gettimeofday(&time2a, NULL);
	float timeALL=time2a.tv_sec+time2a.tv_usec/1e6-(time1a.tv_sec+time1a.tv_usec/1e6);
    
    upsampleDeformationsCL(ux,vx,wx,u1,v1,w1,m,n,o,m1,n1,o1);
	
    float* flow=new float[sz1*3];
	for(int i=0;i<sz1;i++){
		flow[i]=u1[i]; flow[i+sz1]=v1[i]; flow[i+sz1*2]=w1[i];
        //flow[i+sz1*3]=u1i[i]; flow[i+sz1*4]=v1i[i]; flow[i+sz1*5]=w1i[i];
		
	}
    
    //WRITE OUTPUT DISPLACEMENT FIELD AND IMAGE
    writeOutput(flow,outputflow.c_str(),sz1*3);
    warpAffine(warped1,im1,im1b,X,ux,vx,wx);
	
    for(int i=0;i<sz;i++){
        warped1[i]+=thresholdM;
    }
    
    gzWriteNifti(outputfile,warped1,header,m,n,o,1);

    cout<<"SSD before registration: "<<SSD0<<" and after "<<SSD1<<"\n";
    
    // if SEGMENTATION of moving image is provided APPLY SAME TRANSFORM
    if(args.segment){
        short* seg2;
        readNifti(args.moving_seg_file,seg2,header,M,N,O,P);
        
        short* segw=new short[sz];
        fill(segw,segw+sz,0);

        warpAffineS(segw,seg2,X,ux,vx,wx);
        
        
        string outputseg;
        outputseg.append(args.output_stem);
        outputseg.append("_deformed_seg.nii.gz");
        cout<<"outputseg "<<outputseg<<"\n";
        

        
        gzWriteSegment(outputseg,segw,header,m,n,o,1);
    }
    
	cout<<"Finished. Total time: "<<timeALL<<" sec. ("<<timeDataSmooth<<" sec. for MIND+data+reg+trans)\n";
	
	
	return 0;
}
