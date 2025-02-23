
void boxfilter(float* input,float* temp1,float* temp2,int hw,int m,int n,int o){
    
    int sz=m*n*o;
    for(int i=0;i<sz;i++){
        temp1[i]=input[i];
    }
    
    for(int k=0;k<o;k++){
        for(int j=0;j<n;j++){
            for(int i=1;i<m;i++){
                temp1[i+j*m+k*m*n]+=temp1[(i-1)+j*m+k*m*n];
            }
        }
    }
    
    for(int k=0;k<o;k++){
        for(int j=0;j<n;j++){
            for(int i=0;i<(hw+1);i++){
                temp2[i+j*m+k*m*n]=temp1[(i+hw)+j*m+k*m*n];
            }
            for(int i=(hw+1);i<(m-hw);i++){
                temp2[i+j*m+k*m*n]=temp1[(i+hw)+j*m+k*m*n]-temp1[(i-hw-1)+j*m+k*m*n];
            }
            for(int i=(m-hw);i<m;i++){
                temp2[i+j*m+k*m*n]=temp1[(m-1)+j*m+k*m*n]-temp1[(i-hw-1)+j*m+k*m*n];
            }
        }
    }
    
    for(int k=0;k<o;k++){
        for(int j=1;j<n;j++){
            for(int i=0;i<m;i++){
                temp2[i+j*m+k*m*n]+=temp2[i+(j-1)*m+k*m*n];
            }
        }
    }
    
    for(int k=0;k<o;k++){
        for(int i=0;i<m;i++){
            for(int j=0;j<(hw+1);j++){
                temp1[i+j*m+k*m*n]=temp2[i+(j+hw)*m+k*m*n];
            }
            for(int j=(hw+1);j<(n-hw);j++){
                temp1[i+j*m+k*m*n]=temp2[i+(j+hw)*m+k*m*n]-temp2[i+(j-hw-1)*m+k*m*n];
            }
            for(int j=(n-hw);j<n;j++){
                temp1[i+j*m+k*m*n]=temp2[i+(n-1)*m+k*m*n]-temp2[i+(j-hw-1)*m+k*m*n];
            }
        }
    }
    
    for(int k=1;k<o;k++){
        for(int j=0;j<n;j++){
            for(int i=0;i<m;i++){
                temp1[i+j*m+k*m*n]+=temp1[i+j*m+(k-1)*m*n];
            }
        }
    }
    
    for(int j=0;j<n;j++){
        for(int i=0;i<m;i++){
            for(int k=0;k<(hw+1);k++){
                input[i+j*m+k*m*n]=temp1[i+j*m+(k+hw)*m*n];
            }
            for(int k=(hw+1);k<(o-hw);k++){
                input[i+j*m+k*m*n]=temp1[i+j*m+(k+hw)*m*n]-temp1[i+j*m+(k-hw-1)*m*n];
            }
            for(int k=(o-hw);k<o;k++){
                input[i+j*m+k*m*n]=temp1[i+j*m+(o-1)*m*n]-temp1[i+j*m+(k-hw-1)*m*n];
            }
        }
    }
    
    
}


void imshift(float* input,float* output,int dx,int dy,int dz,int m,int n,int o){
    for(int k=0;k<o;k++){
        for(int j=0;j<n;j++){
            for(int i=0;i<m;i++){
                if(i+dy>=0&&i+dy<m&&j+dx>=0&&j+dx<n&&k+dz>=0&&k+dz<o)
                output[i+j*m+k*m*n]=input[i+dy+(j+dx)*m+(k+dz)*m*n];
                else
                output[i+j*m+k*m*n]=input[i+j*m+k*m*n];
            }
        }
    }
}

/*void *distances(void *threadarg)
{
    struct mind_data *my_data;
    my_data = (struct mind_data *) threadarg;
    float* im1=my_data->im1;
    float* d1=my_data->d1;
    int qs=my_data->qs;
    int ind_d1=my_data->ind_d1;
    int m=image_m;
    int n=image_n;
    int o=image_o;*/

void distances(float* im1,float* d1,int m,int n,int o,int qs,int l){
    int sz1=m*n*o;
    float* w1=new float[sz1];
    int len1=6;

    float* temp1=new float[sz1]; float* temp2=new float[sz1];
    int dx[6]={+qs,+qs,-qs,+0,+qs,+0};
    int dy[6]={+qs,-qs,+0,-qs,+0,+qs};
    int dz[6]={0,+0,+qs,+qs,+qs,+qs};
    
        imshift(im1,w1,dx[l],dy[l],dz[l],m,n,o);
        for(int i=0;i<sz1;i++){
            w1[i]=(w1[i]-im1[i])*(w1[i]-im1[i]);
        }
        boxfilter(w1,temp1,temp2,qs,m,n,o);
        for(int i=0;i<sz1;i++){
            d1[i+l*sz1]=w1[i];
        }
    
    delete temp1; delete temp2; delete w1;
}

//__builtin_popcountll(left[i]^right[i]); absolute hamming distances
void descriptor(uint64_t* mindq,float* im1,int m,int n,int o,int qs){
    timeval time1,time2;
    
    //MIND with self-similarity context
    
    int dx[6]={+qs,+qs,-qs,+0,+qs,+0};
    int dy[6]={+qs,-qs,+0,-qs,+0,+qs};
    int dz[6]={0,+0,+qs,+qs,+qs,+qs};
    
    int sx[12]={-qs,+0,-qs,+0,+0,+qs,+0,+0,+0,-qs,+0,+0};
    int sy[12]={+0,-qs,+0,+qs,+0,+0,+0,+qs,+0,+0,+0,-qs};
    int sz[12]={+0,+0,+0,+0,-qs,+0,-qs,+0,-qs,+0,-qs,+0};
    
    int index[12]={0,0,1,1,2,2,3,3,4,4,5,5};
    
    float sigma=0.75;//1.0;//0.75;//1.5;
    int rho=ceil(sigma*1.5)*2+1;
    
    int len1=6;
    const int len2=12;
    
    image_d=12;
    int d=12;
    int sz1=m*n*o;

    pthread_t thread1, thread2, thread3;

    
    //============== DISTANCES USING BOXFILTER ===================
    float* d1=new float[sz1*len1];
    gettimeofday(&time1, NULL);

#pragma omp parallel for
    for(int l=0;l<len1;l++){
        distances(im1,d1,m,n,o,qs,l);
    }

    gettimeofday(&time2, NULL);
    float timeMIND1=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
    gettimeofday(&time1, NULL);

    //quantisation table
    const int val=6;
    
    const unsigned long long power=32;
    
    
#pragma omp parallel for
    for(int k=0;k<o;k++){
        unsigned int tablei[6]={0,1,3,7,15,31};
        float compare[val-1];
        for(int i=0;i<val-1;i++){
            compare[i]=-log((i+1.5f)/val);
        }
        float mind1[12];
        for(int j=0;j<n;j++){
            for(int i=0;i<m;i++){
                for(int l=0;l<len2;l++){
                    if(i+sy[l]>=0&&i+sy[l]<m&&j+sx[l]>=0&&j+sx[l]<n&&k+sz[l]>=0&&k+sz[l]<o){
                        mind1[l]=d1[i+sy[l]+(j+sx[l])*m+(k+sz[l])*m*n+index[l]*sz1];
                    }
                    else{
                        mind1[l]=d1[i+j*m+k*m*n+index[l]*sz1];
                    }
                }
                float minval=*min_element(mind1,mind1+len2);
                float sumnoise=0.0f;
                for(int l=0;l<len2;l++){
                    mind1[l]-=minval;
                    sumnoise+=mind1[l];
                }
                float noise1=max(sumnoise/(float)len2,1e-6f);
                for(int l=0;l<len2;l++){
                    mind1[l]/=noise1;
                }
                unsigned long long accum=0;
                unsigned long long tabled1=1;
                
                for(int l=0;l<len2;l++){
                    //mind1[l]=exp(-mind1[l]);
                    int mind1val=0;
                    for(int c=0;c<val-1;c++){
                        mind1val+=compare[c]>mind1[l]?1:0;
                    }
                    //int mind1val=min(max((int)(mind1[l]*val-0.5f),0),val-1);
                    accum+=tablei[mind1val]*tabled1;
                    tabled1*=power;
                    
                    
                }
                mindq[i+j*m+k*m*n]=accum;

                
                
            }
        }

    }
    

    gettimeofday(&time2, NULL);
    float timeMIND2=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
    delete d1;
    

}


///////////////////////////////////// GPU /////////////////////////////////

__global__ void boxfilter_gpu_2(float* temp1,int m,int n,int o){

    int j = threadIdx.x + blockIdx.x*blockDim.x;
    int k = threadIdx.y + blockIdx.y*blockDim.y;
    if(j<n && k<o){
    for(int i=1;i<m;i++){
        temp1[i+j*m+k*m*n]+=temp1[(i-1)+j*m+k*m*n];
    }
    }
    
}

__global__ void boxfilter_gpu_3(float* temp1,float* temp2,int hw,int m,int n,int o){

    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int j = threadIdx.y + blockIdx.y*blockDim.y;
    int k = threadIdx.z + blockIdx.z*blockDim.z;
    if(i<m && j<n && k<o){
    if(i<(hw+1)){
        temp2[i+j*m+k*m*n]=temp1[(i+hw)+j*m+k*m*n];
    }
    if(i>=(hw+1) && i<(m-hw)){
        temp2[i+j*m+k*m*n]=temp1[(i+hw)+j*m+k*m*n]-temp1[(i-hw-1)+j*m+k*m*n];
    }
    if(i>=(m-hw)){
        temp2[i+j*m+k*m*n]=temp1[(m-1)+j*m+k*m*n]-temp1[(i-hw-1)+j*m+k*m*n];
    }
    }

}

__global__ void boxfilter_gpu_4(float* temp2,int m,int n,int o){

    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int k = threadIdx.y + blockIdx.y*blockDim.y;
    if(i<m && k<o){
    for(int j=1;j<n;j++){
        temp2[i+j*m+k*m*n]+=temp2[i+(j-1)*m+k*m*n];
    }
}
    
}

__global__ void boxfilter_gpu_5(float* temp1,float* temp2,int hw,int m,int n,int o){

    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int j = threadIdx.y + blockIdx.y*blockDim.y;
    int k = threadIdx.z + blockIdx.z*blockDim.z;
    if(i<m && j<n && k<o){
    if(j<(hw+1)){
        temp1[i+j*m+k*m*n]=temp2[i+(j+hw)*m+k*m*n];
    }
    if(j>=(hw+1) && j<(n-hw)){
        temp1[i+j*m+k*m*n]=temp2[i+(j+hw)*m+k*m*n]-temp2[i+(j-hw-1)*m+k*m*n];
    }
    if(j>=(n-hw)){
        temp1[i+j*m+k*m*n]=temp2[i+(n-1)*m+k*m*n]-temp2[i+(j-hw-1)*m+k*m*n];
    }
}
    
}

__global__ void boxfilter_gpu_6(float* temp1,int m,int n,int o){

    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int j = threadIdx.y + blockIdx.y*blockDim.y;

    if(i<m && j<n){
    for(int k=1;k<o;k++){
    
        temp1[i+j*m+k*m*n]+=temp1[i+j*m+(k-1)*m*n];
    
    }
}
}

__global__ void boxfilter_gpu_7(float* input,float* temp1,int hw,int m,int n,int o){

    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int j = threadIdx.y + blockIdx.y*blockDim.y;
    int k = threadIdx.z + blockIdx.z*blockDim.z;

    //for(int j=0;j<n;j++){
        //for(int i=0;i<m;i++){
        if(i<m && j<n && k<o){
            if(k<(hw+1)){
                input[i+j*m+k*m*n]=temp1[i+j*m+(k+hw)*m*n];
            }
            if(k>=(hw+1) && k<(o-hw)){
                input[i+j*m+k*m*n]=temp1[i+j*m+(k+hw)*m*n]-temp1[i+j*m+(k-hw-1)*m*n];
            }
            if(k>=(o-hw)){
                input[i+j*m+k*m*n]=temp1[i+j*m+(o-1)*m*n]-temp1[i+j*m+(k-hw-1)*m*n];
            }
        }
    //}
}
void boxfilter_gpu(float* input,float* temp1,float* temp2,int hw,int m,int n,int o){

    dim3 grid2(floor(n/8),floor(o/4));
    boxfilter_gpu_2<<<grid2,dim3(8,4)>>>(temp1,m,n,o);
    //cudaDeviceSynchronize();

    boxfilter_gpu_3<<<dim3(floor(m/8),floor(n/8),floor(o/4)),dim3(8,8,4)>>>(temp1,temp2,hw,m,n,o);
    //cudaDeviceSynchronize();

    dim3 grid4(floor(m/8),floor(o/4));
    boxfilter_gpu_4<<<grid4,dim3(8,4)>>>(temp2,m,n,o);
    //cudaDeviceSynchronize();

    boxfilter_gpu_5<<<dim3(floor(m/8),floor(n/8),floor(o/4)),dim3(8,8,4)>>>(temp1,temp2,hw,m,n,o);
    //cudaDeviceSynchronize();

    dim3 grid6(floor(m/8),floor(n/8));
    boxfilter_gpu_6<<<grid6,dim3(8,8)>>>(temp1,m,n,o);
    //cudaDeviceSynchronize();

    boxfilter_gpu_7<<<dim3(floor(m/8),floor(n/8),floor(o/4)),dim3(8,8,4)>>>(input,temp1,hw,m,n,o);
    //cudaDeviceSynchronize();


}

__global__ void imshift_gpu(float* input,float* output,int dx,int dy,int dz,int m,int n,int o){
    
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int j = threadIdx.y + blockIdx.y*blockDim.y;
    int k = threadIdx.z + blockIdx.z*blockDim.z;
    
    if(i+dy>=0&&i+dy<m&&j+dx>=0&&j+dx<n&&k+dz>=0&&k+dz<o)
        output[i+j*m+k*m*n]=input[i+dy+(j+dx)*m+(k+dz)*m*n];
    else
        output[i+j*m+k*m*n]=input[i+j*m+k*m*n];
    
}

/*
void imshift(float* input,float* output,int dx,int dy,int dz,int m,int n,int o){
    for(int k=0;k<o;k++){
        for(int j=0;j<n;j++){
            for(int i=0;i<m;i++){
                if(i+dy>=0&&i+dy<m&&j+dx>=0&&j+dx<n&&k+dz>=0&&k+dz<o)
                output[i+j*m+k*m*n]=input[i+dy+(j+dx)*m+(k+dz)*m*n];
                else
                output[i+j*m+k*m*n]=input[i+j*m+k*m*n];
            }
        }
    }
}
*/


__global__ void distance1_gpu(float* w1,float* im1,int sz1){

    int i = threadIdx.x + blockIdx.x*blockDim.x;
    if(i<sz1){
    w1[i]=(w1[i]-im1[i])*(w1[i]-im1[i]);
    }
}

__global__ void distance2_gpu(float* d1,float* w1,int sz1,int l){

    int i = threadIdx.x + blockIdx.x*blockDim.x;
    if(i<sz1){
    d1[i+l*sz1]=w1[i];
    }
}

void distances_gpu(float* d_im1,float* d1,float* d_d1,int m,int n,int o,int qs,int l){
    int sz1=m*n*o;
    float* w1=new float[sz1];
    float* d_w1;
    float* d_temp1; 
    float* d_temp2;
    cudaMalloc((void**)&d_w1,sz1*sizeof(float));
    cudaMalloc((void**)&d_temp1,sz1*sizeof(float));
    cudaMalloc((void**)&d_temp2,sz1*sizeof(float));

    int len1=6;

    float* temp1=new float[sz1]; float* temp2=new float[sz1];
    int dx[6]={+qs,+qs,-qs,+0,+qs,+0};
    int dy[6]={+qs,-qs,+0,-qs,+0,+qs};
    int dz[6]={0,+0,+qs,+qs,+qs,+qs};
    
        dim3 grid(floor_div(m,8),floor_div(n,8),floor_div(o,4));
        imshift_gpu<<<grid,dim3(8,8,4)>>>(d_im1,d_w1,dx[l],dy[l],dz[l],m,n,o);
        cudaDeviceSynchronize();
        
        distance1_gpu<<<m*n*o,1>>>(d_w1,d_im1,sz1);
        cudaDeviceSynchronize();
        
        cudaMemcpy(d_temp1,d_w1,sz1*sizeof(float),cudaMemcpyDeviceToDevice);

        boxfilter_gpu(d_w1,d_temp1,d_temp2,qs,m,n,o);
        
        distance2_gpu<<<sz1,1>>>(d_d1,d_w1,sz1,l);
        cudaDeviceSynchronize();
    
    delete temp1; delete temp2; delete w1;
    cudaFree(d_w1);
    cudaFree(d_temp1);
    cudaFree(d_temp2);
}


__global__ void descriptor_1(uint64_t* mindq,float* compare, float* d1, int qs, int len2, int m, int n,int o){

    const int val=6;
    
    const unsigned long long power=32;
    
    
    int dy[6]={+qs,-qs,+0,-qs,+0,+qs};
    int dz[6]={0,+0,+qs,+qs,+qs,+qs};
    
    int sx[12]={-qs,+0,-qs,+0,+0,+qs,+0,+0,+0,-qs,+0,+0};
    int sy[12]={+0,-qs,+0,+qs,+0,+0,+0,+qs,+0,+0,+0,-qs};
    int sz[12]={+0,+0,+0,+0,-qs,+0,-qs,+0,-qs,+0,-qs,+0};
    
    int index[12]={0,0,1,1,2,2,3,3,4,4,5,5};
    unsigned int tablei[6]={0,1,3,7,15,31};    
    

    int sz1 = m*n*o;
    
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    
    if(i<m && j<n && k<o){

    float mind1[12];
    for(int l=0;l<len2;l++){
        if((i+sy[l])>=0&&(i+sy[l])<m&&(j+sx[l])>=0&&(j+sx[l])<n&&(k+sz[l])>=0&&(k+sz[l])<o){
            mind1[l]=d1[i+sy[l]+(j+sx[l])*m+(k+sz[l])*m*n+index[l]*sz1];
        }
        else{
            mind1[l]=d1[i+j*m+k*m*n+index[l]*sz1];
        }
    }
    
    float minval=mind1[0];//*min_element(mind1,mind1+len2);
    for(int l=1;l<len2;l++){
        
        minval = min(minval,mind1[l]);
    }
    float sumnoise=0.0f;
    for(int l=0;l<len2;l++){
        mind1[l]-=minval;
        sumnoise+=mind1[l];
    }
    
    float noise1=max(sumnoise/(float)len2,1e-6f);
    
    for(int l=0;l<len2;l++){
        mind1[l]/=noise1;
    }
    
    unsigned long accum=0;
    unsigned long tabled1=1;
                
    for(int l=0;l<len2;l++){
        int mind1val=0;
        for(int c=0;c<val-1;c++){
            if(compare[c]>mind1[l])
                mind1val++;    
            //mind1val+=compare[c]>mind1[l]?1:0;
        }
                    
        accum+=tablei[mind1val]*tabled1;
        tabled1*=power;
    }
    
    mindq[i+j*m+k*m*n]=accum;
}

}

void descriptor_GPU(uint64_t* d_mindq,float* d_im1,int m,int n,int o,int qs){
    timeval time1,time2;
    
    //MIND with self-similarity context
    
    int dx[6]={+qs,+qs,-qs,+0,+qs,+0};
    int dy[6]={+qs,-qs,+0,-qs,+0,+qs};
    int dz[6]={0,+0,+qs,+qs,+qs,+qs};
    
    int sx[12]={-qs,+0,-qs,+0,+0,+qs,+0,+0,+0,-qs,+0,+0};
    int sy[12]={+0,-qs,+0,+qs,+0,+0,+0,+qs,+0,+0,+0,-qs};
    int sz[12]={+0,+0,+0,+0,-qs,+0,-qs,+0,-qs,+0,-qs,+0};
    
    int index[12]={0,0,1,1,2,2,3,3,4,4,5,5};
    
    float sigma=0.75;//1.0;//0.75;//1.5;
    int rho=ceil(sigma*1.5)*2+1;
    
    int len1=6;
    const int len2=12;
    
    image_d=12;
    int d=12;
    int sz1=m*n*o;

    

    
    //============== DISTANCES USING BOXFILTER ===================
    float* d1=new float[sz1*len1];
    float* d_d1;
    cudaMalloc((void**)&d_d1,sz1*len1*sizeof(float));

    gettimeofday(&time1, NULL);
    
    for(int l=0;l<len1;l++){
        distances_gpu(d_im1,d1,d_d1,m,n,o,qs,l);
    }


    gettimeofday(&time2, NULL);
    float timeMIND1=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
    gettimeofday(&time1, NULL);

    //quantisation table
    const int val=6;
    
    const unsigned long long power=32;
    unsigned int tablei[6]={0,1,3,7,15,31};
    float* compare = new float[val-1];
    float* d_compare;
    cudaMalloc((void**)&d_compare,(val-1)*sizeof(float));
    
    for(int i=0;i<val-1;i++){
        compare[i]=-log((i+1.5f)/val);
    }

    cudaMemcpy(d_compare,compare,(val-1)*sizeof(float),cudaMemcpyHostToDevice);
    
    dim3 grid(floor_div(m,4),floor_div(n,4),floor_div(o,2));
    dim3 block(8,8,4);
    descriptor_1<<<grid,block>>>(d_mindq,d_compare,d_d1,qs,len2,m,n,o);
    cudaDeviceSynchronize();
    
    gettimeofday(&time2, NULL);
    float timeMIND2=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
    delete d1;
    cudaFree(d_d1);
    cudaFree(d_compare);
    

}


