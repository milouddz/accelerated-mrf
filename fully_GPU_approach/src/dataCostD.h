void interp3xyz(float* datai,float* data,float* datax,float* datay,int len1,int len2){
    //x-interp
    for(int k=0;k<len1;k++){
        for(int j=0;j<len2;j++){
            int j2=(j+1)/2;
            if(j%2==1){
                for(int i=0;i<len1;i++){
                    datax[i+j*len1+k*len1*len2]=data[i+j2*len1+k*len1*len1];
                }
            }
            else
            for(int i=0;i<len1;i++){
                datax[i+j*len1+k*len1*len2]=0.5*(data[i+j2*len1+k*len1*len1]+data[i+(j2+1)*len1+k*len1*len1]);
            }
        }
    }
    
    
    //y-interp
    for(int k=0;k<len1;k++){
        for(int j=0;j<len2;j++){
            for(int i=0;i<len2;i++){
                int i2=(i+1)/2;
                if(i%2==1)
                datay[i+j*len2+k*len2*len2]=datax[i2+j*len1+k*len1*len2];
                else
                datay[i+j*len2+k*len2*len2]=0.5*(datax[i2+j*len1+k*len1*len2]+datax[i2+1+j*len1+k*len1*len2]);
            }
        }
    }
    
    //z-interp
    for(int k=0;k<len2;k++){
        int k2=(k+1)/2;
        if(k%2==1){
            for(int j=0;j<len2;j++){
                for(int i=0;i<len2;i++){
                    datai[i+j*len2+k*len2*len2]=datay[i+j*len2+k2*len2*len2];
                }
            }
            
        }
        else{
            for(int j=0;j<len2;j++){
                for(int i=0;i<len2;i++){
                    datai[i+j*len2+k*len2*len2]=0.5*(datay[i+j*len2+k2*len2*len2]+datay[i+j*len2+(k2+1)*len2*len2]);
                }
            }
        }
    }
    
}

void interp3xyzB(float* datai,float* data,float* datax,float* datay,int len1,int len2){
    //x-interp
    for(int k=0;k<len1;k++){
        for(int j=0;j<len2;j++){
            int j2=(j+1)/2;
            if(j%2==0){
                for(int i=0;i<len1;i++){
                    datax[i+j*len1+k*len1*len2]=data[i+j2*len1+k*len1*len1];
                }
            }
            else
            for(int i=0;i<len1;i++){
                datax[i+j*len1+k*len1*len2]=0.5*(data[i+j2*len1+k*len1*len1]+data[i+(j2-1)*len1+k*len1*len1]);
            }
        }
    }
    
    
    //y-interp
    for(int k=0;k<len1;k++){
        for(int j=0;j<len2;j++){
            for(int i=0;i<len2;i++){
                int i2=(i+1)/2;
                if(i%2==0)
                datay[i+j*len2+k*len2*len2]=datax[i2+j*len1+k*len1*len2];
                else
                datay[i+j*len2+k*len2*len2]=0.5*(datax[i2+j*len1+k*len1*len2]+datax[i2-1+j*len1+k*len1*len2]);
            }
        }
    }
    
    //z-interp
    for(int k=0;k<len2;k++){
        int k2=(k+1)/2;
        if(k%2==0){
            for(int j=0;j<len2;j++){
                for(int i=0;i<len2;i++){
                    datai[i+j*len2+k*len2*len2]=datay[i+j*len2+k2*len2*len2];
                }
            }
            
        }
        else{
            for(int j=0;j<len2;j++){
                for(int i=0;i<len2;i++){
                    datai[i+j*len2+k*len2*len2]=0.5*(datay[i+j*len2+k2*len2*len2]+datay[i+j*len2+(k2-1)*len2*len2]);
                }
            }
        }
    }
    
}


void dataCostCL(uint64_t* data,uint64_t* data2,float* results,int m,int n,int o,int len2,int step1,int hw,float quant,float alpha,int randnum){
    cout<<"d"<<flush;
    
    int len=hw*2+1;
    len2=pow(hw*2+1,3);
    
    int sz=m*n*o;
    int m1=m/step1; int n1=n/step1; int o1=o/step1;
    int sz1=m1*n1*o1;
    
    //cout<<"len2: "<<len2<<" sz1= "<<sz1<<"\n";
    
    
    
    int quant2=quant;
    
    //const int hw2=hw*quant2; == pad1
    
    int pad1=quant2*hw; int pad2=pad1*2;
    
    int mp=m+pad2; int np=n+pad2; int op=o+pad2;
    int szp=mp*np*op;
    uint64_t* data2p=new uint64_t[szp];
    
    for(int k=0;k<op;k++){
        for(int j=0;j<np;j++){
            for(int i=0;i<mp;i++){
                data2p[i+j*mp+k*mp*np]=data2[max(min(i-pad1,m-1),0)+max(min(j-pad1,n-1),0)*m+max(min(k-pad1,o-1),0)*m*n];
            }
        }
    }
    
    
    int skipz=1; int skipx=1; int skipy=1;
    if(step1>4){
        if(randnum>0){
            skipz=2; skipx=2;
        }
        if(randnum>1){
            skipy=2;
        }
    }
    if(randnum>1&step1>7){
        skipz=3; skipx=3; skipy=3;
    }
    if(step1==4&randnum>1)
    skipz=2;
    
    
    float maxsamp=ceil((float)step1/(float)skipx)*ceil((float)step1/(float)skipz)*ceil((float)step1/(float)skipy);
    //printf("randnum: %d, maxsamp: %d ",randnum,(int)maxsamp);
    
    
    float alphai=(float)step1/(alpha*(float)quant);
    
    float alpha1=0.5*alphai/(float)(maxsamp);
    
    //unsigned long buffer[1000];
    
#pragma omp parallel for
    for(int z=0;z<o1;z++){
        for(int x=0;x<n1;x++){
            for(int y=0;y<m1;y++){
                int z1=z*step1; int x1=x*step1; int y1=y*step1;
                /*for(int k=0;k<step1;k++){
                    for(int j=0;j<step1;j++){
                        for(int i=0;i<step1;i++){
                            buffer[i+j*step1+k*step1*step1]=data[i+y1+(j+x1)*m+(k+z1)*m*n];
                        }
                    }
                }*/
                
                for(int l=0;l<len2;l++){
                    int out1=0;
                    int zs=l/(len*len); int xs=(l-zs*len*len)/len; int ys=l-zs*len*len-xs*len;
                    zs*=quant; xs*=quant; ys*=quant;
                    int x2=xs+x1; int z2=zs+z1; int y2=ys+y1;
                    for(int k=0;k<step1;k+=skipz){
                        for(int j=0;j<step1;j+=skipx){
                            for(int i=0;i<step1;i+=skipy){
                                //unsigned int t=buffer[i+j*STEP+k*STEP*STEP]^buf2p[i+j*mp+k*mp*np];
                                //out1+=(wordbits[t&0xFFFF]+wordbits[t>>16]);
                                uint64_t t1=data[i+y1+(j+x1)*m+(k+z1)*m*n];//buffer[i+j*step1+k*step1*step1];
                                uint64_t t2=data2p[i+j*mp+k*mp*np+(y2+x2*mp+z2*mp*np)];
                                out1+=__builtin_popcountll(t1^t2);
                            }
                        }
                    }
                    results[(y+x*m1+z*m1*n1)*len2+l]=out1*alpha1;
                    
                }
                
            }
        }
    }
    
    
    delete data2p;
    
    return;
    
    
}

__global__ void dataCostCL_1(uint64_t* data2,uint64_t* data2p,int m,int n,int o,int mp,int np,int op,int pad1){

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    //for(int k=0;k<op;k++){
        //for(int j=0;j<np;j++){
            //for(int i=0;i<mp;i++){
            if(i<mp && j<np && k<op)
                data2p[i+j*mp+k*mp*np]=data2[max(min(i-pad1,m-1),0)+max(min(j-pad1,n-1),0)*m+max(min(k-pad1,o-1),0)*m*n];
            //}
        //}
    //}
    
}


__global__ void dataCostCL_2(uint64_t* data,
                          uint64_t* data2p,
                          float* results,int step1,int len2,int len,float quant,float alpha1,int mp,int np,
                          int m1,int n1,int o1,int m,int n,int z,int skipx,int skipy,int skipz) {
   
    int x = threadIdx.x + blockIdx.x * blockDim.x;
    int y = threadIdx.y + blockIdx.y * blockDim.y;
    int l = threadIdx.z + blockIdx.z * blockDim.z;
    
    int z1=z*step1; int x1=x*step1; int y1=y*step1;
                
         
                //for(int l=0;l<len2;l++){
                if(x<m1 && y<n1 && l<len2){
                    int out1=0;
                    int zs=l/(len*len); int xs=(l-zs*len*len)/len; int ys=l-zs*len*len-xs*len;
                    zs*=quant; xs*=quant; ys*=quant;
                    int x2=xs+x1; int z2=zs+z1; int y2=ys+y1;
                    
                    for(int k=0;k<step1;k+=skipz){
                        for(int j=0;j<step1;j+=skipx){
                            for(int i=0;i<step1;i+=skipy){
                                out1+=__popcll(data[i+y1+(j+x1)*m+(k+z1)*m*n]^data2p[i+j*mp+k*mp*np+(y2+x2*mp+z2*mp*np)]);
                            }
                        }
                    }
                    results[(y+x*m1+z*m1*n1)*len2+l]=out1*alpha1;
                }
    
}


void dataCostCL_GPU(uint64_t* d_data,uint64_t* d_data2,float* d_results,int m,int n,int o,int len2,int step1,int hw,float quant,float alpha,int randnum){
    cout<<"d"<<flush;
    
    int len=hw*2+1;
    len2=pow(hw*2+1,3);
    
    int sz=m*n*o;
    int m1=m/step1; int n1=n/step1; int o1=o/step1;
    int sz1=m1*n1*o1;
    
    //cout<<"len2: "<<len2<<" sz1= "<<sz1<<"\n";
    
    
    
    int quant2=quant;
    
    //const int hw2=hw*quant2; == pad1
    
    int pad1=quant2*hw; int pad2=pad1*2;
    
    int mp=m+pad2; int np=n+pad2; int op=o+pad2;
    int szp=mp*np*op;
    unsigned long* data2p=new unsigned long[szp];
    uint64_t* d_data2p;
    cudaMalloc((void**)&d_data2p,szp*sizeof(uint64_t));

    /*
    for(int k=0;k<op;k++){
        for(int j=0;j<np;j++){
            for(int i=0;i<mp;i++){
                data2p[i+j*mp+k*mp*np]=data2[max(min(i-pad1,m-1),0)+max(min(j-pad1,n-1),0)*m+max(min(k-pad1,o-1),0)*m*n];
            }
        }
    }
    */
    dim3 grid(mp,np,op);
    dim3 block(1);
    dataCostCL_1<<<grid,block>>>(d_data2,d_data2p,m,n,o,mp,np,op,pad1);
    cudaDeviceSynchronize();
    cudaMemcpy(data2p,d_data2p,szp*sizeof(uint64_t),cudaMemcpyDeviceToHost);

    
    int skipz=1; int skipx=1; int skipy=1;
    if(step1>4){
        if(randnum>0){
            skipz=2; skipx=2;
        }
        if(randnum>1){
            skipy=2;
        }
    }
    if(randnum>1&step1>7){
        skipz=3; skipx=3; skipy=3;
    }
    if(step1==4&randnum>1)
    skipz=2;
    
    
    float maxsamp=ceil((float)step1/(float)skipx)*ceil((float)step1/(float)skipz)*ceil((float)step1/(float)skipy);
    //printf("randnum: %d, maxsamp: %d ",randnum,(int)maxsamp);
    
    
    float alphai=(float)step1/(alpha*(float)quant);
    
    float alpha1=0.5*alphai/(float)(maxsamp);
    
    //unsigned long buffer[1000];
    /*
#pragma omp parallel for
    for(int z=0;z<o1;z++){
        for(int x=0;x<n1;x++){
            for(int y=0;y<m1;y++){
                int z1=z*step1; int x1=x*step1; int y1=y*step1;
                
                for(int l=0;l<len2;l++){
                    int out1=0;
                    int zs=l/(len*len); int xs=(l-zs*len*len)/len; int ys=l-zs*len*len-xs*len;
                    zs*=quant; xs*=quant; ys*=quant;
                    int x2=xs+x1; int z2=zs+z1; int y2=ys+y1;
                    for(int k=0;k<step1;k+=skipz){
                        for(int j=0;j<step1;j+=skipx){
                            for(int i=0;i<step1;i+=skipy){
                                //unsigned int t=buffer[i+j*STEP+k*STEP*STEP]^buf2p[i+j*mp+k*mp*np];
                                //out1+=(wordbits[t&0xFFFF]+wordbits[t>>16]);
                                unsigned long t1=data[i+y1+(j+x1)*m+(k+z1)*m*n];//buffer[i+j*step1+k*step1*step1];
                                unsigned long t2=data2p[i+j*mp+k*mp*np+(y2+x2*mp+z2*mp*np)];
                                out1+=__builtin_popcountll(t1^t2);
                            }
                        }
                    }
                    results[(y+x*m1+z*m1*n1)*len2+l]=out1*alpha1;
                    
                }
                
            }
        }
    }
    */
    /*
    int t;

    //cout << "m1+ n1 " << ((n1+8)/8) << endl;
    if(m1<=64)t=64;
        

    if(m1<=32)t=32;

    dim3 grid1(t/8,t/8,len2/len);
    dim3 block1(8,8,len);
    */
    cout << "m1 : " << m1 << endl;
    dim3 grid1(m1,n1,len2/len);
    dim3 block1(1,1,len);
    
    for(int z=0;z<o1;z++){
    dataCostCL_2<<<grid1,block1>>>(d_data,d_data2p,d_results,step1,len2,len,quant,alpha1,mp,np,m1,
        n1,o1,m,n,z,skipx,skipy,skipz);
    //cudaDeviceSynchronize();
    }


    delete data2p;
    
    return;
    
    
}

void warpImageCL(float* warped,float* im1,float* im1b,float* u1,float* v1,float* w1){
    int m=image_m;
    int n=image_n;
    int o=image_o;
    int sz=m*n*o;
    
    float ssd=0;
    float ssd0=0;
    float ssd2=0;
    
    interp3(warped,im1,u1,v1,w1,m,n,o,m,n,o,true);
    
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            for(int k=0;k<o;k++){
                ssd+=pow(im1b[i+j*m+k*m*n]-warped[i+j*m+k*m*n],2);
                ssd0+=pow(im1b[i+j*m+k*m*n]-im1[i+j*m+k*m*n],2);
            }
        }
    }
    
    ssd/=m*n*o;
    ssd0/=m*n*o;
    SSD0=ssd0;
    SSD1=ssd;
    
}

void warpAffineS(short* warped,short* input,float* X,float* u1,float* v1,float* w1){
    int m=image_m;
    int n=image_n;
    int o=image_o;
    int sz=m*n*o;
    for(int k=0;k<o;k++){
        for(int j=0;j<n;j++){
            for(int i=0;i<m;i++){
                
                float y1=(float)i*X[0]+(float)j*X[1]+(float)k*X[2]+(float)X[3]+v1[i+j*m+k*m*n];
                float x1=(float)i*X[4]+(float)j*X[5]+(float)k*X[6]+(float)X[7]+u1[i+j*m+k*m*n];
                float z1=(float)i*X[8]+(float)j*X[9]+(float)k*X[10]+(float)X[11]+w1[i+j*m+k*m*n];
                int x=round(x1); int y=round(y1);  int z=round(z1);
                
                //if(y>=0&x>=0&z>=0&y<m&x<n&z<o){
                    warped[i+j*m+k*m*n]=input[min(max(y,0),m-1)+min(max(x,0),n-1)*m+min(max(z,0),o-1)*m*n];
                //}
                //else{
                //    warped[i+j*m+k*m*n]=0;
                //}
            }
        }
    }
    
    
}
void warpAffine(float* warped,float* input,float* im1b,float* X,float* u1,float* v1,float* w1){
    int m=image_m;
    int n=image_n;
    int o=image_o;
    int sz=m*n*o;
    
    float ssd=0;
    float ssd0=0;
    float ssd2=0;
    
    for(int k=0;k<o;k++){
        for(int j=0;j<n;j++){
            for(int i=0;i<m;i++){
                
                float y1=(float)i*X[0]+(float)j*X[1]+(float)k*X[2]+(float)X[3]+v1[i+j*m+k*m*n];
                float x1=(float)i*X[4]+(float)j*X[5]+(float)k*X[6]+(float)X[7]+u1[i+j*m+k*m*n];
                float z1=(float)i*X[8]+(float)j*X[9]+(float)k*X[10]+(float)X[11]+w1[i+j*m+k*m*n];
                int x=floor(x1); int y=floor(y1);  int z=floor(z1);
                float dx=x1-x; float dy=y1-y; float dz=z1-z;
                
                
                warped[i+j*m+k*m*n]=(1.0-dx)*(1.0-dy)*(1.0-dz)*input[min(max(y,0),m-1)+min(max(x,0),n-1)*m+min(max(z,0),o-1)*m*n]+
                (1.0-dx)*dy*(1.0-dz)*input[min(max(y+1,0),m-1)+min(max(x,0),n-1)*m+min(max(z,0),o-1)*m*n]+
                dx*(1.0-dy)*(1.0-dz)*input[min(max(y,0),m-1)+min(max(x+1,0),n-1)*m+min(max(z,0),o-1)*m*n]+
                (1.0-dx)*(1.0-dy)*dz*input[min(max(y,0),m-1)+min(max(x,0),n-1)*m+min(max(z+1,0),o-1)*m*n]+
                dx*dy*(1.0-dz)*input[min(max(y+1,0),m-1)+min(max(x+1,0),n-1)*m+min(max(z,0),o-1)*m*n]+
                (1.0-dx)*dy*dz*input[min(max(y+1,0),m-1)+min(max(x,0),n-1)*m+min(max(z+1,0),o-1)*m*n]+
                dx*(1.0-dy)*dz*input[min(max(y,0),m-1)+min(max(x+1,0),n-1)*m+min(max(z+1,0),o-1)*m*n]+
                dx*dy*dz*input[min(max(y+1,0),m-1)+min(max(x+1,0),n-1)*m+min(max(z+1,0),o-1)*m*n];
            }
        }
    }
    
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            for(int k=0;k<o;k++){
                ssd+=pow(im1b[i+j*m+k*m*n]-warped[i+j*m+k*m*n],2);
                ssd0+=pow(im1b[i+j*m+k*m*n]-input[i+j*m+k*m*n],2);
            }
        }
    }
    
    ssd/=m*n*o;
    ssd0/=m*n*o;
    SSD0=ssd0;
    SSD1=ssd;
    
    
}

__global__ void warpAffine_GPU(float* warped,float* input,float* im1b,float* X,float* u1,float* v1,float* w1,
                int m, int n, int o,float* d_SSD0,float* d_SSD1){
    //int m=image_m;
    //int n=image_n;
    //int o=image_o;
    int sz = m*n*o;
    
    //float temp1[sz];
    //float ssd0=0;
    //float ssd2=0;
    
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    
    if(i < m && j < n && k < o){            
                float y1=(float)i*X[0]+(float)j*X[1]+(float)k*X[2]+(float)X[3]+v1[i+j*m+k*m*n];
                float x1=(float)i*X[4]+(float)j*X[5]+(float)k*X[6]+(float)X[7]+u1[i+j*m+k*m*n];
                float z1=(float)i*X[8]+(float)j*X[9]+(float)k*X[10]+(float)X[11]+w1[i+j*m+k*m*n];
                int x=floor(x1); int y=floor(y1);  int z=floor(z1);
                float dx=x1-x; float dy=y1-y; float dz=z1-z;
                
                //printf("i : %i, j : %i, k : %i", i , j, k);
                warped[i+j*m+k*m*n]=(1.0-dx)*(1.0-dy)*(1.0-dz)*input[min(max(y,0),m-1)+min(max(x,0),n-1)*m+min(max(z,0),o-1)*m*n]+
                (1.0-dx)*dy*(1.0-dz)*input[min(max(y+1,0),m-1)+min(max(x,0),n-1)*m+min(max(z,0),o-1)*m*n]+
                dx*(1.0-dy)*(1.0-dz)*input[min(max(y,0),m-1)+min(max(x+1,0),n-1)*m+min(max(z,0),o-1)*m*n]+
                (1.0-dx)*(1.0-dy)*dz*input[min(max(y,0),m-1)+min(max(x,0),n-1)*m+min(max(z+1,0),o-1)*m*n]+
                dx*dy*(1.0-dz)*input[min(max(y+1,0),m-1)+min(max(x+1,0),n-1)*m+min(max(z,0),o-1)*m*n]+
                (1.0-dx)*dy*dz*input[min(max(y+1,0),m-1)+min(max(x,0),n-1)*m+min(max(z+1,0),o-1)*m*n]+
                dx*(1.0-dy)*dz*input[min(max(y,0),m-1)+min(max(x+1,0),n-1)*m+min(max(z+1,0),o-1)*m*n]+
                dx*dy*dz*input[min(max(y+1,0),m-1)+min(max(x+1,0),n-1)*m+min(max(z+1,0),o-1)*m*n];

                atomicAdd(d_SSD1,powf(im1b[i+j*m+k*m*n]-warped[i+j*m+k*m*n],2)/sz);
                atomicAdd(d_SSD0,powf(im1b[i+j*m+k*m*n]-input[i+j*m+k*m*n],2)/sz);
                
                
    }
    
}



__global__ void ssd_gpu(float* d_im1,float* d_im1b,float* d_warped1,int m,int n,int o,float* d_SSD0,float* d_SSD1){

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;
    int sz = m*n*o;
    //*d_SSD0 = 0;
    //*d_SSD1 = 0;

    //dim3 grid(m,n,o);
    //dim3 block(1);

    //for(int i=0;i<m;i++){
        //for(int j=0;j<n;j++){
            //for(int k=0;k<o;k++){
            if(i<m && j<n && k<o){
                atomicAdd(d_SSD1,powf(d_im1b[i+j*m+k*m*n]-d_warped1[i+j*m+k*m*n],2)/sz);
                atomicAdd(d_SSD0,powf(d_im1b[i+j*m+k*m*n]-d_im1[i+j*m+k*m*n],2)/sz);
            }
            //}
        //}
    //}
    
    

}

__global__ void interp3gpu(float* interp,float* input,float* x1,float* y1,float* z1,int m,int n,int o,int m2,int n2,int o2,bool flag){
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    //calxyz(x1, y1, z1, scale_m, scale_n, scale_o, m, n, o);

    if(i < m && j < n && k < o){
                int x=floor(x1[i+j*m+k*m*n]); int y=floor(y1[i+j*m+k*m*n]);  int z=floor(z1[i+j*m+k*m*n]); 
                float dx=x1[i+j*m+k*m*n]-x; float dy=y1[i+j*m+k*m*n]-y; float dz=z1[i+j*m+k*m*n]-z;
                if(flag){
                    x+=j; y+=i; z+=k;
                }
                interp[i+j*m+k*m*n]=(1.0-dx)*(1.0-dy)*(1.0-dz)*input[min(max(y,0),m2-1)+min(max(x,0),n2-1)*m2+min(max(z,0),o2-1)*m2*n2]+
                (1.0-dx)*dy*(1.0-dz)*input[min(max(y+1,0),m2-1)+min(max(x,0),n2-1)*m2+min(max(z,0),o2-1)*m2*n2]+
                dx*(1.0-dy)*(1.0-dz)*input[min(max(y,0),m2-1)+min(max(x+1,0),n2-1)*m2+min(max(z,0),o2-1)*m2*n2]+
                (1.0-dx)*(1.0-dy)*dz*input[min(max(y,0),m2-1)+min(max(x,0),n2-1)*m2+min(max(z+1,0),o2-1)*m2*n2]+
                dx*dy*(1.0-dz)*input[min(max(y+1,0),m2-1)+min(max(x+1,0),n2-1)*m2+min(max(z,0),o2-1)*m2*n2]+
                (1.0-dx)*dy*dz*input[min(max(y+1,0),m2-1)+min(max(x,0),n2-1)*m2+min(max(z+1,0),o2-1)*m2*n2]+
                dx*(1.0-dy)*dz*input[min(max(y,0),m2-1)+min(max(x+1,0),n2-1)*m2+min(max(z+1,0),o2-1)*m2*n2]+
                dx*dy*dz*input[min(max(y+1,0),m2-1)+min(max(x+1,0),n2-1)*m2+min(max(z+1,0),o2-1)*m2*n2];
            }
}
    


void warpImageCL_GPU(float* d_warped,float* d_im1,float* d_im1b,float* d_u1,float* d_v1,float* d_w1,int m,int n,int o,float* d_SSD0,float*d_SSD1){
    
    int sz=m*n*o;
    
    dim3 grid(floor_div(m,4),floor_div(n,4),floor_div(o,2));
    dim3 block(8,8,4);
    interp3gpu<<<grid,block>>>(d_warped,d_im1,d_u1,d_v1,d_w1,m,n,o,m,n,o,true);
    cudaDeviceSynchronize();
    /*
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            for(int k=0;k<o;k++){
                *d_SSD1+=powf(d_im1b[i+j*m+k*m*n]-d_warped[i+j*m+k*m*n],2);
                *d_SSD0+=powf(d_im1b[i+j*m+k*m*n]-d_im1[i+j*m+k*m*n],2);
            }
        }
    }
    */
    //ssd/=m*n*o;
    //ssd0/=m*n*o;
    
    ssd_gpu<<<grid,block>>>(d_im1,d_im1b,d_warped,m,n,o,d_SSD0,d_SSD1);
    cudaDeviceSynchronize();
    //*d_SSD0/=m*n*o;
    //*d_SSD1/=m*n*o;
    
    
}
