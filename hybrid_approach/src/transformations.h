/* several functions to interpolate and symmetrise deformations
 calculates Jacobian and harmonic Energy */


void interp3(float* interp,float* input,float* x1,float* y1,float* z1,int m,int n,int o,int m2,int n2,int o2,bool flag){
	for(int k=0;k<o;k++){
		for(int j=0;j<n;j++){
			for(int i=0;i<m;i++){
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
	}
}



void filter1(float* imagein,float* imageout,int m,int n,int o,float* filter,int length,int dim){
	int i,j,k,f;
	int i1,j1,k1;
	int hw=(length-1)/2;
	
	for(i=0;i<(m*n*o);i++){
		imageout[i]=0.0;
	}
	
	for(k=0;k<o;k++){
		for(j=0;j<n;j++){
			for(i=0;i<m;i++){
				for(f=0;f<length;f++){
					//replicate-padding
					if(dim==1)
						imageout[i+j*m+k*m*n]+=filter[f]*imagein[max(min(i+f-hw,m-1),0)+j*m+k*m*n]; 
					if(dim==2)
						imageout[i+j*m+k*m*n]+=filter[f]*imagein[i+max(min(j+f-hw,n-1),0)*m+k*m*n]; 
					if(dim==3)
						imageout[i+j*m+k*m*n]+=filter[f]*imagein[i+j*m+max(min(k+f-hw,o-1),0)*m*n]; 
				}
			}
		}
	}
}

void volfilter(float* imagein,int m,int n,int o,int length,float sigma){
	
	int hw=(length-1)/2;
	int i,j,f;
	float hsum=0;
	float* filter=new float[length];
	for(i=0;i<length;i++){
		filter[i]=exp(-pow((i-hw),2)/(2*pow(sigma,2)));
		hsum=hsum+filter[i];
	}
	for(i=0;i<length;i++){
		filter[i]=filter[i]/hsum;
	}
	float* image1=new float[m*n*o];
    for(i=0;i<m*n*o;i++){
        image1[i]=imagein[i];
    }
    filter1(image1,imagein,m,n,o,filter,length,1);
	filter1(imagein,image1,m,n,o,filter,length,2);
	filter1(image1,imagein,m,n,o,filter,length,3);
	
	delete image1;
	delete filter;	
	
}



float jacobian(float* u1,float* v1,float* w1,int m,int n,int o,int factor){
	
	float factor1=1.0/(float)factor;
	float jmean=0.0;
	float jstd=0.0;
	int i;
	float grad[3]={-0.5,0.0,0.5};
	float* Jac=new float[m*n*o];
	
	float* J11=new float[m*n*o];
	float* J12=new float[m*n*o];
	float* J13=new float[m*n*o];
	float* J21=new float[m*n*o];
	float* J22=new float[m*n*o];
	float* J23=new float[m*n*o];
	float* J31=new float[m*n*o];
	float* J32=new float[m*n*o];
	float* J33=new float[m*n*o];
	
	for(i=0;i<(m*n*o);i++){
		J11[i]=0.0;
		J12[i]=0.0;
		J13[i]=0.0;
		J21[i]=0.0;
		J22[i]=0.0;
		J23[i]=0.0;
		J31[i]=0.0;
		J32[i]=0.0;
		J33[i]=0.0;
	}
	
	float neg=0; float Jmin=1; float Jmax=1; float J;
	float count=0; float frac;
	
	filter1(u1,J11,m,n,o,grad,3,2);
	filter1(u1,J12,m,n,o,grad,3,1);
	filter1(u1,J13,m,n,o,grad,3,3);
	
	filter1(v1,J21,m,n,o,grad,3,2);
	filter1(v1,J22,m,n,o,grad,3,1);
	filter1(v1,J23,m,n,o,grad,3,3);
	
	filter1(w1,J31,m,n,o,grad,3,2);
	filter1(w1,J32,m,n,o,grad,3,1);
	filter1(w1,J33,m,n,o,grad,3,3);
	
	for(i=0;i<(m*n*o);i++){
		J11[i]*=factor1;
		J12[i]*=factor1;
		J13[i]*=factor1;
		J21[i]*=factor1;
		J22[i]*=factor1;
		J23[i]*=factor1;
		J31[i]*=factor1;
		J32[i]*=factor1;
		J33[i]*=factor1;
	}
	
	for(i=0;i<(m*n*o);i++){
		J11[i]+=1.0;
		J22[i]+=1.0;
		J33[i]+=1.0;
	}
	for(i=0;i<(m*n*o);i++){
		J=J11[i]*(J22[i]*J33[i]-J23[i]*J32[i])-J21[i]*(J12[i]*J33[i]-J13[i]*J32[i])+J31[i]*(J12[i]*J23[i]-J13[i]*J22[i]);
		jmean+=J;
		if(J>Jmax)
			Jmax=J;
		if(J<Jmin)
			Jmin=J;
		if(J<0)
			neg++;
		count++;
		Jac[i]=J;
	}
	jmean/=(m*n*o);
	for(int i=0;i<m*n*o;i++){
		jstd+=pow(Jac[i]-jmean,2.0);
	}
	jstd/=(m*n*o-1);
	jstd=sqrt(jstd);
	frac=neg/count;
	cout<<"std(J)="<<round(jstd*100)/100.0;
	//cout<<"Range: ["<<Jmin<<", "<<Jmax<<"]round(jmean*100)/100.0<<
    cout<<" (J<0)="<<round(frac*1e7)/100.0<<"e-7  ";
	delete []Jac;
	
	
	delete []J11;
	delete []J12;
	delete []J13;
	delete []J21;
	delete []J22;
	delete []J23;
	delete []J31;
	delete []J32;
	delete []J33;
	
	return jstd;
	
	
}



void consistentMappingCL(float* u,float* v,float* w,float* u2,float* v2,float* w2,int m,int n,int o,int factor){
    float factor1=1.0/(float)factor;
    float* us=new float[m*n*o];
    float* vs=new float[m*n*o];
    float* ws=new float[m*n*o];
    float* us2=new float[m*n*o];
    float* vs2=new float[m*n*o];
    float* ws2=new float[m*n*o];
    
    for(int i=0;i<m*n*o;i++){
        us[i]=u[i]*factor1; vs[i]=v[i]*factor1;	ws[i]=w[i]*factor1;
        us2[i]=u2[i]*factor1; vs2[i]=v2[i]*factor1;	ws2[i]=w2[i]*factor1;
    }
    
    for(int it=0;it<10;it++){
        interp3(u,us2,us,vs,ws,m,n,o,m,n,o,true);
        interp3(v,vs2,us,vs,ws,m,n,o,m,n,o,true);
        interp3(w,ws2,us,vs,ws,m,n,o,m,n,o,true);
        for(int i=0;i<m*n*o;i++){
            u[i]=0.5*us[i]-0.5*u[i];
            v[i]=0.5*vs[i]-0.5*v[i];
            w[i]=0.5*ws[i]-0.5*w[i];
            
        }
        interp3(u2,us,us2,vs2,ws2,m,n,o,m,n,o,true);
        interp3(v2,vs,us2,vs2,ws2,m,n,o,m,n,o,true);
        interp3(w2,ws,us2,vs2,ws2,m,n,o,m,n,o,true);
        for(int i=0;i<m*n*o;i++){
            u2[i]=0.5*us2[i]-0.5*u2[i];
            v2[i]=0.5*vs2[i]-0.5*v2[i];
            w2[i]=0.5*ws2[i]-0.5*w2[i];
        }
        
        for(int i=0;i<m*n*o;i++){
            us[i]=u[i]; vs[i]=v[i]; ws[i]=w[i];
            us2[i]=u2[i]; vs2[i]=v2[i]; ws2[i]=w2[i];
        }
        
    }
    
    
    for(int i=0;i<m*n*o;i++){
        u[i]*=(float)factor;
        v[i]*=(float)factor;
        w[i]*=(float)factor;
        u2[i]*=(float)factor;
        v2[i]*=(float)factor;
        w2[i]*=(float)factor;
    }
    
    
    delete us; delete vs; delete ws;
    delete us2; delete vs2; delete ws2;
}

__global__ void interp3_gpu(float* interp,float* input,float* x1,float* y1,float* z1,int m,int n,int o,int m2,int n2,int o2,bool flag){
	
	int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;
    if(i<m && j<n && k<o){
//	for(int k=0;k<o;k++){
//		for(int j=0;j<n;j++){
//			for(int i=0;i<m;i++){
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
//		}
//	}
}


__global__ void consistentMappingCL_gpu_1(float* u,float*v,float*w,float* u2,float*v2,float*w2,float* us,float*vs,float*ws,float* us2,float*vs2,float*ws2,float factor1){

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    
    
    us[i]=u[i]*factor1; vs[i]=v[i]*factor1; ws[i]=w[i]*factor1;
    us2[i]=u2[i]*factor1; vs2[i]=v2[i]*factor1; ws2[i]=w2[i]*factor1;
    
    
}

__global__ void consistentMappingCL_gpu_2(float*u,float*v,float*w,float*us,float*vs,float*ws,int m,int n,int o){

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    float tempu;
    float tempv;
    float tempw;
    tempu = u[i];
    tempv = v[i];
    tempw = w[i];
    
    
    u[i]=0.5*us[i]-0.5*tempu;
    v[i]=0.5*vs[i]-0.5*tempv;
    w[i]=0.5*ws[i]-0.5*tempw;
            
    
}

__global__ void consistentMappingCL_gpu_3(float*u,float*v,float*w,float*u2,float*v2,float*w2,float*us,float*vs,float*ws,float*us2,float*vs2,float*ws2){

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    
    us[i]=u[i]; vs[i]=v[i]; ws[i]=w[i];
    us2[i]=u2[i]; vs2[i]=v2[i]; ws2[i]=w2[i];
    

}

__global__ void consistentMappingCL_gpu_4(float*u,float*v,float*w,float*u2,float*v2,float*w2,int factor){

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    
    u[i]*=(float)factor;
    v[i]*=(float)factor;
    w[i]*=(float)factor;
    u2[i]*=(float)factor;
    v2[i]*=(float)factor;
    w2[i]*=(float)factor;
    
}

void consistentMappingCL_GPU(float* u,float* v,float* w,float* u2,float* v2,float* w2,int m,int n,int o,int factor){
    float factor1=1.0/(float)factor;
    float* us;
    float* vs;
    float* ws;
    float* us2;
    float* vs2;
    float* ws2;

    cudaMalloc((void**)&us,m*n*o*sizeof(float));
    cudaMalloc((void**)&vs,m*n*o*sizeof(float));
    cudaMalloc((void**)&ws,m*n*o*sizeof(float));
    cudaMalloc((void**)&us2,m*n*o*sizeof(float));
    cudaMalloc((void**)&vs2,m*n*o*sizeof(float));
    cudaMalloc((void**)&ws2,m*n*o*sizeof(float));

    
    consistentMappingCL_gpu_1<<<(m*n*o),1>>>(u,v,w,u2,v2,w2,us,vs,ws,us2,vs2,ws2,factor1);
    cudaDeviceSynchronize();

    dim3 grid(m,n,o);
    dim3 block(1);

    for(int it=0;it<10;it++){
        
        interp3_gpu<<<grid,block>>>(u,us2,us,vs,ws,m,n,o,m,n,o,true);
        cudaDeviceSynchronize();
        
        interp3_gpu<<<grid,block>>>(v,vs2,us,vs,ws,m,n,o,m,n,o,true);
        cudaDeviceSynchronize();
        
        interp3_gpu<<<grid,block>>>(w,ws2,us,vs,ws,m,n,o,m,n,o,true);
        cudaDeviceSynchronize();

        consistentMappingCL_gpu_2<<<m*n*o,1>>>(u,v,w,us,vs,ws,m,n,o);
        cudaDeviceSynchronize();
        
        interp3_gpu<<<grid,block>>>(u2,us,us2,vs2,ws2,m,n,o,m,n,o,true);
        cudaDeviceSynchronize();
        
        interp3_gpu<<<grid,block>>>(v2,vs,us2,vs2,ws2,m,n,o,m,n,o,true);
        cudaDeviceSynchronize();
        
        interp3_gpu<<<grid,block>>>(w2,ws,us2,vs2,ws2,m,n,o,m,n,o,true);
        cudaDeviceSynchronize();

        consistentMappingCL_gpu_2<<<m*n*o,1>>>(u2,v2,w2,us2,vs2,ws2,m,n,o);
        cudaDeviceSynchronize();
        consistentMappingCL_gpu_3<<<(m*n*o),1>>>(u,v,w,u2,v2,w2,us,vs,ws,us2,vs2,ws2);
        cudaDeviceSynchronize();
        
    }
    
    consistentMappingCL_gpu_4<<<(m*n*o),1>>>(u,v,w,u2,v2,w2,factor);
    cudaDeviceSynchronize();

    
    
    cudaFree(us);cudaFree(vs);cudaFree(ws);
    cudaFree(us2);cudaFree(vs2);cudaFree(ws2);

}


void upsampleDeformationsCL(float* u1,float* v1,float* w1,float* u0,float* v0,float* w0,int m,int n,int o,int m2,int n2,int o2){
    
    
    float scale_m=(float)m/(float)m2;
    float scale_n=(float)n/(float)n2;
    float scale_o=(float)o/(float)o2;
    
    float* x1=new float[m*n*o];
    float* y1=new float[m*n*o];
    float* z1=new float[m*n*o];
    for(int k=0;k<o;k++){
        for(int j=0;j<n;j++){
            for(int i=0;i<m;i++){
                x1[i+j*m+k*m*n]=j/scale_n;
                y1[i+j*m+k*m*n]=i/scale_m;
                z1[i+j*m+k*m*n]=k/scale_o;
            }
        }
    }
    
    interp3(u1,u0,x1,y1,z1,m,n,o,m2,n2,o2,false);
    interp3(v1,v0,x1,y1,z1,m,n,o,m2,n2,o2,false);
    interp3(w1,w0,x1,y1,z1,m,n,o,m2,n2,o2,false);
    
    delete []x1;
    delete []y1;
    delete []z1;
    
    //for(int i=0;i<m2*n2*o2;i++){
    //	u2[i]*=scale_n;
    //	v2[i]*=scale_m;
    //	w2[i]*=scale_o;
    //}
    
}

__global__ void upsampleDeformationsCLgpu(float* x1,float* y1,float* z1,int m,int n,int o,int m2,int n2,int o2){
    
    
    float scale_m=(float)m/(float)m2;
    float scale_n=(float)n/(float)n2;
    float scale_o=(float)o/(float)o2;
    
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;


    if(i < m && j < n && k < o){
        x1[i+j*m+k*m*n]=j/scale_n;
        y1[i+j*m+k*m*n]=i/scale_m;
        z1[i+j*m+k*m*n]=k/scale_o;
    }
        
    
}

__global__ void interp33gpu(float* ux,float* u0,float* vx,float* v0,float* wx,float* w0,float* x1,float* y1,float* z1,int m,int n,int o,int m2,int n2,int o2){
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    //calxyz(x1, y1, z1, scale_m, scale_n, scale_o, m, n, o);

    if(i < m && j < n && k < o){
                int x=floor(x1[i+j*m+k*m*n]); int y=floor(y1[i+j*m+k*m*n]);  int z=floor(z1[i+j*m+k*m*n]); 
                float dx=x1[i+j*m+k*m*n]-x; float dy=y1[i+j*m+k*m*n]-y; float dz=z1[i+j*m+k*m*n]-z;
                
                int ii = min(max(y,0),m2-1);
                int iii = min(max(y+1,0),m2-1);
                int jj = min(max(x,0),n2-1);
                int jjj = min(max(x+1,0),n2-1);
                int kk = min(max(z,0),o2-1);
                int kkk = min(max(z+1,0),o2-1);

                ux[i+j*m+k*m*n]=(1.0-dx)*(1.0-dy)*(1.0-dz)*u0[ii+jj*m2+kk*m2*n2]+(1.0-dx)*dy*(1.0-dz)*u0[iii+jj*m2+kk*m2*n2]+
                dx*(1.0-dy)*(1.0-dz)*u0[ii+jjj*m2+kk*m2*n2]+(1.0-dx)*(1.0-dy)*dz*u0[ii+jj*m2+kkk*m2*n2]+
                dx*dy*(1.0-dz)*u0[iii+jjj*m2+kk*m2*n2]+(1.0-dx)*dy*dz*u0[iii+jj*m2+kkk*m2*n2]+
                dx*(1.0-dy)*dz*u0[ii+jjj*m2+kkk*m2*n2]+dx*dy*dz*u0[iii+jjj*m2+kkk*m2*n2];

                vx[i+j*m+k*m*n]=(1.0-dx)*(1.0-dy)*(1.0-dz)*v0[ii+jj*m2+kk*m2*n2]+(1.0-dx)*dy*(1.0-dz)*v0[iii+jj*m2+kk*m2*n2]+
                dx*(1.0-dy)*(1.0-dz)*v0[ii+jjj*m2+kk*m2*n2]+(1.0-dx)*(1.0-dy)*dz*v0[ii+jj*m2+kkk*m2*n2]+
                dx*dy*(1.0-dz)*v0[iii+jjj*m2+kk*m2*n2]+(1.0-dx)*dy*dz*v0[iii+jj*m2+kkk*m2*n2]+
                dx*(1.0-dy)*dz*v0[ii+jjj*m2+kkk*m2*n2]+dx*dy*dz*v0[iii+jjj*m2+kkk*m2*n2];

                wx[i+j*m+k*m*n]=(1.0-dx)*(1.0-dy)*(1.0-dz)*w0[ii+jj*m2+kk*m2*n2]+(1.0-dx)*dy*(1.0-dz)*w0[iii+jj*m2+kk*m2*n2]+
                dx*(1.0-dy)*(1.0-dz)*w0[ii+jjj*m2+kk*m2*n2]+(1.0-dx)*(1.0-dy)*dz*w0[ii+jj*m2+kkk*m2*n2]+
                dx*dy*(1.0-dz)*w0[iii+jjj*m2+kk*m2*n2]+(1.0-dx)*dy*dz*w0[iii+jj*m2+kkk*m2*n2]+
                dx*(1.0-dy)*dz*w0[ii+jjj*m2+kkk*m2*n2]+dx*dy*dz*w0[iii+jjj*m2+kkk*m2*n2];

            }
}

int floor_div(int num, int denom) {
	if(num%denom==0)
    return num/denom;
	else return num/denom+1;
    
}

void upsampleDeformationsCL_GPU(float* d_u1,float* d_v1,float* d_w1,float* d_u0,float* d_v0,float* d_w0,int m,int n,int o,int m2,int n2,int o2){

	float* d_x1;
    float* d_y1;
    float* d_z1;
    cudaMalloc((void**)&d_x1,m*n*o*sizeof(float));
    cudaMalloc((void**)&d_y1,m*n*o*sizeof(float));
    cudaMalloc((void**)&d_z1,m*n*o*sizeof(float));


	dim3 grid4(floor_div(m,4),floor_div(n,4),floor_div(o,2));
    dim3 block4(4,4,2);
    cout << "m n o : " << m << " " << n << " " << o << endl;
    cout << "m n o floor : " << floor_div(m,4) <<" " << floor_div(n,4) <<" " << floor_div(o,2) << endl;
	upsampleDeformationsCLgpu<<<grid4,block4>>>(d_x1,d_y1,d_z1,m,n,o,m2,n2,o2);
    cudaDeviceSynchronize();

    interp33gpu<<<grid4,block4>>>(d_u1,d_u0,d_v1,d_v0,d_w1,d_w0,d_x1,d_y1,d_z1,m,n,o,m2,n2,o2);
    cudaDeviceSynchronize();

    

}