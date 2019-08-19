#include <stdio.h>
#include <stdlib.h>
#include <time.h>

__device__ float logsumexp(float a, float b)
{
	if(a > b)
	{
		return a + log(1.0+exp(b-a));
	}
	else
	{
		return b + log(1.0+exp(a-b));
	}
}

/*
__global__ void felsensteinfast(const int alphabet, const int numcols, const int numnodes, const int startnode, const int* left, const int* right, const float* transprobs, const float* freqs, float* data, float* collogliks)
{
	__shared__ float liks [8192];
	__shared__ int indices[512];
	__shared__ int revindices[512];
	__shared__ int index = 0;
	for(int i = 0 ; i < 512 ; i++)
	{
		indices[i] = -1;
		revindices[i] = -1;
	}
	__syncthreads();

	int col = blockIdx.x*blockDim.x + threadIdx.x; // column
	int nodeindex = 0;
	int lindex = 0;
	int rindex = 0;
	int leftindex = 0;
	int rightindex = 0;
	float m = 0.0;
	float bsum = 0.0;
	float csum = 0.0;
	float v = 0.0;
	int alphabetplus1 = alphabet+1;
	
	if(col < numcols)
	{		
		for(int node = startnode ; node < numnodes ; node += 1) // post-order tree traversal, calculate starting at tips and ending at root. 
		{	
			nodeindex = threadIdx.x*numnodes*alphabetplus1 + node*alphabetplus1;	
			lindex = left[node]; // left child
			rindex = right[node]; // right child				
			leftindex = threadIdx.x*numnodes*alphabetplus1 + lindex*alphabetplus1;
			rightindex = threadIdx.x*numnodes*alphabetplus1 + rindex*alphabetplus1;
			m = 0.0;
			int li = indices[leftindex];
			if(li == -1)
			{
				for(int a = 0 ; a < alphabet ; a++)
				{
					liks[index*al
				}
				index++;
			}
			int ri = indices[rightindex];

			for(int a = 0 ; a < alphabet ; a++)
			{
				bsum = 0.0;
				csum = 0.0;								
				for(int d = 0 ; d < alphabet ; d++)
				{
					bsum += transprobs[lindex*alphabet*alphabet + a*alphabet + d]*liks[leftindex+d];
					csum += transprobs[rindex*alphabet*alphabet + a*alphabet + d]*liks[rightindex+d];
				}
				v = bsum*csum;
				liks[nodeindex+a] = v;
				if(v > m)
				{
					m = v;	
				}				
			}

			for(int a = 0 ; a < alphabet ; a++)
			{
				liks[nodeindex+a] /= m;
			}				
			liks[nodeindex+alphabet] = log(m) + liks[leftindex+alphabet] + liks[rightindex+alphabet];
		}
	
		float logm = liks[threadIdx.x*numnodes*alphabetplus1 + (numnodes-1)*alphabetplus1 + alphabet];		
		collogliks[col] = 0.0;
		for(int a = 0 ; a < alphabet ; a++)
		{
			collogliks[col] += freqs[a]*liks[threadIdx.x*numnodes*alphabetplus1 + (numnodes-1)*alphabetplus1 + a]; 
		}
		collogliks[col] = log(collogliks[col]) + logm;
	}
}*/

__global__ void felsensteinleaves(const int alphabet, const int numcols, const int numpairs, const int* x, const int* y, float* data, float* res)
{
	int alphabetplus1 = alphabet+1;
	int a = 0;
	int e = 0;
	int f = 0;
	int thread = blockIdx.x*blockDim.x + threadIdx.x;
	if (thread < numpairs)
	{
		for(a = 0 ; a < alphabet ; a++)
		{
			e = a / 4;
			f = a % 4;
			res[thread*alphabetplus1+a] = data[x[thread]*4 + e]*data[y[thread]*4 + f];
		}
		res[thread*alphabetplus1+alphabet] = 0.0;
	}
}

__global__ void felsensteinleaves16(const int numcols, const int numpairs, const int* x, const int* y, float* data, float* res)
{
	int a = 0;
	int thread = blockIdx.x*blockDim.x + threadIdx.x;
	int startcol = 0;
	if (thread < numpairs)
	{
		int xindex = x[thread]*4;
		int yindex = y[thread]*4;
		#pragma unroll
		for(a = 0 ; a < 16 ; a++)
		{
			res[thread*17+a] = data[xindex + (a / 4)]*data[yindex + (a % 4)];
		}
		res[thread*17+16] = 0.0;
	}
}

__global__ void storearr(const int numpairs, const float logw, float* dest, float* src)
{
	int thread = blockIdx.x*blockDim.x + threadIdx.x;
	if (thread < numpairs)
	{
		dest[thread] = logw+src[thread];
	}
}

__global__ void logsumexparr(const int numpairs, const float logw1, float* dest, const float logw2, float* src)
{
	int thread = blockIdx.x*blockDim.x + threadIdx.x;
	if (thread < numpairs)
	{
		dest[thread] = logsumexp(logw1+dest[thread], logw2+src[thread]);
	}
}

__global__ void felsensteinhelper16(const int numcols, const int numpairs, const int numnodes, const int node, const int lindex, const int rindex, const int* leftindices, const int* rightindices, const float* transprobs, const float* left, const float* right, float* res)
{
	float m = 0.0;
	int thread = blockIdx.x*blockDim.x + threadIdx.x;
	float bsum = 0.0;
	float csum = 0.0;
	int lefttransstride = lindex*16*16;
	int lefttrans = 0;
	int righttransstride = rindex*16*16;
	int righttrans = 0;
	float v = 0.0;
	int a = 0;
	int d = 0;	
	int leftindex = leftindices[thread]*17;
	int rightindex = rightindices[thread]*17;
	int nodeindex = thread*17;
	if (thread < numpairs)
	{
		#pragma unroll
		for(a = 0 ; a < 16 ; a++)
		{
			bsum = 0.0;
			csum = 0.0;	
			lefttrans = lefttransstride + a*16;							
			righttrans = righttransstride + a*16;

			#pragma unroll
			for(d = 0 ; d < 16 ; d++)
			{
				bsum += transprobs[lefttrans+d]*left[leftindex+d];
				csum += transprobs[righttrans+d]*right[rightindex+d];
			}
			v = bsum*csum;
			res[nodeindex+a] = v;
			if(v > m)
			{
				m = v;	
			}						
		}
	
		#pragma unroll
		for(a = 0 ; a < 16 ; a++)
		{
			res[nodeindex+a] /= m;
		}
		res[nodeindex+16] = log(m) + left[leftindex+16] + right[rightindex+16];			
	}
}

__global__ void felsensteinleaves16paired(const int numcols, int numratecats, const int numpairs, const int* x, const int* y, float* data, float* res)
{
	int a = 0;
	int thread = blockIdx.x*blockDim.x + threadIdx.x;
	int pair = thread % numpairs;
	if (thread < numratecats*numpairs)
	{
		int xindex = x[pair]*4;
		int yindex = y[pair]*4;
		#pragma unroll
		for(a = 0 ; a < 16 ; a++)
		{
			res[thread*17+a] = data[xindex + (a / 4)]*data[yindex + (a % 4)];
		}
		res[thread*17+16] = 0.0;
	}
}

__global__ void felsensteinhelper16paired(const int numcols, int numratecats, const int numpairs, const int leftnumpairs, const int rightnumpairs, const int numnodes, const int node, const int lindex, const int rindex, const int* leftindices, const int* rightindices, const float* transprobs, const float* left, const float* right, float* res)
{
	float m = 0.0;
	int thread = blockIdx.x*blockDim.x + threadIdx.x;
	float bsum = 0.0;
	float csum = 0.0;	
	int lefttrans = 0;
	int righttrans = 0;
	float v = 0.0;
	int a = 0;
	int d = 0;	
	int ratecat = thread / numpairs;
	int pair = thread % numpairs;
	int leftindex = ratecat*leftnumpairs*17 + leftindices[pair]*17;
	int rightindex = ratecat*rightnumpairs*17 + rightindices[pair]*17;
	int nodeindex = thread*17;
	int lefttranstride = ratecat*numnodes*16*16 + lindex*16*16;
	int righttranstride = ratecat*numnodes*16*16 + rindex*16*16;

	if (thread < numratecats*numpairs)
	{
		#pragma unroll
		for(a = 0 ; a < 16 ; a++)
		{			
			lefttrans =  lefttranstride + a*16;							
			righttrans = righttranstride + a*16;

			bsum = 0.0;
			csum = 0.0;
			#pragma unroll
			for(d = 0 ; d < 16 ; d++)
			{
				bsum += transprobs[lefttrans+d]*left[leftindex+d];
				csum += transprobs[righttrans+d]*right[rightindex+d];
			}
			v = bsum*csum;
			res[nodeindex+a] = v;
			if(v > m)
			{
				m = v;	
			}						
		}
	
		#pragma unroll
		for(a = 0 ; a < 16 ; a++)
		{
			res[nodeindex+a] /= m;
		}
		res[nodeindex+16] = log(m) + left[leftindex+16] + right[rightindex+16];			
	}
}

__global__ void sumfinalpaired(const int alphabet, const int numratecats, const int numpairs, float* freqs, float* rootliks, float* finallogliks)
{
	int alphabetplus1 = alphabet+1;
	int thread = blockIdx.x*blockDim.x + threadIdx.x; // column
	int freqindex = (thread / numpairs)*alphabet;
	int a = 0;
	if(thread < numratecats*numpairs)
	{
		finallogliks[thread] = 0.0;
		for(a = 0 ; a < alphabet ; a++)
		{
			finallogliks[thread] += freqs[freqindex+a]*rootliks[thread*alphabetplus1+a];
		}
		finallogliks[thread] = log(finallogliks[thread]) + rootliks[thread*alphabetplus1+alphabet];
	}
}

__global__ void sumcats(const int numratecats, const int numpairs, float* logweights, float* logliks, float* logfinalliks)
{
	int thread = blockIdx.x*blockDim.x + threadIdx.x; // column
	if(thread < numpairs)
	{
		logfinalliks[thread] = logweights[0] + logliks[0*numpairs + thread];
		for(int r = 1 ; r < numratecats ; r++)
		{
			logfinalliks[thread] = logsumexp(logfinalliks[thread], logweights[r] + logliks[r*numpairs + thread]);
		}		
	}
}

__global__ void felsensteinhelper(const int alphabet, const int numcols, const int numpairs, const int numnodes, const int node, const int lindex, const int rindex, const int* leftindices, const int* rightindices, const float* transprobs, const float* left, const float* right, float* res)
{
	int alphabetplus1 = alphabet+1;
	float m = 0.0;
	int thread = blockIdx.x*blockDim.x + threadIdx.x;
	float bsum = 0.0;
	float csum = 0.0;
	int lefttrans = 0;
	int righttrans = 0;
	float v = 0.0;
	int a = 0;
	int d = 0;	
	int leftindex = 0;
	int rightindex = 0;
	int nodeindex = 0;
	if (thread < numpairs)
	{
		leftindex = leftindices[thread]*alphabetplus1;
		rightindex = rightindices[thread]*alphabetplus1;
		nodeindex = thread*alphabetplus1;
		for(a = 0 ; a < alphabet ; a++)
		{
			bsum = 0.0;
			csum = 0.0;	
			lefttrans = lindex*alphabet*alphabet + a*alphabet;							
			righttrans = rindex*alphabet*alphabet + a*alphabet;
			for(d = 0 ; d < alphabet ; d++)
			{
				bsum += transprobs[lefttrans+d]*left[leftindex+d];
				csum += transprobs[righttrans+d]*right[rightindex+d];
			}
			v = bsum*csum;
			res[nodeindex+a] = v;
			if(v > m)
			{
				m = v;	
			}						
		}
	
		for(a = 0 ; a < alphabet ; a++)
		{
			res[nodeindex+a] /= m;
		}
		res[nodeindex+alphabet] = log(m) + left[leftindex+alphabet] + right[rightindex+alphabet];			
	}
}

__global__ void felsensteindinucleotide(const int alphabet, const int numcols, const int numcategories, const int numpairs, const int numnodes, const int numleaves, const int* left, const int* right, const float* transprobs, const float* freqs, int* pairs, float* data, float* liks, float* collogliks)
{
	int thread = blockIdx.x*blockDim.x + threadIdx.x;	
	int nodeindex = 0;
	int lindex = 0;
	int rindex = 0;
	int leftindex = 0;
	int rightindex = 0;
	float m = 0.0;
	float bsum = 0.0;
	float csum = 0.0;
	float v = 0.0;
	int alphabetplus1 = alphabet+1;
	int node = 0;
	int a = 0;
	int e = 0;
	int f = 0;
	if(thread < numpairs*numcategories)
	{		
		int cat = thread / numpairs;
		int col = thread % numpairs;
		int pair = pairs[col];				
		int col1 = pair / numcols;
		int col2 = pair % numcols;
		int stride = cat*numpairs*numnodes*alphabetplus1 + col*numnodes*alphabetplus1;
		for(node = 0 ; node < numleaves ; node += 1)
		{
			nodeindex = stride + node*alphabetplus1;	
			for(a = 0 ; a < alphabet ; a++)
			{
				e = a / 4;
				f = a % 4;
				liks[nodeindex+a] = data[node*numcols*4 + col1*4 + e]*data[node*numcols*4 + col2*4+f];
			}
			liks[nodeindex+alphabet] = 0.0;
		}
		int stridetrans = cat*numnodes*alphabet*alphabet;
		for(node = numleaves ; node < numnodes ; node += 1)
		{	
			nodeindex = stride + node*alphabetplus1;
			lindex = left[node];
			rindex = right[node];					
			leftindex = stride + lindex*alphabetplus1;
			rightindex = stride + rindex*alphabetplus1;
			m = 0.0;				
			for(a = 0 ; a < alphabet ; a++)
			{
				bsum = 0.0;
				csum = 0.0;								
				for(int d = 0 ; d < alphabet ; d++)
				{
					bsum += transprobs[stridetrans + lindex*alphabet*alphabet + a*alphabet + d]*liks[leftindex+d];
					csum += transprobs[stridetrans + rindex*alphabet*alphabet + a*alphabet + d]*liks[rightindex+d];
				}
				v = bsum*csum;
				liks[nodeindex+a] = v;
				if(v > m)
				{
					m = v;	
				}				
			}

			for(a = 0 ; a < alphabet ; a++)
			{
				liks[nodeindex+a] /= m;
			}				
			liks[nodeindex+alphabet] = log(m) + liks[leftindex+alphabet] + liks[rightindex+alphabet];
		}

		float logm = liks[stride + (numnodes-1)*alphabetplus1 + alphabet];		
		collogliks[thread] = 0.0;
		for(a = 0 ; a < alphabet ; a++)
		{
			collogliks[thread] += freqs[cat*alphabet+a]*liks[stride + (numnodes-1)*alphabetplus1 + a]; 
		}
		collogliks[thread] = log(collogliks[thread]) + logm;
	}
}

__global__ void sumfinal(const int alphabet, const int numpairs, float* freqs, float* rootliks, float* finallogliks)
{
	int alphabetplus1 = alphabet+1;
	int thread = blockIdx.x*blockDim.x + threadIdx.x; // column
	int a = 0;
	if(thread < numpairs)
	{
		finallogliks[thread] = 0.0;
		for(a = 0 ; a < alphabet ; a++)
		{
			finallogliks[thread] += freqs[a]*rootliks[thread*alphabetplus1+a];
		}
		finallogliks[thread] = log(finallogliks[thread]) + rootliks[thread*alphabetplus1+alphabet];
	}
}

__global__ void sumcategories(const int numcategories, const int numcols, const float* catlogprobs, const float* collogliks, float* finallogliks)
{
	int col = blockIdx.x*blockDim.x + threadIdx.x; // column
	if(col < numcols)
	{
		finallogliks[col] = catlogprobs[0] + collogliks[col];
		for(int cat = 1 ; cat < numcategories ; cat++)
		{
			finallogliks[col] = logsumexp(finallogliks[col], catlogprobs[cat] + collogliks[cat*numcols + col]);
		}
	}
}


__global__ void felsensteinfast(const int alphabet, const int numcols, const int numnodes, const int startnode, const int* left, const int* right, const float* transprobs, const float* freqs, float* liks, float* collogliks)
{
	int col = blockIdx.x*blockDim.x + threadIdx.x; // column
	int nodeindex = 0;
	int lindex = 0;
	int rindex = 0;
	int leftindex = 0;
	int rightindex = 0;
	float m = 0.0;
	float bsum = 0.0;
	float csum = 0.0;
	float v = 0.0;
	int alphabetplus1 = alphabet+1;
	if(col < numcols)
	{
		for(int node = startnode ; node < numnodes ; node += 1) // post-order tree traversal, calculate starting at tips and ending at root. 
		{	
			nodeindex = col*numnodes*alphabetplus1 + node*alphabetplus1;	
			lindex = left[node]; // left child
			rindex = right[node]; // right child					
			leftindex = col*numnodes*alphabetplus1 + lindex*alphabetplus1;
			rightindex = col*numnodes*alphabetplus1 + rindex*alphabetplus1;
			m = 0.0;
			for(int a = 0 ; a < alphabet ; a++)
			{
				bsum = 0.0;
				csum = 0.0;								
				for(int d = 0 ; d < alphabet ; d++)
				{
					bsum += transprobs[lindex*alphabet*alphabet + a*alphabet + d]*liks[leftindex+d];
					csum += transprobs[rindex*alphabet*alphabet + a*alphabet + d]*liks[rightindex+d];
				}
				v = bsum*csum;
				liks[nodeindex+a] = v;
				if(v > m)
				{
					m = v;	
				}				
			}

			for(int a = 0 ; a < alphabet ; a++)
			{
				liks[nodeindex+a] /= m;
			}				
			liks[nodeindex+alphabet] = log(m) + liks[leftindex+alphabet] + liks[rightindex+alphabet];
		}
	
		float logm = liks[(col*numnodes*alphabetplus1) + (numnodes-1)*alphabetplus1 + alphabet];		
		collogliks[col] = 0.0;
		for(int a = 0 ; a < alphabet ; a++)
		{
			collogliks[col] += freqs[a]*liks[(col*numnodes*alphabetplus1) + (numnodes-1)*alphabetplus1 + a]; 
		}
		collogliks[col] = log(collogliks[col]) + logm;
	}
}

/*
__global__ void felsenstein(int numcols, int numnodes, int* left, int* right, float* logtransprobs, float* logfreqs, float* logliks, float* collogliks)
{
	int col = blockIdx.x*blockDim.x + threadIdx.x; // column
	if(col < numcols)
	{
		for(int node = 0 ; node < numnodes ; node += 1) // post-order tree traversal, calculate starting at tips and ending at root. 
		{	
			int lognodeindex = col*numnodes*alphabet + node*alphabet;	
			int lindex = left[node]; // left child
			int rindex = right[node]; // right child
			if(lindex == -1 && rindex == -1) // if 'node' is leaf node
			{

			}
			else // if 'node' is internal node
			{			
				int logleftindex = col*numnodes*alphabet + lindex*alphabet;
				int logrightindex = col*numnodes*alphabet + rindex*alphabet;
				for(int a = 0 ; a < alphabet ; a++)
				{
					float bsum = -1e10;
					float csum = -1e10;								
					for(int d = 0 ; d < alphabet ; d++)
					{
						int transindex = node*alphabet*alphabet + a*alphabet + d;
						float logtransprob = logtransprobs[transindex]; // transition probability
						bsum = logsumexp(bsum, logtransprob+logliks[logleftindex+d]);
						csum = logsumexp(csum, logtransprob+logliks[logrightindex+d]);
					}
					logliks[lognodeindex+a] = bsum+csum;
				}
			}
		}
	
		collogliks[col] = -1e10;
		for(int a = 0 ; a < alphabet ; a++)
		{
			collogliks[col] = logsumexp(collogliks[col], logfreqs[a] + logliks[col*numnodes*alphabet + (numnodes-1)*alphabet + a]); 
		}
	}
}

float randfloat()
{
	return ((float)rand()/(float)(RAND_MAX)) * 1.0;
}

int main(void)
{
	int code = 0;

	int numcols = 10000;
	int numnodes = 250;
	
	int *left, *right, *d_left, *d_right;
	left = (int*)malloc(numnodes*sizeof(int));
	right = (int*)malloc(numnodes*sizeof(int));
	for(int node = 0 ; node < numnodes ; node += 1)
	{
		//left[node] = rand() % numnodes;
		//right[node] = rand() % numnodes;
	}
	code = cudaMalloc(&d_left, numnodes*sizeof(int));
	printf("A %d\n", code);
	code = cudaMalloc(&d_right, numnodes*sizeof(int));
	printf("B %d\n", code);
 	code = cudaMemcpy(d_left, left, numnodes*sizeof(int), cudaMemcpyHostToDevice);
	printf("C %d\n", code);
	code = cudaMemcpy(d_right, right, numnodes*sizeof(int), cudaMemcpyHostToDevice);
	printf("D %d\n", code);

	float *logtransprobs, *d_logtransprobs;
	logtransprobs = (float*)malloc(numnodes*alphabet*alphabet*sizeof(float));
	for(int node = 0 ; node < numnodes ; node += 1)
	{
		for(int a = 0 ; a < alphabet ; a++)
		{
			float sum = 0.0;
			for(int b = 0 ; b < alphabet ; b++)
			{
				logtransprobs[node*alphabet*alphabet+a*alphabet+b] = randfloat();
				sum += logtransprobs[node*alphabet*alphabet+a*alphabet+b];
			}
			for(int b = 0 ; b < alphabet ; b++)
			{
				logtransprobs[node*alphabet*alphabet+a*alphabet+b] = log(logtransprobs[node*alphabet*alphabet+a*alphabet+b]/sum);
			}
			
		}
	}

	float *logfreqs, *d_logfreqs;
	logfreqs = (float*)malloc(alphabet*sizeof(float));
	for(int a = 0 ; a < alphabet ; a++)
	{
		logfreqs[a] = log(randfloat()); 
	}
	cudaMalloc(&d_logfreqs, alphabet*sizeof(float));
	cudaMemcpy(d_logfreqs, logfreqs, alphabet*sizeof(float), cudaMemcpyHostToDevice);

	code = cudaMalloc(&d_logtransprobs, numnodes*alphabet*alphabet*sizeof(float));
	printf("E %d\n", code);
	code = cudaMemcpy(d_logtransprobs, logtransprobs, numnodes*alphabet*alphabet*sizeof(float), cudaMemcpyHostToDevice);
	printf("F %d\n", code);
	
	float *logliks, *d_logliks;
	logliks = (float*)malloc(numcols*numnodes*alphabet*sizeof(float));
	for(int col = 0 ; col < numcols ; col += 1)
	{
		for(int node = 0 ; node < numnodes ; node += 1)
		{
			int lognodeindex = col*numnodes*alphabet + node*alphabet;
			for(int a = 0 ; a < alphabet ; a += 1) 
			{
				logliks[lognodeindex+a] = -1e10;
			}
			logliks[lognodeindex + (rand() % alphabet)] = 0.0;
		}
	}
	cudaMalloc(&d_logliks, numcols*numnodes*alphabet*sizeof(float));
	cudaMemcpy(d_logliks, logliks, numcols*numnodes*alphabet*sizeof(float), cudaMemcpyHostToDevice);

	float *collogliks, *d_collogliks;
	collogliks = (float*)malloc(numcols*sizeof(float));
	for(int col = 0 ; col < numcols ; col++)
	{
		collogliks[col] = 20.0;
	}
	cudaMalloc(&d_collogliks, numcols*sizeof(float));
	cudaMemcpy(d_collogliks, collogliks, numcols*sizeof(float), cudaMemcpyHostToDevice);
	felsenstein<<<(numcols+255)/256, 256>>>(numcols, numnodes, d_left, d_right, d_logtransprobs, d_logfreqs, d_logliks, d_collogliks);
	code = cudaMemcpy(collogliks,  d_collogliks, numcols*sizeof(float), cudaMemcpyDeviceToHost);
	printf("finished %d\n", code);
	
	for(int col = 0 ; col < numcols ; col++)
	{
		printf("%d\t%lf\n",col,collogliks[col]);
	}
	

	cudaFree(d_left);
	cudaFree(d_right);
	cudaFree(d_logtransprobs);
	cudaFree(d_logliks);
	cudaFree(d_collogliks);
	free(left);
	free(right);
	free(logtransprobs);
	free(logliks);
	free(collogliks);
	return 0;
}*/
