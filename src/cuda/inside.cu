#include <stdio.h>
#include <stdlib.h>
#include <time.h>

__device__ float logsumexp(float a, float b)
{

	if(a <= -1e20f)
	{
		return b;
	}
	else
	if(b <= -1e20f)
	{
		return a;
	}

	/*float diff = a-b;
	if (diff < -20.0f)
	{
		return b;
	}
	else
	if (diff > 20.0f)
	{
		return a;
	}*/
	

	if(a > b)
	{
		return a + log(1.0f+exp(b-a));
	}
	else
	{
		return b + log(1.0f+exp(a-b));
	}
}

__device__ float safeadd(float a, float b)
{
	if(a <= -1e20f)
	{
		return b;
	}
	else
	if(b <= -1e20f)
	{
		return a;
	}
	return a+b;
}

/*
push!(rules, Rule('S', "LS",0.868534, 3))
    push!(rules, Rule('S', "s",0.117609877998, 1))
    push!(rules, Rule('S', "dFd",0.013856122002, 2))
    push!(rules, Rule('F', "dFd",0.787640, 2))
    push!(rules, Rule('F', "LS",0.21236, 3))
    push!(rules, Rule('L', "s",0.894603, 1))
    push!(rules, Rule('L', "dFd",0.105397, 2))
    type1rules = Rule[rules[2],rules[6]]
    type2rules = Rule[rules[3],rules[4],rules[7]]
    type3rules = Rule[rules[1],rules[5]]
    ruleindex = Dict('S' => 1, 'F' => 2, 'L' => 3)
*/

__constant__ int S = 0;
__constant__ int F = 1;
__constant__ int L = 2;
	
__constant__ float r1logprob = -0.14094854611f;
__constant__ float r2logprob = -2.14038225046f;
__constant__ float r3logprob = -4.27902812221f;
__constant__ float r4logprob = -0.2387141463f;
__constant__ float r5logprob = -1.549472331f;
__constant__ float r6logprob = -0.11137523453f;
__constant__ float r7logprob = -2.25002110628f;

/*
__global__ void initialiseinside(float *inside, const float* unpairedlogprobs, int len)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i < len)
	{
		for(int j=0 ; j < len ; j++)
		{
			if(i == j)
			{
				inside[S*len*len + i*len + j] = unpairedlogprobs[j] + r2logprob;
				inside[L*len*len + i*len + j] = unpairedlogprobs[j] + r6logprob;
				inside[F*len*len + i*len + j] = -1e20f;
			}
			else
			{
				inside[S*len*len + i*len + j] = -1e20f;
				inside[L*len*len + i*len + j] = -1e20f;
				inside[F*len*len + i*len + j] = -1e20f;
			}
		}
	}
}

__global__ void insidealgorithm(float* inside, const float* pairedlogprobs, const float* unpairedlogprobs, const int b, const int len, const float BT)
{	
	int j = blockIdx.x*blockDim.x + threadIdx.x;
	if (j < len-b)
	{		
		int index = 0;
		// type 3 rules

		// rule 1
		float tmp = -1e20f;
		for(int h=j ; h < j+b ; h++)
		{			
			tmp = logsumexp(tmp, inside[L*len*len + j*len + h] + inside[S*len*len + (h+1)*len + (j+b)]);
		}
		index = S*len*len + j*len + j+b;
		inside[index] = logsumexp(inside[index], r1logprob + tmp);

		// rule 5
		index = F*len*len + j*len + j+b;
		inside[index] = logsumexp(inside[index], r5logprob + tmp);

		// type 2 rules

		float v = pairedlogprobs[j*len+j+b]*BT + inside[F*len*len+(j+1)*len+ (j+b-1)];

		// rule 3
		index = S*len*len + j*len + j+b;
		inside[index] = logsumexp(inside[index], r3logprob + v);

		// rule 4
		index = F*len*len + j*len + j+b;
		inside[index] = logsumexp(inside[index], r4logprob + v);

		// rule 7
		index = L*len*len + j*len + j+b;
		inside[index] = logsumexp(inside[index], r7logprob + v);  
 	}
	
}

__global__ void insidez(const float* inside, float* Z, const int len)
{
	int j = blockIdx.x*blockDim.x + threadIdx.x;
	if(j == 0)
	{
		Z[j] = inside[len-1];
	}
}*/


/*
__global__ void initialiseinside(float *inside, const float* unpairedlogprobs, int len)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i < len)
	{
		for(int j=0 ; j < len ; j++)
		{
			if(i == j)
			{
				inside[i*len*3 + j*3 + S] = unpairedlogprobs[j] + r2logprob;
				inside[i*len*3 + j*3 + L] = unpairedlogprobs[j] + r6logprob;
				inside[i*len*3 + j*3 + F] = -1e20f;
			}
			else
			{
				inside[i*len*3 + j*3 + S] = -1e20f;
				inside[i*len*3 + j*3 + L] = -1e20f;
				inside[i*len*3 + j*3 + F] = -1e20f;
			}
		}
	}
}

__global__ void insidealgorithm(float* inside, const float* pairedlogprobs, const float* unpairedlogprobs, const int b, const int len, const float BT)
{	
	int j = blockIdx.x*blockDim.x + threadIdx.x;
	if (j < len-b)
	{		
		int index = 0;
		// type 3 rules

		// rule 1
		float tmp = -1e20f;
		for(int h=j ; h < j+b ; h++)
		{			
			tmp = logsumexp(tmp, inside[j*len*3 + h*3 + L] + inside[(h+1)*len*3 + (j+b)*3 + S]);
		}
		index = j*len*3 + (j+b)*3 + S;
		inside[index] = logsumexp(inside[index], r1logprob + tmp);

		// rule 5
		index = j*len*3 + (j+b)*3 + F;
		inside[index] = logsumexp(inside[index], r5logprob + tmp);

		// type 2 rules

		float v = pairedlogprobs[j*len+j+b] + inside[(j+1)*len*3 + (j+b-1)*3 + F];

		// rule 3
		index = j*len*3 + (j+b)*3 + S;
		inside[index] = logsumexp(inside[index], r3logprob + v);

		// rule 4
		index = j*len*3 + (j+b)*3 + F;
		inside[index] = logsumexp(inside[index], r4logprob + v);

		// rule 7
		index = j*len*3 + (j+b)*3 + L;
		inside[index] = logsumexp(inside[index], r7logprob + v);  
 	}
	
}

__global__ void insidez(const float* inside, float* Z, const int len)
{
	int j = blockIdx.x*blockDim.x + threadIdx.x;
	if(j == 0)
	{
		Z[j] = inside[(len-1)*3];
	}
}*/


__global__ void initialiseinside(float* insideS, float* insideL, float* insideF, const float* unpairedlogprobs, int len)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i < len)
	{
		for(int j=0 ; j < len ; j++)
		{
			if(i == j)
			{
				insideS[i*len + j] = unpairedlogprobs[j] + r2logprob;
				insideL[i*len + j] = unpairedlogprobs[j] + r6logprob;
				insideF[i*len + j] = -1e20f;
			}
			else
			{
				insideS[i*len + j] = -1e20f;
				insideL[i*len + j] = -1e20f;
				insideF[i*len + j] = -1e20f;
			}
		}
	}
}

__global__ void insidealgorithm(float* insideS, float* insideL, float* insideF, const float* pairedlogprobs, const float* unpairedlogprobs, const int b, const int len, const float BT)
{	
	int j = blockIdx.x*blockDim.x + threadIdx.x;
	if (j < len-b)
	{		
		int index = j*len + j+b;
		// type 3 rules

		// rule 1
		float tmp = -1e20f;
		for(int h=j ; h < j+b ; h++)
		{			
			tmp = logsumexp(tmp, insideL[j*len + h] + insideS[(h+1)*len + (j+b)]);
		}
		insideS[index] = logsumexp(insideS[index], r1logprob + tmp);

		// rule 5
		insideF[index] = logsumexp(insideF[index], r5logprob + tmp);

		// type 2 rules

		float v = pairedlogprobs[index]*BT + insideF[(j+1)*len+ (j+b-1)];

		// rule 3
		insideS[index] = logsumexp(insideS[index], r3logprob + v);

		// rule 4
		insideF[index] = logsumexp(insideF[index], r4logprob + v);

		// rule 7
		insideL[index] = logsumexp(insideL[index], r7logprob + v);  
 	}
	
}

__global__ void posteriordecoding(float* ematrix, int* smatrix,  const float* pairprobs, const float* singleprobs, const int datalen, const int diag, const float alpha)
{	
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < datalen-diag)
	{		
		int j = i + diag;
		float e1 = singleprobs[i] + ematrix[(i+1)*datalen + j];
		float e2 = alpha*pairprobs[i*datalen+j] + ematrix[(i+1)*datalen + j-1];
		float maxe3 = -1e10;
		int maxk = 0;
		for(int k=i+1 ; k <= j-1 ; k++)
		{
			float v = alpha*pairprobs[i*datalen + k] + ematrix[(i+1)*datalen + k-1] + ematrix[(k+1)*datalen + j];
			if(v > maxe3)
			{
				maxe3 = v;
				maxk = k;
			}
		}
		
		float maxval = e1;
		smatrix[i*datalen + j] = -1;
		if(e2 > maxval)
		{
			maxval = e2;
			smatrix[i*datalen + j] = -2;
		}
		if(maxe3 > maxval)
		{
			maxval = maxe3;
			smatrix[i*datalen + j]  = maxk+1;
		}		
		ematrix[i*datalen + j] = maxval;
 	}
	
}

__global__ void insidez(const float* insideS, float* Z, const int len)
{
	int j = blockIdx.x*blockDim.x + threadIdx.x;
	if(j == 0)
	{
		Z[j] = insideS[len-1];
	}
}

/*
__global__ void initialiseinside(float* insideS, float* insideL, float* insideF, const float* unpairedlogprobs, int len, const int stride)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i < len)
	{
		for(int j=0 ; j < len ; j++)
		{
			if(i == j)
			{
				insideS[i*stride + j] = unpairedlogprobs[j] + r2logprob;
				insideL[i*stride + j] = unpairedlogprobs[j] + r6logprob;
				insideF[i*stride + j] = -1e20f;
			}
			else
			{
				insideS[i*stride + j] = -1e20f;
				insideL[i*stride + j] = -1e20f;
				insideF[i*stride + j] = -1e20f;
			}
		}
	}
}

__global__ void insidealgorithm(float* insideS, float* insideL, float* insideF, const float* pairedlogprobs, const float* unpairedlogprobs, const int b, const int len, const int stride, const float BT)
{	
	int j = blockIdx.x*blockDim.x + threadIdx.x;
	if (j < len-b)
	{		
		int index = j*stride + j+b;
		// type 3 rules

		// rule 1
		float tmp = -1e20f;
		for(int h=j ; h < j+b ; h++)
		{			
			tmp = logsumexp(tmp, insideL[j*stride + h] + insideS[(h+1)*stride + (j+b)]);
		}
		insideS[index] = logsumexp(insideS[index], r1logprob + tmp);

		// rule 5
		insideF[index] = logsumexp(insideF[index], r5logprob + tmp);

		// type 2 rules

		float v = pairedlogprobs[j*len + j+b]*BT + insideF[(j+1)*stride + (j+b-1)];

		// rule 3
		insideS[index] = logsumexp(insideS[index], r3logprob + v);

		// rule 4
		insideF[index] = logsumexp(insideF[index], r4logprob + v);

		// rule 7
		insideL[index] = logsumexp(insideL[index], r7logprob + v);  
 	}
	
}

__global__ void insidez(const float* insideS, float* Z, const int len)
{
	int j = blockIdx.x*blockDim.x + threadIdx.x;
	if(j == 0)
	{
		Z[j] = insideS[len-1];
	}
}*/


/*
push!(rules, Rule('S', "LS",0.868534, 3))
    push!(rules, Rule('S', "s",0.117609877998, 1))
    push!(rules, Rule('S', "dFd",0.013856122002, 2))
    push!(rules, Rule('F', "dFd",0.787640, 2))
    push!(rules, Rule('F', "LS",0.21236, 3))
    push!(rules, Rule('L', "s",0.894603, 1))
    push!(rules, Rule('L', "dFd",0.105397, 2))
    type1rules = Rule[rules[2],rules[6]]
    type2rules = Rule[rules[3],rules[4],rules[7]]
    type3rules = Rule[rules[1],rules[5]]
    ruleindex = Dict('S' => 1, 'F' => 2, 'L' => 3)
*/
__global__ void outsidealgorithm(float* outside, const float* inside, const float* pairedlogprobs, const float* unpairedlogprobs, const int b, const int len, const float BT)
{	
	int j = blockIdx.x*blockDim.x + threadIdx.x;
	if (j < len - b)
	{
		int index = 0;
		// type 3 rules

		// rule 1 Rule('S', "LS",0.868534, 3))
		float tmp = -1e20f;
		for (int k = j + b + 1; k < len; k++)
		{
			tmp = logsumexp(tmp, outside[S*len*len + j*len + k] + inside[S*len*len + (j+b+1)*len + k]);
		}
		index = L*len*len + j*len + j+b;
		outside[index] = logsumexp(outside[index], r1logprob*BT + tmp);

		tmp = -1e20f;
		for (int k = 0 ; k < j ; k++)
		{
			tmp = logsumexp(tmp, outside[S*len*len + k*len + j+b] + inside[L*len*len + k*len + j-1]);
		}
		index = S*len*len + j*len + j+b;
		outside[index] = logsumexp(outside[index], r1logprob*BT + tmp);

		// rule 5 Rule('F', "LS",0.21236, 3)
		tmp = -1e20f;
		for (int k = j + b + 1; k < len; k++)
		{
			tmp = logsumexp(tmp, outside[F*len*len + j*len + k] + inside[S*len*len + (j+b+1)*len + k]);
		}
		index = L*len*len + j*len + j+b;
		outside[index] = logsumexp(outside[index], r5logprob*BT + tmp);

		tmp = -1e20f;
		for (int k = 0 ; k < j ; k++)
		{
			tmp = logsumexp(tmp, outside[F*len*len + k*len + j+b] + inside[L*len*len + k*len + j-1]);
		}
		index = S*len*len + j*len + j+b;
		outside[index] = logsumexp(outside[index], r5logprob*BT + tmp);
			
		// type 2 rules
		if ((j>=1) && (j+b+1<len))
		{
			float v = pairedlogprobs[(j-1)*len+(j+b+1)]*BT;
			index = F*len*len + j*len + j+b;
			// rule 3 Rule('S', "dFd",0.013856122002, 2)
			outside[index] = logsumexp(outside[index], r3logprob*BT + outside[S*len*len + (j-1)*len + j+b+1] + v);

			// rule 4 Rule('F', "dFd",0.787640, 2)
			outside[index] = logsumexp(outside[index], r4logprob*BT + outside[F*len*len + (j-1)*len + j+b+1] + v);

			// rule 7 Rule('L', "dFd",0.105397, 2)
			outside[index] = logsumexp(outside[index], r7logprob*BT + outside[L*len*len + (j-1)*len + j+b+1] + v);
		}
 	}
	
}
/*
int main()
{
  int len = 8000;
  int N = len;
  float* inside = (float*)malloc(3*len*len*sizeof(float));
  float* pairedlogprobs = (float*)malloc(len*len*sizeof(float));
  float* unpairedlogprobs = (float*)malloc(len*sizeof(float));

  float* d_inside;
  float* d_pairedlogprobs;
  float* d_unpairedlogprobs;
  cudaMalloc(&d_inside, 3*len*len*sizeof(float)); 
  cudaMalloc(&d_pairedlogprobs, len*len*sizeof(float)); 
  cudaMalloc(&d_unpairedlogprobs, len*len*sizeof(float));


  cudaMemcpy(d_inside, inside, 3*len*len*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_pairedlogprobs, pairedlogprobs, len*len*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_unpairedlogprobs, unpairedlogprobs, len*sizeof(float), cudaMemcpyHostToDevice);

  for(int b=1 ; b < len ; b++)
  {	
  	insidealgorithm<<<(N+511)/512, 512>>>(d_inside, d_pairedlogprobs, d_unpairedlogprobs, b, len, 1.0);
  }

  int code = cudaMemcpy(inside, d_inside, 3*len*len*sizeof(float), cudaMemcpyDeviceToHost);
  printf("exitcode %d\n", code);
  printf("Z %f\n", inside[len-1]);
}*/

