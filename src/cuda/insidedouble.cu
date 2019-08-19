#include <stdio.h>
#include <stdlib.h>
#include <time.h>

__device__ double logsumexp(double a, double b)
{

	if(a <= -1e20)
	{
		return b;
	}
	else
	if(b <= -1e20)
	{
		return a;
	}

	/*double diff = a-b;
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
		return a + log(1.0+exp(b-a));
	}
	else
	{
		return b + log(1.0+exp(a-b));
	}
}

__device__ double safeadd(double a, double b)
{
	if(a <= -1e20)
	{
		return b;
	}
	else
	if(b <= -1e20)
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
	
__constant__ double r1logprob = -0.14094854611;
__constant__ double r2logprob = -2.14038225046;
__constant__ double r3logprob = -4.27902812221;
__constant__ double r4logprob = -0.2387141463;
__constant__ double r5logprob = -1.549472331;
__constant__ double r6logprob = -0.11137523453;
__constant__ double r7logprob = -2.25002110628;

__global__ void insidealgorithm(double* inside, const double* pairedlogprobs, const double* unpairedlogprobs, const int b, const int len, const double BT)
{	
	int j = blockIdx.x*blockDim.x + threadIdx.x;
	if (j < len-b)
	{		
		int index = 0;
		// type 3 rules

		// rule 1
		double tmp = -1e20;
		for(int h=j ; h < j+b ; h++)
		{
			tmp = logsumexp(tmp, inside[L*len*len + j*len + h] + inside[S*len*len + (h+1)*len + (j+b)]);
		}
		index = S*len*len + j*len + j+b;
		inside[index] = logsumexp(inside[index], r1logprob*BT + tmp);

		// rule 5
		/*tmp = -1e20f;
		for(int h=j ; h < j+b ; h++)
		{			
			double prob1 = inside[L*len*len + j*len + h];
			double prob2 = inside[S*len*len + (h+1)*len + j+b];
			tmp = logsumexp(tmp, prob1 + prob2);
		}*/
		index = F*len*len + j*len + j+b;
		inside[index] = logsumexp(inside[index], r5logprob*BT + tmp);

		// type 2 rules

		double v = pairedlogprobs[j*len+j+b]*BT + inside[F*len*len+(j+1)*len+ (j+b-1)];

		// rule 3
		index = S*len*len + j*len + j+b;
		inside[index] = logsumexp(inside[index], r3logprob*BT + v);

		// rule 4
		index = F*len*len + j*len + j+b;
		inside[index] = logsumexp(inside[index], r4logprob*BT + v);

		// rule 7
		index = L*len*len + j*len + j+b;
		inside[index] = logsumexp(inside[index], r7logprob*BT + v);  
 	}
	
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
__global__ void outsidealgorithm(double* outside, const double* inside, const double* pairedlogprobs, const double* unpairedlogprobs, const int b, const int len, const double BT)
{	
	int j = blockIdx.x*blockDim.x + threadIdx.x;
	if (j < len - b)
	{
		int index = 0;
		// type 3 rules

		// rule 1 Rule('S', "LS",0.868534, 3))
		double tmp = -1e20;
		for (int k = j + b + 1; k < len; k++)
		{
			tmp = logsumexp(tmp, outside[S*len*len + j*len + k] + inside[S*len*len + (j+b+1)*len + k]);
		}
		index = L*len*len + j*len + j+b;
		outside[index] = logsumexp(outside[index], r1logprob*BT + tmp);

		tmp = -1e20;
		for (int k = 0 ; k < j ; k++)
		{
			tmp = logsumexp(tmp, outside[S*len*len + k*len + j+b] + inside[L*len*len + k*len + j-1]);
		}
		index = S*len*len + j*len + j+b;
		outside[index] = logsumexp(outside[index], r1logprob*BT + tmp);

		// rule 5 Rule('F', "LS",0.21236, 3)
		tmp = -1e20;
		for (int k = j + b + 1; k < len; k++)
		{
			tmp = logsumexp(tmp, outside[F*len*len + j*len + k] + inside[S*len*len + (j+b+1)*len + k]);
		}
		index = L*len*len + j*len + j+b;
		outside[index] = logsumexp(outside[index], r5logprob*BT + tmp);

		tmp = -1e20;
		for (int k = 0 ; k < j ; k++)
		{
			tmp = logsumexp(tmp, outside[F*len*len + k*len + j+b] + inside[L*len*len + k*len + j-1]);
		}
		index = S*len*len + j*len + j+b;
		outside[index] = logsumexp(outside[index], r5logprob*BT + tmp);
			
		// type 2 rules
		if ((j>=1) && (j+b+1<len))
		{
			double v = pairedlogprobs[(j-1)*len+(j+b+1)]*BT;
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
