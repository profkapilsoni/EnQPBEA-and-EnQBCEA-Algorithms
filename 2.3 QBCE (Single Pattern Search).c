/*
*********************************************************************************************************
*	About: This code implements existing QBCE Algorithm as restricted singleton pattern set |P|=1 as P	*
*	Usage: Run in command prompt "QBCE.exe <input_file>"												*
*	1. 	<input_file> contains character subset of biological gene text sequence of SARS-CoV-2 for Human	*
* 		A file size {128, 256, 512} characters are intentioally prepared for feasible QuEST Simulation 	*
*	2.	The searching of existing QBCE Algorithm is performed for the fixed pattern P of the length 3	*
*		Quantum algorithmic circuit is designed for a searching pattern P = "A T G" within the text		*
*********************************************************************************************************
*/

// Instead of applying unitary circuit we use their corresponding unitary matrices
// Some realization of quantum circuits are done classically, instead of quantum
// NOTE: This is done to make the QuEST simulation as memory efficient and effective

#include <stdio.h>
#include <conio.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "QuEST.h"

//The classical function gen_Psym is not optimized and runs in brute-force 
int gen_Psym(char P[], int size, char **Psym)
{
	*Psym = (char *)malloc(sizeof(char) * size);
	int count = 0;
	for(int i = 0; i < size; i++)
	{
		int found = 0;
		for(int j = 0; j < count; j++)
		{
			if((*Psym)[j] == P[i]) found = 1;
		}
		if(!found) (*Psym)[count++] = P[i];
	}
	return count;
}

//Used in building Algebraic Normal Form (ANF)
int check_term(int t, int a, int width)
{
	for (int i = 0; i < width; i++)
	{
		if ((t & (1 << i)) == 0 && (a & (1 << i)) != 0) return 0;
	}
	return 1;
}

//Build or Construct Algebraic Normal Form (ANF): Realized for QMEM, Quantum Adder Operation of QAF Filtering, and QEM operation
void anf_calc(int f[], int n, int ***anf, int **anf_size)
{
	*anf = (int **)malloc(sizeof(int *) * n);
	for(int i = 0; i < n; i++)
		(*anf)[i] = (int *)malloc(sizeof(int) * (int)pow(2, n));
	*anf_size = (int *)malloc(sizeof(int) * n);
	for(int i = 0; i < n; i++)
		(*anf_size)[i] = 0;

	for (int i = n - 1; i >= 0; i--)
	{
		for (int term = 0; term < (int)pow(2, n); term++)
		{
			int count = 0;
			for (int fi = 0; fi < (int)pow(2, n); fi++)
			{
				if (check_term(term, fi, n) && (f[fi] & (1 << i)) != 0) count++;
			}
			if (count % 2 == 1) (*anf)[i][(*anf_size)[i]++] = term;
		}
	}
}

//Return 2's complement of an integer of "n" bits
int _2s_complement(int num, int bits)
{
	for(int i = 0; i < bits; i++)
		num ^= (1 << i);
	return (num + 1);
}

//Binary Adder
int bin_add(int a, int b, int width)
{
	int carry = 0, sum = 0;
	for(int i = 0; i < width; i++)
	{
		int _a = ((a & (1 << i)) >> i);
		int _b = ((b & (1 << i)) >> i);
		sum ^= ((_a ^ _b ^ carry) << i);
		if((_a + _b + carry) >= 2)
			carry = 1;
		else carry = 0;
	}
	return sum;
}

int is_at(char *T, int T_size, int idx, char *P, int P_size)
{
	int found = 1;
	int temp_i = 0;
	if((idx + P_size) > T_size) return 0;
	for(int i = idx; i < idx + P_size; i++)
		if(T[i] != P[temp_i++]) found = 0;
	return found;
}

int main (int narg, char *varg[])
{
	if(narg != 2)
	{
		printf("Usage: %s <input_file>", varg[0]);
		return 0;
	}
	
	//STARTS - Reading and Storing the Input File in the text array "T"
	//The code assumes that the number of elements in the input file is in power of 2
	//for the purpose of easier processing and binary equivalent matching

	FILE *fp = fopen(varg[1], "r");
	int count = 0;
	char ch;
	while((ch = fgetc(fp)) != EOF)
	{
		if(tolower(ch) == 'a' || tolower(ch) == 't' || tolower(ch) == 'g' || tolower(ch) == 'c')
			count++;
	}
	fclose(fp);
	
	char *T = (char *)malloc(sizeof(char) * count);
	int idx = 0;
	fp = fopen(varg[1], "r");
	while((ch = fgetc(fp)) != EOF)
	{
		if(tolower(ch) == 'a' || tolower(ch) == 't' || tolower(ch) == 'g' || tolower(ch) == 'c')
			T[idx++] = ch;
	}
	fclose(fp);

	//ENDS - Reading and Storing Input File in the array "T"

	//STARTS - Finds "Psym" and "i_j" corresponding to the pattern "P"
	//An arbitrary pattern P = "ATG" used for searching within the text "T"

	const int P_size = 3;
	char P[] = { 'a', 't', 'g' };
	char *Psym;
	int Psym_size = gen_Psym(P, P_size, &Psym); //Generates Psym corresponding to P

	int *i_j = (int *)malloc(sizeof(int) * Psym_size);
	for(int i = 0; i < Psym_size; i++)
	{
		for(int j = 0; j < P_size; j++)
		{
			if(Psym[i] == P[j])
			{
				i_j[i] = j;
				break;
			}
		}
	}

	printf("\nDNA sequence is:\n");
	for(int i = 1; i <= count; i++)
	{
		printf("%c", T[i - 1]);
		if(i % 32 == 0) printf("\n");
	}
	printf("Pattern to search: ");
	for(int i = 0; i < P_size; i++)
		printf("%c", P[i]);
	printf("\nPsym is: ");
	for(int i = 0; i < Psym_size; i++)
		printf("%c", Psym[i]);
	printf("\ni_j is: ");
	for(int i = 0; i < Psym_size; i++)
		printf("%d ", i_j[i]);
	printf("\n");

	//ENDS - Finds "Psym" and "i_j" corresponding to the pattern "P"

	//STARTS - Preliminaries for quantum circuit corresponding to U_Loc and U_Sub
	//NOTE: This is not the efficient approach, however, can be improved for unrestricted simulation
	
	int n = (int)(log(count)/log(2));

	int **anf, *anf_size;
    int *f = (int *)malloc(sizeof(int) * (int)pow(2, n));
    for(int i = 0; i < count; i++) //count = 2^n
    {
    	idx = 0;
    	for(int j = 0; j < Psym_size; j++)
    	{	
    		if(tolower(T[i]) == tolower(Psym[j])) 
    		{
    			idx = i_j[j];
    			break;
    		}
    	}
    	f[i] = bin_add(i, _2s_complement(idx, n), n);
    }
    anf_calc(f, n, &anf, &anf_size);

	//ENDS - Preliminaries for quantum circuit corresponding to U_Loc

	//STARTS - Quantum Approximate Filtering (QAF) Call through Quantum Adder Operation
	//Quantum Environment (env) Realizing Single Quantum System (qubits) = Simulates One Quantum Core over Multithreaded Classical Core

    QuESTEnv env = createQuESTEnv();

    Qureg qubits = createQureg(2 * n, env);
    initZeroState(qubits);
    printf("\nQuantum Parameters for Filtering Part:\n");
    reportQuregParams(qubits);
    reportQuESTEnv(env);

    //Single qubit unitary for XOR gate
    ComplexMatrix2 ux = {
    	.real={{0,1},{1,0}},
		.imag={{0,0},{0,0}}
    };

    for(int i = n; i < 2 * n; i++)
    	hadamard(qubits, i);

    //Evaluation of U_Loc and U_Sub by below quantum circuit
    for(int idx = 0; idx < n; idx++)
    {
    	for(int i = 0; i < anf_size[idx]; i++)
    	{
    		if(anf[idx][i] == 0)
    		{
    			pauliX(qubits, idx);
    			continue;
    		}
    		int *ctrls = (int *)malloc(sizeof(int) * n);
    		int ctrl_size = 0;
    		int term = anf[idx][i];
    		int qb = 0;
    		while(term)
    		{
    			if(term & 1) ctrls[ctrl_size++] = (qb + n);
    			term >>= 1;
    			qb++;
    		}
    		multiControlledUnitary(qubits, ctrls, ctrl_size, idx, ux);
    		free(ctrls);
    	}
    }
	
    for(int i = n; i < 2 * n; i++)
    	hadamard(qubits, i);

    //Prints the quantum state before applying Grover's algorithm
    /*
	printf("Required Quantum State is constructed.\n");
	printf("Do you want to view the constructed quantum state?(y/n):");
	ch = getch();
	if(ch == 'y')
	{
		printf("\nConstructed quantum state before applying Grover's search is:");
	    qreal prob;
	    int mask = (int)pow(2, n) - 1;
	    for(int i = 0; i < (int)pow(2, 2 * n); i++)
	    {
	    	prob = getProbAmp(qubits, i);
	    	if(prob != 0.0)
	    		printf("\ni = %d, i-ij = %d has prob %f", (i >> n), (i & mask), prob);
	   	}
	}
	*/

    //Reading all the highest probability indexes corresponding to the threshold
    int d = Psym_size - 1; //Threshold of the approximate search
    double req_prob = (double)d / (double)pow(2, n);
    req_prob *= req_prob;
    //printf("Threshold probability is %f\n", req_prob);
    //Counting and storing likely indexes
    int res_count = 0;
    int realloc_size = 1;
    int *Loc = (int *)malloc(sizeof(int) * realloc_size);
    qreal prob;
    for(int i = 0; i < (int)pow(2, n); i++)
    {
    	prob = getProbAmp(qubits, i);
    	if((int)(prob * 1000000) >= (int)(req_prob * 1000000))
    	{
    		if(realloc_size == res_count)
    		{
    			realloc_size *= 2;
    			Loc = (int *)realloc(Loc, sizeof(int) * realloc_size);
    		}
    		Loc[res_count++] = i;
    	}
    }

	printf("\nLikely indexes are:\n");
	for(int i = 0; i < res_count; i++)
		printf("Index %d\n", Loc[i]);

	destroyQureg(qubits, env);
    destroyQuESTEnv(env);

	//ENDS - Quantum Approximate Filtering (QAF) Call

    if(res_count < 2)
    {
    	printf("Grover's Searching is not required since only %d index is found!\n", res_count);
    	if(res_count == 1 && is_at(T, count, Loc[0], P, P_size)) printf("Classical checked that the pattern is at index %d.\n", Loc[0]);
    	else printf("The pattern is not at index %d.\n", Loc[0]);
    	return 0;
    }

	//STARTS - Existing QBCE ALgorithm for Single Pattern Matching through Quantum Exact Match (QEM) & Grover's Search Operator (GSO) Call
	//Quantum Environment (env) Realizing Single Quantum System (qubits) = Simulates One Quantum Core over Multithreaded Classical Core

    env = createQuESTEnv();

    int req_qbs = ceil(log(res_count)/log(2));
    qubits = createQureg(req_qbs, env);
    initZeroState(qubits);
    printf("\nQuantum Parameters for Grover's Part:\n");
    reportQuregParams(qubits);
    reportQuESTEnv(env);

	int T_size = count;
	count = 0;
	ComplexMatrixN e = createComplexMatrixN(req_qbs);
    for(int i = 0; i < (int)pow(2, req_qbs); i++)
    {
    	if(i < res_count && is_at(T, T_size, Loc[i], P, P_size))
    	{
    		e.real[i][i] = -1;
    		count++;
    	}
    	else e.real[i][i] = 1;
    }

	int *targs = (int *)malloc(sizeof(int) * n);
	for(int i = 0; i < req_qbs; i++)
		targs[i] = i;
    int times =  ceil(3.14 * (pow(2, req_qbs / 2) / sqrt(count)) / 4);
    printf("\nRunning Grover's Algorithm..\n");
    for(int i = 0; i < req_qbs; i++)
        hadamard(qubits, i);
    for(int gi = 0; gi < times; gi++)
    {
    	//Marking
    	multiQubitUnitary(qubits, targs, req_qbs, e);
    	
        //Diffusion
        for(int i = 0; i < req_qbs; i++)
            hadamard(qubits, i);
        for(int i = 0; i < req_qbs; i++)
            pauliX(qubits, i);
        multiControlledPhaseFlip(qubits, targs, req_qbs);
        for(int i = 0; i < req_qbs; i++)
            pauliX(qubits, i);
        for(int i = 0; i < req_qbs; i++)
            hadamard(qubits, i);
    }

    qreal max = 0.0;
    prob = 0.0;
    for(int i = 0; i < (int)pow(2, req_qbs); i++)
    {
    	prob = getProbAmp(qubits, i);
    	if(max <= prob) max = prob;
    }
    for(int i = 0; i < (int)pow(2, req_qbs); i++)
    {
    	prob = getProbAmp(qubits, i);
    	if(fabs(max - prob) <= 0.0000001) printf("Correct Index %d\n", Loc[i]);
    }

    destroyQureg(qubits, env);
	destroyComplexMatrixN(e);
    destroyQuESTEnv(env);

	//ENDS - Existing QBCE ALgorithm for Single Pattern Matching

    return 0;
}