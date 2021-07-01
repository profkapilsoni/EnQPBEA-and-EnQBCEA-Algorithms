/*
*************************************************************************************************************
*	About: The code implements Enhanced quantum based combined exact algorithm (EnQBCEA-MPM) for Equal Size	*
*	Usage: Run in command prompt "EnQBCEA.exe <input_file>"													*
*	1. 	<input_file> contains character subset of biological gene text sequence of SARS-CoV-2 for Human.	*
* 		A file sizes {128, 256, 512} characters are intentioally prepared for feasible QuEST Simulation 	*
*	2.	EnQBCEA-MPM matching is performed for patterns set P = {P1, P2, P3} of equal lengths {2, 3, 4}		*
		Multiple Quantum Core (QCore = 3) realization based EnQBCEA-MPM Algorithm simulation.				*
*		Quantum circuits are designed to search P1 = "T A A", P2 = "T A G", P3 = "T G A" within text T		*
*************************************************************************************************************
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

	printf("\nDNA sequence is:\n");
	for(int i = 1; i <= count; i++)
	{
		printf("%c", T[i - 1]);
		if(i % 32 == 0) printf("\n");
	}

	//ENDS - Reading and Storing Input File in the array "T"

	//STARTS - Finds "Psym" and "i_j" corresponding to the pattern "P1"
	//An arbitrary pattern P1 = "TAA" used for searching within the text "T"
	
	const int P_size1 = 3;
	char P1[] =  { 't', 'a', 'a' };
	char *Psym1;
	int Psym_size1 = gen_Psym(P1, P_size1, &Psym1); //Generates Psym corresponding to P

	int *i_j1 = (int *)malloc(sizeof(int) * Psym_size1);
	for(int i = 0; i < Psym_size1; i++)
	{
		for(int j = 0; j < P_size1; j++)
		{
			if(Psym1[i] == P1[j])
			{
				i_j1[i] = j;
				break;
			}
		}
	}
	printf("\n\tCORE 1:\n");
	printf("------------------------------------------------\n");
	printf("Pattern to search: ");
	for(int i = 0; i < P_size1; i++)
		printf("%c", P1[i]);
	printf("\nPsym is: ");
	for(int i = 0; i < Psym_size1; i++)
		printf("%c", Psym1[i]);
	printf("\ni_j is: ");
	for(int i = 0; i < Psym_size1; i++)
		printf("%d ", i_j1[i]);
	printf("\n------------------------------------------------\n\n");

	//ENDS - Finds "Psym" and "i_j" corresponding to the pattern "P1"


	//STARTS - Finds "Psym" and "i_j" corresponding to the pattern "P2"
	//An arbitrary pattern P2 = "TAG" used for searching within the text "T"

	const int P_size2 = 3;
	char P2[] =  { 't', 'a', 'g' };
	char *Psym2;
	int Psym_size2 = gen_Psym(P2, P_size2, &Psym2); //Generates Psym corresponding to P

	int *i_j2 = (int *)malloc(sizeof(int) * Psym_size2);
	for(int i = 0; i < Psym_size2; i++)
	{
		for(int j = 0; j < P_size2; j++)
		{
			if(Psym2[i] == P2[j])
			{
				i_j2[i] = j;
				break;
			}
		}
	}
	printf("\n\tCORE 2:\n");
	printf("------------------------------------------------\n");
	printf("Pattern to search: ");
	for(int i = 0; i < P_size2; i++)
		printf("%c", P2[i]);
	printf("\nPsym is: ");
	for(int i = 0; i < Psym_size2; i++)
		printf("%c", Psym2[i]);
	printf("\ni_j is: ");
	for(int i = 0; i < Psym_size2; i++)
		printf("%d ", i_j2[i]);
	printf("\n------------------------------------------------\n\n");

	//ENDS - Finds "Psym" and "i_j" corresponding to the pattern "P2"


	//STARTS - Finds "Psym" and "i_j" corresponding to the pattern "P3"
	//An arbitrary pattern P3 = "TGA" used for searching within the text "T"

	const int P_size3 = 3;
	char P3[] =  { 't', 'g', 'a' };
	char *Psym3;
	int Psym_size3 = gen_Psym(P3, P_size3, &Psym3); //Generates Psym corresponding to P

	int *i_j3 = (int *)malloc(sizeof(int) * Psym_size3);
	for(int i = 0; i < Psym_size3; i++)
	{
		for(int j = 0; j < P_size3; j++)
		{
			if(Psym3[i] == P3[j])
			{
				i_j3[i] = j;
				break;
			}
		}
	}
	printf("\n\tCORE 3:\n");
	printf("------------------------------------------------\n");
	printf("Pattern to search: ");
	for(int i = 0; i < P_size3; i++)
		printf("%c", P3[i]);
	printf("\nPsym is: ");
	for(int i = 0; i < Psym_size3; i++)
		printf("%c", Psym3[i]);
	printf("\ni_j is: ");
	for(int i = 0; i < Psym_size3; i++)
		printf("%d ", i_j3[i]);
	printf("\n------------------------------------------------\n\n");

	//ENDS - Finds "Psym" and "i_j" corresponding to the pattern "P3"


	//STARTS - Preliminaries for quantum circuit that will be realized seperately on each quantum core
	//NOTE: This is not the efficient approach, however, can be improved for unrestricted simulation
	
	int n1 = (int)(log(count)/log(2));
	int **anf1, *anf_size1;
    int *f1 = (int *)malloc(sizeof(int) * (int)pow(2, n1));
    for(int i = 0; i < count; i++) //count = 2^n
    {
    	idx = 0;
    	for(int j = 0; j < Psym_size1; j++)
    	{	
    		if(tolower(T[i]) == tolower(Psym1[j])) 
    		{
    			idx = i_j1[j];
    			break;
    		}
    	}
    	f1[i] = bin_add(i, _2s_complement(idx, n1), n1);
    }
    anf_calc(f1, n1, &anf1, &anf_size1);


	int n2 = (int)(log(count)/log(2));
	int **anf2, *anf_size2;
    int *f2 = (int *)malloc(sizeof(int) * (int)pow(2, n2));
    for(int i = 0; i < count; i++) //count = 2^n
    {
    	idx = 0;
    	for(int j = 0; j < Psym_size2; j++)
    	{	
    		if(tolower(T[i]) == tolower(Psym2[j])) 
    		{
    			idx = i_j2[j];
    			break;
    		}
    	}
    	f2[i] = bin_add(i, _2s_complement(idx, n2), n2);
    }
    anf_calc(f2, n2, &anf2, &anf_size2);


	int n3 = (int)(log(count)/log(2));
	int **anf3, *anf_size3;
    int *f3 = (int *)malloc(sizeof(int) * (int)pow(2, n3));
    for(int i = 0; i < count; i++) //count = 2^n
    {
    	idx = 0;
    	for(int j = 0; j < Psym_size3; j++)
    	{	
    		if(tolower(T[i]) == tolower(Psym3[j])) 
    		{
    			idx = i_j3[j];
    			break;
    		}
    	}
    	f3[i] = bin_add(i, _2s_complement(idx, n3), n3);
    }
    anf_calc(f3, n3, &anf3, &anf_size3);

	//END - Preliminaries for quantum circuit that will be realized seperately on each quantum core


	//STARTS - Quantum Approximate Filtering (QAF) on the Multiple Quantum Core (QCore = 3) through Quantum Adder Operation
	//Quantum Environment (env) Realizing Three Quantum System (qubits1, qubits2, qubits3) = Simulates 3 Quantum Core over Multithreaded Classical Core

	ComplexMatrix2 ux = {
    	.real={{0,1},{1,0}},
		.imag={{0,0},{0,0}}
    };


    QuESTEnv env = createQuESTEnv();


	printf("\n\tCORE 1:\n");
	printf("------------------------------------------------\n");
    Qureg qubits1 = createQureg(2 * n1, env);
    initZeroState(qubits1);
    printf("\nQuantum Parameters for Filtering Part:\n");
    reportQuregParams(qubits1);
    reportQuESTEnv(env);
	printf("\n------------------------------------------------\n\n");


	printf("\n\tCORE 2:\n");
	printf("------------------------------------------------\n");
    Qureg qubits2 = createQureg(2 * n2, env);
    initZeroState(qubits2);
    printf("\nQuantum Parameters for Filtering Part:\n");
    reportQuregParams(qubits2);
    reportQuESTEnv(env);
	printf("\n------------------------------------------------\n\n");


	printf("\n\tCORE 3:\n");
	printf("------------------------------------------------\n");
    Qureg qubits3 = createQureg(2 * n3, env);
    initZeroState(qubits3);
    printf("\nQuantum Parameters for Filtering Part:\n");
    reportQuregParams(qubits3);
    reportQuESTEnv(env);
	printf("\n------------------------------------------------\n\n");


	//STARTS - Quantum Approximate Filtering (QAF) Call for Simulated QCore 1

	printf("\n\tCORE 1:\n");
	printf("------------------------------------------------\n");
    for(int i = n1; i < 2 * n1; i++)
    	hadamard(qubits1, i);
    for(int idx = 0; idx < n1; idx++)
    {
    	for(int i = 0; i < anf_size1[idx]; i++)
    	{
    		if(anf1[idx][i] == 0)
    		{
    			pauliX(qubits1, idx);
    			continue;
    		}
    		int *ctrls = (int *)malloc(sizeof(int) * n1);
    		int ctrl_size = 0;
    		int term = anf1[idx][i];
    		int qb = 0;
    		while(term)
    		{
    			if(term & 1) ctrls[ctrl_size++] = (qb + n1);
    			term >>= 1;
    			qb++;
    		}
    		multiControlledUnitary(qubits1, ctrls, ctrl_size, idx, ux);
    		free(ctrls);
    	}
    }
    for(int i = n1; i < 2 * n1; i++)
    	hadamard(qubits1, i);
    int d1 = Psym_size1 - 1; //Threshold of the approximate search
    double req_prob1 = (double)d1 / (double)pow(2, n1);
    req_prob1 *= req_prob1;
    int res_count1 = 0;
    int realloc_size1 = 1;
    int *Loc1 = (int *)malloc(sizeof(int) * realloc_size1);
    qreal prob1;
    for(int i = 0; i < (int)pow(2, n1); i++)
    {
    	prob1 = getProbAmp(qubits1, i);
    	if(prob1 >= req_prob1)
    	{
    		if(realloc_size1 == res_count1)
    		{
    			realloc_size1 *= 2;
    			Loc1 = (int *)realloc(Loc1, sizeof(int) * realloc_size1);
    		}
    		Loc1[res_count1++] = i;
    	}
    }
	printf("\nLikely indexes are:\n");
	for(int i = 0; i < res_count1; i++)
		printf("Index %d\n", Loc1[i]);
	destroyQureg(qubits1, env);
	printf("\n------------------------------------------------\n\n");

	//END - Quantum Approximate Filtering (QAF) Call for Simulated QCore 1


	//STARTS - Quantum Approximate Filtering (QAF) Call for Simulated QCore 2

	printf("\n\tCORE 2:\n");
	printf("------------------------------------------------\n");
    for(int i = n2; i < 2 * n2; i++)
    	hadamard(qubits2, i);
    for(int idx = 0; idx < n2; idx++)
    {
    	for(int i = 0; i < anf_size2[idx]; i++)
    	{
    		if(anf2[idx][i] == 0)
    		{
    			pauliX(qubits2, idx);
    			continue;
    		}
    		int *ctrls = (int *)malloc(sizeof(int) * n2);
    		int ctrl_size = 0;
    		int term = anf2[idx][i];
    		int qb = 0;
    		while(term)
    		{
    			if(term & 1) ctrls[ctrl_size++] = (qb + n2);
    			term >>= 1;
    			qb++;
    		}
    		multiControlledUnitary(qubits2, ctrls, ctrl_size, idx, ux);
    		free(ctrls);
    	}
    }
    for(int i = n2; i < 2 * n2; i++)
    	hadamard(qubits2, i);
    int d2 = Psym_size2 - 1; //Threshold of the approximate search
    double req_prob2 = (double)d2 / (double)pow(2, n2);
    req_prob2 *= req_prob2;
    int res_count2 = 0;
    int realloc_size2 = 1;
    int *Loc2 = (int *)malloc(sizeof(int) * realloc_size2);
    qreal prob2;
    for(int i = 0; i < (int)pow(2, n2); i++)
    {
    	prob2 = getProbAmp(qubits2, i);
    	if(prob2 >= req_prob2)
    	{
    		if(realloc_size2 == res_count2)
    		{
    			realloc_size2 *= 2;
    			Loc2 = (int *)realloc(Loc2, sizeof(int) * realloc_size2);
    		}
    		Loc2[res_count2++] = i;
    	}
    }
	printf("\nLikely indexes are:\n");
	for(int i = 0; i < res_count2; i++)
		printf("Index %d\n", Loc2[i]);
	destroyQureg(qubits2, env);
	printf("\n------------------------------------------------\n\n");

	//END - Quantum Approximate Filtering (QAF) Call for Simulated QCore 2


	//START - Quantum Approximate Filtering (QAF) Call for Simulated QCore 3

	printf("\n\tCORE 3:\n");
	printf("------------------------------------------------\n");
    for(int i = n3; i < 2 * n3; i++)
    	hadamard(qubits3, i);
    for(int idx = 0; idx < n3; idx++)
    {
    	for(int i = 0; i < anf_size3[idx]; i++)
    	{
    		if(anf3[idx][i] == 0)
    		{
    			pauliX(qubits3, idx);
    			continue;
    		}
    		int *ctrls = (int *)malloc(sizeof(int) * n3);
    		int ctrl_size = 0;
    		int term = anf3[idx][i];
    		int qb = 0;
    		while(term)
    		{
    			if(term & 1) ctrls[ctrl_size++] = (qb + n3);
    			term >>= 1;
    			qb++;
    		}
    		multiControlledUnitary(qubits3, ctrls, ctrl_size, idx, ux);
    		free(ctrls);
    	}
    }
    for(int i = n3; i < 2 * n3; i++)
    	hadamard(qubits3, i);
    int d3 = Psym_size3 - 1; //Threshold of the approximate search
    double req_prob3 = (double)d3 / (double)pow(2, n3);
    req_prob3 *= req_prob3;
    int res_count3 = 0;
    int realloc_size3 = 1;
    int *Loc3 = (int *)malloc(sizeof(int) * realloc_size3);
    qreal prob3;
    for(int i = 0; i < (int)pow(2, n3); i++)
    {
    	prob3 = getProbAmp(qubits3, i);
    	if(prob3 >= req_prob3)
    	{
    		if(realloc_size3 == res_count3)
    		{
    			realloc_size3 *= 2;
    			Loc3 = (int *)realloc(Loc3, sizeof(int) * realloc_size3);
    		}
    		Loc3[res_count3++] = i;
    	}
    }
	printf("\nLikely indexes are:\n");
	for(int i = 0; i < res_count3; i++)
		printf("Index %d\n", Loc3[i]);
	destroyQureg(qubits3, env);
	printf("\n------------------------------------------------\n\n");

	//END - Quantum Approximate Filtering (QAF) Call for Simulated QCore 3

	destroyQuESTEnv(env);

	//Destroy Quantum Environment (env) Realized Three Quantum System (qubits1, qubits2, qubits3)
	//END - Quantum Approximate Filtering (QAF) on the Multiple Quantum Core (QCore = 3)


	//STARTS - EnQBCEA-MPM Algorithm based on the Multiple Quantum Core (QCore = 3) through Quantum Exact Match (QEM) & Grover's Search Operator (GSO)
	//Quantum Environment (env) Realizing Three Quantum System (qubits1, qubits2, qubits3) = Simulates 3 Quantum Core over Multithreaded Classical Core


	//START - Checking on individual core for the single pattern occurrence to be verified at dedicated text index

	int flag1 = 1, flag2 = 1, flag3 = 1;

	
    if(res_count1 < 2)
    {
		printf("\n\tCORE 1:\n");
		printf("------------------------------------------------\n");
    	printf("Grover's Searching is not required since only %d index is found!\n", res_count1);
    	if(res_count1 == 1 && is_at(T, count, Loc1[0], P1, P_size1)) printf("Classical checked that the pattern is at index %d.\n", Loc1[0]);
    	else printf("The pattern is not at index %d.\n", Loc1[0]);
		printf("\n------------------------------------------------\n\n");
    	flag1 = 0;
    }
	


	if(res_count2 < 2)
    {
		printf("\n\tCORE 2:\n");
		printf("------------------------------------------------\n");
    	printf("Grover's Searching is not required since only %d index is found!\n", res_count2);
    	if(res_count2 == 1 && is_at(T, count, Loc2[0], P2, P_size2)) printf("Classical checked that the pattern is at index %d.\n", Loc2[0]);
    	else printf("The pattern is not at index %d.\n", Loc2[0]);
		printf("\n------------------------------------------------\n\n");
    	flag2 = 0;
    }


	
	if(res_count3 < 2)
    {
		printf("\n\tCORE 3:\n");
		printf("------------------------------------------------\n");
    	printf("Grover's Searching is not required since only %d index is found!\n", res_count3);
    	if(res_count3 == 1 && is_at(T, count, Loc3[0], P3, P_size3)) printf("Classical checked that the pattern is at index %d.\n", Loc3[0]);
    	else printf("The pattern is not at index %d.\n", Loc3[0]);
		printf("\n------------------------------------------------\n\n");
    	flag3 = 0;
    }
	
	//END - Checking on individual core for the single pattern occurrence to be verified at dedicated text index


	
	//START - Checking on individual core for the all "t" exact pattern occurrence to be verified at dedicated text index


	//STARTS - EnQBCEA Multiple Pattern Matching Algorithm Call for Simulated QCore 1

    env = createQuESTEnv();
	int T_size = count;

	if(flag1)
	{
		printf("\n\tCORE 1:\n");
		printf("------------------------------------------------\n");
		int req_qbs1 = ceil(log(res_count1)/log(2));
		qubits1 = createQureg(req_qbs1, env);
		initZeroState(qubits1);
		printf("\nQuantum Parameters for Grover's Part:\n");
		reportQuregParams(qubits1);
		reportQuESTEnv(env);
		count = 0;
		ComplexMatrixN e1 = createComplexMatrixN(req_qbs1);
		for(int i = 0; i < (int)pow(2, req_qbs1); i++)
		{
			if(i < res_count1 && is_at(T, T_size, Loc1[i], P1, P_size1))
			{
				e1.real[i][i] = -1;
				count++;
			}
			else e1.real[i][i] = 1;
		}
		int *targs1 = (int *)malloc(sizeof(int) * n1);
		for(int i = 0; i < req_qbs1; i++)
			targs1[i] = i;
		int times1 =  ceil(3.14 * (pow(2, req_qbs1 / 2) / sqrt(count)) / 4);
		printf("\nRunning Grover's Algorithm..\n");
		for(int i = 0; i < req_qbs1; i++)
			hadamard(qubits1, i);
		for(int gi = 0; gi < times1; gi++)
		{
			//Marking
			multiQubitUnitary(qubits1, targs1, req_qbs1, e1);
			
			//Diffusion
			for(int i = 0; i < req_qbs1; i++)
				hadamard(qubits1, i);
			for(int i = 0; i < req_qbs1; i++)
				pauliX(qubits1, i);
			multiControlledPhaseFlip(qubits1, targs1, req_qbs1);
			for(int i = 0; i < req_qbs1; i++)
				pauliX(qubits1, i);
			for(int i = 0; i < req_qbs1; i++)
				hadamard(qubits1, i);
		}
		qreal max1 = 0.0;
		prob1 = 0.0;
		for(int i = 0; i < (int)pow(2, req_qbs1); i++)
		{
			prob1 = getProbAmp(qubits1, i);
			if(max1 <= prob1) max1 = prob1;
		}
		for(int i = 0; i < (int)pow(2, req_qbs1); i++)
		{
			prob1 = getProbAmp(qubits1, i);
			if(fabs(max1 - prob1) <= 0.0000001) printf("Correct Index %d\n", Loc1[i]);
		}
		destroyQureg(qubits1, env);
		destroyComplexMatrixN(e1);
		printf("\n------------------------------------------------\n\n");
	}
	
	//END - EnQBCEA Multiple Pattern Matching Algorithm Call for Simulated QCore 1


	//STARTS - EnQBCEA Multiple Pattern Matching Algorithm Call for Simulated QCore 2

	if(flag2)
	{
		printf("\n\tCORE 2:\n");
		printf("------------------------------------------------\n");
		int req_qbs2 = ceil(log(res_count2)/log(2));
		qubits2 = createQureg(req_qbs2, env);
		initZeroState(qubits2);
		printf("\nQuantum Parameters for Grover's Part:\n");
		reportQuregParams(qubits2);
		reportQuESTEnv(env);
		count = 0;
		ComplexMatrixN e2 = createComplexMatrixN(req_qbs2);
		for(int i = 0; i < (int)pow(2, req_qbs2); i++)
		{
			if(i < res_count2 && is_at(T, T_size, Loc2[i], P2, P_size2))
			{
				e2.real[i][i] = -1;
				count++;
			}
			else e2.real[i][i] = 1;
		}
		int *targs2 = (int *)malloc(sizeof(int) * n2);
		for(int i = 0; i < req_qbs2; i++)
			targs2[i] = i;
		int times2 =  ceil(3.14 * (pow(2, req_qbs2 / 2) / sqrt(count)) / 4);
		printf("\nRunning Grover's Algorithm..\n");
		for(int i = 0; i < req_qbs2; i++)
			hadamard(qubits2, i);
		for(int gi = 0; gi < times2; gi++)
		{
			//Marking
			multiQubitUnitary(qubits2, targs2, req_qbs2, e2);
			
			//Diffusion
			for(int i = 0; i < req_qbs2; i++)
				hadamard(qubits2, i);
			for(int i = 0; i < req_qbs2; i++)
				pauliX(qubits2, i);
			multiControlledPhaseFlip(qubits2, targs2, req_qbs2);
			for(int i = 0; i < req_qbs2; i++)
				pauliX(qubits2, i);
			for(int i = 0; i < req_qbs2; i++)
				hadamard(qubits2, i);
		}
		qreal max2 = 0.0;
		prob2 = 0.0;
		for(int i = 0; i < (int)pow(2, req_qbs2); i++)
		{
			prob2 = getProbAmp(qubits2, i);
			if(max2 <= prob2) max2 = prob2;
		}
		for(int i = 0; i < (int)pow(2, req_qbs2); i++)
		{
			prob2 = getProbAmp(qubits2, i);
			if(fabs(max2 - prob2) <= 0.0000001) printf("Correct Index %d\n", Loc2[i]);
		}
		destroyQureg(qubits2, env);
		destroyComplexMatrixN(e2);
		printf("\n------------------------------------------------\n\n");
	}

	//END - EnQBCEA Multiple Pattern Matching Algorithm Call for Simulated QCore 2


	//STARTS - EnQBCEA Multiple Pattern Matching Algorithm Call for Simulated QCore 3

	if(flag3)
	{
		printf("\n\tCORE 3:\n");
		printf("------------------------------------------------\n");
		int req_qbs3 = ceil(log(res_count3)/log(2));
		qubits3 = createQureg(req_qbs3, env);
		initZeroState(qubits3);
		printf("\nQuantum Parameters for Grover's Part:\n");
		reportQuregParams(qubits3);
		reportQuESTEnv(env);
		count = 0;
		ComplexMatrixN e3 = createComplexMatrixN(req_qbs3);
		for(int i = 0; i < (int)pow(2, req_qbs3); i++)
		{
			if(i < res_count3 && is_at(T, T_size, Loc3[i], P3, P_size3))
			{
				e3.real[i][i] = -1;
				count++;
			}
			else e3.real[i][i] = 1;
		}
		int *targs3 = (int *)malloc(sizeof(int) * n3);
		for(int i = 0; i < req_qbs3; i++)
			targs3[i] = i;
		int times3 =  ceil(3.14 * (pow(2, req_qbs3 / 2) / sqrt(count)) / 4);
		printf("\nRunning Grover's Algorithm..\n");
		for(int i = 0; i < req_qbs3; i++)
			hadamard(qubits3, i);
		for(int gi = 0; gi < times3; gi++)
		{
			//Marking
			multiQubitUnitary(qubits3, targs3, req_qbs3, e3);
			
			//Diffusion
			for(int i = 0; i < req_qbs3; i++)
				hadamard(qubits3, i);
			for(int i = 0; i < req_qbs3; i++)
				pauliX(qubits3, i);
			multiControlledPhaseFlip(qubits3, targs3, req_qbs3);
			for(int i = 0; i < req_qbs3; i++)
				pauliX(qubits3, i);
			for(int i = 0; i < req_qbs3; i++)
				hadamard(qubits3, i);
		}
		qreal max3 = 0.0;
		prob3 = 0.0;
		for(int i = 0; i < (int)pow(2, req_qbs3); i++)
		{
			prob3 = getProbAmp(qubits3, i);
			if(max3 <= prob3) max3 = prob3;
		}
		for(int i = 0; i < (int)pow(2, req_qbs3); i++)
		{
			prob3 = getProbAmp(qubits3, i);
			if(fabs(max3 - prob3) <= 0.0000001) printf("Correct Index %d\n", Loc3[i]);
		}
		destroyQureg(qubits3, env);
		destroyComplexMatrixN(e3);
		printf("\n------------------------------------------------\n\n");
	}

	//END - EnQBCEA Multiple Pattern Matching Algorithm Call for Simulated QCore 3


	//END - Checking on individual core for the all "t" exact pattern occurrence to be verified at dedicated text index


	destroyQuESTEnv(env);

	//Destroy Quantum Environment (env) Realized Three Quantum System (qubits1, qubits2, qubits3)
	//END - EnQBCEA-MPM Algorithm based on the Multiple Quantum Core (QCore = 3)

	getch();

    return 0;
}