/*
*************************************************************************************************************
*	About: The code implements Enhanced QMEM processing based exact algorithm (EnQPBEA-MPM) for Unequal Size*
*	Usage: Run in command prompt "EnQPBEA.exe <input_file>"													*
*	1. 	<input_file> contains character subset of biological gene text sequence of SARS-CoV-2 for Human.	*
* 		A file sizes {128, 256, 512} characters are intentioally prepared for feasible QuEST Simulation 	*
*	2.	EnQPBEA-MPM matching is performed for patterns set P = {P1, P2, P3} of equal lengths {2, 3, 4}		*
		Multiple Quantum Core (QCore = 3) realization based EnQPBEA-MPM Algorithm simulation.				*
*		Quantum circuits are designed to search P1 = "T A", P2 = "T A G", P3 = "T G A C" within text T		*
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

//Used in building Algebraic Normal Form (ANF)
int check_term(int t, int a, int width)
{
	for (int i = 0; i < width; i++)
	{
		if ((t & (1 << i)) == 0 && (a & (1 << i)) != 0) return 0;
	}
	return 1;
}

//Build or Construct Algebraic Normal Form (ANF): Realized for Quantum Memory (QMEM) and Quantum Exact Match (QEM) operation
void anf_calc(int f[], int m, int n, int ***anf, int **anf_size)
{
	*anf = (int **)malloc(sizeof(int *) * n);
	for(int i = 0; i < n; i++)
		(*anf)[i] = (int *)malloc(sizeof(int) * (int)pow(2, m));
	*anf_size = (int *)malloc(sizeof(int) * n);
	for(int i = 0; i < n; i++)
		(*anf_size)[i] = 0;

	for (int i = n - 1; i >= 0; i--)
	{
		for (int term = 0; term < (int)pow(2, m); term++)
		{
			int count = 0;
			for (int fi = 0; fi < (int)pow(2, m); fi++)
			{
				if (check_term(term, fi, m) && (f[fi] & (1 << i)) != 0) count++;
			}
			if (count % 2 == 1) (*anf)[n - 1 - i][(*anf_size)[n - 1 - i]++] = term;
		}
	}
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
	int T_size = count;
	count = 0;
	fp = fopen(varg[1], "r");
	while((ch = fgetc(fp)) != EOF)
	{
		if(tolower(ch) == 'a' || tolower(ch) == 't' || tolower(ch) == 'g' || tolower(ch) == 'c')
			T[count++] = ch;
	}
	fclose(fp);
	printf("\nDNA sequence is:\n");
	for(int i = 1; i <= T_size; i++)
	{
		printf("%c", T[i - 1]);
		if(i % 32 == 0) printf("\n");
	}

	//ENDS - Reading and Storing Input File in the array "T"
	
	//The arbitrary patterns P1 = "T A", P2 = "T A G", P3 = "T G A" are used for searching within the text sequence database "T"
	//START - Printing individual search pattern on each core that belongs to patterns set P = {P1, P2, P3} of equal lengths {2, 3, 4}
	
	printf("\n\tCORE 1:\n");
	printf("------------------------------------------------\n");
	char P1[] =  { 't', 'a' };
	const int P_size1 = 2;
	printf("\nPattern to search: ");
	for(int i = 0; i < P_size1; i++)
		printf("%c", P1[i]);
	printf("\n------------------------------------------------\n\n");


	printf("\n\tCORE 2:\n");
	printf("------------------------------------------------\n");
	char P2[] =  { 't', 'a', 'g' };
	const int P_size2 = 3;
	printf("\nPattern to search: ");
	for(int i = 0; i < P_size2; i++)
		printf("%c", P2[i]);
	printf("\n------------------------------------------------\n\n");


	printf("\n\tCORE 3:\n");
	printf("------------------------------------------------\n");
	char P3[] =  { 't', 'g', 'a', 'c' };
	const int P_size3 = 4;
	printf("\nPattern to search: ");
	for(int i = 0; i < P_size3; i++)
		printf("%c", P3[i]);
	printf("\n------------------------------------------------\n\n");

	//END - Printing individual search pattern on each core that belongs to patterns set P = {P1, P2, P3} of equal lengths {2, 3, 4}

	//STARTS - Preliminaries for quantum circuit that will be realized seperately on each quantum core
	//NOTE: This is not the efficient approach, however, can be improved for unrestricted simulation

	//printf("\n\tCORE 1:\n");
	//printf("------------------------------------------------\n");
	int n1 = (int)(log(count)/log(2));
	int **anf1, *anf_size1;
    int *f1 = (int *)malloc(sizeof(int) * (int)pow(2, n1));
    for(int i = 0; i < T_size; i++) //T_size = 2^n
    {
    	if(tolower(T[i]) == 'a') f1[i] = 0;
    	else if(tolower(T[i]) == 't') f1[i] = 1;
    	else if(tolower(T[i]) == 'c') f1[i] = 2;
    	else if(tolower(T[i]) == 'g') f1[i] = 3;
    }
    anf_calc(f1, n1, 2, &anf1, &anf_size1);
	//printf("\n------------------------------------------------\n\n");



	//printf("\n\tCORE 2:\n");
	//printf("------------------------------------------------\n");
	int n2 = (int)(log(count)/log(2));
	int **anf2, *anf_size2;
    int *f2 = (int *)malloc(sizeof(int) * (int)pow(2, n2));
    for(int i = 0; i < T_size; i++) //T_size = 2^n
    {
    	if(tolower(T[i]) == 'a') f2[i] = 0;
    	else if(tolower(T[i]) == 't') f2[i] = 1;
    	else if(tolower(T[i]) == 'c') f2[i] = 2;
    	else if(tolower(T[i]) == 'g') f2[i] = 3;
    }
    anf_calc(f2, n2, 2, &anf2, &anf_size2);
	//printf("\n------------------------------------------------\n\n");



	//printf("\n\tCORE 3:\n");
	//printf("------------------------------------------------\n");
	int n3 = (int)(log(count)/log(2));
	int **anf3, *anf_size3;
    int *f3 = (int *)malloc(sizeof(int) * (int)pow(2, n3));
    for(int i = 0; i < T_size; i++) //T_size = 2^n
    {
    	if(tolower(T[i]) == 'a') f3[i] = 0;
    	else if(tolower(T[i]) == 't') f3[i] = 1;
    	else if(tolower(T[i]) == 'c') f3[i] = 2;
    	else if(tolower(T[i]) == 'g') f3[i] = 3;
    }
    anf_calc(f3, n3, 2, &anf3, &anf_size3);
	//printf("\n------------------------------------------------\n\n");

	//END - Preliminaries for quantum circuit that will be realized seperately on each quantum core

	//STARTS - EnQPBEA-MPM Algorithm based on the Multiple Quantum Core (QCore = 3) through Quantum Exact Match (QEM) & Grover's Search Operator (GSO)
	//Quantum Environment (env) Realizing Three Quantum System (qubits1, qubits2, qubits3) = Simulates 3 Quantum Core over Multithreaded Classical Core

    QuESTEnv env = createQuESTEnv();

	ComplexMatrix2 ux = {
    	.real={{0,1},{1,0}},
		.imag={{0,0},{0,0}}
    };

	//STARTS - EnQPBEA Multiple Pattern Matching Algorithm Call for Simulated QCore 1

	printf("\n\tCORE 1:\n");
	printf("------------------------------------------------\n");
	int req_qubits1 = (2 * P_size1) + n1 + 2 + 2;
    Qureg qubits1 = createQureg(req_qubits1, env);
    initZeroState(qubits1);
    printf("\nQuantum Parameters for QPBE:\n");
    reportQuregParams(qubits1);
    reportQuESTEnv(env);
    printf("\nConstructing required quantum state..\n");
    for(int i = 0; i < n1; i++)
    	hadamard(qubits1, i);
    int set_qb1 = (n1 + 1);
    for(int idx = 0; idx < 2; idx++)
    {
    	for(int i = 0; i < anf_size1[idx]; i++)
    	{
    		if(anf1[idx][i] == 0)
    		{
    			pauliX(qubits1, set_qb1);
    			continue;
    		}
    		int *ctrls = (int *)malloc(sizeof(int) * n1);
    		int ctrl_size = 0;
    		int term = anf1[idx][i];
    		int qb = 0;
    		while(term)
    		{
    			if(term & 1) ctrls[ctrl_size++] = qb;
    			term >>= 1;
    			qb++;
    		}
    		multiControlledUnitary(qubits1, ctrls, ctrl_size, set_qb1, ux);
    		free(ctrls);
    	}
    	set_qb1--;
    }
    set_qb1 = (2*P_size1 + n1 + 2);
    for(int i = 0; i < P_size1; i++)
    {
    	int set = 0;
    	if(tolower(P1[i]) == 'a') set = 0;
    	else if(tolower(P1[i]) == 't') set = 1;
    	else if(tolower(P1[i]) == 'c') set = 2;
    	else if(tolower(P1[i]) == 'g') set = 3;
    	if(set & 2) pauliX(qubits1, set_qb1 - 1);
    	if(set & 1) pauliX(qubits1, set_qb1 - 2);
    	set_qb1 -= 2;
    }
	set_qb1 = (2*P_size1 + n1 + 3); //Ancilla qubits index (2*P_size + n + 3) & (2*P_size + n + 2)
	int eval_form1[][2] = { { 0, 0 }, { 1, 1 } };
	int ctrls_P1[2] = { set_qb1 - 2, set_qb1 - 3 }; //Qubits that store P1
	int ctrls_d_i1[2] = { n1 + 1, n1 }; //Qubits that store d_i
	for(int i = 0; i < 2; i++)
	{
		int ctrls[] = { ctrls_P1[i], ctrls_d_i1[i] };
		multiStateControlledUnitary(qubits1, ctrls, eval_form1[0], 2, set_qb1, ux);
		multiStateControlledUnitary(qubits1, ctrls, eval_form1[1], 2, set_qb1, ux);
		set_qb1--;
	}
	count = 0;
	ComplexMatrixN e1 = createComplexMatrixN(n1);
    for(int i = 0; i < (int)pow(2, n1); i++)
    {
    	if(is_at(T, T_size, i, P1, P_size1))
    	{
    		e1.real[i][i] = -1;
    		count++;
    	}
    	else e1.real[i][i] = 1;
    }
	int *targs1 = (int *)malloc(sizeof(int) * n1);
	for(int i = 0; i < n1; i++)
		targs1[i] = i;
	set_qb1 = (2*P_size1 + n1 + 3);
	int ctrls1[] = { set_qb1, set_qb1 - 1 }; //Using ancilla as control qubits
    printf("\nRunning Grover's Algorithm..\n");
    if(count != 0)
    {
	    int times =  (int)(3.14 * (pow(2, n1 / 2) / sqrt(count)) / 4);
	    for(int gi = 0; gi < times; gi++)
	    {
	    	//Marking
	    	multiControlledMultiQubitUnitary(qubits1, ctrls1, 2, targs1, n1, e1);
	    	
	        //Diffusion
	        for(int i = 0; i < n1; i++)
	            hadamard(qubits1, i);
	        for(int i = 0; i < n1; i++)
	            pauliX(qubits1, i);
	        multiControlledPhaseFlip(qubits1, targs1, n1);
	        for(int i = 0; i < n1; i++)
	            pauliX(qubits1, i);
	        for(int i = 0; i < n1; i++)
	            hadamard(qubits1, i);
	    }
	    qreal prob, max = 0.0;
	    for(int i = 0; i < (int)pow(2, (2 * P_size1) + n1 + 2 + 2); i++)
	    {
	    	prob = getProbAmp(qubits1, i);
	    	if(max <= prob) max = prob;
	    }
	    int mask = (int)pow(2, n1) - 1;
	    for(int i = 0; i < (int)pow(2, (2 * P_size1) + n1 + 2 + 2); i++)
	    {
	    	prob = getProbAmp(qubits1, i);
	    	if(fabs(max - prob) <= 0.0000001) printf("Correct Index %d with prob %f\n", (i & mask), prob);
	    }
	}
    else printf("No index exist containing the pattern.\n");
	destroyQureg(qubits1, env);
	destroyComplexMatrixN(e1);
	printf("\n------------------------------------------------\n\n");


	//END - EnQPBEA Multiple Pattern Matching Algorithm Call for Simulated QCore 1


	//START - EnQPBEA Multiple Pattern Matching Algorithm Call for Simulated QCore 2


	printf("\n\tCORE 2:\n");
	printf("------------------------------------------------\n");
	int req_qubits2 = (2 * P_size2) + n2 + 2 + 2;
    Qureg qubits2 = createQureg(req_qubits2, env);
    initZeroState(qubits2);
    printf("\nQuantum Parameters for QPBE:\n");
    reportQuregParams(qubits2);
    reportQuESTEnv(env);
    printf("\nConstructing required quantum state..\n");
    for(int i = 0; i < n2; i++)
    	hadamard(qubits2, i);
    int set_qb2 = (n2 + 1);
    for(int idx = 0; idx < 2; idx++)
    {
    	for(int i = 0; i < anf_size2[idx]; i++)
    	{
    		if(anf2[idx][i] == 0)
    		{
    			pauliX(qubits2, set_qb2);
    			continue;
    		}
    		int *ctrls = (int *)malloc(sizeof(int) * n2);
    		int ctrl_size = 0;
    		int term = anf2[idx][i];
    		int qb = 0;
    		while(term)
    		{
    			if(term & 1) ctrls[ctrl_size++] = qb;
    			term >>= 1;
    			qb++;
    		}
    		multiControlledUnitary(qubits2, ctrls, ctrl_size, set_qb2, ux);
    		free(ctrls);
    	}
    	set_qb2--;
    }
    set_qb2 = (2*P_size2 + n2 + 2);
    for(int i = 0; i < P_size2; i++)
    {
    	int set = 0;
    	if(tolower(P2[i]) == 'a') set = 0;
    	else if(tolower(P2[i]) == 't') set = 1;
    	else if(tolower(P2[i]) == 'c') set = 2;
    	else if(tolower(P2[i]) == 'g') set = 3;
    	if(set & 2) pauliX(qubits2, set_qb2 - 1);
    	if(set & 1) pauliX(qubits2, set_qb2 - 2);
    	set_qb2 -= 2;
    }
	set_qb2 = (2*P_size2 + n2 + 3); //Ancilla qubits index (2*P_size + n + 3) & (2*P_size + n + 2)
	int eval_form2[][2] = { { 0, 0 }, { 1, 1 } };
	int ctrls_P2[2] = { set_qb2 - 2, set_qb2 - 3 }; //Qubits that store P1
	int ctrls_d_i2[2] = { n2 + 1, n2 }; //Qubits that store d_i
	for(int i = 0; i < 2; i++)
	{
		int ctrls[] = { ctrls_P2[i], ctrls_d_i2[i] };
		multiStateControlledUnitary(qubits2, ctrls, eval_form2[0], 2, set_qb2, ux);
		multiStateControlledUnitary(qubits2, ctrls, eval_form2[1], 2, set_qb2, ux);
		set_qb2--;
	}
	count = 0;
	ComplexMatrixN e2 = createComplexMatrixN(n2);
    for(int i = 0; i < (int)pow(2, n2); i++)
    {
    	if(is_at(T, T_size, i, P2, P_size2))
    	{
    		e2.real[i][i] = -1;
    		count++;
    	}
    	else e2.real[i][i] = 1;
    }
	int *targs2 = (int *)malloc(sizeof(int) * n2);
	for(int i = 0; i < n2; i++)
		targs2[i] = i;
	set_qb2 = (2*P_size2 + n2 + 3);
	int ctrls2[] = { set_qb2, set_qb2 - 1 }; //Using ancilla as control qubits
    printf("\nRunning Grover's Algorithm..\n");
    if(count != 0)
    {
	    int times =  (int)(3.14 * (pow(2, n2 / 2) / sqrt(count)) / 4);
	    for(int gi = 0; gi < times; gi++)
	    {
	    	//Marking
	    	multiControlledMultiQubitUnitary(qubits2, ctrls2, 2, targs2, n2, e2);
	    	
	        //Diffusion
	        for(int i = 0; i < n2; i++)
	            hadamard(qubits2, i);
	        for(int i = 0; i < n2; i++)
	            pauliX(qubits2, i);
	        multiControlledPhaseFlip(qubits2, targs2, n2);
	        for(int i = 0; i < n2; i++)
	            pauliX(qubits2, i);
	        for(int i = 0; i < n2; i++)
	            hadamard(qubits2, i);
	    }
	    qreal prob, max = 0.0;
	    for(int i = 0; i < (int)pow(2, (2 * P_size2) + n2 + 2 + 2); i++)
	    {
	    	prob = getProbAmp(qubits2, i);
	    	if(max <= prob) max = prob;
	    }
	    int mask = (int)pow(2, n2) - 1;
	    for(int i = 0; i < (int)pow(2, (2 * P_size2) + n2 + 2 + 2); i++)
	    {
	    	prob = getProbAmp(qubits2, i);
	    	if(fabs(max - prob) <= 0.0000001) printf("Correct Index %d with prob %f\n", (i & mask), prob);
	    }
	}
    else printf("No index exist containing the pattern.\n");
	destroyQureg(qubits2, env);
	destroyComplexMatrixN(e2);
	printf("\n------------------------------------------------\n\n");


	//END - EnQPBEA Multiple Pattern Matching Algorithm Call for Simulated QCore 2


	//START - EnQPBEA Multiple Pattern Matching Algorithm Call for Simulated QCore 3


	printf("\n\tCORE 3:\n");
	printf("------------------------------------------------\n");
	int req_qubits3 = (2 * P_size3) + n3 + 2 + 2;
    Qureg qubits3 = createQureg(req_qubits3, env);
    initZeroState(qubits3);
    printf("\nQuantum Parameters for QPBE:\n");
    reportQuregParams(qubits3);
    reportQuESTEnv(env);
    printf("\nConstructing required quantum state..\n");
    for(int i = 0; i < n3; i++)
    	hadamard(qubits3, i);
    int set_qb3 = (n3 + 1);
    for(int idx = 0; idx < 2; idx++)
    {
    	for(int i = 0; i < anf_size3[idx]; i++)
    	{
    		if(anf3[idx][i] == 0)
    		{
    			pauliX(qubits3, set_qb3);
    			continue;
    		}
    		int *ctrls = (int *)malloc(sizeof(int) * n3);
    		int ctrl_size = 0;
    		int term = anf3[idx][i];
    		int qb = 0;
    		while(term)
    		{
    			if(term & 1) ctrls[ctrl_size++] = qb;
    			term >>= 1;
    			qb++;
    		}
    		multiControlledUnitary(qubits3, ctrls, ctrl_size, set_qb3, ux);
    		free(ctrls);
    	}
    	set_qb3--;
    }
    set_qb3 = (2*P_size3 + n3 + 2);
    for(int i = 0; i < P_size3; i++)
    {
    	int set = 0;
    	if(tolower(P3[i]) == 'a') set = 0;
    	else if(tolower(P3[i]) == 't') set = 1;
    	else if(tolower(P3[i]) == 'c') set = 2;
    	else if(tolower(P3[i]) == 'g') set = 3;
    	if(set & 2) pauliX(qubits3, set_qb3 - 1);
    	if(set & 1) pauliX(qubits3, set_qb3 - 2);
    	set_qb3 -= 2;
    }
	set_qb3 = (2*P_size3 + n3 + 3); //Ancilla qubits index (2*P_size + n + 3) & (2*P_size + n + 2)
	int eval_form3[][2] = { { 0, 0 }, { 1, 1 } };
	int ctrls_P3[2] = { set_qb3 - 2, set_qb3 - 3 }; //Qubits that store P1
	int ctrls_d_i3[2] = { n3 + 1, n3 }; //Qubits that store d_i
	for(int i = 0; i < 2; i++)
	{
		int ctrls[] = { ctrls_P3[i], ctrls_d_i3[i] };
		multiStateControlledUnitary(qubits3, ctrls, eval_form3[0], 2, set_qb3, ux);
		multiStateControlledUnitary(qubits3, ctrls, eval_form3[1], 2, set_qb3, ux);
		set_qb3--;
	}
	count = 0;
	ComplexMatrixN e3 = createComplexMatrixN(n3);
    for(int i = 0; i < (int)pow(2, n3); i++)
    {
    	if(is_at(T, T_size, i, P3, P_size3))
    	{
    		e3.real[i][i] = -1;
    		count++;
    	}
    	else e3.real[i][i] = 1;
    }
	int *targs3 = (int *)malloc(sizeof(int) * n3);
	for(int i = 0; i < n3; i++)
		targs3[i] = i;
	set_qb3 = (2*P_size3 + n3 + 3);
	int ctrls3[] = { set_qb3, set_qb3 - 1 }; //Using ancilla as control qubits
    printf("\nRunning Grover's Algorithm..\n");
    if(count != 0)
    {
	    int times =  (int)(3.14 * (pow(2, n3 / 2) / sqrt(count)) / 4);
	    for(int gi = 0; gi < times; gi++)
	    {
	    	//Marking
	    	multiControlledMultiQubitUnitary(qubits3, ctrls3, 2, targs3, n3, e3);
	    	
	        //Diffusion
	        for(int i = 0; i < n3; i++)
	            hadamard(qubits3, i);
	        for(int i = 0; i < n3; i++)
	            pauliX(qubits3, i);
	        multiControlledPhaseFlip(qubits3, targs3, n3);
	        for(int i = 0; i < n3; i++)
	            pauliX(qubits3, i);
	        for(int i = 0; i < n3; i++)
	            hadamard(qubits3, i);
	    }
	    qreal prob, max = 0.0;
	    for(int i = 0; i < (int)pow(2, (2 * P_size3) + n3 + 2 + 2); i++)
	    {
	    	prob = getProbAmp(qubits3, i);
	    	if(max <= prob) max = prob;
	    }
	    int mask = (int)pow(2, n3) - 1;
	    for(int i = 0; i < (int)pow(2, (2 * P_size3) + n3 + 2 + 2); i++)
	    {
	    	prob = getProbAmp(qubits3, i);
	    	if(fabs(max - prob) <= 0.0000001) printf("Correct Index %d with prob %f\n", (i & mask), prob);
	    }
	}
    else printf("No index exist containing the pattern.\n");
	destroyQureg(qubits3, env);
	destroyComplexMatrixN(e3);
	printf("\n------------------------------------------------\n\n");


	//END - EnQPBEA Multiple Pattern Matching Algorithm Call for Simulated QCore 3

	
    destroyQuESTEnv(env);

	
	//Destroy Quantum Environment (env) Realized Three Quantum System (qubits1, qubits2, qubits3)
	//END - EnQPBEA-MPM Algorithm based on the Multiple Quantum Core (QCore = 3)


	getch();

    return 0;
}