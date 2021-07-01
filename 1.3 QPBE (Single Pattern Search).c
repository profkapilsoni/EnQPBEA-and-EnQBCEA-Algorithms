/*
*********************************************************************************************************
*	About: This code implements existing QPBE Algorithm as restricted singleton pattern set |P|=1 as P	*
*	Usage: Run in command prompt "QPBE.exe <input_file>"												*
*	1. 	<input_file> contains character subset of biological gene text sequence of SARS-CoV-2 for Human	*
* 		A file size {128, 256, 512} characters are intentioally prepared for feasible QuEST Simulation 	*
*	2.	The searching of existing QPBE Algorithm is performed for the fixed pattern P of the length 3	*
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

	//ENDS - Reading and Storing Input File in the array "T"
	
	//An arbitrary pattern P = "ATG" used for searching within the text "T"
	char P[] =  { 'a', 't', 'g' };
	const int P_size = 3;

	//Printing all the inputs
	printf("\nDNA sequence is:\n");
	for(int i = 1; i <= T_size; i++)
	{
		printf("%c", T[i - 1]);
		if(i % 32 == 0) printf("\n");
	}
	printf("\nPattern to search: ");
	for(int i = 0; i < P_size; i++)
		printf("%c", P[i]);
	printf("\n");

	//Display Warning if necessary
	
	int n = (int)(log(count)/log(2));

	ch = 'y';
	int req_qubits = (2 * P_size) + n + 2 + 2;
	if(req_qubits > 25)
	{
		printf("\nWARNING: We do not recommend to run the QPBE Algorithm since it requires %d qubits!", req_qubits);
		printf("\nEnter your choice if you want to continue (y/n): ");
		ch = getch();
		ch = tolower(ch);
	}
	else printf("\nThe QPBE Algorithm requires %d qubits.\n", req_qubits);

	if(ch != 'y')
	{
		printf("\nExiting!\n");
		return 0;
	}
	else if(ch == 'y' && req_qubits > 25)
	{
		printf("\n\nYou chose to ignore our warning and continue!\nGOOD LUCK!\n");
	}

	//STARTS - Preliminaries for quantum circuit
	//NOTE: This is not the efficient approach, however, can be improved for unrestricted simulation

	int **anf, *anf_size;
    int *f = (int *)malloc(sizeof(int) * (int)pow(2, n));
    for(int i = 0; i < T_size; i++) //T_size = 2^n
    {
    	if(tolower(T[i]) == 'a') f[i] = 0;
    	else if(tolower(T[i]) == 't') f[i] = 1;
    	else if(tolower(T[i]) == 'c') f[i] = 2;
    	else if(tolower(T[i]) == 'g') f[i] = 3;
    }
    anf_calc(f, n, 2, &anf, &anf_size);

	//ENDS - Preliminaries for quantum circuit

	//STARTS - Existing QPBE ALgorithm for Single Pattern Matching through Quantum Exact Match (QEM) & Grover's Search Operator (GSO) Call
	//Quantum Environment (env) Realizing Single Quantum System (qubits) = Simulates One Quantum Core over Multithreaded Classical Core

    QuESTEnv env = createQuESTEnv();

    //"n" qubits to store indexes "i", 2 qubits to store "d_i", 2*P_size qubits to store "P"
    //A = 00, T = 01, C = 10, G = 11
    //The state is stored in following sequence: |an>|P>|d_i>|i>
    // where |an> represents ancilla qubits (2-qubit)
    Qureg qubits = createQureg(req_qubits, env);
    initZeroState(qubits);
    printf("\nQuantum Parameters for QPBE:\n");
    reportQuregParams(qubits);
    reportQuESTEnv(env);

    //Single qubit unitary for XOR gate
    ComplexMatrix2 ux = {
    	.real={{0,1},{1,0}},
		.imag={{0,0},{0,0}}
    };

    //Print status of quantum state construction
    printf("\nConstructing required quantum state..\n");

    //Create superposition of all indexes from qubits indexed from 0 to (n - 1)
    for(int i = 0; i < n; i++)
    	hadamard(qubits, i);

    //Store "d_i" corresponding to each "i"
    //"d_i" are stored in qubits indexed from (n + 1) to (n)
    int set_qb = (n + 1);
    for(int idx = 0; idx < 2; idx++)
    {
    	for(int i = 0; i < anf_size[idx]; i++)
    	{
    		if(anf[idx][i] == 0)
    		{
    			pauliX(qubits, set_qb);
    			continue;
    		}
    		int *ctrls = (int *)malloc(sizeof(int) * n);
    		int ctrl_size = 0;
    		int term = anf[idx][i];
    		int qb = 0;
    		while(term)
    		{
    			if(term & 1) ctrls[ctrl_size++] = qb;
    			term >>= 1;
    			qb++;
    		}
    		multiControlledUnitary(qubits, ctrls, ctrl_size, set_qb, ux);
    		free(ctrls);
    	}
    	set_qb--;
    }

    //Store the pattern "P" in qubits indexed from "2*P_size n + 2 - 1" to "n + 2"
    //If "P" = (P1, P2, ..., Pm) then
    //"P1" is stored in (2*P_size n + 2 - 1)-th and (2*P_size n + 2 - 2)-th qubit and
    // so on to "Pm" in (n + 3)-th and (n + 2)-th qubit
    //A = 00, T = 01, C = 10, G = 11
    set_qb = (2*P_size + n + 2);
    for(int i = 0; i < P_size; i++)
    {
    	int set = 0;
    	if(tolower(P[i]) == 'a') set = 0;
    	else if(tolower(P[i]) == 't') set = 1;
    	else if(tolower(P[i]) == 'c') set = 2;
    	else if(tolower(P[i]) == 'g') set = 3;
    	if(set & 2) pauliX(qubits, set_qb - 1);
    	if(set & 1) pauliX(qubits, set_qb - 2);
    	set_qb -= 2;
    }

/*
    qreal prob;
    int mask = (int)pow(2, n) - 1;
    for(int i = 0; i < (int)pow(2, (2 * P_size) + n + 2); i++)
    {
    	prob = getProbAmp(qubits, i);
    	if(prob != 0.0)
    		printf("i = %d d = %X P = %X has prob %f\n", 
    			(i & mask),
    			(i & ((1 << ((n + 1)) ^ (1 << n)))) >> n,
    			(i >> (n + 2)),
    			prob);
    }
*/

	//Evaluate d_i == P_1
	set_qb = (2*P_size + n + 3); //Ancilla qubits index (2*P_size + n + 3) & (2*P_size + n + 2)
	int eval_form[][2] = { { 0, 0 }, { 1, 1 } };
	int ctrls_P1[2] = { set_qb - 2, set_qb - 3 }; //Qubits that store P1
	int ctrls_d_i[2] = { n + 1, n }; //Qubits that store d_i
	for(int i = 0; i < 2; i++)
	{
		int ctrls[] = { ctrls_P1[i], ctrls_d_i[i] };
		multiStateControlledUnitary(qubits, ctrls, eval_form[0], 2, set_qb, ux);
		multiStateControlledUnitary(qubits, ctrls, eval_form[1], 2, set_qb, ux);
		set_qb--;
	}


	//Prints the quantum state before applying Grover's algorithm
/*	printf("Required Quantum State is constructed.\n");
	printf("Do you want to view the constructed quantum state?(y/n):");
	ch = getch();
	if(ch == 'y')
	{
		printf("\nConstructed quantum state before applying Grover's search:\n");
		qreal prob;
	    int mask = (int)pow(2, n) - 1;
	    for(int i = 0; i < (int)pow(2, (2 * P_size) + n + 2 + 2); i++)
	    {
	    	prob = getProbAmp(qubits, i);
	    	if(prob != 0.0)
	    		printf("i = %d d = %X P1 = %X (d_i & P1) = %X has prob %f\n", 
	    			(i & mask),
	    			(i & ((1 << ((n + 1)) ^ (1 << n)))) >> n,
	    			(i >> (2*(P_size - 1) + n + 2)) & ((int)pow(2, 2) - 1),
	    			(i >> (2*P_size + n + 2)),
	    			prob);
	    }
	}
*/
	count = 0;
	ComplexMatrixN e = createComplexMatrixN(n);
    for(int i = 0; i < (int)pow(2, n); i++)
    {
    	if(is_at(T, T_size, i, P, P_size))
    	{
    		e.real[i][i] = -1;
    		count++;
    	}
    	else e.real[i][i] = 1;
    }

	int *targs = (int *)malloc(sizeof(int) * n);
	for(int i = 0; i < n; i++)
		targs[i] = i;
	set_qb = (2*P_size + n + 3);
	int ctrls[] = { set_qb, set_qb - 1 }; //Using ancilla as control qubits
    printf("\nRunning Grover's Algorithm..\n");
    if(count != 0)
    {
	    int times =  (int)(3.14 * (pow(2, n / 2) / sqrt(count)) / 4);
	    for(int gi = 0; gi < times; gi++)
	    {
	    	//Marking
	    	multiControlledMultiQubitUnitary(qubits, ctrls, 2, targs, n, e);
	    	
	        //Diffusion
	        for(int i = 0; i < n; i++)
	            hadamard(qubits, i);
	        for(int i = 0; i < n; i++)
	            pauliX(qubits, i);
	        multiControlledPhaseFlip(qubits, targs, n);
	        for(int i = 0; i < n; i++)
	            pauliX(qubits, i);
	        for(int i = 0; i < n; i++)
	            hadamard(qubits, i);
	    }

	    qreal prob, max = 0.0;
	    for(int i = 0; i < (int)pow(2, (2 * P_size) + n + 2 + 2); i++)
	    {
	    	prob = getProbAmp(qubits, i);
	    	if(max <= prob) max = prob;
	    }
	    int mask = (int)pow(2, n) - 1;
	    for(int i = 0; i < (int)pow(2, (2 * P_size) + n + 2 + 2); i++)
	    {
	    	prob = getProbAmp(qubits, i);
	    	if(fabs(max - prob) <= 0.0000001) printf("Correct Index %d with prob %f\n", (i & mask), prob);
	    }
	}
    else printf("No index exist containing the pattern.\n");

    //ENDS - Existing QPBE ALgorithm for Single Pattern Matching

	destroyQureg(qubits, env);
	destroyComplexMatrixN(e);
    destroyQuESTEnv(env);

	getch();

    return 0;
}