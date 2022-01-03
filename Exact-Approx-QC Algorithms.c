/*
*************************************************************************************************************************
*	About: This code implements Quantum Counting Algorithms Exact-QC and Approx.-QC.				*
*	Usage: Run in command prompt "Exact-Approx-QC Algorithms.exe <input_file>"					*					
*	1. 	<input_file> contains character subset of biological gene text sequence of SARS-CoV-2 for Human.	*
* 		A file sizes {128, 256, 512} characters are intentioally prepared for feasible QuEST Simulation 	*
*	2.	Exact-QC is implemented with log_N and log_t precision qubits for EnQPBEA-MPM and EnQBCEA-MPM. 		*
*	3.	Approx.-QC is implemented with log_N-2 and log_t-2 precision qubits for EnQPBEA-MPM and EnQBCEA-MPM.	*
*	4.	No. of pattern occurrences "t" is being identified for equal and unequal sized patterns.  		*
*	5.	The possible cases (m <= C) and (m > C) are considered in our experimentations os simulation. 		*											*
*	6.	Exact-QC and Approx.-QC are repeated 10 times on individual quantum core for each pattern.		*	
*************************************************************************************************************************
*/


#include <stdio.h>
#include <conio.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "QuEST.h"

//Value of Pi
const double pi = 3.14;

//Logic realized for Quantum Exact Match (QEM) operation
int is_at(char *T, int T_size, int idx, char *P, int P_size)
{
	int found = 1;
	int temp_i = 0;
	if((idx + P_size) >= T_size) return 0;
	for(int i = idx; i < idx + P_size; i++)
		if(T[i] != P[temp_i++]) found = 0;
	return found;
}


//Quantum counting circuit is implemented here for Exact-QC and Approx.-QC.
//Quantum environment is created as to run on individual quantum cores.
double count_est_amp(char *T, int count, char *Pat, int p_len, int l)
{
	int n = (int)(log(count)/log(2));
	int N = pow(2, n);


    QuESTEnv env = createQuESTEnv();

    //"p" represents the precision qubits
    //"P = 2^p"
    int p = l, P = pow(2, p);
    Qureg qubits = createQureg((p + n), env);
    initZeroState(qubits);
    //reportQuregParams(qubits);
    //reportQuESTEnv(env);

    //Unitary matrix creation
	ComplexMatrixN e = createComplexMatrixN(n);
    for(int i = 0; i < N; i++)
    {
    	if(is_at(T, count, i, Pat, p_len)) e.real[i][i] = -1;
    	else e.real[i][i] = 1;
    }

    //Apply Hadamard on all the (p + n)-qubits to create superposition
    //Here each element in last n-qubit superposition is considered to represent "x"
    for(int i = 0; i < (p + n); i++)
        hadamard(qubits, i);
    
    //Target qubits on which grover's search must be performed to count the number of correct result
    int targs[n];
    for(int i = 0; i < n; i++)
    	targs[i] = p + i;

    //Applying controlled grover's operator
    int iterations = 1;
    for(int gi = p - 1; gi >= 0; gi--)
    {
    	int ctrls[] = { gi };
    	for(int j = 0; j < iterations; j++)
    	{
    		//Marking
	    	multiControlledMultiQubitUnitary(qubits, ctrls, 1, targs, n, e);
	    	
	        //Diffusion
	        for(int i = p; i < (n + p); i++)
	            hadamard(qubits, i);
	        for(int i = p; i < (n + p); i++)
	            pauliX(qubits, i);
	        multiControlledPhaseFlip(qubits, targs, n);
	        for(int i = p; i < (n + p); i++)
	            pauliX(qubits, i);
	        for(int i = p; i < (n + p); i++)
	            hadamard(qubits, i);
    	}
    	iterations *= 2;
    }
	
    //Inverse QFT on the first p-qubits
    for(int i = 0; i < p; i++)
    {
    	for(int j = 0; j < i; j++)
    	{
    		qreal angle = -2.0 * pi / (pow(2, i - j + 1));
    		controlledPhaseShift(qubits, j, i, angle);
    	}
    	hadamard(qubits, i);
    }
    
    //Measures the qubits to get the value of first p-qubits
    qreal prob = 0, max_prob = 0;
    int p_val = 0;
    for(int i = 0; i < P; i++)
    {
    	prob = getProbAmp(qubits, i);
    	if(prob > max_prob)
    	{
    		max_prob = prob;
    		p_val = i;
    	}
    }

    destroyQureg(qubits, env);
    destroyComplexMatrixN(e);
    destroyQuESTEnv(env);

	if(p_val > (P / 2)) p_val = P - p_val;
	return (N * pow(sin(p_val * pi / P), 2));
}


int main (int narg, char *varg[])
{
	if(narg != 3)
	{
		printf("Usage: %s <input_file> <pattern_string>", varg[0]);
		return 0;
	}
	
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

	int p_len = 0;
	while(varg[2][p_len] != '\0')
		p_len++;

	/*int real_count = 0;
	for(int i = 0; i < count - p_len + 1; i++)
		if(is_at(T, count, i, varg[2], p_len)) real_count++;
	printf("\nReal Count: %d", real_count);*/

	printf("\nExact Count = %.0f", round(count_est_amp(T, count, varg[2], p_len, (int)(log(count)/log(2)))));
	printf("\nApprox Count = %.0f\n", round(count_est_amp(T, count, varg[2], p_len, (int)(log(count)/log(2) - 2))));
	//printf("\n%f\n", Basic_Approx_Count(T, count, varg[2], p_len, 0.1));

    return 0;
}