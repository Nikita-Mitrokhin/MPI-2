#include "stdio.h"
#include "mpi.h"
#include "fstream"
#include "ctime"
#include "iostream"
double startT, stopT;


int* mergeArrays(int *pMas1, int size1, int *pMas2, int size2, int *result)

{

	int i = 0, j = 0, k = 0;

	while (i < size1 && j < size2)

		if (pMas1[i] < pMas2[j])

			result[k++] = pMas1[i++];

		else

			result[k++] = pMas2[j++];

	if (i == size1)

		while (j < size2)

			result[k++] = pMas2[j++];

	if (j == size2)

		while (i < size1)

			result[k++] = pMas1[i++];

	return result;

}

void swap(int* v, int i, int j)
{
	int t;
	t = v[i];
	v[i] = v[j];
	v[j] = t;
}

void sort(int* v, int n)
{
	int i, j;

	for (i = n - 2; i >= 0; i--)
		for (j = 0; j <= i; j++)
			if (v[j] > v[j + 1])
				swap(v, j, j + 1);
}

using namespace std;

void main(int argc, char* argv[])
{
	int* data;
	int* resultant_array;
	int* sub;

	int m, n;
	int Rank, Size;
	int i, s;
	int step;

	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	MPI_Comm_size(MPI_COMM_WORLD, &Size);
	n = 4;

	if (argc >= 2) {
		n = atoi(argv[1]);
	}

	s = n / Size;

	data = new int[n];

	srand(unsigned int(MPI_Wtime()));

	for (i = 0; i < n; i++) {
		data[i] = rand() % 10000;
		}

	if (Rank == 0) {

		

		startT = MPI_Wtime();
		MPI_Bcast(&s, 1, MPI_INT, 0, MPI_COMM_WORLD);
		resultant_array = new int[s];
		MPI_Scatter(data, s, MPI_INT, resultant_array, s, MPI_INT, 0, MPI_COMM_WORLD);
		sort(resultant_array, s);
	}
	else {

		MPI_Bcast(&s, 1, MPI_INT, 0, MPI_COMM_WORLD);
		resultant_array = new int[s];
		MPI_Scatter(data, s, MPI_INT, resultant_array, s, MPI_INT, 0, MPI_COMM_WORLD);
		sort(resultant_array, s);
	}

	step = 1;

	while (step < Size) {
		if (Rank % (2 * step) == 0) {
			if (Rank + step < Size) {
				MPI_Recv(&m, 1, MPI_INT, Rank + step, 0, MPI_COMM_WORLD, &status);
				sub = new int[m];
				MPI_Recv(sub, m, MPI_INT, Rank + step, 0, MPI_COMM_WORLD, &status);
				int *tmp = new int[s + m];
				resultant_array = mergeArrays(resultant_array, s, sub, m, tmp);
				s = s + m;
			}
		}
		else {
			int near = Rank - step;
			MPI_Send(&s, 1, MPI_INT, near, 0, MPI_COMM_WORLD);
			MPI_Send(resultant_array, s, MPI_INT, near, 0, MPI_COMM_WORLD);
			break;
		}

		step = step * 2;
	}


	if (Rank == 0) {
		stopT = MPI_Wtime();
		double parallelTime = stopT - startT;

		for (i = 0; i < n; i++)
		{
			cout << resultant_array[i] << " ";
		}

		printf("\n\n\n Parallel Time: %f", parallelTime);

		startT = MPI_Wtime();
		sort(data, n);
		stopT = MPI_Wtime();

		double poslTime = stopT - startT;
		printf("\n Not Parallel Time: %f \n\t", stopT - startT);
		printf("\n SpeedUp: %f \n\n\n", poslTime / parallelTime);

		
	}

	MPI_Finalize();
}