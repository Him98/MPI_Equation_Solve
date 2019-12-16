#include <bits/stdc++.h>
#include <mpi.h>
#include <sys/time.h>
using namespace std;
typedef long long int ll;
#define pb push_back
#define mp make_pair
#define f first
#define s second
#define sc(n) scanf("%d",&n)
#define scl(n) scanf("%lld",&n)
#define pr(n) printf("%d",n)
#define prl(n) printf("%lld",n)
#define nl printf("\n")
#define fast_io ios_base::sync_with_stdio(false);cin.tie(nULL)
void swap(double *b,double *c)
{
	double tmp = *b;
	*b = *c;
	*c = tmp;
}

int n,m;
typedef struct{
		double item;
		int row;
	} double_int_pair;
// void msg_passing(int fl, double (*arr),int irow, int pivot, MPI_Comm comm, MPI_Status status,int curr_id)
// {
// 	if(fl)
// 	{
// 		MPI_Recv(&arr[irow], n+1, MPI_DOUBLE, curr_id, 77, comm, &status);
// 		MPI_Send(&arr[pivot], n+1, MPI_DOUBLE, curr_id, 88, comm);
// 	}
// 	return;
// }
int main(int argc, char **argv)
{
	int curr_id, total, curr_id_from, curr_id_to;
	// Initializing MPI
	MPI_Init(&argc, &argv);
	MPI_Comm comm;
	MPI_Status status;
	comm = MPI_COMM_WORLD;

	// Get current id or rank
	MPI_Comm_rank(comm, &curr_id);
	if(curr_id == 0)
		scanf("%d %d\n",&n,&m);
	MPI_Barrier(comm);
	MPI_Bcast(&n, 1, MPI_INT, 0, comm);
	MPI_Bcast(&m, 1, MPI_INT, 0, comm);
	// Total number of processors
	MPI_Comm_size(comm, &total);
	timeval start_time, end_time;
	int irow,jrow, i,j, pivot;
	double t, amul, item, item_max, time_taken = 0.0;
	double aug_mat[n][n+1];
	double x[n+1];


	if(curr_id==0)
	{
		srand(time(NULL));
		double val;
		for(i=0;i<n;++i)
		{
			for(j=0;j<n;++j)
			{
				scanf("%lf",&val);
				aug_mat[i][j]=val;
			}
			// b value stored along with a, denoted as augmented matrix
			scanf("%lf",&val);
			aug_mat[i][j]=(-1*val);
		}
		gettimeofday(&start_time, NULL);
	}
	MPI_Barrier(comm);
	MPI_Bcast(aug_mat, n*(n+1), MPI_DOUBLE, 0, comm);

	for(irow=0;irow<n-1;irow++)
	{
		double_int_pair local, global;
		item_max=fabs(aug_mat[irow][irow]);
		local.item = item_max;
		//the thread who has the next row to eliminate
		// Broadcasting this row to all others
		curr_id_from = n-1;
		curr_id_to = n-1;
		curr_id_from = (curr_id_from-irow)%total;
		MPI_Bcast(&item_max, 1, MPI_DOUBLE, curr_id_from, comm);
		
		//pivoting step
		pivot=irow;
		local.row = pivot;
		jrow=n-1-curr_id;
		for(;jrow>=irow+1;jrow-=total)
			if(max(item_max,fabs(aug_mat[jrow][irow]))>item_max)
			{
				local.item=max(item_max,fabs(aug_mat[jrow][irow]));
				local.row=jrow;
				item_max = local.item;
				pivot = local.row;
			}
		//Interchange the rows, depending on pivot

		MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
		
		// The rank of the child processes with whom pivot will be exchanged
		// curr_id_to=n-1;
		curr_id_to=(curr_id_to-(global.row))%total;
		pivot = global.row;
		if(curr_id_from==curr_id_to)
		{ 
			//replacement at the same slave no communication
			if(curr_id==curr_id_from && pivot!=irow)
				for(j=0;j<n+1;j++)
					swap(&aug_mat[irow][j],&aug_mat[pivot][j]);
		}
		else if(curr_id_from != curr_id_to)
		{ //communication is needed for pivot exchange
			for(j=n;j>=irow;j--)
				if(curr_id==curr_id_from)
					aug_mat[pivot][j]=aug_mat[irow][j];
				else if(curr_id==curr_id_to)
					aug_mat[irow][j]=aug_mat[pivot][j];
			if(curr_id==curr_id_from)
			{
				// msg_passing(1,(*aug_mat),irow,pivot,comm,status,curr_id_to);
				MPI_Recv(&aug_mat[irow], n+1, MPI_DOUBLE, curr_id_to, 77, comm, &status);
				MPI_Send(&aug_mat[pivot], n+1, MPI_DOUBLE, curr_id_to, 88, comm);
			}
			if(curr_id==curr_id_to)
			{
				// for(j=n;j>=irow;j--)
				// 	aug_mat[irow][j]=aug_mat[pivot][j];
				MPI_Send(&aug_mat[irow], n+1, MPI_DOUBLE, curr_id_from, 77, comm);
				MPI_Recv(&aug_mat[pivot], n+1, MPI_DOUBLE, curr_id_from, 88, comm, &status);
			}
		}

		MPI_Bcast(&aug_mat[irow], n+1, MPI_DOUBLE, curr_id_from, comm);

		// Eliminating the x variable accordingly
		t=1.0/aug_mat[irow][irow];
		for(jrow=n-1-curr_id;jrow>irow;)
		{
			amul = aug_mat[jrow][irow] * t;
			for(j=n;j>=0;j--)
				aug_mat[jrow][j] -= (amul * aug_mat[irow][j]);
			jrow -= total;
		}
	}

	if(curr_id==0)
	{
		// cout<<"The upper triangular matrix"<<endl;
		// PrintMatrix();

		//Back Substitution
		for(i=n-1;i>=0;i--)
		{
			x[i] = (-1.0 * aug_mat[i][n]);
			for(j=i-1;j>=0;j--)
			{
				aug_mat[j][n] += ((aug_mat[j][i]*x[i])/aug_mat[i][i]);
				aug_mat[j][i] = 0;
			}
			x[i] /= aug_mat[i][i];
			
		}
		gettimeofday(&end_time, NULL);
		for(i=0;i<n;i++)
			printf("%.5lf ", x[i]);
		cout<<endl;
		time_taken = (end_time.tv_sec - start_time.tv_sec);
		time_taken *= 1000.0;
		double micro = (end_time.tv_usec - start_time.tv_usec); 
		micro /= 1000.0;
		time_taken += micro;
		printf("Time taken: %lf milliseconds\n", time_taken);
	}

	MPI_Finalize();
}