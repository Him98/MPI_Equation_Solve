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
#define fast_io ios_base::sync_with_stdio(false);cin.tie(NULL)
MPI_Comm comm;
MPI_Status status;
void swap(double *b,double *c)
{
	double tmp = *b;
	*b = *c;
	*c = tmp;
}
const double eps = 1e-4;
int curr_id, total, tt, start, end, size, size2, tag;
int k1 = 100000, k2 = 200000, k3 = 300000;
int M, n;
int col;
double *t1, *t2, *x;
double ans = 0.0, summation = 0.0;
double *nr, *res;
int k = 1,i,j,m;
double* incr(double *r)
{
	for(i = 0; i < n; i++)
		r[i] += res[i];
	return r;
}

void proc_mul()
{
	if(++tt == 1)
		MPI_Recv(t1, col * M, MPI_DOUBLE, 0, k1 + k * total + curr_id, comm, &status);
	summation = 0.0;
	memset(res,0,n*sizeof(double));
	MPI_Recv(t2, col, MPI_DOUBLE, 0, k2 + k * total + curr_id, comm, &status);
	for(i = 0;i < n;i++)
	{
		summation = 0.0;
		for (j = 0;j < col;j++)
			summation += t1[j * n + i] * t2[j];
		res[i] += summation;
	}
	MPI_Send(res, M, MPI_DOUBLE, 0, k3 + k * total + curr_id, comm);
}

void vector_dot(double *a, double *b)
{
	ans = 0.0;
	i = n-1;
	for (;i>=0;)
	{
		double var = a[i] * b[i];
		ans += var;
		i--;
	}
}

void mat_multiply(double **a, double *b, double *r) {
	//Dividing matrix, sending to different processes
	if(++tt == 1)
	{
		for(m = 1; m < total; m++)
		{
			size = 0, start = (m-1), end = m;
			start *= col;
			j = start;
			end *= col;
			for(;j < end; j++)
				for(i = 0; i < M; )
					t1[size++] = a[i++][j];
			MPI_Send(t1, size, MPI_DOUBLE, m, k1 + k * total + m, comm);
		}
	}
	for(m = 1;m < total;m++)
	{
		size2 = 0, start = (m-1), end = m;
		start *= col;
		j = start;
		end *= col;
		for(;j < end;)
			t2[size2++] = b[j++];
		MPI_Send(t2, size2, MPI_DOUBLE, m, k2 + k * total + m, comm);
	}
	memset(r,0,n*sizeof(double));
	//Receiving answers from other processes
	for(m = 1; m < total; m++)
	{
		tag = k * total;
		tag += k3;
		tag += m;
		MPI_Recv(res, M, MPI_DOUBLE, m, tag, comm, &status);
		r = incr(r);
	}
	//Returning value in r
	start = (total - 1) * col;
	for(i = 0;i < n; i++)
		for(j = start; j < n; j++)
			r[i] += a[i][j] * b[j];
}

void solve2(double **a, double *b, double *r, double *p)
{
	int flag;
	proc_mul();
	int con = 0;
	while(!con)
	{
		MPI_Recv(&flag, 1, MPI_INT, 0, (k+1) * total + curr_id, comm, &status);
		k++;
		if(flag == 1)
			break;
		proc_mul();
	}
}
bool check(double *r)
{
	double res = 0;
	for (i = 0; i < n; i++)
	{
		res += fabs(r[i]);
		if(i == n-1 && res>=eps)
			return false;
		else if(i == n-1)
			return true;	
	}
}
// Main function to run the original root process
void solve(double **a, double *b, double *r, double *p)
{
	double *tmp;
	tmp = new double[n];
	memset(x,(rand()%100)+1,n*sizeof(double));
	nr = new double[n];
	mat_multiply(a, x, tmp);
	for (i = 0; i < M; i++)
	{
		r[i] = b[i] - tmp[i];
		p[i] = r[i];
	}
	int flag,g = 1;
	while(g == 1)
	{
		k++;
		flag = 0;
		for(m = 1; m < total; m++)
			MPI_Send(&flag, 1, MPI_INT, m, k * total + m, comm);
		ans = 0.0;
		vector_dot(r, r);
		double sqr_r = ans;
		mat_multiply(a, p, tmp);
		vector_dot(p, tmp);
		double alpha = sqr_r/ans;
		for(i = 0; i < n; i++)
		{
			x[i] += alpha * p[i];
			nr[i] = r[i] - alpha * tmp[i];
		}
		if(check(nr))
			break;
		vector_dot(nr, nr);
		double sqr_nr = ans;
		vector_dot(r, r);
		double beta = sqr_nr/ans;
		for(i = 0; i < n; i++)
		{
			p[i] = nr[i] + beta * p[i];
			r[i] = nr[i];
		}
	}
	k++;
	if(!flag)
		flag = 1;
	for (m = 1; m < total; m++)
		MPI_Send(&flag, 1, MPI_INT, m, k * total + m, comm);
}
int main(int argc, char **argv)
{
	double time_taken = 0.0;
	double **a;
	double *b,*r, *p;
	
	MPI_Init(&argc, &argv);
	
	comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &curr_id);
    if(curr_id == 0)
		scanf("%d %d\n",&n,&m);
	MPI_Barrier(comm);
	MPI_Bcast(&n, 1, MPI_INT, 0, comm);
	M = n;
	MPI_Comm_size(comm, &total);
    col = n/total;
    a = new double *[n];
	for(i = 0; i < n; i++)
		a[i] = new double [n];
    t1 = new double[n * col];
    timeval start_time, end_time;
	
	t2 = new double[col];
	res = new double[n];
	// memset(b,0,n*sizeof(double));
	x = new double [n];
	memset(x,0,n*sizeof(double));
	r = new double [n];
	// memset(r,0,n*sizeof(double));
	p = new double [n];
	b = new double [n];
	if (curr_id == 0)
	{
		for(i = 0; i < n; i++)
		{
			for(j = 0; j < n; j++)
				scanf("%lf",&a[i][j]);
			scanf("%lf",&b[i]);
		}
		
		gettimeofday(&start_time, NULL);
		solve(a, b, r, p);
		gettimeofday(&end_time, NULL);
		// Output x values
		for (i = 0; i < n; i++)
			cout<<x[i]<<endl;
		time_taken = (end_time.tv_sec - start_time.tv_sec);
		time_taken *= 1000.0;
		double micro = (end_time.tv_usec - start_time.tv_usec); 
		micro /= 1000.0;
		time_taken += micro;
		printf("Time taken: %lf milliseconds\n", time_taken);		
	}
	else
		solve2(a, b, r, p);
	MPI_Finalize();
	return 0;
}