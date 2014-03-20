#include<glpk.h>
#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<climits>
#include<sstream>

#include<algorithm>
#include<string>
#include<vector>
#include<set>
#include<map>
#include<utility>
#include<stack>
#include<queue>
#include<deque>
#include<list>
#include<bitset>

using namespace std;

typedef vector<int> vi; 
typedef vector<vi> vvi; 
typedef vector<string> vs; 
typedef pair<int,int> ii;
typedef long long int LLI;
typedef unsigned long long int ULLI;

#define sz(a)                        int((a).size()) 
#define pb                           push_back 
#define mp                           make_pair
#define F                            first
#define S                            second
#define present(c,x)                 ((c).find(x) != (c).end()) 
#define cpresent(c,x)                (find(all(c),x) != (c).end())
#define tr(c,i)                      for(typeof((c).begin()) i = (c).begin(); i != (c).end(); i++)
#define all(c)                       (c).begin(),(c).end()
#define si(n)                        scanf("%d",&n)
#define sl(n)                        scanf("%lld",&n)
#define sf(n)                        scanf("%f",&n)
#define sd(n)                        scanf("%lf",&n)
#define ss(n)                        scanf("%s",n)

#define abs(x)                       ((x)<0?-(x):(x))
#define fill(a,v)                    memset((a),(v),sizeof (a))
#define INF                          INT_MAX
#define LINF                         (long long)1e18
#define EPS                          1e-9
#define MAX 100010
#define MAX2 110

struct point
{
	int i;
	int j;
	int k;

	point()
	{
		i = j = k = 0;
	}

	point(int ii, int jj, int kk)
	{
		i = ii;
		j = jj;
		k = kk;
	}
};

point inp[MAX];

int ia[MAX];
int ja[MAX];
int X[MAX2][MAX2];
int n;
int already_filled;
double ar[MAX];
glp_prob *lp;

int getpos(int i, int j, int k)
{
	--i, --j;
	return (i*n*n + j*n + k);
}


void set_rows()
{
	int constraint_cnt = 4*n*n + already_filled;
	glp_add_rows(lp,constraint_cnt);

	for(int i=1; i<=constraint_cnt; i++)
	{
		glp_set_row_bnds(lp,i,GLP_FX,1,1);
	}

}

void set_cols()
{
	int No_of_variables = n*n*n;
	glp_add_cols(lp,No_of_variables);

	for(int i=1; i<=No_of_variables; i++)
	{
		glp_set_col_kind(lp, i, GLP_BV);
		glp_set_col_bnds(lp, i, GLP_DB, 0, 1);
		glp_set_obj_coef(lp, i, 0);

	}
}

void set_matrix()
{
	int q,m;

	//For the first set of constraints - only one k in each column
	int c=1, r=1;
	for(int j=1; j<=n; j++)
	{
		for(int k=1; k<=n; k++)
		{
			for(int i=1; i<=n; i++)
			{
				ia[c] = r;
				ja[c] = getpos(i,j,k);
				ar[c] = 1;
				c++;
			}
			r++;
		}
	}

	// For the second set of constraints - only one k in each row
	for(int i=1; i<=n; i++)
	{
		for(int k=1; k<=n; k++)
		{
			for(int j=1; j<=n; j++)
			{
				ia[c] = r;
				ja[c] = getpos(i,j,k);
				ar[c] = 1;
				c++;
			}
			r++;
		}
	}

	//For the third set of constraints - only one k in each block
	m = sqrt(n);
	for(int k=1; k<=n; k++)
	{
		for(int p=1; p<=m; p++)
		{
			for(q=1; q<=m; q++)
			{
				for(int i=m*p-m+1; i<=m*p; i++)
				{
					for(int j=m*q-m+1; j<=m*q; j++)
					{
						ia[c] = r;
						ja[c] = getpos(i,j,k);
						ar[c] = 1;
						c++;
					}
				}
				r++;
			}
		}
	}


	//For the fourth set of constraints - every element of the grid must be filled

	for(int i=1; i<=n; i++)
	{
		for(int j=1; j<=n; j++)
		{
			for(int k=1; k<=n; k++)
			{
				ia[c] = r;
				ja[c] = getpos(i,j,k);
				ar[c] = 1;
				c++;
			}
			r++;
		}
	}

	//Finally, the already inp elements are added as constraints

	for(int i=0; i<already_filled; i++)
	{

		ia[c] = r;
		ja[c] = getpos(inp[i].i,inp[i].j,inp[i].k);
		ar[c] = 1;
		c++;
		r++;       
	}
	c--;
	printf("Number of constraints: %d\n",c);

	//Solving
	glp_load_matrix(lp, c, ia, ja, ar);
}

int main(int argc, const char *argv[])
{

	lp = glp_create_prob();
	glp_set_prob_name(lp,"Sudoku solved");
	glp_set_obj_dir(lp, GLP_MIN);

	printf("Enter N = ");
	scanf("%d",&n);

	printf("Enter the number of given elements: ");
	scanf("%d",&already_filled);

	if(already_filled!=0)
	{
		printf("Enter i j k suct that A[i][j]=k:\n");
		for(int i=0; i<already_filled; i++)
			scanf("%d%d%d",&inp[i].i,&inp[i].j,&inp[i].k);
	}

	//Solving
	set_rows();
	set_cols();
	set_matrix();

	glp_simplex(lp, NULL);
	glp_intopt(lp,NULL); 


	// Display results
	int No_of_variables = n*n*n;
	for(int idx=1; idx<=No_of_variables; idx++)
	{
		int x = glp_mip_col_val(lp, idx);
		int tmp=idx;
		if(x==1)
		{
			int i,j,k;
			--idx;
			i = idx/(n*n)+1;
			idx = idx%(n*n);
			j = idx/n + 1;
			idx = idx%n;
			k = idx+1;

			X[i][j] = k;
		}
		idx = tmp;

	}

	for(int i=1; i<=n; i++)
	{
		for(int j=1; j<=n; j++)
			printf("%d ",X[i][j]);
		printf("\n");
	}

	/* housekeeping */
	glp_delete_prob(lp);
	glp_free_env();
	return 0;
}
