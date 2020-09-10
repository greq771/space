#include <iostream>
#include <vector>
#include <cstring>
#include <iomanip>
using namespace std;

//square matrix//

template<class T> using SquareMatrix = vector< vector<T> >;

template<class T>
void Initialize( SquareMatrix<T> &vec,int n)
{
	int i;
	vector<T> vi;
	for(i=0;i<n;i++)
		vi.push_back(0);
	for(i=0;i<n;i++)
		vec.push_back(vi);
	vi.clear();
}


//auxiliary functions//
int mod_numer(int k,int n) 
{
	if(k <= n)
		return n-1;
	else
		return n;
}

int minus_number(int n, int m)
{
	if ((n+m)%2==0)
		return 1;
	else
		return -1;
}

//main functions//

template<class T>
void Minor( int k, int m, SquareMatrix<T>  &mat, SquareMatrix<T>  &mino)
{
	int len = mat.size();
	for(int i=0;i<len;i++)
	{
		for(int j=0;j<len;j++)
		{
			if(i!=k && j!=m)
				mino[ mod_numer(k,i) ][ mod_numer(m,j) ] = mat[i][j];
			else
				continue;
		}
	}
}

template<class T>
void Transpose( SquareMatrix<T>  &mat, SquareMatrix<T>  &trans)
{
	int len = mat.size();

	for(int i=0;i<len;i++)
		for(int j=0;j<len;j++)
			trans[i][j] = mat[j][i];
}

template<class T>
T Determinant(SquareMatrix<T>  &mat)
{
	int len = mat.size();

	if(len<=2)
	{
		if(len==2)
			return(mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0]);
		else
			return mat[0][0];
	}
	else
	{	
		SquareMatrix<T> min;
		Initialize<T>(min,len-1);
		T sum=0;
		for(int i=0;i<len;i++)
		{
			Minor(i,0,mat,min);
			sum += minus_number(i,0)*mat[i][0]*Determinant( min );
		}
		return sum;
	}
}

template<class T>
void Inverse( SquareMatrix<T>  &mat, SquareMatrix<T>  &inv)
{

	T det;
	det = Determinant(mat);
	if(det==0)
		cout << "Non-invertible matrix\n";
	else
	{
		int len = mat.size();
		Initialize<T>(inv,len);

		SquareMatrix<T> min;
		Initialize<T>(min,len-1);

		SquareMatrix<T> adj;
		Initialize<T>(adj,len);

		T sum=0;
		int i,j;
		for(i=0;i<len;i++)
		{
			for(j=0;j<len;j++)
			{
				Minor(i,j,mat,min);
				adj[i][j] = minus_number(i,j)*Determinant(min);
			}
		}
		Transpose(adj,inv);

		for(i=0;i<len;i++)
			for(j=0;j<len;j++)
				inv[i][j] /= det;

	}
}
//printing//

template<class T>
int c_length(T num)
{
	string tmp;
	tmp = to_string(num);
	return( tmp.length() );
}

template<class T>
void Print(SquareMatrix<T>  &mat)
{//Printing matrices//
	int len = mat.size();
	T tmp;
	int i,j;
	int max=0;
	for(i=0;i<len;i++)
	{
		for(j=0;j<len;j++)
		{
			if( c_length( mat[i][j] ) >= max )
				max = c_length( mat[i][j] );
		}
	}
	
	for(i=0;i<len;i++)
	{
		for(j=0;j<len;j++)
		{ 
			cout << right << setw(max+3) << ((mat[i][j]==0)? 0 : mat[i][j]) << " ";
		}
		cout << endl;

	}
}



int main()
{
	cout << "\nLet m be a square matrix:\n";

	SquareMatrix<double> m;
	m = {{-1,1,2},{-2,2,12},{3,21,3}};

	/* 
	Alternatively: 

	Initialize<double>(m,3);
	m[0][0]=-1, m[0][1]=1, m[0][2]=2,
	m[1][0] = -2,...
	*/

	Print(m);

	cout << "Determinant = " << Determinant(m) << endl;

	cout << "Inverse is " << endl;
	
	SquareMatrix<double> im; //inverse of m//
	Inverse(m,im);
	Print(im);





	return 0;
}
