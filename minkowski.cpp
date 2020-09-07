#include <iostream>
#include <vector>
#include <initializer_list>

using namespace std;

// c++11 or later required//

/*This program allows defining vectors, their duals (one-forms) and supports 
basic operations in non-Euclidean Minkowski space in arbitrary dimensions using 
Cartesian coordinates. This inlcudes summation, multiplication, raising/lowering
operations and calculating the scalar product.

Convention: metric signature (-,+,+,+,...)
Vectors are labelled by additional label +1
one-forms (covectors) are lebelled by -1
 */



template<class T> T** a_matrix(int rows, int cols, T **mat)
{	// allacates a matrix of a type T //
    mat = new T*[rows];
    for (int i = 0; i < rows; i++)
        mat[i] = new T[cols];
    return mat;
}

double KroneckerDelta(int i,int j) { 
	if (i==j)
		return 1;
	else
		return 0;
 }//Kronecker delta function//


double** metric = NULL; // covariant Minkowski metric (matrix) //
int dimension; //spacetime dimension//

class Minkowski {
public:
	Minkowski(){}
	Minkowski(int dim);
	~Minkowski() {}
	
	int GetDimension() const;
	virtual void Print() const; 
};



class Vector:public Minkowski {
private:
	vector<double> vec;
	int covariance; // positive for vectors, zero or negative //
	//for covectors (1-forms) //
public:
	Vector();
	Vector(initializer_list<double> init,int cov) : vec(init),
	covariance( (cov<=0)? -1:1 ) {}//covariant or contravariant vector//
	Vector(int cov) : covariance( (cov<=0)? -1:1 ) {} // vector without specifying coordinates //
	Vector(const Vector & rhs);
	~Vector() {}

	void Print() const; 
	double GetCoord(int crd) const;
	double operator() (int crd) {return vec[crd];} // convenient way to get the coordinate //
	int GetCov() const;

	void SetVector();
	void SetVectorCoord(int crd,double val);
	void SetVectorCov(int cov);

//raising and lowering labels://
	void lower();
	void raise();
	//the norm//
	double norm();

	// Operators //

	friend Vector operator+ (const Vector& v1, const Vector& v2);
	friend Vector operator- (const Vector& v1, const Vector& v2);

	friend Vector operator* (double x, const Vector& v);
	friend Vector operator* (const Vector& v, double x);
	friend double operator* (const Vector& v1, const Vector& v2); //scalar product//

};



int main()
{

	//class Minkowski - spacetime dimension and the metrix (matrix) //
	cout << "\nOne starts defining spacetime dimension by introducing a member of Minkowski class.\n";

	Minkowski *minkowski = new Minkowski(4); //dimension d=4//
	
	cout << "\nThe metric reads:\n";
	minkowski->Print(); // d=4 metric //

	cout << "\nDefining vectors requires specifying the covariance label:\n";
	cout << "Convention: +1 for (contravariant) vectors and -1 for (covariant) one-form.\n";
	cout << "\nLet u be a contravariant vector = ";
	Vector u({2,2,0,0},1); //contravariant vector//
	u.Print();
	cout << "corresponds the moment in time t = " << u(0) << ".\n";
	cout << "This is a zero-norm vactor: norm = u*u = " << u.norm() << ".\n";
	cout << "\nLet v be a covariant one-form, v = ";
  Vector v({-1,1,0,2},-1); //covariant vector//
	v.Print();
	cout << "Rising labels gives\n";
	v.raise();
	v.Print();
	cout << "\nLet w be a combination: w = 2 u+ 3 v. One finds w = ";
  Vector w(1);  //contravariant vector without specifying coordinates // 
  
	w=2*u+3*v; 
	w.Print();
	cout << "\nThis operation is defined to be only possible for vectors of the same covariance\n";
	cout << "(two contravariant vectors or two one-forms).\n";
	cout << "\nThe scalar product u*w reads: \n" << v*w << endl; 


	
	

	return 0;
}

	Minkowski::Minkowski(int dim) {
		metric = a_matrix(dim,dim,metric);
		dimension = dim;
		for(int i=0;i<dim;i++)
		{
			for(int j=0;j<dim;j++)
			{
				if(i==0&&j==0)
					metric[i][j]=-1;
				else
					metric[i][j]=KroneckerDelta(i,j);
			}
		}
	}

	int Minkowski::GetDimension() const {return dimension;}

void Minkowski::Print() const {
	for(int i=0;i<dimension;i++)
	{
		for(int j=0;j<dimension;j++)
		{
			if(j==0&&i!=0)
				cout << " ";
			cout <<metric[i][j] << " ";
		}
		cout << endl;
	}
}

	Vector::Vector() {
		for(int i=0;i<dimension;i++)
			vec.push_back(0); // zero vector //
		covariance=1; // contravariant by default //
	} 

	Vector::Vector(const Vector& rhs)
	{
		for(int i=0;i<dimension;i++)
			vec[i]=rhs.vec[i];
		covariance = rhs.covariance;
	}

	void Vector::Print() const
	{
		cout << "( ";
		for(int i=0;i<dimension;i++)
			cout << vec[i] << " ";
		cout << ")" << endl;
	}

	void Vector::SetVector() 
	{ // set vector by coordinates //
		vec.clear();
		double tmp;
		int cov;
		for(int i=0;i<dimension;i++)
		{
			cout << "Coordinate " << i+1 << " of " << dimension << ": ";
			cin >> tmp;
			vec.push_back(tmp);
		}
		cout << "Covariance: positive for vectors, negative or zero for covectors:\n";
		cin >> cov;
		covariance = (  (cov<=0)? -1:1  );
	
	}

	void Vector::SetVectorCoord(int crd,double val)
	{
		vec[crd] = val;
	}

	void Vector::SetVectorCov(int cov)
	{
		covariance = (  (cov<=0)? -1:1  );
	}

	double Vector::GetCoord(int crd) const {return vec[crd];}
	int Vector::GetCov() const {return covariance; }

	double Vector::norm() { //norm of the vector//
		double x = 0;
		for(int i=0;i<dimension;i++)
				for(int j=0;j<dimension;j++)
				 	x += metric[i][j] * vec[i] * vec[j];
		return x;
	}

	void Vector::raise() //making contravariant//
	{
		if(covariance>0)
			cout << "The vector is already contravariant!\n";
		else
		{
			for(int i=0;i<dimension;i++)
			{
				double tmp=0;
				for(int j=0;j<dimension;j++)
				{
					tmp+=metric[i][j]*vec[j]; //summation//
				}
				vec[i]=tmp;
			}
			covariance = 1;
		}
	}

	void Vector::lower() //making covariant//
	{
		if(covariance<=0)
			cout << "The vector is already covariant!\n";
		else
		{
			for(int i=0;i<dimension;i++)
			{
				double tmp=0;
				for(int j=0;j<dimension;j++)
				{
					tmp+=metric[i][j]*vec[j]; //summation//
				}
				vec[i]=tmp;
			}
			covariance = -1;
		}
	}

	//Operators//


	Vector operator + (const Vector &v1,const Vector &v2) //sum//
	{
		Vector vr;
		if(  (v1.covariance)*(v2.covariance) < 0 )
			{
				cout << "Cannot combine covariant with contravariant!\n";
				vr.Print();
				return vr;
			}
		else
			{
				for(int i=0;i<dimension;i++)
					vr.vec[i] = v1.vec[i]+v2.vec[i];
				vr.covariance = v1.covariance;
				return vr;
			}
	}

	Vector operator - (const Vector &v1,const Vector &v2) //sum//
	{
		Vector vr;
		if(  (v1.covariance)*(v2.covariance) < 0 )
			{
				cout << "Cannot combine covariant with contravariant!\n";
				vr.Print();
				return vr;
			}
		else
			{
				for(int i=0;i<dimension;i++)
					vr.vec[i] = v1.vec[i]-v2.vec[i];
				vr.covariance = v1.covariance;
				return vr;
			}
	}


	Vector operator* (double x, const Vector& v)
	{
		Vector vr;
		for(int i=0;i<dimension;i++)
			vr.vec[i] = v.vec[i] * x;
		vr.covariance=(v.covariance);
		return vr;
	}

	Vector operator* (const Vector& v, double x)
	{
		Vector vr;
		for(int i=0;i<dimension;i++)
			vr.vec[i] = v.vec[i] * x;
		vr.covariance=(v.covariance);
		return vr;
	}	

	double operator* (const Vector& v1, const Vector& v2)
	{ //scalar product: using metric if necessarily //
		double x=0;
		if(v1.covariance * v2.covariance < 0)
		{
			for(int i=0;i<dimension;i++)
				x += v1.vec[i] * v2.vec[i];
		}
		else
		{
			for(int i=0;i<dimension;i++)
				for(int j=0;j<dimension;j++)
				 	x += metric[i][j] * v1.vec[i] * v2.vec[j];
		}

		return x;
	}
