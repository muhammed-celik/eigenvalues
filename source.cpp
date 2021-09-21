#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <tuple>
using std::tuple;
using namespace std;

class Matrix {
    int m,n;  // mxn matrix
    vector<vector<double> > my_matrix;  // defining matrix itself as 2-d vector
    
  public:
    Matrix (int,int,double);
	Matrix (const char*);  // constructor for reading matrix from txt file and copying its contents to our private 2-d vector
	~Matrix(){} 
    
	Matrix transpose();  // matrix transpose method
	Matrix operator + (Matrix&);  // matrix addition operator " + "
    Matrix operator - (Matrix&);
	Matrix operator * (Matrix&);  // matrix multiplication operator " * "
	Matrix operator * (double); // scalar multiplication operator " * "
    Matrix operator / (double); // scalar division operator " / "

	double operator () (const int&, const int&);  // returns the element of the matrix when asked in the form matrix(x,y)
    int getRows() const {return m;}
    int getCols() const {return n;}
    void setValues(int,int,double);
	
	
};
// assistant member function definitons
void Matrix::setValues(int a, int b, double c){
	my_matrix[a][b] = c;
}
Matrix Matrix::transpose(){
	Matrix result(n,m,0);
	int i,j;
	for(i=0;i<m;i++){
		for(j=0;j<n;j++){
			result.setValues(i,j,my_matrix[j][i]);
		}
	}
	return result;
}
// constructor definitions
Matrix::Matrix(int rowNum, int colNum, double initialVal){
    int i;
	m = rowNum;
    n = colNum;
	my_matrix.resize(rowNum);
    
    for (i=0;i<m;i++){
    	my_matrix[i].resize(n,initialVal);
	}
}
Matrix::Matrix(const char* filename){
	
	ifstream file;   // file for input
	file.open(filename); // opening file for taking matrix A input. First command line argument goes there
	string temp;
    int dim,i;

    while(getline(file,temp)){ 
		file >> temp;  // cursor is going through the file
		dim++;  // dim counter counts the number of loops so that matrix dimension is calculated
	}
	
	n = dim; 
    m = dim;
    
	file.clear(); // cursor is at the end of the file
	file.seekg(0,ios::beg); // setting cursor to beginning of the file

    double** matrixA = new double*[dim]; // firstly 1-d dynamic array is created
	
	for (int i=0; i<dim; i++) {
		matrixA[i] = new double[dim];  // second dimension of the dynamic array is created
	}
	
	
	// copying entries from txt file to the array
	for(int i=0; i<dim; i++){
		for(int j=0; j<dim; j++){
			file >> matrixA[i][j];  //looping through file and copying matrix entries to our dynamic array
		}
	}
	file.close(); // closing file

	my_matrix.resize(dim);
    
    for (i=0;i<dim;i++){
    	my_matrix[i].resize(dim,0);
	}
	
	for(int i=0; i<dim; i++){
		for(int j=0; j<dim; j++){
			my_matrix[i][j] = matrixA[i][j];  
		}
	}
	
	for (int i= 0; i < dim; i++) { 
		delete[] matrixA[i]; // deleting second dimension of matrix A
	}
	delete[] matrixA; // deleting matrix A
	
}
// operator overloads
double Matrix::operator () (const int & a, const int & b){
	return my_matrix[a][b];
}
Matrix Matrix::operator * (Matrix& S){
	Matrix result(m,S.getCols(),0);
	int i,j,k;
	double temp;
	if (n != S.getRows()){
		cout << "Multiplication Error! Matrix dimensions doesn't match." << endl;
	}
	else{
		for(i=0;i<m;i++){
			for(j=0;j<S.getCols();j++){
				temp = 0;
				for(k=0;k<n;k++){
					temp += my_matrix[i][k] * S(k,j);
				}
				result.setValues(i,j,temp);
			}
		}
	}
	return result;
	
}
Matrix Matrix::operator * (double s){
	int i,j;
	Matrix result(m,n,0);
	
	for(i=0;i<m;i++){
		for(j=0;j<n;j++){
			result.setValues(i,j,my_matrix[i][j] * s);
		}
	}
	
	return result;
} 
Matrix Matrix::operator / (double s){
	int i,j;
	Matrix result(m,n,0);
	
	for(i=0;i<m;i++){
		for(j=0;j<n;j++){
			result.setValues(i,j,my_matrix[i][j] / s);
		}
	}
	
	return result;
}
Matrix Matrix::operator + (Matrix& S){
	int i,j;
	Matrix result(m,n,0);
	if (m != S.getRows() || n != S.getCols()){
		cout << "Addition Error! Dimensions of the two matrices does not match." << endl;
	}
	else{
		for(i=0;i<m;i++){
			for(j=0;j<n;j++){
				result.setValues(i,j,my_matrix[i][j] + S(i,j));
			}
		}
	}
	return result;
} 
Matrix Matrix::operator - (Matrix& S){
	int i,j;
	Matrix result(m,n,0);
	if (m != S.getRows() || n != S.getCols()){
		cout << "Substraction Error! Dimensions of the two matrices does not match." << endl;
	}
	else{
		for(i=0;i<m;i++){
			for(j=0;j<n;j++){
				result.setValues(i,j,(my_matrix[i][j] - S(i,j)));
			}
		}
	}
	return result;
} 

tuple<Matrix, double, double, int> powerIteration(Matrix A,double tolerance){
	
	int i,m,n,location;
	
	m = A.getRows();
	n = A.getCols();
	Matrix old_eig_vec(n,1,1);
	Matrix new_eig_vec(n,1,1);
	double old_eig_val = 0;
	double new_eig_val = 1;
	double inf_norm,nonzero;
	int sign = 1;
	
	

	while (abs(new_eig_val-old_eig_val) > tolerance){
		
		old_eig_val = new_eig_val;
		
		new_eig_vec = A * old_eig_vec;
		
		inf_norm = abs(new_eig_vec(0,0));
		
		
		for(i=0;i<m;i++){
			if (abs(new_eig_vec(i,0)) > inf_norm){
				inf_norm = abs(new_eig_vec(i,0));
			}		
		}
		
		
		new_eig_val = inf_norm;
		new_eig_vec = new_eig_vec / inf_norm;
		
		if (old_eig_vec(0,0) * new_eig_vec(0,0) < 0){
			sign = -1;
		}
		old_eig_vec = new_eig_vec;
		
	}
	if (sign == -1){
		new_eig_val = -1 * new_eig_val;
	}
	for(i=0;i<m;i++){
		if(new_eig_vec(i,0) != 0){
			nonzero = new_eig_vec(i,0);
			location = i;
		}
	}
	
	
	return make_tuple(new_eig_vec,new_eig_val,nonzero,location);
}
Matrix deflation(Matrix A,Matrix eig_vector,double nonzero_val,int loc){
	int i;
	Matrix nonzero_vec(1,eig_vector.getRows(),1);
	for(i=0;i<eig_vector.getRows();i++){
		nonzero_vec.setValues(0,i,A(loc,i));
	}
	Matrix inter = eig_vector * nonzero_vec;
	Matrix inter2 = inter / nonzero_val; 
	Matrix B = A - inter2;
	return B;
}

int main(){
	
	Matrix sample("A.txt");
	Matrix eigVector1(sample.getRows(),1,1);
	Matrix eigVector2(sample.getRows(),1,1);
	double precision = 1e-6;
	double eigenvalue1,eigenvalue2,nonzero;
	
	
	
	int i,j,location;
	
	for(i=0;i<sample.getRows();i++){
		for(j=0;j<sample.getCols();j++){
			cout << sample(i,j) << " ";
		}
		cout << endl;
	}
	
	tie(eigVector1, eigenvalue1,nonzero,location) = powerIteration(sample,precision);
	
	
	
	Matrix B = deflation(sample,eigVector1,nonzero,location);
	tie(eigVector2,eigenvalue2,nonzero,location) = powerIteration(B,precision);
	

	ofstream myfile; 
  	myfile.open ("output.txt");  // creating a file named "x.txt" to write the solutions
  	
	myfile << "Eigenvalue#1: " << eigenvalue1 << endl;
	for(i=0;i<sample.getRows();i++){
		myfile << eigVector1(i,0) << endl; // writing the solutions
	}
	myfile << "Eigenvalue#2: " << eigenvalue2 << endl; 
	myfile.close();
	
	return 0;
}