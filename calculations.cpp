#include <random>
#include <NXOpen/Point.hxx>
#include <NXOpen/PointCollection.hxx>
#include <chrono>
#include<iostream> 
#include<cmath>
#include <vector>
#include <algorithm>
#include <sstream>
#define _USE_MATH_DEFINES
#include <numeric> 
#include <algorithm>
#include <functional>
#include <limits>
#include <cfloat>
#include <math.h>
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif
#ifndef M_PIl
#define M_PIl (3.14159265358979323846264338327950288)
#endif
class CombinationsIndexArray {
	std::vector<int> index_array;
	int last_index;
public:
	CombinationsIndexArray(int number_of_things_to_choose_from, int number_of_things_to_choose_in_one_combination) {
		last_index = number_of_things_to_choose_from - 1;
		for (int i = 0; i < number_of_things_to_choose_in_one_combination; i++) {
			index_array.push_back(i);
		}
	}
	int operator[](int i) {
		return index_array[i];
	}
	int size() {
		return index_array.size();
	}
	bool advance() {
		int i = index_array.size() - 1;
		if (index_array[i] < last_index) {
			index_array[i]++;
			return true;
		}
		else {
			while (i > 0 && index_array[i - 1] == index_array[i] - 1) {
				i--;
			}
			if (i == 0) {
				return false;
			}
			else {
				index_array[i - 1]++;
				while (i < index_array.size()) {
					index_array[i] = index_array[i - 1] + 1;
					i++;
				}
				return true;
			}
		}
	}
};
class calculation {
public:
	double calculation::randomInRange(double lowerLimit, double upperLimit) {
		double f = (double)rand() / RAND_MAX;
		return lowerLimit + f * (upperLimit - lowerLimit);
	}
	double calculation::circleArea(double diameter) {
		return ((diameter*diameter) / 4)*M_PIl;
	}
	bool collinear(int x1, int y1, int x2, int y2, int x3, int y3)
	{
		// Calculation the area of  
		// triangle. We have skipped  
		// multiplication with 0.5  
		// to avoid floating point  
		// computations  
		double a = x1 * (y2 - y3) +
			x2 * (y3 - y1) +
			x3 * (y1 - y2);
		if (a == 0)
			return true;
		else
			return false;
	}
	double distance(double x1, double y1, double x2, double y2) {
		return sqrt(((x1 - x2)*(x1 - x2)) + ((y1 - y2)*(y1 - y2)));
	}
	int compute_rank(std::vector<std::vector<double>> A) {
		const double EPS = 1E-9;
		int n = A.size();
		int m = A[0].size();
		int rank = 0;
		std::vector<bool> row_selected(n, false);
		for (int i = 0; i < m; ++i) {
			int j;
			for (j = 0; j < n; ++j) {
				if (!row_selected[j] && abs(A[j][i]) > EPS)
					break;
			}
			if (j != n) {
				++rank;
				row_selected[j] = true;
				for (int p = i + 1; p < m; ++p)
					A[j][p] /= A[j][i];
				for (int k = 0; k < n; ++k) {
					if (k != j && abs(A[k][i]) > EPS) {
						for (int p = i + 1; p < m; ++p)
							A[k][p] -= A[j][p] * A[k][i];
					}
				}
			}
		}
		return rank;
	}
	double getDeterminant(const std::vector<std::vector<double>> vect) {
		if(vect.size() != vect[0].size()) {
			throw std::runtime_error("Matrix is not quadratic");
		} 
		int dimension = vect.size();
		if(dimension == 0) {
			return 1;
		}
		if(dimension == 1) {
			return vect[0][0];
		}
		//Formula for 2x2-matrix
		if(dimension == 2) {
			return vect[0][0] * vect[1][1] - vect[0][1] * vect[1][0];
		}
		double result = 0;
		int sign = 1;
		for(int i = 0; i < dimension; i++) {
			//Submatrix
			std::vector<std::vector<double>> subVect(dimension - 1, std::vector<double> (dimension - 1));
			for(int m = 1; m < dimension; m++) {
				int z = 0;
				for(int n = 0; n < dimension; n++) {
					if(n != i) {
						subVect[m-1][z] = vect[m][n];
						z++;
					}
				}
			}
			//recursive call
			result = result + sign * vect[0][i] * getDeterminant(subVect);
			sign = -sign;
		}
		return result;
	}
	std::vector<std::vector<double>> getTranspose(const std::vector<std::vector<double>> matrix1) {
		//Transpose-matrix: height = width(matrix), width = height(matrix)
		std::vector<std::vector<double>> solution(matrix1[0].size(), std::vector<double> (matrix1.size()));
		//Filling solution-matrix
		for(size_t i = 0; i < matrix1.size(); i++) {
			for(size_t j = 0; j < matrix1[0].size(); j++) {
				solution[j][i] = matrix1[i][j];
			}
		}
		return solution;
	}
	std::vector<std::vector<double>> getCofactor(const std::vector<std::vector<double>> vect) {
		if(vect.size() != vect[0].size()) {
			throw std::runtime_error("Matrix is not quadratic");
		} 
		std::vector<std::vector<double>> solution(vect.size(), std::vector<double> (vect.size()));
		std::vector<std::vector<double>> subVect(vect.size() - 1, std::vector<double> (vect.size() - 1));
		for(std::size_t i = 0; i < vect.size(); i++) {
			for(std::size_t j = 0; j < vect[0].size(); j++) {
				int p = 0;
				for(size_t x = 0; x < vect.size(); x++) {
					if(x == i) {
						continue;
					}
					int q = 0;
					for(size_t y = 0; y < vect.size(); y++) {
						if(y == j) {
							continue;
						}
						subVect[p][q] = vect[x][y];
						q++;
					}
					p++;
				}
				solution[i][j] = pow(-1, i + j) * getDeterminant(subVect);
			}
		}
		return solution;
	}
	std::vector<std::vector<double>> getInverse(const std::vector<std::vector<double>> vect) {
		if(getDeterminant(vect) == 0) {
			throw std::runtime_error("Determinant is 0");
		} 
		double d = 1.0/getDeterminant(vect);
		std::vector<std::vector<double>> solution(vect.size(), std::vector<double> (vect.size()));
		for(size_t i = 0; i < vect.size(); i++) {
			for(size_t j = 0; j < vect.size(); j++) {
				solution[i][j] = vect[i][j]; 
			}
		}
		solution = getTranspose(getCofactor(solution));
		for(size_t i = 0; i < vect.size(); i++) {
			for(size_t j = 0; j < vect.size(); j++) {
				solution[i][j] *= d;
			}
		}
		return solution;
	}
	std::vector<std::vector<double>> multiply(std::vector<std::vector<double>> A,std::vector<std::vector<double>> B)
	{
		//	If A is an m × n matrix and B is an n × p matrix,
		//then the matrix product C = AB (denoted without multiplication signs or dots) is defined to be the m × p matrix
		int m_,n_,p_;
		m_=A.size();
		n_=B.size();	
		p_=B[0].size();
		std::vector<std::vector<double>>  C(m_);
		for(int m=0;m<m_;m++)
		{
			for(int p=0;p<p_;p++)
			{
				double c=0;
				for(int n=0;n<n_;n++)
				{
					c=c+ (A[m][n]*B[n][p]); 
				}
				C[m].push_back(c);
			}
		}
		return C;
	}  
	std::vector<std::vector<double>> negative(const std::vector<std::vector<double>> vect) {
		std::vector<std::vector<double>> result(vect.size());
		for (int i=0;i<vect.size();i++)
		{
			for(int j=0;j<vect[0].size();j++)
			{
				result[i].push_back( -1*	vect[i][j]);
			}
		}
		return result;
	}
	double standardDeviation(std::vector<double> v) {
		double sum = std::accumulate(v.begin(), v.end(), 0.0);
		double mean = sum / v.size();
		std::vector<double> diff(v.size());
		std::transform(v.begin(), v.end(), diff.begin(),
			std::bind2nd(std::minus<double>(), mean));
		double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
		double stdev = std::sqrt(sq_sum / v.size());
		return stdev;
	}
	double rand_gen() {
		// return a uniformly distributed random value
		return ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
	}
	double normalGenerator(double mean, double stdDev)
	{
		// std::default_random_engine generator(time(0)); // We initialize a generator
		//std::normal_distribution<double> distribution(mean,stdDev); // And a distribution
		// return distribution(generator); // distribution operator() always take the same generator with the same state
		/* Create random engine with the help of seed */
		unsigned seed =std:: chrono::steady_clock::now().time_since_epoch().count(); 
		std::default_random_engine e (seed); 
		/* declaring normal distribution object 'distN' and initializing its mean and standard deviation fields. */
		/* Mean and standard deviation are distribution parameters of Normal distribution. Here, we have used mean=5, and standard deviation=2. You can take mean and standard deviation as per your choice */
		std::normal_distribution<double> distN(mean,stdDev); 
		return distN(e);
	}
	double randomNormalDistribution(int lowerLimit,int upperLimit,double range,double mean,double segma)
	{
		double result;
		double x;
		for(int i=0;i<100;i++)
		{
			x = normalGenerator(0,segma);
			result=(x*range)+mean;
			if(result<lowerLimit|| result>upperLimit)
			{ 
				result=mean; continue;
			}else
			{
				break;
			}
		}
		return result; 
	}
	std::string matrixToString(	std::vector<std::vector< double>> &matrix)
	{
		std::string data;
		for(int i=0;i<matrix.size();i++)
		{
			for(int j=0;j<matrix[0].size();j++)
			{
				if(j!=0)
				{
					data=data +",";
				}
				data=data + doubleToString(matrix[i][j]);
			}
			data=data +" \n";
		}
		data=data +" \n\n\n"; 
		return data;
	}
	std::string matrixToString(	int ** matrix,int v,int h)
	{
		std::string data;
		for (int i = 0; i < v; i++) {
			for (int j = 0; j < h; j++) 
			{
				if(j!=0)
				{
					data=data +",";
				}
				data=data + std::to_string(matrix[j][i]);
			}
			data=data +" \n";
		}
		data=data +" \n\n\n"; 
		return data;
	}
	std::string matrixToString(	std::vector< double> &matrix)
	{
		std::string data;
		for(int j=0;j<matrix.size();j++)
		{
			data=data +"\n"+ doubleToString(matrix[j]);
		}
		data=data +" \n\n\n"; 
		return data;
	}
	std::string matrixToOneLine(	std::vector< double> &matrix)
	{
		std::string data;
		for(int j=0;j<matrix.size();j++)
		{
			data=data +" , "+ doubleToString(matrix[j]);
		}
		return data;
	}
	std::string matrixToOneLine(	std::vector< int> &matrix)
	{
		std::string data;
		for(int j=0;j<matrix.size();j++)
		{
			data=data +" , "+ std::to_string(matrix[j]);
		}
		return data;
	}
	std::string matrixToOneLine(	double matrix[],int size)
	{
		std::string data;
		for(int j=0;j<size;j++)
		{
			data=data +" , "+ std::to_string(matrix[j]);
		}
		return data;
	}
	std::string matrixToString(	std::vector< std::string> &matrix)
	{
		std::string data;
		for(int j=0;j<matrix.size();j++)
		{
			data=data +"\n"+ matrix[j];
		}
		data=data +" \n\n\n"; 
		return data;
	}
	std::string matrixToString(	std::vector< int> &matrix)
	{
		std::string data;
		for(int j=0;j<matrix.size();j++)
		{
			data=data +"\n"+ std::to_string(matrix[j]);
		}
		data=data +" \n\n\n"; 
		return data;
	}
	std::string matrixToString(	 int matrix[], int arraySize)
	{
		std::string data;
		for(int j=0;j<arraySize;j++)
		{
			data=data +"\n"+ std::to_string(matrix[j]);
		}
		data=data +" \n\n\n"; 
		return data;
	}
	std::string matrixToString(	 double matrix[], int arraySize)
	{
		std::string data;
		for(int j=0;j<arraySize;j++)
		{
			data=data +"\n"+ std::to_string(matrix[j]);
		}
		data=data +" \n\n\n"; 
		return data;
	}
	double randomOfCustomProbabilty1(double a,double c,double z,double l)
	{
		double e=2.71828182845904523536;
		double result=0.0;
		//generate random in positive range
		double x=0.0;
		double temp_range=(		(z*c)		+		(	(z/l)		*		(1	-	pow(e,	l*	(c-a)	)	)	)			)	;
		double	y=randMToN(0,temp_range);
		if(y>0 && y<=z*c)
		{
			x=y/z;
		}else if(y>z*c	 &&		y<temp_range	)
		{
			//log her calculate for base euler (e) which means ln
			x=c - (		(1/l)		 *		 (	 log( 1	-	(	(	l*	(	y		-	(z*c)	)	)	/z		)	) ));
		}
		//	x=c - (		(1/l)		 *		 (	 log( 1	-	(	(	l*	(	y		-	(z*c)	)	)	/z		)	) ));
		if(rand() % 2==0)
		{
			result=x;
		}else
		{
			result=-x;
		}
		return result;
	}
	double randomOfCustomProbabilty2(double a,double c)
	{
		double z=100;
		if(c-a==0)
		{
			a+=0.0000001;
		}
		double l=(log(1/z))/	(c-a);
		return randomOfCustomProbabilty1(a,c,z,l);
	}
	double randMToN(double M, double N)
	{
		return M + (rand() / ( RAND_MAX / (N-M) ) ) ;  
	}
	void smoothTheDervative(std::vector<		std::vector <double>>&pathPointsUnitNorm_x,std::vector<		std::vector <double>>&pathPointsUnitNorm_y)
	{
		//dervative smoothing
		std::vector<double>temp_newSmoothedValues_x;
		std::vector<double>temp_newSmoothedValues_y;
		for(int i=0;i<pathPointsUnitNorm_x[0].size();i++)
		{
			double previusValue_x,nextValue_x;
			double previusValue_y,nextValue_y;
			if(i==0)
			{
				previusValue_x=pathPointsUnitNorm_x[0][pathPointsUnitNorm_x[0].size()-1];
				previusValue_y=pathPointsUnitNorm_y[0][pathPointsUnitNorm_y[0].size()-1];
			}else
			{
				previusValue_x=pathPointsUnitNorm_x[0][i-1];
				previusValue_y=pathPointsUnitNorm_y[0][i-1];
			}
			if(i==pathPointsUnitNorm_x[0].size()-1)
			{
				nextValue_x=pathPointsUnitNorm_x[0][0];
				nextValue_y=pathPointsUnitNorm_y[0][0];
			}else
			{
				nextValue_x=pathPointsUnitNorm_x[0][i+1];
				nextValue_y=pathPointsUnitNorm_y[0][i+1];
			}
			double temp_x=(previusValue_x*0.25)+(pathPointsUnitNorm_x[0][i]*0.5)+(nextValue_x*0.25);
			double temp_y=(previusValue_y*0.25)+(pathPointsUnitNorm_y[0][i]*0.5)+(nextValue_y*0.25);
			temp_newSmoothedValues_x.push_back(temp_x/sqrt((temp_x*temp_x)+(temp_y*temp_y)));
			temp_newSmoothedValues_y.push_back(temp_y/sqrt((temp_x*temp_x)+(temp_y*temp_y)));
		}
		pathPointsUnitNorm_x[0]=temp_newSmoothedValues_x;
		pathPointsUnitNorm_y[0]=temp_newSmoothedValues_y;
	}
	void thetaProportionalSmooth(std::vector<		std::vector <double>>&pathPointsUnitNorm_x,std::vector<		std::vector <double>>&pathPointsUnitNorm_y,double theta)
	{
		double e=2.71828182845904523536;
		//dervative smoothing
		std::vector<double>temp_newSmoothedValues_x;
		std::vector<double>temp_newSmoothedValues_y;
		int numberOfSmoothingPairs=50 *pow(e,-theta);
		for(int i=0;i<pathPointsUnitNorm_x[0].size();i++)
		{
			double sum_x=0,sum_y=0;
			for(int j=1;j<=numberOfSmoothingPairs;j++)
			{
				//next value
				sum_x+=pathPointsUnitNorm_x[0][(i+j)%pathPointsUnitNorm_x[0].size()];
				sum_y+=pathPointsUnitNorm_y[0][(i+j)%pathPointsUnitNorm_y[0].size()];
				//previous value
				int position_temp=i-j;
				if(position_temp>=0)
				{
					sum_x+=pathPointsUnitNorm_x[0][(position_temp)];
					sum_y+=pathPointsUnitNorm_y[0][(position_temp)];
				}else
				{									
					sum_x+=pathPointsUnitNorm_x[0][pathPointsUnitNorm_x[0].size()+position_temp];
					sum_y+=pathPointsUnitNorm_y[0][pathPointsUnitNorm_y[0].size()+position_temp];
				}
			}
			double temp_x=(sum_x+pathPointsUnitNorm_x[0][i])/(numberOfSmoothingPairs+1);
			double temp_y=(sum_y+pathPointsUnitNorm_y[0][i])/(numberOfSmoothingPairs+1);
			temp_newSmoothedValues_x.push_back(temp_x/sqrt((temp_x*temp_x)+(temp_y*temp_y)));
			temp_newSmoothedValues_y.push_back(temp_y/sqrt((temp_x*temp_x)+(temp_y*temp_y)));
		}
		pathPointsUnitNorm_x[0]=temp_newSmoothedValues_x;
		pathPointsUnitNorm_y[0]=temp_newSmoothedValues_y;
	}
	void thetaProportionalSmooth(	std::vector <double>&values,double alpha)
	{ 
		double e=2.71828182845904523536;
		//dervative smoothing
		std::vector<double>temp_newSmoothedValues;
		int numberOfSmoothingPairs=50 *pow(e,-alpha);
		for(int i=0;i<values.size();i++)
		{
			double sum=0;
			for(int j=1;j<=numberOfSmoothingPairs;j++)
			{
				//next value
				sum+=values[(i+j)%values.size()];
				//previous value
				int position_temp=i-j;
				if(position_temp>=0)
				{
					sum+=values[(position_temp)];
				}else
				{							
					sum+=values[values.size()+position_temp];
				}
			}
			double temp=(sum+values[i])/(numberOfSmoothingPairs+1);
			temp_newSmoothedValues.push_back(temp);
		}
		values=temp_newSmoothedValues;
	}
	void extractTwoSubsets(std::vector<		std::vector <double>>&pathPointsUnitNorm_x,std::vector<		std::vector <double>>&pathPointsUnitNorm_y,std::vector<	int>& mapToSubSetOne,std::vector<	int>& mapToSubSetTwo)
	{
		for(int i=0;i<	pathPointsUnitNorm_x[0].size();i++)
		{
			//map the two subsets
			if(abs(pathPointsUnitNorm_x[0][i]!=0 && abs(pathPointsUnitNorm_y[0][i])/abs(pathPointsUnitNorm_x[0][i])<1))
			{
				//first subset
				mapToSubSetOne.push_back(i);
			}else
			{
				//second subset
				mapToSubSetTwo.push_back(i);
			}
			;
		}
	}
	std::vector<std::vector	<	double>> getNt(	std::vector <double>*pathPointsUnitNorm_x,	std::vector <double>*pathPointsUnitNorm_y,int locatorsPointsLocation[3])
	{
		std::vector<std::vector	<	double>> Nt(3);
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<6;j++)
			{
				if(j==i*2)
				{
					Nt[i].push_back((*pathPointsUnitNorm_x)[locatorsPointsLocation[i]]);
				}else if(j==i*2+1)
				{
					Nt[i].push_back((*pathPointsUnitNorm_y)[locatorsPointsLocation[i]]);
				}else {
					Nt[i].push_back(0);
				}
			}
		}
		return Nt;
	}
	double poolStandardDeviation(std::vector<	double> standardDev,int samplesSize)
	{
		double PooledStandardDeviation;
		//calculate pooled standard deviation
		double x_sum=0;
		for(int z=0;z<standardDev.size();z++)
		{
			x_sum+=(samplesSize-1)*standardDev[z]*standardDev[z];
		}
		return		PooledStandardDeviation=sqrt(x_sum/((samplesSize*standardDev.size())-standardDev.size()));
	}
	void squareAllMatrixMembers(std::vector<std::vector	<	double>> &matrix)
	{
		std::vector<std::vector	<	double>> resultMatrix;
		for(int i=0;i<matrix.size();i++)
		{
			for(int j=0;j<matrix[0].size();j++)
			{
				double theSquare=matrix[i][j]*matrix[i][j];
				matrix[i][j]=theSquare;
			}
		}
	}
	std::vector <double> curvatureCalculation(std::vector <double> FirstDervitive,std::vector <double> SecondDervitive)
	{
		std::vector <double>curvature_k;
		for(int i=0;i<FirstDervitive.size();i++)
		{double numerator,denominator ;
		numerator=SecondDervitive[i];
		denominator=pow( 1+(FirstDervitive[i]*FirstDervitive[i]),3/2);
		curvature_k.push_back(numerator/denominator);
		}
		return curvature_k;
	}
	std::vector <double>	normBased_firstDervative(std::vector <double>UnitNorm_x,	std::vector <double>UnitNorm_y,double ebsilon)
	{	std::vector <double>FirstDervitive;
	for(int i=0;i<UnitNorm_x.size();i++)
	{
		double x,y;
		if(UnitNorm_y[i]==0)
		{
			x=ebsilon; //to avoid singularity
		}else
		{
			x=UnitNorm_y[i];
		}
		y=-UnitNorm_x[i];
		FirstDervitive.push_back(y/x);
	} return FirstDervitive;
	}
	double slope( std::vector<double>& x_points,  std::vector<double>& y_points,double ebsilon,int pivotPoint)
	{
		size_t n = x_points.size();
		double xDif_sm=0,yDif_sum=0;
		for(int i=0; i<n; ++i)
		{
			if(i==pivotPoint)
			{
				continue;
			}
			yDif_sum+=(n-i)*(y_points[i]-y_points[pivotPoint]);
			xDif_sm+=(n-i)*(x_points[i]-x_points[pivotPoint]);
		}
		if(xDif_sm == 0.0)
		{
			xDif_sm=ebsilon;
			xDif_sm=xDif_sm+( (n-1)*x_points[0] )-((n-1)*ebsilon) ;
		}
		return (  yDif_sum / xDif_sm);
	}
	double determinant(std::vector< std::vector<double> > A )
	{	const double SMALL = 1.0E-30;
	int n = A.size();
	double det = 1;
	// Row operations for i = 0, ,,,, n - 2 (n-1 not needed)
	for ( int i = 0; i < n - 1; i++ )
	{
		// Partial pivot: find row r below with largest element in column i
		int r = i;
		double maxA = abs( A[i][i] );
		for ( int k = i + 1; k < n; k++ )
		{
			double val = abs( A[k][i] );
			if ( val > maxA )
			{
				r = k;
				maxA = val;
			}
		}
		if ( r != i )
		{
			for ( int j = i; j < n; j++ ) std::swap( A[i][j], A[r][j] );
			det = -det;
		}
		// Row operations to make upper-triangular
		double pivot = A[i][i];
		if ( abs( pivot ) < SMALL ) return 0.0;              // Singular matrix
		for ( int r = i + 1; r < n; r++ )                    // On lower rows
		{
			double multiple = A[r][i] / pivot;                // Multiple of row i to clear element in ith column
			for ( int j = i; j < n; j++ ) A[r][j] -= multiple * A[i][j];
		}
		det *= pivot;                                        // Determinant is product of diagonal
	}
	det *= A[n-1][n-1];
	return det;
	}
	int newPointGenerator(double k,int mean, int rangeSize)
	{
		double a_min=6;
		double c_min=1;
		int newPoint;
		double	a_max=(rangeSize-1)/2;
		double c_max=a_max/2;
		double	a=a_min+k*(a_max-a_min);
		double c=	c_min+k*k*(c_max-c_min);
		int	temp=randomOfCustomProbabilty2(a,c);
		newPoint=mean+temp;
		if(newPoint<0 || newPoint>=rangeSize)
		{
			newPoint=mean-temp;
		}
		return newPoint;
	}
	double newMetrapolisCriterian(int iterations,int currentIteration,double cost,double minCost,double maxcost,double T_min,double T_max)
	{
		double e=2.71828182845904523536;
		double T;
		//simple cooling function
		double k_cooling=-log(T_min/T_max)/((double)iterations-1);
		T=T_max*pow(e,-(k_cooling*(double)currentIteration));
		double	k_B;
		double		k_B_min;
		double		k_B_max=0.45;
		if(minCost==0)
		{
			k_B_min=0.000006;
		}else
		{
			k_B_min=-(minCost)/(T_min*log(0.2));
		}
		if(maxcost==0)
		{
			k_B=k_B_min;
		}else
		{
			k_B_max=-(maxcost)/(T_max*log(0.8));
			k_B=k_B_min+		(	pow(e,-(T_max-T))	* (k_B_max-k_B_min));
		}
		return  	(double)pow(e,-(cost/(T*k_B)));
	}
	double average(std::vector<double> v)
	{
		auto n = v.size(); 
		float average = 0.0f;
		if ( n != 0) {
			average =std:: accumulate( v.begin(), v.end(), 0.0) / n; 
		}
		return (double)average;
	}
	double enyropy(std::vector <double>vector_values)
	{
		std::vector <double>curveture2_temp=vector_values;
		double total_size=curveture2_temp.size();
		std::vector<	int> sets_size;
		while (!curveture2_temp.empty())
		{
			int temp_size=0;
			double temp_value=curveture2_temp[0];
			//detect phase
			std::vector <double>vector_temp;
			for(int i=0;i<curveture2_temp.size();i++)
			{
				if(curveture2_temp[i]==temp_value)
				{
					temp_size+=1;
				}else
				{
					vector_temp.push_back(curveture2_temp[i]);
				}
			}
			curveture2_temp.clear();
			if(!vector_temp.empty())
			{
				curveture2_temp=vector_temp;
			}
			sets_size.push_back(temp_size);
		}
		double entropy_sum=0;
		for(int i=0;i<sets_size.size();i++)
		{
			double propapilty=(double)sets_size[i]/total_size;
			entropy_sum+=-propapilty*(log(propapilty)/log(2));
		}
		return entropy_sum;
	}
	void linearNormilaization(	std::vector <double>& vec)
	{
		double sum=0;
		for(int i=0;i<vec.size();i++)
		{
			sum+=vec[i];
		}
		for(int i=0;i<vec.size();i++)
		{	
			vec[i]=(double)vec[i]/sum;
		}
	}
	std::string doubleToString(double number)
	{
		std::ostringstream streamObj;
		streamObj << number;
		return streamObj.str();
	}
};