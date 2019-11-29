#include <ctime>
#include <string>
#include <vector>
#include "gauss.h"
//#include "tests.h"
using namespace std;

void Reflecsions(double *a, int n,double Norma);
void ReflMatr(double * U, double * x , int n, int k);
void print_matrix_t(double* m, int size);
void LeftMult(double *a, double *b , double *x,  int n, int k);
void RightMult(double *a, double *b , double *x,  int n, int k);
void EigenVal(double *a, double *eig, int n);
void matr_sub_sdvig(double *a ,double sdvig, int n, int k1, int k2);
double  CheckNorma(double * a , int n);
void RLmult(double *a , int n, int k1, int k2);
void LRdeform(double *a, int n, int k1, int k2);
double Trace(double *a, int n);
void EigenValImpl(double* a, int size, double* eig, int begin, int end, double Norma);
void EigValDim2(double *a , int size, double *eig, int begin, int end);
int EigValDim2(double *a , int size, int begin, int end);
double finds_k(double *a, int size, int begin,int end, double Norma);
void print_matrix_t(double* m, int size,int begin , int end);
	
	
	
static int COMPLEXVAL=0;


int main(int argc, char** argv)
{
	srand(time(0));
	double *a,*ev;
	int n = 0 ;
	string name = "";
	
	
	
	
	if (!((argc == 4) || (argc == 3)) || ((n = stoi(argv[1])) <= 0) || ((EPS = stod(argv[2])) <=0) ) {
	    printf("usage: %s n EPS [name] \n ", argv[0]);
	    return -1;
	}
	
	try{a = new double[n * n];}
	catch(...){
	    printf("Not enough memory a ");
	    return -2;
	}
	
	
	if (argc == 4) {
	    name = argv[3];
	    if (init_matrix_file(a, n, name) < 0) {
                delete [] a;
                return 0;
	 
	 };
	}

	if (argc == 3) {
	    formula_matr(a, n);
	}

	try{ev = new double[n];}
	catch(...){
	    printf("Not enough memory x ");
		delete []a;
	    return -2;
	}

	printf("Init. Matrix \n");
	print_matrix_t(a, n);
	
	double t = clock();
	EigenVal(a,ev,n);
	t = (clock() - t) / CLOCKS_PER_SEC;
	
	cout << "TIME: " << t << endl;

	delete[] a;
	delete[] ev;

    return 0;
}

void EigenVal(double *a, double *eig, int n){
	double FirstCheck = Trace(a,n), SecondCheck=CheckNorma(a,n);
	LOG(EPS);
	double Norma = norma_matr(a,n);
    Norma = min(Norma, 1.0);
	if(n>=3){
		Reflecsions(a,n,Norma);
	}
	
	
	//LOG("do rec");
	EigenValImpl(a, n, eig, 0, n , Norma);
	
	if(COMPLEXVAL!=-1){
		
		printf("\nVector Eigen Value\n");
		print_vector(eig,n);
	
		double firstRes = 0, secondRes = 0;
		for(int i = 0; i < n ; i++){
			firstRes+=eig[i];
			secondRes+=eig[i]*eig[i];
		}
	
        cout << "\nI check: nado " << FirstCheck << ", est' " << firstRes << "\n";
        cout << "II check: nado " << SecondCheck << ", est' " << secondRes << "\n";
	}else{
		printf("EXIST COMPLEX EIGEN VAL\n");
	}
}

void EigenValImpl(double* a, int size, double* eig, int begin, int end, double Norma) {
	if (end - begin == 1){
		eig[begin] = a[begin * size + begin];
		return;
	}
//	LOG("posle 1 if");
	if (end - begin == 2){
		EigValDim2(a,size,eig,begin,end);
		return;
	}
	
	double s_k = a[(end - 2) * size + (end - 1)];
	
	for (int i = begin + 1; i < end; i++) {
		if (fabs(a[(i - 1) * size + i]) < EPS * Norma) {
			EigenValImpl(a, size, eig, begin, i, Norma);
			EigenValImpl(a, size, eig, i, end, Norma);
			return;
		}
	}
	
    int iterCount = 0;
	while (fabs(a[(end - 2) * size + (end - 1)]) > EPS * Norma) {
	//	s_k = finds_k(a,size,begin, end,Norma);
		//LOG(s_k);
		s_k = a[(end - 1) * size + (end - 1)];
		// LOG(s_k);
		matr_sub_sdvig(a, s_k, size, begin, end);
		LRdeform(a, size, begin, end);
		RLmult(a, size, begin, end);
		matr_sub_sdvig(a, -s_k, size, begin, end);
        iterCount++;
	}
  //  printf("Na %d znachenie %d iteratcii\n", end - 1, iterCount);
	eig[end - 1] = a[(end - 1) * size + (end - 1)];
	
	EigenValImpl(a, size, eig, begin, end - 1, Norma);
}

double finds_k(double *a, int size, int begin,int end, double Norma){
	double val1=EigValDim2(a,size,begin,end);
	double val2=0.5*a[(end - 2) * size + (end - 1)] + a[(end - 1) * size + (end - 1)];
	double val3=a[(end - 1) * size + (end - 1)];
	if(fabs(val1) > EPS*Norma){ 
		return val1;
	}else if(fabs(val2) > EPS*Norma){	
		return  val2;
	}else if(fabs(val3) > EPS*Norma){	
		return  val3;
	}
    return -1;
}


void Reflecsions(double *a, int n,double Norma){
	
	double sk=0, colNorm=0, xNorma=0;
	double *p, *x , *b;
	

	b = new double [n*n];
	
	memcpy(b,a,n*n*sizeof(double));
	
	x = new double [n];
	for(int i = 0;i<n;i++){
		x[i]=0;
	}
	
	for(int k = 0; k < n-2;k++){

		sk=0; colNorm=0; xNorma=0;
		p = a + k*n;
		for(int i= k + 2; i < n; i++){
			sk+= p[i]*p[i]; 
		}

		colNorm = sqrt(p[k + 1 ]*p[k + 1 ] + sk);

		x[k + 1] = p[(k + 1)] - colNorm;
		memcpy(x + k+2 , a + k*n + (k + 2) , (n - (k + 2))*sizeof(double)); 

		xNorma= sqrt(x[k+1]*x[k+1] + sk	);

		if(xNorma < EPS*Norma ){
			//cout<<"ERROR: vector x is 0 - method REFLECTION INAPPL. \n"<<endl;
			continue;
		}
		for(int i=k+1 ; i < n; i++){
			x[i]/=xNorma;
		}

		memcpy(b,a,n*n*sizeof(double));
		LeftMult(a,b,x,n,k);
		
		memcpy(a,b,n*n*sizeof(double));

		RightMult(b,a,x,n,k);

	}
	
	delete [] x;
	delete [] b;
}

void print_matrix_t(double* m, int size) {
        for (int y = 0; y < ((size < 10) ? size : 10) ; y++) {
                for (int x = 0; x < ((size < 10) ? size : 10); x++) {
                                cout << m[x* size + y] << "\t";
                        }
                cout << endl;
        }
}

void print_matrix_t(double* m, int size,int begin , int end) {
	int l = ((end - begin) < 10 ? end : begin + 10)  ;
	
        for (int y = begin; y <l ; y++) {
                for (int x = begin; x < l; x++) {
                                cout << m[x* size + y] << "\t";
                        }
                cout << endl;
        }
}

void LeftMult(double *a, double *b , double *x,  int n, int k){
	double scalar = 0;
	double *p1, *p2;
	int l = min((n - (k+1)) % 6 , n);
//	LOG(l);
	for( int i = 0; i< n;i++){
		scalar=0; p1 = a + i*n; p2 = b + i*n;

		for(int j = k + 1 ; j < k + 1 + l  ; j++){
			scalar+= x[j]*p1[j];
		}
		
		for ( int j = k + 1 + l; j < n ;j+=6){
			scalar+= x[j]*p1[j];
			scalar+= x[j + 1]*p1[j + 1];
			scalar+= x[j + 2]*p1[j + 2];
			scalar+= x[j + 3]*p1[j + 3];
			scalar+= x[j + 4]*p1[j + 4];
			scalar+= x[j + 5]*p1[j + 5];
		}
		for(int j = k + 1 ; j < k + 1 + l  ; j++){
			p2[(j)] = p1[(j)] - 2*x[(j)]*scalar;
		}
		
		for ( int j = k + 1 + l; j < n ;j+=6){
			p2[j] = p1[j] - 2*x[j]*scalar;
			p2[j + 1] = p1[j + 1] - 2*x[j + 1]*scalar;
			p2[j + 2] = p1[j + 2] - 2*x[j + 2]*scalar;
			p2[j + 3] = p1[j + 3] - 2*x[j + 3]*scalar;
			p2[j + 4] = p1[j + 4] - 2*x[j + 4]*scalar;
			p2[j + 5] = p1[j + 5] - 2*x[j + 5]*scalar;
			
		}
		
	}

}	

void RightMult(double *a, double *b , double *x,  int n, int k){
	double scalar = 0;
	int l = min((n - (k+1) ) % 6 , n);
	for( int i = 0; i< n;i++){
		scalar=0; 
		// for ( int j = k + 1; j < n;j++){
			// scalar+= x[j]*a[j*n + i];
		// }
		
		
		for(int j = k + 1 ; j < k + 1 + l  ; j++){
			scalar+= x[j]*a[j*n + i];
		}
		
		for ( int j = k + 1 + l; j < n ;j+=6){
			scalar+= x[j] *  a[j*n + i];
			scalar+= x[j + 1]*a[(j + 1)*n + i];
			scalar+= x[j + 2]*a[(j + 2)*n + i];
			scalar+= x[j + 3]*a[(j + 3)*n + i];
			scalar+= x[j + 4]*a[(j + 4)*n + i];
			scalar+= x[j + 5]*a[(j + 5)*n + i];
		}
		
		
		
		for(int j = k + 1 ; j < k + 1 + l  ; j++){
			b[(j)*n + i] = a[(j)*n + i] - 2*x[(j)]*scalar;
		}
		
		
		for ( int j = k + 1 + l; j < n ;j+=6){
			b[j*n + i] = a[j*n + i] - 2*x[j]*scalar;
			b[(j + 1)*n + i] = a[(j + 1)*n + i] - 2*x[j + 1]*scalar;
			b[(j + 2)*n + i] = a[(j + 2)*n + i] - 2*x[j + 2 ]*scalar;
			b[(j + 3)*n + i] = a[(j + 3)*n + i] - 2*x[(j + 3)]*scalar;
			b[(j + 4)*n + i] = a[(j + 4)*n + i] - 2*x[j + 4]*scalar;
			b[(j + 5)*n + i] = a[(j + 5)*n + i] - 2*x[j + 5]*scalar;
		}
		
		
		
	}

}	

void EigValDim2(double *a , int size, double *eig, int begin, int end){
	double trace=a[begin*size + begin] + a[(begin + 1)*size + (begin + 1)];
	double det= a[begin*size + begin]*a[(begin + 1)*size + (begin + 1)] - a[(begin + 1)*size + begin]*a[(begin)*size + (begin + 1)] ;
	double D = trace*trace - 4*det;
	if(D < 0) {
		printf("Eigen value Eig%d and Eig%d is COMPLEX \n",begin,begin + 1 );
		COMPLEXVAL=-1;
		return;
	}else{
		eig[begin] = (trace - sqrt(D))/2;
		eig[begin + 1] = (trace + sqrt(D))/2;
	}

}	

int EigValDim2(double *a , int size, int begin, int end){
	begin = end - 2;
	double trace=a[begin*size + begin] + a[(begin + 1)*size + (begin + 1)];
	double det= a[begin*size + begin]*a[(begin + 1)*size + (begin + 1)] - a[(begin + 1)*size + begin]*a[(begin)*size + (begin + 1)] ;
	double D = trace*trace - 4*det;
	//LOG(D);
	if(D < 0) {
		return 0 ;
	}else{
		
		return max((trace - sqrt(D))/2,(trace + sqrt(D))/2);
	}

}	


void matr_sub_sdvig(double *a ,double sdvig, int n, int k1 ,int  k2){
	for(int i = k1; i < k2; i ++ ){
		a[i*n + i]-=sdvig;
	}
}

double  CheckNorma(double * a , int n){
	double sum = 0;
	
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			sum+=a[i*n + j] * a[j*n + i];
		}
	}
	return sum;
}


void RLmult(double *a , int n, int k1, int k2){
	
	
	for(int l = k1 ; l < k2 -1 ; l ++){ //stolci
			a[l*n + k1]+=a[(l+1)*n + k1]*a[(l)*n + l + 1];
	}
	
	for(int i = k1 + 1 ; i < k2 ; i++){ // stroki
	    double templ = a[i*n + i];
		for(int l = i ; l < k2 -1 ; l ++){ //stolci
			a[l*n + i]+=a[(l+1)*n + i]*a[(l)*n + l + 1];
		}
		a[(i-1)*n + i]*=templ;
	}
	
	
}


void LRdeform(double *a, int n, int k1, int k2 ){
	double li = 0;
	for( int i = k1 + 1 ; i < k2 ;i++){
		//LOG(i);
		a[(i-1)*n + i ]/=a[(i-1)*n + (i-1)];
		li = a[(i-1)*n + i ];
		for ( int l = i; l < k2 ; l++){
			a[l*n + i]-=li*a[l*n + (i-1)];
		}
		
	}
}

double Trace(double *a, int n){
	
	double sum=0;
	for(int i = 0; i < n; i++){
		sum+=a[i*n + i];
	}
	return sum;
}



