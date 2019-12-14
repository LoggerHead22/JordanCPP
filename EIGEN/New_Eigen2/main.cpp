#include <ctime>
#include <string>
#include <vector>
#include <algorithm>
#include "gauss.h"
//#include "tests.h"
using namespace std;

void Reflecsions(double* a, int n, double Norma);
void ReflMatr(double* U, double* x, int n, int k);
void print_matrix_t(double* m, int size);
void LeftMult(double* a, double* b, double* x, int n, int k);
void RightMult(double* a, double* b, double* x, int n, int k);
void EigenVal(double* a, double* eig, int n);
void matr_sub_sdvig(double* a, double sdvig, int n, int k1, int k2);
double  CheckNorma(double* a, int n);
void RLmult(double* a, int n, int k1, int k2);
int LRdeform(double* a, int n, int k1, int k2, double Norma);
double Trace(double* a, int n);
void EigenValImpl(double* a, int size, double* eig, int begin, int end, double Norma, double* save);
void EigValDim2(double* a, int size, double* eig, int begin, int end);
void EigValDim2Sdvig(double* a, int size, int begin, int end, double* val1, double* val2);
double finds_k(double* a, int size, int begin, int end, double Norma);
void print_matrix_t(double* m, int size, int begin, int end);
int LRdeform_simm(double* a, int n, int k1, int k2, double Norma);
void RLmult_simm(double* a, int n, int k1, int k2);


static int COMPLEXVAL = 0;
static int iterCount = 0;
static int simmetric = 0;

int main(int argc, char** argv)
{
	srand(time(0));
	double* a, * ev;
	int n = 0;
	string name = "";




	if (!((argc == 4) || (argc == 3)) || ((n = stoi(argv[1])) <= 0) || ((EPS = stod(argv[2])) <= 0)) {
		printf("usage: %s n EPS [name] \n ", argv[0]);
		return -1;
	}

	try { a = new double[n * n]; }
	catch (...) {
		printf("Not enough memory a ");
		return -2;
	}

	LOG(1);
	if (argc == 4) {
		name = argv[3];
		LOG(name);
		if (init_matrix_file(a, n, name) < 0) {
			delete[] a;
			return 0;

		};
	}
	LOG(2);

	if (argc == 3) {
		formula_matr(a, n);
	}

	try { ev = new double[n]; }
	catch (...) {
		printf("Not enough memory x ");
		delete[]a;
		return -2;
	}

	printf("Init. Matrix \n");
	print_matrix_t(a, n);

	//double t = clock();
	EigenVal(a, ev, n);
	//t = (clock() - t) / CLOCKS_PER_SEC;

	//cout << "TIME: " << t << endl;

	delete[] a;
	delete[] ev;

	return 0;
}

void Transpose(double* a, int n) {

	for (int i = 1; i < n; i++) {
		swap(a[(i - 1) * n + i], a[i * n + (i - 1)]);
	}

	for (int i = 0; i < n - 2; i++) {
		for (int j = i + 2; j < n; j++) {
			a[i * n + j] = a[j * n + i];
			a[j * n + i] = 0;
		}
	}

}

void Transpose_refl(double* a, int n, int k) {

	
	for (int i = k; i < n; i++) {
		for (int j = 0; j < i; j++) {
			swap(a[i * n + j],a[j * n + i]);
		}
	}

}

void isSimmetric(double* a, int n, int Norma) {

	
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			if(fabs(a[i * n + j] - a[j * n + i]) > EPS*Norma) {
				simmetric=-1;
				return;
			}
		}
	}
	cout<<"Simmetric matrix"<<endl;
}



void pickUpGoodSk(double& s_k, double* a, int n, int begin, int end, double Norma, double* save) {
	//  LOG("pickUpSk");

	std::copy(a, a + n * n, save);
	double possibleValues[4] = { 0, 0, 0, 0 };

	double val11 = 0, val12 = 0;
	EigValDim2Sdvig(a, n, end - 2, end, &val11, &val12);
	if(fabs(a[(end - 1) * n + (end - 1)] - val11) < fabs(a[(end - 1) * n + (end - 1)] - val12) ){
			val12 = val11;
	}
	double val2 = 0.5 * a[(end - 1) * n + (end - 2)] + a[(end - 1) * n + (end - 1)];
	double val3 = a[(end - 1) * n + (end - 1)];
	if (fabs(val11) > EPS* Norma) {
		possibleValues[0] = (fabs(val11) > fabs (val12)?val11:val12);
	}
	if (fabs(val12) > EPS* Norma) {
		possibleValues[1] = val12;
	}
	if (fabs(val2) > EPS* Norma) {
		possibleValues[3] = val2;
	}
	if (fabs(val3) > EPS* Norma) {
		possibleValues[2] = val3;
	}
	int flag;
	//	cout<<val11<<" "<<val12<<" "<<val2<<" "<<val3<<endl;
	//	cout<<possibleValues[0]<<" "<<possibleValues[1]<<" "<<possibleValues[2]<<" "<<possibleValues[3]<<endl;
	for (int k = 1; k < 4 ; k++) {
		s_k = possibleValues[k];
		//LOG(s_k);
		if (fabs(s_k) < EPS * Norma) continue;
		//	LOG(s_k);
		matr_sub_sdvig(a, s_k, n, begin, end);
		// if(simmetric==0){
			// flag = LRdeform_simm(a, n, begin, end, Norma);
		// }else{
			// flag = LRdeform(a, n, begin, end, Norma);
		// }
		flag = LRdeform(a, n, begin, end, Norma);
		if (flag == 0) {

			return;
		}
		else {
			copy(save, save + n * n, a);
		}
	}

}

void EigenVal(double* a, double* eig, int n) {
	double FirstCheck = Trace(a, n), SecondCheck = CheckNorma(a, n);
	LOG(EPS);
	double Norma = norma_matr(a, n);
	Norma = min(Norma, 1.0);
//	LOG(Norma);
	double time =clock(), ttime = clock();
	
	if (n >= 3) {
		isSimmetric(a,n,Norma);
		Reflecsions(a, n, Norma);
	}

	//print_matrix_t(a,n);

	double* save = new double[n * n];
	if(simmetric!=0){
		Transpose(a, n);
	}
	
	ttime = clock() - ttime;
	cout << "REFLECTION  TIME: " << ttime / CLOCKS_PER_SEC << endl;
	LN;
	//print_matrix(a,n);
//	return;
	Norma = norma_matr(a, n);
	Norma = min(Norma, 1.0);
	LOG(Norma);
	EigenValImpl(a, n, eig, 0, n, Norma, save);

	if (COMPLEXVAL != -1) {

		printf("\nVector Eigen Value\n");
		print_vector(eig, n);

		double firstRes = 0, secondRes = 0;
		for (int i = 0; i < n; i++) {
			firstRes += eig[i];
			secondRes += eig[i] * eig[i];
		}

		cout << "\nI check: nado " << FirstCheck << ", est' " << firstRes << ", differ: " << firstRes - FirstCheck << endl;
		cout << "II check: nado " << SecondCheck << ", est' " << secondRes << ", differ: " << secondRes - SecondCheck << endl;
		cout << "ITER COUNT : " << iterCount << endl;
	}
	else {
		printf("EXIST COMPLEX EIGEN VAL\n");
	}
	delete[] save;
	time=clock() - time;
	printf("REFL.TIME %f \nLR TIME %f \nFULL TIME %f\n",ttime/CLOCKS_PER_SEC,(time -ttime)/CLOCKS_PER_SEC,time/CLOCKS_PER_SEC);
}

void EigenValImpl(double* a, int size, double* eig, int begin, int end, double Norma, double* save) {
	//printf("Recurent func with begin: %d and end: %d\n",begin ,end);
	if (end - begin == 1) {
		eig[begin] = a[begin * size + begin];
		return;
	}

	if (end - begin == 2) {
		EigValDim2(a, size, eig, begin, end);
		return;
	}

	double s_k = 0;

	for (int i = begin + 1; i < end; i++) {
		if (fabs(a[(i - 1) * size + i]) < EPS * Norma) {
		//	printf("Reccurent func with (%d,%d) and (%d,%d)\n",begin,i,i,end);
			EigenValImpl(a, size, eig, begin, i, Norma, save);
			EigenValImpl(a, size, eig, i, end, Norma, save);
			return;
		}
	}


	while (fabs(a[(end - 1) * size + (end - 2)]) > EPS* Norma) {

		s_k = 0;

		pickUpGoodSk(s_k, a, size, begin, end, Norma, save);

		// if(simmetric==0){
			// RLmult_simm(a, size, begin, end);
		// }else{
			// RLmult(a, size, begin, end);
		// }
		RLmult(a, size, begin, end);
		matr_sub_sdvig(a, -s_k, size, begin, end);
		iterCount++;
		// for (int i = begin + 1; i < end; i++) {
			// if (fabs(a[(i - 1) * size + i]) < EPS * Norma) {
				// printf("Reccurent func with (%d,%d) and (%d,%d)\n",begin,i,i,end);
				// EigenValImpl(a, size, eig, begin, i, Norma, save);
				// EigenValImpl(a, size, eig, i, end, Norma, save);
				// return;
			// }
		// }
	}
	//  printf("Na %d znachenie %d iteratcii\n", end - 1, iterCount);
	eig[end - 1] = a[(end - 1) * size + (end - 1)];

	EigenValImpl(a, size, eig, begin, end - 1, Norma, save);
}

// double finds_k(double *a, int size, int begin,int end, double Norma) {
	// double val1=EigValDim2(a,size,begin,end);
	// double val2=0.5*a[(end - 2) * size + (end - 1)] + a[(end - 1) * size + (end - 1)];
	// double val3=a[(end - 1) * size + (end - 1)];
	// if(fabs(val1) > EPS*Norma){ 
		// return val1;
	// }else if(fabs(val2) > EPS*Norma){	
		// return  val2;
	// }else if(fabs(val3) > EPS*Norma){	
		// return  val3;
	// }
	// return -1;
// }


void Reflecsions(double* a, int n, double Norma) {

	double sk = 0, colNorm = 0, xNorma = 0;
	double* p, * x;


	//	b = new double [n*n];

		//memcpy(b,a,n*n*sizeof(double));

	x = new double[n];
	for (int i = 0; i < n; i++) {
		x[i] = 0;
	}

	for (int k = 0; k < n - 2; k++) {

		sk = 0; colNorm = 0; xNorma = 0;
		p = a + k * n;
		for (int i = k + 2; i < n; i++) {
			sk += p[i] * p[i];
		}
		if (sk < EPS * Norma ) {
			//cout<<"ERROR: vector x is 0 - method REFLECTION INAPPL. \n"<<endl;
			continue;
		}
		colNorm = sqrt(p[k + 1] * p[k + 1] + sk);
		if (colNorm < EPS *  Norma ) {
			//cout<<"ERROR: vector x is 0 - method REFLECTION INAPPL. \n"<<endl;
			continue;
		}
		x[k + 1] = p[(k + 1)] - colNorm;
		memcpy(x + k + 2, a + k * n + (k + 2), (n - (k + 2)) * sizeof(double));

		xNorma = sqrt(x[k + 1] * x[k + 1] + sk);

		if (xNorma < EPS * Norma ) {
			//cout<<"ERROR: vector x is 0 - method REFLECTION INAPPL. \n"<<endl;
			continue;
		}
		for (int i = k + 1; i < n; i++) {
			x[i] /= xNorma;
		}

		//	memcpy(b,a,n*n*sizeof(double));
		LeftMult(a, a, x, n, k);
		Transpose_refl(a,n,k);
		LeftMult(a, a, x, n, k);
		if(simmetric != 0){
			Transpose_refl(a,n,k);
		}
		//		memcpy(a,b,n*n*sizeof(double));

		//RightMult(a, a, x, n, k);

	}

	delete[] x;
	//	delete [] b;
}

void print_matrix_t(double* m, int size) {
	for (int y = 0; y < ((size < 10) ? size : 10); y++) {
		for (int x = 0; x < ((size < 10) ? size : 10); x++) {
			cout << m[x * size + y] << "\t";
		}
		cout << endl;
	}
}

void print_matrix_t(double* m, int size, int begin, int end) {
	int l = ((end - begin) < 10 ? end : begin + 10);

	for (int y = begin; y < l; y++) {
		for (int x = begin; x < l; x++) {
			cout << m[x * size + y] << "\t";
		}
		cout << endl;
	}
}

void LeftMult(double* a, double* b, double* x, int n, int k) {
	double scalar = 0;
	double* p1, * p2;
	int l = min((n - (k + 1)) % 6, n);
	//	LOG(l);
	for (int i = 0; i < n; i++) {
		scalar = 0; p1 = a + i * n; p2 = b + i * n;

		for (int j = k + 1; j < k + 1 + l; j++) {
			scalar += x[j] * p1[j];
		}

		for (int j = k + 1 + l; j < n; j += 6) {
			scalar += x[j] * p1[j];
			scalar += x[j + 1] * p1[j + 1];
			scalar += x[j + 2] * p1[j + 2];
			scalar += x[j + 3] * p1[j + 3];
			scalar += x[j + 4] * p1[j + 4];
			scalar += x[j + 5] * p1[j + 5];
		}
		for (int j = k + 1; j < k + 1 + l; j++) {
			p2[(j)] = p1[(j)] - 2 * x[(j)] * scalar;
		}

		for (int j = k + 1 + l; j < n; j += 6) {
			p2[j] = p1[j] - 2 * x[j] * scalar;
			p2[j + 1] = p1[j + 1] - 2 * x[j + 1] * scalar;
			p2[j + 2] = p1[j + 2] - 2 * x[j + 2] * scalar;
			p2[j + 3] = p1[j + 3] - 2 * x[j + 3] * scalar;
			p2[j + 4] = p1[j + 4] - 2 * x[j + 4] * scalar;
			p2[j + 5] = p1[j + 5] - 2 * x[j + 5] * scalar;

		}

	}

}

void RightMult(double* a, double* b, double* x, int n, int k) {
	double scalar = 0;
	int l = min((n - (k + 1)) % 6, n);
	for (int i = 0; i < n; i++) {
		scalar = 0;
		// for ( int j = k + 1; j < n;j++){
			// scalar+= x[j]*a[j*n + i];
		// }


		for (int j = k + 1; j < k + 1 + l; j++) {
			scalar += x[j] * a[j * n + i];
		}

		for (int j = k + 1 + l; j < n; j += 6) {
			scalar += x[j] * a[j * n + i];
			scalar += x[j + 1] * a[(j + 1) * n + i];
			scalar += x[j + 2] * a[(j + 2) * n + i];
			scalar += x[j + 3] * a[(j + 3) * n + i];
			scalar += x[j + 4] * a[(j + 4) * n + i];
			scalar += x[j + 5] * a[(j + 5) * n + i];
		}



		for (int j = k + 1; j < k + 1 + l; j++) {
			b[(j)*n + i] = a[(j)*n + i] - 2 * x[(j)] * scalar;
		}


		for (int j = k + 1 + l; j < n; j += 6) {
			b[j * n + i] = a[j * n + i] - 2 * x[j] * scalar;
			b[(j + 1) * n + i] = a[(j + 1) * n + i] - 2 * x[j + 1] * scalar;
			b[(j + 2) * n + i] = a[(j + 2) * n + i] - 2 * x[j + 2] * scalar;
			b[(j + 3) * n + i] = a[(j + 3) * n + i] - 2 * x[(j + 3)] * scalar;
			b[(j + 4) * n + i] = a[(j + 4) * n + i] - 2 * x[j + 4] * scalar;
			b[(j + 5) * n + i] = a[(j + 5) * n + i] - 2 * x[j + 5] * scalar;
		}



	}

}

void EigValDim2(double* a, int size, double* eig, int begin, int end) {
	double trace = a[begin * size + begin] + a[(begin + 1) * size + (begin + 1)];
	double det = a[begin * size + begin] * a[(begin + 1) * size + (begin + 1)] - a[(begin + 1) * size + begin] * a[(begin)*size + (begin + 1)];
	double D = trace * trace - 4 * det;
	double LocalNorma = max(fabs(a[begin * size + begin]) + fabs(a[(begin ) * size + begin  + 1 ]), fabs(a[ (begin + 1) * size + begin ]) + fabs(a[(begin + 1) * size + begin + 1]));
	LocalNorma=min(LocalNorma,1.0);
	LOG(LocalNorma);
	LOG(D);
	if (D < 0 && fabs(D)>EPS*10* LocalNorma) {
		printf("Eigen value Eig%d and Eig%d is COMPLEX with D: %f \n", begin, begin + 1, D);
		COMPLEXVAL = -1;
		return;
	}
	else {
		D = fabs(D);
		if(fabs(D) < EPS*10* LocalNorma) D = 0;
		if(trace > 0){
               eig[begin]  = (trace - sqrt(D)) / 2;
				   if(fabs(eig[begin]) > EPS*LocalNorma){
						eig[begin + 1]  = det / (eig[begin]);
				   }else{ eig[begin + 1] = (trace);
				   }
        }else {
				   
                  eig[begin]  = (trace + sqrt(D)) / 2;
                  if(fabs(eig[begin]) > EPS*LocalNorma){
						eig[begin + 1]  = det / (eig[begin]);
				  }else{ eig[begin + 1] = (trace);
				  }
        }
	}

}

void EigValDim2Sdvig(double* a, int size, int begin, int end, double* val1, double* val2) {
	begin = end - 2;
	double trace = a[begin * size + begin] + a[(begin + 1) * size + (begin + 1)];
	double det = a[begin * size + begin] * a[(begin + 1) * size + (begin + 1)] - a[(begin + 1) * size + begin] * a[(begin)*size + (begin + 1)];
	double D = trace * trace - 4 * det;
	double LocalNorma = max(fabs(a[begin * size + begin]) + fabs(a[(begin ) * size + begin  + 1 ]), fabs(a[ (begin + 1) * size + begin ]) + fabs(a[(begin + 1) * size + begin + 1]));
	LocalNorma=min(LocalNorma,1.0);
	//LOG(D);
	if (D < 0 && fabs(D)>EPS* LocalNorma) {
		return;
	}
	else {
		D = fabs(D);
		if(fabs(D) < EPS*10* LocalNorma) D = 0;
		   if(trace > 0){
               *val1 = (trace - sqrt(D)) / 2;
					if(fabs(*val1) > EPS * LocalNorma){
						*val2 = det / (*val1);
					}else{
						*val2 = trace;
					}
               }else {
                   *val1 = (trace + sqrt(D)) / 2;
                    if(fabs(*val1) > EPS * LocalNorma){
						*val2 = det / (*val1);
					}else{
						*val2 = trace;
					}
        }
	}

}


void matr_sub_sdvig(double* a, double sdvig, int n, int k1, int  k2) {
	for (int i = k1; i < k2; i++) {
		a[i * n + i] -= sdvig;
	}
}

double  CheckNorma(double* a, int n) {
	double sum = 0;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			sum += a[i * n + j] * a[j * n + i];
		}
	}
	return sum;
}


void RLmult(double* a, int n, int k1, int k2) {


	for (int l = k1; l < k2 - 1; l++) { //stolci
		a[k1 * n + l] += a[k1 * n + (l + 1)] * a[(l + 1) * n + l];
	}

	for (int i = k1 + 1; i < k2; i++) { // stroki
		double templ = a[i * n + i];
		int krai = min( (k2 -1 -i)%4,k2-1);
		
		for (int l = i; l < i + krai ; l++) { //stolci
			a[i * n + l ] += a[i * n + (l + 1)] * a[(l + 1) * n + l ];
		}
		
		for (int l = i + krai; l < k2 - 1; l+=4) { //stolci
			a[i * n + l ]    += a[i * n + (l + 1)] * a[(l + 1) * n + l ];
			a[i * n + l + 1] += a[i * n + (l + 2)] * a[(l + 2) * n + l + 1];
			a[i * n + l + 2] += a[i * n + (l + 3)] * a[(l + 3) * n + l + 2];
			a[i * n + l + 3] += a[i * n + (l + 4)] * a[(l + 4) * n + l + 3];
		}
		
		a[(i)*n + i - 1] *= templ;
	}


}


int LRdeform(double* a, int n, int k1, int k2, double Norma) {
	double li = 0;
	//LOG("IN LRdeform ");
	for (int i = k1 + 1; i < k2; i++) {
		
		int krai = min( (k2 - i)%4,k2);
		
		if (fabs(a[(i - 1) * n + (i - 1)]) < EPS * Norma) return -1;
		
		a[(i)*n + i - 1] /= a[(i - 1) * n + (i - 1)];
		li = a[(i)*n + i - 1];
		for (int l = i ; l < i + krai; l++) {
			a[i * n + l] -= li * a[(i - 1) * n + l];
		}
		for (int l = i + krai ; l < k2; l+=4) {
			a[i * n + l] -= li * a[(i - 1) * n + l];
			a[i * n + l + 1] -= li * a[(i - 1) * n + l + 1];
			a[i * n + l + 2] -= li * a[(i - 1) * n + l + 2];
			a[i * n + l + 3] -= li * a[(i - 1) * n + l + 3];
		}
		

	}
	return 0;
}

int LRdeform_simm(double* a, int n, int k1, int k2, double Norma) {

	//LOG("IN simm LRdeform ");
	for (int i = k1 + 1; i < k2; i++) {
		
		if (fabs(a[(i - 1) * n + (i - 1)]) < EPS * Norma) return -1;
		a[i*n + i-1] = a[i*n + i-1] / a[(i-1)*n + i - 1];
		
		a[i*n + i] = a[i*n + i] - a[i*n + i-1]*a[(i-1)*n + i];
		if(i < k2 - 1){
			a[i*n + i + 1] = a[i*n + i + 1] - a[i*n + i-1]*a[(i-1)*n + i + 1];
		}

	}
	return 0;
}

void RLmult_simm(double* a, int n, int k1, int k2) {


	for(int i = k1; i < k2 -1 ; i++){
        a[i*n + i] = a[i*n + i] + a[i*n + i + 1] * a[(i+1)*n + i];

        a[(i+1)*n + i] = a[(i+1)*n + i + 1]*a[(i+1)*n + i];
       
    }


}


double Trace(double* a, int n) {

	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += a[i * n + i];
	}
	return sum;
}


