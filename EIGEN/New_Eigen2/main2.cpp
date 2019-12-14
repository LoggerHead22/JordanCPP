#include <pthread.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <cassert>
#include <cstring>
#include <sched.h>
//#include <sys/time.h>
//#include <sys/sysinfo.h>
//#include <sys/resource.h>
#include <float.h>
using namespace std;

#define ERR_CANNOT_OPEN -1
#define ERR_READ -2 
#define ERR_DEG_MATRIX -3   // РІС‹СЂРѕР¶РґРµРЅРЅР°СЏ РјР°С‚СЂРёС†Р°
#define EPS 5e-15

#define LOG(...) std::cout<< #__VA_ARGS__<< " = " <<(__VA_ARGS__)<< "\n"
#define LN       std::cout << "\n";






struct Args{
    char *filename;
    double *a,*b, thr_time;
    int p,n1,n2;
    int num;
    int error = 0;
};

struct krai{
    double first,last;
    int beg,end;
};

static pthread_mutex_t total_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_barrier_t barrier;




void *f(void * args);
void print_matrix(double* m, int size1,int size2);
void formula_matr(double* a, int n1,int n2);
int read_matrix(double* a, int n, const string& name);
int init_matrix_file(double* a, int n, const string& name);
double get_time();
double get_full_time();




void print_copy(double *b, int n2, int thr_num, int p, int beg, int len){
	pthread_mutex_lock(&total_mutex);

		printf("Beg: %d and len %d of %d potoka\n",beg,len,thr_num);
		printf("Copy of %d potoka\n",thr_num);
		for(int i = 0;i< 2* n2;i++){
			if( i % n2 ==0 && i!=0) LN;
			cout<<b[i]<<" ";
			
		}
		LN;
		printf("Down of %d potoka\n",thr_num);
		for(int i = 0;i< 2 * n2;i++){
			cout<<(b + 2*n2)[i]<<" ";
			if( i % n2 ==0 && i!=0) LN;
		}
		
		LN;

	pthread_mutex_unlock(&total_mutex);
}

void *f(void * args){
    Args& data= *((Args *) args);
    double *a=data.a;
    double *b=data.b;
    int n1=data.n1, n2=data.n2, thr_num=data.num, p=data.p;
	
    //data.thr_time=get_time();
    if(n1 > 4 && n2 > 4){
        int len=(n1%p > 0 && (thr_num >= p - n1%p)? int(n1/p) + 1: n1/p );
        int beg=((thr_num >= p - n1%p)? thr_num * int(n1/p) + thr_num - p  + n1%p  :thr_num*(n1/p));

        // if(thr_num!=0){
			// memcpy(b,a +  ((beg - 2)<0?0 : beg - 2)*n2,2*n2*sizeof(double));
        // }else {
            // beg+=2;
            // memcpy(b,a +  (beg - 2)*n2,2*n2*sizeof(double));
            // beg-=2;
        // }
		 memcpy(b,a +  ((beg - 2)<0?0 : beg - 2)*n2,2*n2*sizeof(double));
		 memcpy(b + 2*n2 ,a + ((beg + len + 2)>=n1?n1 -2:beg + len )*n2,2*n2*sizeof(double));
        // if(thr_num!=p-1){

            // memcpy(b + 2*n2 ,a + ((beg + len)>=n1-1?n1 -2:beg + len)*n2,2*n2*sizeof(double));
			   
        // }else {
            // len-=2;
            // memcpy(b + 2*n2 ,a + (beg + len)*n2,2*n2*sizeof(double));
            // len+=2;
        // }
        double buf;

        if(thr_num==0){
            beg+=2;
            len-=2;
        }
		
        pthread_barrier_wait(&barrier);
		
		//print_copy(b,n2, thr_num, p,  beg,len);
		
		
        for(int i=beg;i<beg + len - 2;i++){

            double *copy=b + ((i - beg)%2)*n2;
            copy[0] = a[i*n2];
            copy[1] = a[i*n2  + 1];
            for(int k = 2; k < n2-2;k++){
      //         cout<<k <<" "<<copy[k -2]<<" "<< a[i*n2 + k + 2]<<" "<<copy[k]<<" "<<a[(i + 2)*n2 + k]<<endl;
                buf =  a[i*n2 + k ];
                a[i*n2 + k ] =  (copy[k -2] +  a[i*n2 + k + 2]  +  copy[k] + a[(i + 2)*n2 + k]) / 4 ;
                copy[k] = buf;
            }
            copy[n2 -2 ] = a[i*n2 + n2 -2 ];
            copy[n2 -1 ] = a[i*n2  + n2 - 1];
        }
          pthread_barrier_wait(&barrier);
        if(thr_num!=p-1){
				if(len !=1){
                    for(int i = beg + len - 2; i< beg + len ; i++){
						if(i < 2){
							
							continue;
						}
						if( p == (n1 + 1)/2 && thr_num==1){
							for(int i = 0; i<n2; i++){
							swap(b[i], b[i + n2]);	
								
							}
						}
					
						
                        double *copy=b + ((i - beg)%2)*n2;

                        double *down = b + 2*n2 + n2*(i - (beg + len - 2));
						// if(thr_num==1){
							// LOG("copy");
							// for(int i = 0;i< n2;i++){
								// cout<<copy[i]<<" ";
								// if( i % n2 ==0 && i!=0) LN;
							// }
							// LN;
							// LOG("down");
							// for(int i = 0;i< n2;i++){
								// cout<<down[i]<<" ";
							// }
							// LN;
						// }
                        copy[0] = a[i*n2];
                        copy[1] = a[i*n2  + 1];
                        for(int k = 2; k < n2-2;k++){
                            // if(thr_num==1){
                                    // cout<<k<<" "<<copy[k -2]<<" "<< a[i*n2 + k + 2]<<" "<<copy[k]<<" "<<down[k]<<endl;
                            // }
                            buf =  a[i*n2 + k ];
                            a[i*n2 + k ] =  (copy[k -2] +  a[i*n2 + k + 2]  +  copy[k] + down[k]) / 4 ;
                            copy[k] = buf;
                        }
                        copy[n2 -2 ] = a[i*n2 + n2 -2 ];
                        copy[n2 -1 ] = a[i*n2  + n2 - 1];
                    }
				}else if(beg > 1 && beg < p -2){
						int i = beg;
						
					
						
						// if( p == n1 && thr_num==p-2){
							// for(int i = 0; i<n2; i++){
							// swap(b[2*n2 + i], b[n2*2 + i + n2]);	
								
							// }
						// }
						double *copy=b + ((i - beg)%2)*n2;

						double *down = b + 2*n2 + n2;
						
					
						copy[0] = a[i*n2];
                        copy[1] = a[i*n2  + 1];
                        for(int k = 2; k < n2-2;k++){
                            // if(thr_num==1){
                                    // cout<<k<<" "<<copy[k -2]<<" "<< a[i*n2 + k + 2]<<" "<<copy[k]<<" "<<down[k]<<endl;
                            // }
                            buf =  a[i*n2 + k ];
                            a[i*n2 + k ] =  (copy[k -2] +  a[i*n2 + k + 2]  +  copy[k] + down[k]) / 4 ;
                            copy[k] = buf;
                        }
                        copy[n2 -2 ] = a[i*n2 + n2 -2 ];
                        copy[n2 -1 ] = a[i*n2  + n2 - 1];
					
					
				}
        }

    }
 //   data.thr_time=get_time() - data.thr_time;

    pthread_barrier_wait(&barrier);
    return 0;


}



void print_matrix(double* m, int size1,int size2) {
        for (int y = 0; y < ((size1 < 30) ? size1 : 30) ; y++) {
                for (int x = 0; x < ((size2 < 30) ? size2 : 30); x++) {
				cout << m[y * size2 + x] << "\t";
			}
		cout << endl;
	}
}



void formula_matr(double* a, int n1, int n2) {
	for (int i = 0; i < n1; i++) {
		double* p = a + i * n2;
		for (int k = 0; k < n2; k++) {
			// p[k] = abs(i - k);
                    //    p[k] = int(-5 + Random() * 10);
                        p[k] = fabs(i - k);
                  //  p[k] = i*n2 + k;
                //        p[k] = (double) 1/( i + k + 1);
			//p[k] = n - max(i,k);
		}
	}
}

int read_matrix(double* a, int n, const string& name) {
	ifstream file(name);
	if (!file.is_open()) {
		printf("read_matrix: can\'t open the file %s \n", name.c_str());
		return ERR_CANNOT_OPEN;
	}

	for (int i = 0; i < n ; i++) {
        //	LOG(i);
		if (!(file >> a[i])) {
                	LOG(i);
			printf("read_matrix: can\'t read the file %s (invalid format)\n", name.c_str());
			return ERR_READ;
		}
	}
	return 0;
}

int init_matrix_file(double* a, int n1,int n2, const string& name) {
    int res =-1;
        if (name != "") {
                int res = read_matrix(a, n1*n2, name);
                if (res < 0) {
                        switch (res) {
                        case ERR_CANNOT_OPEN: {
                                printf("Cannot open %s \n", name.c_str());
                                break;
                        }
                        case ERR_READ: {
                                printf("Invalid data in %s \n", name.c_str());
                                break;
                        }
                        default: {
                                printf("Unknown error %d in %s n", res, name.c_str());
                                break;
                        }
                        }

                }
                return res;
        }
        return res;
}






 // double get_time(){
         // struct rusage buf;
         // getrusage(RUSAGE_THREAD,&buf);
         // return (double)buf.ru_utime.tv_sec+(double)buf.ru_utime.tv_usec/1000000.;
         // return 1;
 // }

 // double get_full_time(){
         // struct timeval buf;
         // gettimeofday(&buf,0);
         // return (double)buf.tv_sec+(double)buf.tv_usec/1000000.;
 // }



int main(int argc , char *argv[]){

    pthread_t *threads;
    Args *args;
    
    if(argc < 4){
        printf ("Invalid input data");
        return 1;
    }
        double *a;
	int n1 = 0 , n2 = 0 , p =0;
	string name = "";
    
	if (!((argc == 4) || (argc == 5)) || ((p = stoi(argv[1])) <= 0) || ((n1  = stoi(argv[2])) <= 0) || ((n2= stoi(argv[3])) <= 0) ) {
            printf("usage: %s p n1 n2 [name] \n ", argv[0]);
	    return -1;
	}

	a = new double[n1 * n2];
	if (!a) {
	    printf("Not enough memory a ");
            delete []a;
	    return -2;
	}
	
    if (argc == 5) {
        name = argv[4];
		LOG(name);
	    if (init_matrix_file(a, n1,n2, name) < 0) {
                delete [] a;
                printf("Problema 4teniya v %s \n",name.c_str());
                return 0;
	    };
	}

    if (argc == 4) {
	    formula_matr(a, n1,n2);
	}
	print_matrix(a, n1,n2);
	
        double *copyrow;
        copyrow = new double[4 * p * n2];
        if (!copyrow) {
	    printf("Not enough memory b ");
            delete []a;

	    return -2;
	}


    args = new Args[p];
    

    for(int i=0;i<p;i++){
        args[i].p =p;
        args[i].n1= n1;
        args[i].n2=n2;
        args[i].a=a;
        args[i].b=copyrow + 4*n2*i;
        args[i].num=i;
    }

    
    threads = new pthread_t[p];
    pthread_barrier_init (&barrier, nullptr, p);


  //  double time=get_full_time();

    for(int i = 0;i<p; i++){
        if(pthread_create(threads+ i, 0 , f , (void *) (args + i))){
            printf("Cant create thread %d\n",i);
            return 2;
        }
    }


    for(int i = 0; i < p;i++){
        if(pthread_join (threads[i],0)){
              printf("Cant wait thread %d\n",i);
        }
    }
       // time=get_full_time() - time;
        for (int i = 0; i < p;i++){
      //      printf("Time of %i potok : %f\n",i,args[i].thr_time);

        }
	//cout<<"TIME : "<<time/CLOCKS_PER_SEC<<endl;
	LN;
	LN;
      print_matrix(a,n1,n2);



   pthread_barrier_destroy(&barrier);

  

    delete [] copyrow;
    delete [] threads;
    delete [] args;
   delete []a;


}
