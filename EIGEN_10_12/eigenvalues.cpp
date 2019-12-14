#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string>
//#include <ctime>
#include <sys/time.h>
#include <sys/resource.h>
//#include<sys/sysinfo.h>

using namespace std;

double getTime(){
    struct timeval t;
    gettimeofday(&t,0);
    return (double)t.tv_sec + (double)t.tv_usec/1000000.;
}

double matrixNorm(double* matrix, int n){
    int count = 0;
    double max = 0;
    double current = 0;
    for(int j = 0; j < n; j++){
        max += fabs(matrix[j]);
    }
    for(int i = 1; i < n; i++){
        current = 0;
        for(int j = 0; j < n; j++){
            current += fabs(matrix[i*n + j]);
        }
        if(max < current){
            max = current;
        }
    }
    return max;
}

void printMatrix(double* matrix, int n){
    int maxPrint = n;
    if(n > 10)
        maxPrint = 12;
    for(int i = 0; i < maxPrint; i++){
        for(int j = 0; j < maxPrint; j++){
            printf("%lf ", matrix[i*n +j]);
        }
        printf("\n");
    }
    printf("\n");
}

double value(int i, int j, int n){
   /* int a = n / 3;
    if((i + j == 2*a - 1) && (i < 2*a*2) && (j < 2*a)) return 7;
    if((i < n-1)&& (j < n-1) && (i >= 2*a)&&(j>= 2*a)&&(i == j))return 1;
    if((i == n-1) && (i >= 2*a) && (j >= 2*a)) return j -2*a + 1;
    if(( j == n - 1) && (j >= 2*a) && ( i >= 2*a)) return i - 2*a + 1;
    return 0;*/
   /* if(i == j) return 2;
    if(abs(i - j) == 1 ) return -1;
    return 0;*/
    if (i == n - 1) return j + 1;
    if (j == n - 1) return i + 1;
    if( j == n - 1 && i == n-1) return n;
    if(i == j) return 1;
    return 0;
    //return 1.0/(i+j + 1);
    
}

int readMatrix(double* matrix,string name, int n){
    FILE *file;
    double current;
    int count = 0;
    if(name != ""){
        file = fopen(name.c_str(), "r");
        if(file == NULL){
            printf("File wasn't found\n");
            return -1;
        }
        while (fscanf(file, "%lf", &current) == 1) {
            matrix[count] = current;
            count++;
        }
        if(count != n*n){
            printf("Wrong data \n");
            return -1;
        }
    } else {
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                matrix[i*n +j] = value(i,j,n);//fabs(i-j);//value(i,j,n);//rand() % 100;//fabs(i-j);
            }
        }
    }
    return 0;
    
}

double findWillShift(double* matrix, int n, int current) {
    double shift,secondShift, firstShift,d;
    d = pow(matrix[(current -1)*n + current - 1] + matrix[current*n + current],2) - 4*( matrix[(current -1)*n + current - 1]*matrix[current*n + current]
                                                                                       - matrix[(current -1)*n + current]*matrix[(current)*n + current - 1]);
    if(d < 0){
        
        // printf("D < 0 \n");
        // printf("matrix[current*n + current - 1] is %le ,   %le\n",matrix[current*n + current - 1], matrix[current*n + current]);
        shift = matrix[(current -1)*n + current - 1] + matrix[current*n + current];
        
    } else {
        d = sqrt(d);
        if(matrix[(current -1)*n + current - 1] + matrix[current*n + current] < 0){
            firstShift = (matrix[(current -1)*n + current - 1] + matrix[current*n + current] - d ) / 2;
        } else {
            firstShift = (matrix[(current -1)*n + current - 1] + matrix[current*n + current] + d ) / 2;
            //shift = ( matrix[(current -1)*n + current - 1]*matrix[current*n + current] - matrix[(current -1)*n + current]*matrix[(current)*n + current - 1]) / shift;
        }
        secondShift = ( matrix[(current -1)*n + current - 1]*matrix[current*n + current] - matrix[(current -1)*n + current]*matrix[(current)*n + current - 1]) / firstShift;
        if(fabs(firstShift - matrix[current*n + current]) < fabs(secondShift - matrix[current*n + current])){
            shift = firstShift;
        } else {
            shift = secondShift;
        }
    }
    return shift;
    
}




int findEigenvalues(double* matrix, int n, double norm,double eps){
    double shift, d, firstRoot, secondRoot, past;
    int count = 0;
    int k = 0;
    int current = n - 1;
    
    while(current > 1){
        
        if(fabs(matrix[current*n + current - 1]) < eps*norm){
            
           // printf(" %d \n",current);
            current--;
         
            if(current == 1){
                break;
            }
            continue;
        }
        past = matrix[current*n + current - 1];
        shift = matrix[current*n + current];
      
    
        //printf("shift is %lf \n", shift);
        for(int i = 0; i <= current; i++){
            matrix[i*n + i] -= shift;
        }
       // printf("befor\n");
       // printMatrix(matrix,n);
       
        // находим L and R
        for(int i = 1; i <= current; i++){
            if(fabs(matrix[(i-1)*n + i - 1]) > eps*norm){
                matrix[i*n + i-1] = matrix[i*n + i-1] / matrix[(i-1)*n + i - 1];
            } else {
                printf("ehh \n");
                //printf("this is %lf , %d\n", matrix[(i-1)*n + i - 1],(i-1)*n + i - 1);
            }
            for(int k = i; k <= current; k++){
                matrix[i*n + k] = matrix[i*n + k] - matrix[i*n + i-1]*matrix[(i-1)*n + k];
            }
        }

      //  printf(" l R \n");
      //  printMatrix(matrix,n);
        
     // умножаем R and L
        for(int i = 0; i <= current; i++){
            for(int k = i; k < current; k++){
                matrix[i*n + k] = matrix[i*n + k] + matrix[i*n + k + 1] * matrix[(k+1)*n + k];
            }
            if(i < current){
                matrix[(i+1)*n + i] = matrix[(i+1)*n + i + 1]*matrix[(i+1)*n + i];
            }
        }
        
     //   printf("R l \n");
      //  printMatrix(matrix,n);
      //  printf("rrr is %le \n", matrix[current*n + current - 1]);
        
        for(int i = 0; i <= current; i++){
            matrix[i*n + i] += shift;
           // printMatrix(matrix,n);
        }
        if(past < matrix[current*n + current - 1]){
            k++;
        }
        
        count++;
        
    }
    if(n > 1){
       
        d = pow(matrix[0] + matrix[n + 1],2) - 4*( matrix[0]*matrix[n+1] - matrix[1]*matrix[n]);
        if(d < 0){
            printf(" D < 0 \n");
            return count;
        }
        d = sqrt(d);
        if(matrix[0] + matrix[n + 1] < 0){
            firstRoot = (matrix[0] + matrix[n + 1] - d ) / 2;
        } else {
            firstRoot = (matrix[0] + matrix[n + 1] + d ) / 2;
        }
        secondRoot = (matrix[0]*matrix[n+1] - matrix[1]*matrix[n] ) / firstRoot;
        matrix[0] = firstRoot;
        matrix[n + 1] = secondRoot;
    }
    
    return count;
    
}



void makeTriangularMatrix(double* matrix, int n, double norm, double eps){
    double cosF = 0;
    double sinF = 0;
    double root = 0;
    double sum1 = 0;
    double sum2 = 0;
    //double eps = 1e-14;
    for(int i = 1; i < n - 1; i++){
       // printf("i is %d \n",i);
        for(int j = i + 1; j < n; j++){
            if(fabs(matrix[j*n + i - 1]) < eps*norm){
               // printf("i is %d \n",i);
                continue;
            }
            //printf("i is %le \n",matrix[j*n + i - 1]);
            root = sqrt(matrix[i*n + i - 1]*matrix[i*n + i - 1] + matrix[j*n + i - 1]*matrix[j*n + i - 1]);
            if(root > eps*norm){
                cosF = matrix[ i*n + i - 1] / root;
                sinF = -matrix[j*n + i - 1] / root;
                matrix[i*n + i - 1] = root;
                matrix[j*n + i - 1] = 0;
            } else {
                continue;
            }
       
            for(int l = i; l < n; l++){
                sum1 = cosF*matrix[i*n + l] - sinF*matrix[j*n + l];
                sum2 = sinF*matrix[i*n + l] + cosF*matrix[j*n + l];
                matrix[i*n + l] = sum1;
                matrix[j*n + l] = sum2;
            }
           
            for(int l = 0; l < n; l++){
                sum1 = cosF*matrix[i + l*n] - sinF*matrix[j + l*n];
                sum2 = sinF*matrix[i + l*n] + cosF*matrix[j + l*n];
                matrix[l*n + i] = sum1;
                matrix[l*n + j] = sum2;
            }
           // printf("ok \n");
        }
    }
}

int isSymmetric(double *matrix, int n, double eps, double norm){
    for(int i = 0; i < n; i++){
        for(int j = i; j < n;j++){
            if(fabs(fabs(matrix[i*n + j]) - fabs(matrix[j*n + i])) > eps*norm){
                return -1;
            }
        }
    }
    return 1;
}





int findEigenvaluesForSymmetric(double *matrix, int n, double norm, double eps){
    double shift, d, firstRoot, secondRoot, past, firstShift, secondShift;
    int count = 0;
    //int k = 0;
    int current = n - 1;
    
    while(current > 1){
        
        if(fabs(matrix[current*n + current - 1]) < eps*norm){
            
           // printf(" %d \n",current);
            //printf("rrr is %le, eps*norm is %le,   %le \n", matrix[current*n + current - 1], eps*norm, fabs(matrix[current*n + current - 1]));
            current--;
            
            if(current == 1){
                
                break;
            }
            continue;
        }
       
       /*if(k == 0){
           shift = matrix[current*n + current - 1] /2  + matrix[current*n + current];
        } else if(k > 0 && k <= 2){
        
            
        }
        if(k > 2){
            shift = matrix[current*n + current - 1] /2  + matrix[current*n + current];
            k=1;
        }*/
        
        if(count != 0){
            if(fabs(findWillShift(matrix,n,current-1)) > eps*norm){
                if(fabs(findWillShift(matrix,n,current)/findWillShift(matrix,n,current-1) - 1) < 2){  // 2 можно юзать для 2 -1 , обычно используем 0.5
                    shift = findWillShift(matrix,n,current);
                    
                } else {//if (fabs(1 - matrix[(current-1)*n + current-1]/matrix[current*n + current]) < 1/3){
                    shift = matrix[current*n + current];
                }
            }
            else {
               
                shift = matrix[current*n + current];
               
            }
        } else {
            shift = matrix[current*n + current - 1] /2  + matrix[current*n + current] ;
        }
        past = matrix[current*n + current - 1];
        
        //printf("shift is %lf \n", shift);
        for(int i = 0; i <= current; i++){
            matrix[i*n + i] -= shift;
        }
        // printf("befor\n");
        // printMatrix(matrix,n);
        
        // находим L and R
        //printf("count is %d\n", count);
        for(int i = 1; i <= current; i++){
            if(fabs(matrix[(i-1)*n + i - 1]) > eps*norm){
                matrix[i*n + i-1] = matrix[i*n + i-1] / matrix[(i-1)*n + i - 1];
            } else {
                printf("ehh \n");
            }
            matrix[i*n + i] = matrix[i*n + i] - matrix[i*n + i-1]*matrix[(i-1)*n + i];
            if(i < current){
                matrix[i*n + i + 1] = matrix[i*n + i + 1] - matrix[i*n + i-1]*matrix[(i-1)*n + i + 1];
            }
            
        }
        
        //  printf(" l R \n");
        //  printMatrix(matrix,n);
        
        // умножаем R and L
        for(int i = 0; i < current; i++){
            matrix[i*n + i] = matrix[i*n + i] + matrix[i*n + i + 1] * matrix[(i+1)*n + i];
            if(i < current){
                matrix[(i+1)*n + i] = matrix[(i+1)*n + i + 1]*matrix[(i+1)*n + i];
            }
        }
        
        //   printf("R l \n");
        //  printMatrix(matrix,n);
        printf("rrr is %le \n", matrix[current*n + current - 1]);
        
        for(int i = 0; i <= current; i++){
            matrix[i*n + i] += shift;
            // printMatrix(matrix,n);
        }
        if(past < matrix[current*n + current - 1]){
          //  k++;
        }
        
        count++;
        //printf("count is %d\n", count);
        
    }
    if(n > 1){
        
        d = pow(matrix[0] + matrix[n + 1],2) - 4*( matrix[0]*matrix[n+1] - matrix[1]*matrix[n]);
        if(d < 0){
            printf(" D < 0 \n");
            return count;
        }
        d = sqrt(d);
        if(matrix[0] + matrix[n + 1] < 0){
            firstRoot = (matrix[0] + matrix[n + 1] - d ) / 2;
        } else {
            firstRoot = (matrix[0] + matrix[n + 1] + d ) / 2;
        }
        secondRoot = (matrix[0]*matrix[n+1] - matrix[1]*matrix[n] ) / firstRoot;
        matrix[0] = firstRoot;
        matrix[n + 1] = secondRoot;
    }
    
    return count;
}


int main(int argc, char *argv[]) {
    srand(time(0));
    int n,count;
    double* matrix, norm, sum = 0, sumOfSquares = 0, time = 0;
    double eps;
    string name = "";
    if(argc < 3){
        printf("Недостаточно данных \n");
        return -1;
    }
    n = atoi(argv[1]);
    eps = atof(argv[2]);
    if(argc == 4){
        name = argv[3];
    }
    
    matrix = new double[n*n];
    if(readMatrix(matrix, name, n) != 0) {
        return -1;
    }
    for(int i = 0; i < n; i++){
        sum += matrix[i*n + i];
        
    }
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            sumOfSquares += matrix[i*n+j] * matrix[j*n +i];
        }
    }
   
    time = getTime();
    printMatrix(matrix,n);
    norm = matrixNorm(matrix,n);
    printf("norm is %lf \n", norm);
    makeTriangularMatrix(matrix,n, norm, eps);
    
    printMatrix(matrix,n);
    if(isSymmetric(matrix,n,eps,norm) == -1) {
        printf("no\n");
        count = findEigenvalues(matrix,n,norm, eps);
        
    } else {
        printf("symmetric \n");
        count = findEigenvaluesForSymmetric(matrix,n,norm,eps);
    }
    //printMatrix(matrix,n);
    
    time = -time + getTime();
    
    
    for(int i = 0; i < n; i++){
        sum -= matrix[i*n + i];
        sumOfSquares -= matrix[i*n + i] * matrix[i*n + i];
        //printf("eigenvalue is %le \n", matrix[i*n + i]);
    }
    
    
    
    printf("Sum is %le\nSum of squares is %le \n",fabs(sum),fabs(sumOfSquares));
    
    printf("Time is %lf \n", time);
    printf("Count of stages %d \n", count);
    
    return 0;
}

