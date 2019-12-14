#include "gauss.h"


double Random() {
        return rand() / double(RAND_MAX);
}

void generate_file(int n) {
        ofstream file("data.txt");
        for (int i = 0; i < n * n; i++) {
                if (i % n == 0 && i != 0) {
                        file << "\n";
                }
                file << 5 + Random() * 10 << "\t";
        }
}

void print_matrix(double* m, int size) {
        for (int y = 0; y < ((size < 10) ? size : 10) ; y++) {
                for (int x = 0; x < ((size < 10) ? size : 10); x++) {
                                cout << m[y * size + x] << "\t";
                        }
                cout << endl;
        }
}

void print_equa(double* a, double* b, int size) {
        for (int y = 0; y < ((size < 10) ? size : 10); y++) {
                for (int x = 0; x < ((size < 10) ? size : 10); x++) {
                                cout << a[y * size + x] << "\t";
                }
                cout << " | " << b[y] << endl;
        }
}

void print_vector(int* m, int size) {
        for (int y = 0; y  <((size < 10) ? size : 10); y++) {
                cout << m[y] << "  ";
        }
        cout << endl;
}

void print_vector(double* m, int size) {
        for (int y = 0; y < ((size < 10) ? size : 10); y++) {
                cout << m[y] << "  ";
        }
        cout << endl;
}

void print_matrix_rect(double* m, int w, int h) {
        for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                        if (abs(m[y * w + x]) < 1e-7) {
                                cout << 0 << "\t";
                        }
                        else {
                                cout << m[y * w + x] << "\t";
                        }
                }
                cout << endl;
        }
}

void matr_to_E(double* m, int size) {
        for (int y = 0; y < size; y++) {
                for (int x = 0; x < size; x++) {
                        m[y * size + x] = x == y;
                }
        }
}

void matr_to_NULL(double* m, int size) {
        for (int y = 0; y < size * size; y++) {
                m[y] = 0;
        }
}

void formula_matr(double* a, int n) {
	   for (int x = 0; x < n; x++) {
        for (int y = 0; y < n; y++) {
              //   a[x * n + y] = int(-5 + Random() * 10);
              // a[x * n + y] = fabs(x - y);
			//	a[x * n + y] = n  - max(x,y);
                               // if (x == y) {
                                        // a[y* n + x] = 2;
                                      // } else if ((y > 0) && (y < n - 1) && (x == y + 1 || x == y - 1)) {
                                        // a[y* n + x] = -1;
                                       // } else if ((y == 0 && x == 1) || (y == n - 1 && x == n - 2)) {
                                        // a[y* n + x] = -1;
                                       // } else {
                                        // a[y* n + x] = 0;
                               // }
				
                               a[x*n + y]= (double) 1/( 1 + x+ y);
                               // if (x == n-1) { a[x*n + y]= y+1;
                               // } else if (y == n-1) { a[x*n + y]= x+1;
                               // } else if (x == y) { a[x*n + y]=1;
                               // } else { a[x*n + y]= 0; }
            }
        }
}

int JordanInv(double* a, int size, double* b, double norma) {
        matr_to_E(b, size);
        //print_matrix(a,size);
       // norma = norma_matr(a,size);
        if(norma > 1){
                norma = 1;
        }

        double *p1,*p2;
        int* columns = new int[size];
        for (int i = 0; i < size; i++) {
                columns[i] = i;
        }
//	LOG("Nachinaen Jordan");
        for (int q = 0; q < size; q++) {
                // Р С‘РЎвЂ°Р ВµР С РЎРЊР В»Р ВµР СР ВµР Р…РЎвЂљ РЎРѓ Р СР В°Р С”РЎРѓР С‘Р СР В°Р В»РЎРЉР Р…Р С•Р в„– Р Р…Р С•РЎР‚Р СР С•Р в„–
                p1 = a + q*size;
                p2 = b + q*size;
                int index = q;
                for (int j = q + 1; j < size; j++) {
                        if (fabs(p1[columns[j]]) > fabs(p1[columns[index]])) {
                                index = j;
                        }
                }
                if (fabs(p1[columns[index]]) < EPS * norma) {
                        delete[] columns;
                        return ERR_DEG_MATRIX;
                }
                swap(columns[q], columns[index]);

                // РЎС“Р СР Р…Р С•Р В¶Р С‘Р С РЎРѓРЎвЂљРЎР‚Р С•Р С”РЎС“ Р Р…Р В° 1 / a_qq
                double k = 1 / a[q * size + columns[q]];
                p1[columns[q]]=1;
                for (int j = q+1; j < size; j++) {
                        p1[columns[j]] *= k;
                }
                for (int j = 0; j < size; j++) {
                        p2[j] *= k;
                }

                // Р Р†РЎвЂ№РЎвЂЎРЎвЂљР ВµР С
                for (int i = 0; i < size; i++) {
                        if (i == q) {
                                continue;
                        }

                        double k = a[i * size + columns[q]];
                        for (int j = q ; j < size; j++) {
                                a[i * size + columns[j]] -= p1[columns[j]] * k;
                        }
                        for (int j = 0; j < size; j++) {
                                b[i * size + j] -=p2[j] * k;
                        }
                }
        }


        for (int i = 0; i < size; i++) {
                p1 = a + columns[i] * size;
                p2 = b + i*size;
                for (int j = 0; j < size; j++) {
                        p1[j] = p2[j];
                }
        }

        delete[] columns;
        return 0;
}

int read_matrix(double* a, int n, const string& name) {

    ifstream file(name);
	
    if (!file.is_open()) {
		
            printf("read_matrix: can\'t open the file %s \n", name.c_str());
            return ERR_CANNOT_OPEN;
    }
	
    for (int y = 0; y < n; y++) {
        for (int x = 0; x < n; x++) {
//				LOG(x);
//				LOG(y);
            if (!(file >> a[x * n + y])) {
				
                printf("read_matrix: can\'t read the file %s (invalid format)\n", name.c_str());
                return ERR_READ;
            }
		
        }
    }
    return 0;
}

int init_matrix_file(double* a, int n, const string& name) {
    int res =-1;
        if (name != "") {
                int res = read_matrix(a, n, name);
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

double norma_matr(double* a, int n) {
        double ret = fabs(a[0]);
        for (int y = 0; y < n; y++) {
                double sum = fabs(a[y * n]);
                for (int x = 1; x < n; x++) {
                        sum += fabs(a[y * n + x]);
                }
                if (sum > ret) {
                        ret = sum;
                }
        }
        return ret;
}

void get_block(double* a, int aSize, int blockSize, int x, int y, double* c3) {
//	printf("Get block %ix%i with size %i\n", y, x, blockSize);
#ifdef DEBUG
        if (aSize <= 0 || blockSize <= 0 || x >= aSize || y >= aSize) {
                std::cout << "get_block: invalid args\n";
                exit(-1);
        }
#endif
        x *= blockSize;
        y *= blockSize;
        double* p;
        for (int i = 0; i < blockSize; i++) {
                p = a + (y + i) * aSize + x;
                for (int k = 0; k < blockSize; k++) {
                        c3[i * blockSize + k] = p[k];
                }
        }
}

void get_block_rect(double* a, int aSize, int blockSize, int blockW, int blockH, int x, int y, double* c3) {
//	printf("Get block %ix%i with size %ix%i (blockSize=%i)\n", y, x, blockW, blockH, blockSize);
        if (aSize <= 0 || blockW <= 0 || blockH <= 0 || x >= aSize || y >= aSize) {
                std::cout << "get_block_rect: invalid args\n";
                exit(-1);
        }
        x *= blockSize;
        y *= blockSize;
        double* p;
        for (int i = 0; i < blockH; i++) {
                p = a + (y + i) * aSize + x;
                for (int k = 0; k < blockW; k++) {
                        c3[i * blockSize + k] = p[k];
                }
        }
}

void push_block(double* a, int aSize, int blockSize, int x, int y, double* c3) {
//	printf("Push block %ix%i with size %i\n", y, x, blockSize);
        if (aSize <= 0 || blockSize <= 0 || x >= aSize || y >= aSize) {
                std::cout << "push_block: invalid args\n";
                exit(-1);
        }
        x *= blockSize;
        y *= blockSize;
        double* p;
        for (int i = 0; i < blockSize; i++) {
                p = a + (y + i) * aSize + x;
                for (int k = 0; k < blockSize; k++) {
                        p[k] = c3[i * blockSize + k];
                }
        }
}

void push_block_rect(double* a, int aSize, int blockSize, int blockW, int blockH, int x, int y, double* c3) {
//	printf("Push block %ix%i with size %ix%i (blockSize=%i)\n", y, x, blockW, blockH, blockSize);
        if (aSize <= 0 || blockW <= 0 || blockH <= 0 || x >= aSize || y >= aSize) {
                std::cout << "get_block_rect: invalid args\n";
                exit(-1);
        }
        x *= blockSize;
        y *= blockSize;
        double* p;
        for (int i = 0; i < blockH; i++) {
                p = a + (y + i) * aSize + x;
                for (int k = 0; k < blockW; k++) {
                        p[k] = c3[i * blockSize + k];
                }
        }
}

void get_b_block(double* b, int blockSize, int getSize, int row, double* c3) {
        if (getSize <= 0) {
                std::cout << "push_block: invalid args\n";
                exit(-1);
        }
        double* p = b + row * blockSize;
        for (int k = 0; k < getSize; k++) {
                c3[k * blockSize] = p[k];
        }
}

void push_b_block(double* b, int blockSize, int pasteSize, int row, double* c3) {
        if (pasteSize <= 0) {
                std::cout << "push_block: invalid args\n";
                exit(-1);
        }
        double* p = b + row * blockSize;
        for (int k = 0; k < pasteSize; k++) {
                p[k] = c3[k * blockSize];
        }
}

void mult_matrix(double *a, double *b, double *c, int n){
        int i,j,k;
        double *pc=c,*pa=a,*pb=b,sum[9];

        for(i=0;i<n;i++)
                for(j=0;j<n;j++) c[i*n+j]=0.;

        for(i=0;i<n;i+=3)
                for(j=0;j<n;j+=3){
                        sum[0]=sum[1]=sum[2]=sum[3]=sum[4]=sum[5]=sum[6]=sum[7]=sum[8]=0.;
                        for(k=0;k<n;k++){
                                pa=a+i*n+k;
                                pb=b+k*n+j;
                                sum[0]+=pa[0]*pb[0];
                                sum[1]+=pa[0]*pb[1];
                                sum[2]+=pa[0]*pb[2];
                                sum[3]+=pa[n]*pb[0];
                                sum[4]+=pa[n]*pb[1];
                                sum[5]+=pa[n]*pb[2];
                                sum[6]+=pa[2*n]*pb[0];
                                sum[7]+=pa[2*n]*pb[1];
                                sum[8]+=pa[2*n]*pb[2];
                        }
                        pc=c+i*n+j;
                        pc[0]			=sum[0];
                        pc[1]			=sum[1];
                        pc[2]			=sum[2];
                        pc[n]			=sum[3];
                        pc[n+1]			=sum[4];
                        pc[n+2]			=sum[5];
                        pc[2*n]			=sum[6];
                        pc[2*n+1]		=sum[7];
                        pc[2*n+2]		=sum[8];
                }
}

void matr_sub_matr(double* a, double* b, double* c, int w, int h) {
        for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                        c[y * w + x] = a[y * w + x] - b[y * w + x];
                }
        }
}

void Matr_mlt(double* a, double* b, double* c, int _i, int _r, int _j) {
    for (int i = 0; i < _i; i++) {
                //LOG(i);
        for (int j = 0; j < _j; j++) {
            c[i * _j + j] = 0;
            for (int r = 0; r < _r; r++) {
                        //	LOG(a[i * _r + r]);
                        //	LOG(b[r * _j + j]);
                c[i * _j + j] += a[i * _r + r] * b[r * _j + j];
                                //LOG(c[i * _j + j]);
            }
        }
    }
}

int Jordan_Solving_System(double* a, double* x, double* b, int size, int blockSize) {
    assert(blockSize > 0 && blockSize % 3 == 0);
    const int blockCount = size / blockSize;
    const int m = blockSize;
    const int p = size % blockSize;
    const int l = p;
    const int k = blockCount;
    int ErrorInv = 0, minNormaCol = -1;

    double* C1 = new double[blockSize * blockSize];
    double* C2 = new double[blockSize * blockSize];
    double* C3 = new double[blockSize * blockSize];
    double aNorma = norma_matr(a, size);

    int* columns = new int[blockCount];
    for (int i = 0; i < blockCount; i++) {
        columns[i] = -1;
    }

    for (int row = 0; row < blockCount; row++) {

        double minNorma = DBL_MAX;
        minNormaCol = -1;
        //print_matrix(C2,m);
        //LOG(minNorma);

        //LOG("Ishem Normu");
        // ###############################################################
        //LOG("NNNOOORMMMMIIII STROKI");
        //LOG("row");
        for (int col = 0; col < blockCount; col++) {
            if(columns[col]!=-1) continue;
            get_block(a, size, blockSize, col, row, C2);
            ErrorInv = JordanInv(C2, blockSize, C1, aNorma);


            if (ErrorInv == 0) {
                double n = norma_matr(C2, blockSize);

                if (n < minNorma) {
                    minNormaCol = col;
                    minNorma = n;
                }
            }
        }
        //LOG(minNorma);
       
        //LOG(minNorma);
        if (minNormaCol < 0) {
            delete[] C1;
            delete[] C2;
            delete[] C3;
            delete[] columns;
            cout << "Matrix is Degenerate" << endl;
            return ERR_DEG_MATRIX;
        }
        // printf("Column with min norm: %i\n", minNormaCol);
        columns[minNormaCol] = row;

        //if (minNormaCol != blockCount - 1) {
        get_block(a, size, blockSize, minNormaCol, row, C2);
        ErrorInv = JordanInv(C2, blockSize, C1, aNorma);
        // }
        // cout << "##################################" << ErrorInv << "##############" << endl;
        //LOG("OBRATNAYA");
        //print_matrix(C2,m);
        // cout << "#########################################" << endl;



        //C2 - obratnaya (row,row)
        //#############################################NORMIRUEM STROKU#############################
        //matr_to_E(C1, blockSize);
        //push_block(a, size, blockSize, columns[row], row, C1);

        for (int col = 0; col < blockCount; col++) {
            if(columns[col]!=-1) continue;
            get_block(a, size, blockSize, col, row, C1);   // C1 = A[mxm]{0, i}
            mult_matrix(C2, C1, C3, blockSize);  // C3 = C2 * C1
            push_block(a, size, blockSize, col, row, C3);  // C3 = A[mxm]{0, i}
        }
        //cout<<size<<" "<<blockSize<<" "<<size % blockSize;
        if(l>0){
            matr_to_NULL(C1,blockSize);
            get_block_rect(a, size, blockSize, size % blockSize, blockSize, blockCount, row, C1);   // C1 = A[mxp]{0, k}
            mult_matrix(C2, C1, C3, blockSize);                           // C3[mxp] = C2 * C1
            push_block_rect(a, size, blockSize, size % blockSize, blockSize, blockCount, row, C3);
        }
        //// print_matrix_rect(C3,l,m);
        //// print_matrix(a,size);

        matr_to_NULL(C1,blockSize);
        get_b_block(b, blockSize, blockSize, row, C1);
        // LOG("C2");
        // print_matrix(C2,m);
        // LOG("C1");
        // print_matrix(C1,m);
        mult_matrix(C2, C1, C3, m);
        // LOG("C3");
        // print_matrix(C3,m);
        push_b_block(b, blockSize, blockSize, row, C3);
        // LOG("NORMIROVKASTROKI");
        // print_equa(a,b,size);
        // push_block_rect(a, size, blockSize, size % blockSize, blockSize, blockCount, row, C3);
        // Matr_mlt(C2, b + row * blockSize, C3, blockSize, blockSize, 1);
        // push_b_block(b, blockSize, blockSize, row, C3);
        //print_equa(a,b,size);
        //######################OBNULAYEM COLONKI##################################################
        for (int i = 0; i < blockCount; i++) {
            if (i == row) {
                continue;
            }
            get_block(a, size, blockSize, minNormaCol, i, C1);         // C1 = A[mxm]{i, 0}
            //matr_to_NULL(C3, blockSize);
            //push_block(a, size, blockSize, columns[row], i, C3);
            for (int j = 0; j < blockCount; j++) {
                if(columns[j]!=-1) continue;
                get_block(a, size, blockSize, j, row, C2);     // C2 = A[mxm]{0, j}
                mult_matrix(C1, C2, C3, blockSize);  // C3 = C2 * C1         // TODO: РЎРѓР В»Р ВµР Р†Р В° Р С‘Р В»Р С‘ РЎРѓР С—РЎР‚Р В°Р Р†Р В°?
                get_block(a, size, blockSize, j, i, C2);       // C2 = A[mxm]{i, j}
                matr_sub_matr(C2, C3, C2, blockSize, blockSize);        // C2 = C2 - C3
                push_block(a, size, blockSize, j, i, C2);
            }

            if(l>0){
                //get_block(a, size, m, columns[row], i, C1);
                matr_to_NULL(C2, blockSize);           // C1 = A[mxm]{i, 0}
                get_block_rect(a, size, blockSize, l, m, k, row, C2);
                mult_matrix(C1, C2, C3, m);

                get_block_rect(a, size, blockSize, l, m, k, i, C2);     // C1 = A[mxl]{i, k}
                matr_sub_matr(C2, C3, C2, m, m);                        // C2 = C1 - C3
                push_block_rect(a, size, blockSize, l, m, k, i, C2);
            }
            // LOG("DO");
            //// print_equa(a,b,size);
            matr_to_NULL(C2, blockSize);
            get_b_block(b, blockSize, blockSize, row, C2);
            mult_matrix(C1, C2, C3, m);
            get_b_block(b, blockSize, blockSize, i, C2);
            matr_sub_matr(C2, C3, C2, m, m);
            push_b_block(b, blockSize, blockSize, i, C2);
            //// print_matrix(a,size);
            // LOG("POSLE");
            //print_equa(a,b,size);
        }
        if(l>0){
            matr_to_NULL(C1, blockSize);
            get_block_rect(a, size, blockSize, m, l, minNormaCol, k, C1); //
            for (int i = 0; i < blockCount; i++) {
                if(columns[i]!=-1) continue;
                get_block(a, size, m, i, row, C2);
                mult_matrix(C1, C2, C3,m);
                get_block_rect(a, size, blockSize, m, l, i, k, C2);
                matr_sub_matr(C2, C3, C2, m, m);
                push_block_rect(a, size, blockSize, m, l, i, k, C2);
            }
            matr_to_NULL(C2, blockSize);
            get_block_rect(a, size, blockSize, l, m, k, row, C2);
            // LOG("C2");
            // print_matrix(C2,m);
            // LOG("C1");
            // print_matrix(C1,m);
            mult_matrix(C2, C1, C3, m);
            // C2 = A[mxp]{0, k}
            mult_matrix(C1, C2, C3, m);
            // LOG("C3");
            // print_matrix(C3,m);                                 // C1 = C2 * C3
            get_block_rect(a, size, blockSize, l, l, k, k, C2);             // C2 = A[mxm]{k, k}
            matr_sub_matr(C2, C3, C2, m, m);                                // C2 = C2 - C1
            push_block_rect(a, size, blockSize, p, p, k, k, C2);
            /////////////////////////////////////////////////////
            matr_to_NULL(C2, blockSize);
            get_b_block(b, blockSize, blockSize, row, C2);
            mult_matrix(C1, C2, C3,m);
            get_b_block(b, blockSize, l, k, C2);
            matr_sub_matr(C2, C3, C2, m, m);
            push_b_block(b, blockSize, l, k, C2);
            //////////////////////////////////////////////////
            //// print_matrix(a,size);
            // LOG("OBNUKIKICOLONKU");

            // LOG("Po novoi");
        }
        // LOG(row);
        // print_equa(a,b,size);

    }

    if(l>0){
        matr_to_E(C2, m);
        //matr_to_E(C1, l);

        get_block_rect(a, size, blockSize, l, l, k, k, C2);
        //LOG("DO");
        //print_matrix(C2,l);
        //print_matrix(C1,m);
        //double n = norma_matr(C2,l);
        ErrorInv = JordanInv(C2, m, C1, aNorma);
        //LOG("POSLE");
        //print_matrix(C2,l);
        if (ErrorInv < 0) {
            LOG("Matrix is Degenerate");
            delete[] C1;
            delete[] C2;
            delete[] C3;
            delete[] columns;
            return ERR_DEG_MATRIX;
        }
        // LOG("Obratnaya");
        // print_matrix(C2,m);
        matr_to_NULL(C1, m);
        get_b_block(b, blockSize, l, k, C1);
        mult_matrix(C2,C1,C3,m);
        // LOG("KUSOK");
        // print_matrix(C3,m);
        push_b_block(b, blockSize, l, k, C3);//C3 - b(k+1)
        // LOG("C3-posl");
        // print_matrix(C3,m);
        for (int i = 0; i < blockCount; i++) {
            get_block_rect(a, size, blockSize, l, m, k, i, C1);
            //		// LOG("C1-posl");
            //		// print_matrix(C1,m);
            //		// LOG("C3-posl");
            //		// print_matrix(C3,m);
            mult_matrix(C1, C3, C2, m);
            //		// LOG("C2-posl");
            //		// print_matrix(C2,m);
            get_b_block(b, blockSize, m, i, C1);
            //		// LOG("C1-posl");
            //		// print_matrix(C1,m);
            matr_sub_matr(C1, C2, C1, m, m);
            //		// LOG("C1-posl");
            //		// print_matrix(C1,m);
            push_b_block(b, m, m, i, C1);
        }
    }
  
    for (int i = 0; i < blockCount; i++) {
        for (int j = 0; j < blockSize; j++) {
            x[i* blockSize + j] = b[columns[i] * blockSize + j];
        }
    }
    if(l>0){
        for (int i = 0; i < l; i++) {
            x[blockSize * blockCount + i] = b[blockSize * blockCount + i];
        }
    }



    // print_vector(x,size);

    //print_vector(b, size);
    delete[] C1;
    delete[] C2;
    delete[] C3;
    delete[] columns;

    return 0;
}

void RightSide(double* a, double* b, int n) {
        for (int i = 0; i < n; i++) {
                double* p = a + i * n;
                b[i] = 0;
                for (int k = 0; k < n; k += 2) {
                        b[i] += p[k];
                }
        }
}
// void RightSide(double* a, double* b, int n) {
        // for (int i = 0; i < n; i++) {
                // double* p = a + i * n;
                // b[i] = 0;
                // for (int k = 0; k < n; k ++) {
                        // b[i] += k*p[k];
                // }
        // }
// }


double norma_vec(double *b,int n){
        double norm = 0;
        for(int i = 0;i < n;i++){
         if(fabs(b[i]) > norm){
                norm = fabs(b[i]);
         }
        }
        return norm;
}

double Residual(double *a, double *b, double *x,int n){
        double *c = new double[n];
        Matr_mlt(a,x,c,n ,n ,1);
        for(int i = 0;i < n;i++){
                c[i]-=b[i];
        }
        //print_vector(c,n);
        double norm = norma_vec(c,n);
        delete []c;
        return norm;
}

// double Residual(double *a, double *b, double *x,int n){


        // double norm;
        // norm = residual_inf(a,x,b,n);
        // cout<<"Norma inf:  "<<norm<<endl;
        // norm = residual_1(a,x,b,n);
        // cout<<"Norma 1:  "<<norm<<endl;
        // norm = residual_2(a,x,b,n);
        // cout<<"Norma 2:  "<<norm<<endl;


        // return norm;
// }

// double residual_inf(double *a, double *x, double *b, int n)
// {
        // double sum, max = 0;
        // for (int i=0;i<n;i++) {
                // sum = 0;
                // for (int j=0;j<n;j++)
                        // sum += a[i*n+j]*x[j];
                // sum = fabs(sum-b[i]);
                // if (sum>max)
                        // max = sum;
        // }
        // return max;
// }

// double residual_1(double *a, double *x, double *b, int n)
// {
        // double sum, res = 0;
        // for (int i=0;i<n;i++) {
                // sum = 0;
                // for (int j=0;j<n;j++)
                        // sum += a[i*n+j]*x[j];
                // res += fabs(sum-b[i]);
        // }
        // return res;
// }

// double residual_2(double *a, double *x, double *b, int n)
// {
        // double sum, res = 0;
        // for (int i=0;i<n;i++) {
                // sum = 0;
                // for (int j=0;j<n;j++)
                        // sum += a[i*n+j]*x[j];
                // res += (sum-b[i])*(sum-b[i]);
        // }
        // return sqrt(res);
// }

double Error( double *x,int n){
        for(int i = 0;i<n;i+=2 ){
                x[i]-=1;
        }
        return norma_vec(x,n);
}
