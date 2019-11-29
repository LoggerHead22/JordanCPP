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

