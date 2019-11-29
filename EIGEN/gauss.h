#pragma once
#include <float.h>
#include <cmath>
#include <string>
#include <fstream>
#include <cassert>
#include <iostream>
#include <cstring>
using namespace std;

#define ERR_CANNOT_OPEN -1
#define ERR_READ -2 
#define ERR_DEG_MATRIX -3   // вырожденная матрица
extern double  EPS;

#define LOG(...) std::cout<< #__VA_ARGS__<< " = " <<(__VA_ARGS__)<< "\n"
#define LN       std::cout << "\n";

//#define DEBUG 



double Random();
void generate_file(int n);
void print_matrix(double* m, int size);
void print_equa(double* a, double* b, int size);
void print_vector(int* m, int size);
void print_vector(double* m, int size);
void print_matrix_rect(double* m, int w, int h);
void matr_to_E(double* m, int size);
void matr_to_NULL(double* m, int size);
void formula_matr(double* a, int n);
int JordanInv(double* a, int size, double* b, double norma=1);
int read_matrix(double* a, int n, const string& name);
int init_matrix_file(double* a, int n, const string& name);
double norma_matr(double* a, int n);
void get_block(double* a, int aSize, int blockSize, int x, int y, double* c3);
void get_block_rect(double* a, int aSize, int blockSize, int blockW, int blockH, int x, int y, double* c3);
void push_block(double* a, int aSize, int blockSize, int x, int y, double* c3);
void push_block_rect(double* a, int aSize, int blockSize, int blockW, int blockH, int x, int y, double* c3);
void get_b_block(double* b, int blockSize, int getSize, int row, double* c3);
void push_b_block(double* b, int blockSize, int pasteSize, int row, double* c3);
void mult_matrix(double *a, double *b, double *c, int n);
void matr_sub_matr(double* a, double* b, double* c, int w, int h);
void Matr_mlt(double* a, double* b, double* c, int _i, int _r, int _j);
// get_block - учесть что size % blockCount не всегда != 0
// учесть что blockCount == 0
int Jordan_Solving_System(double* a, double* x, double* b, int size, int blockSize);
void RightSide(double* a, double* b, int n);
double norma_vec(double *b,int n);
double Residual(double *a, double *b, double *x,int n);
double Error( double *x,int n);
double residual_2(double *a, double *x, double *b, int n);
double residual_1(double *a, double *x, double *b, int n);
double residual_inf(double *a, double *x, double *b, int n);
