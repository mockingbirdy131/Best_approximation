#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <random>
#include <cmath>
#include <Eigen/Dense>
using namespace std;

const double a = -1;       
const double b = 1;   
const int K = 4;         
const int N = 3;   // узлов     
const int L = 10;       
const int M = (N-1)*K+1;
const int M_vis = 1000;

double func(double x);
void print(double value, int size, int space, ofstream &fout);
void mass(double *randd, double *x, double *y);
void deltas (double *abs, double *d, int size, ofstream &fout, const string name);
void denominator(double znam[N], double *x);
double base(double tochka, double *x, double znam[N], int i, int seg);
void block(double *randd, double *x, double znam[N], double mat[][N], int seg);
void Gram(double **gram, double *randd, double *x, double znam[N]);
void Vector_b (double *B, double *randd, double *x, double znam[N]);
void SLAU(double **gram, double *B, double *c);
double approximate(double *c, double *x, double znam[N], double tochka);

double func(double x){
    return pow(x, 3);
}

void print(double value, int size, int space, ofstream &fout){
    fout << scientific << setprecision(size) << setw(space) << right << value;
}
// генерим по L точек внутри каждого отрезка
void mass(double *randd, double *x, double *y){    
    int i, j = 0, k, begin;
    double res, delta = 1e-16;
    random_device rand;   
    mt19937 gen(rand());  
    for (k = 0; k < K; k ++){
        begin = k*(N-1);
        uniform_real_distribution<> seg(x[begin],x[begin+N-1]); 
        for (i = 0; i < L; i ++){
            res = seg(gen);
            if (res-x[begin] > delta && x[begin+N-1]-res > delta){
                randd[j] = res;
                y[j] = func(randd[j]);
                j ++;
            }
        }
    }
}

void deltas (double *abs, double *d, int size, ofstream &fout, const string name){
    int i;
    double tmp;
    double otn_l1 = 0.0, otn_l2 = 0.0, otn_oo = fabs(abs[0]);
    double abs_l1 = 0.0, abs_l2 = 0.0, abs_oo = fabs(abs[0]-d[0]);
    for (i = 0; i < size; i ++){
        tmp = fabs(abs[i]-d[i]);
        abs_l1 += tmp;
        abs_l2 += tmp * tmp;
        if (tmp > abs_oo)
            abs_oo = tmp;
        tmp = fabs(abs[i]);
        otn_l1 += tmp;              // считаем знамнатель относительной погрешности
        otn_l2 += tmp * tmp;
        if (tmp > otn_oo)
            otn_oo = tmp;
    } 
    abs_l2 = sqrt(abs_l2);
    otn_l2 = abs_l2 / sqrt(otn_l2);
    otn_l1 = abs_l1 / otn_l1;
    otn_oo = abs_oo / otn_oo;
    
    fout << "\n----- " << name << " -----\n";
    fout << "otn_l1 = "; print(otn_l1, 3, 9, fout); fout << "      abs_l1 = "; print(abs_l1, 3, 9, fout); fout << "\n"; 
    fout << "otn_l2 = "; print(otn_l2, 3, 9, fout); fout << "      abs_l2 = "; print(abs_l2, 3, 9, fout); fout << "\n"; 
    fout << "otn_oo = "; print(otn_oo, 3, 9, fout); fout << "      abs_oo = "; print(abs_oo, 3, 9, fout); fout << "\n"; 
}

// знаменатели одинаковые на подотрезках
void denominator(double znam[N], double *x){
    int i, j;
    for (i = 0; i < N; i ++){
        znam[i] = 1.0;
        for (j = 0; j < N; j ++)
            if (i != j)
                znam[i] *= (x[i]-x[j]);
    }
}

// i - локальный индекс узла на подотрезке
double base(double tochka, double *x, double znam[N], int i, int seg){
    double res = 1.0;
    int j, begin = seg*(N-1);
    for (j = 0; j < N; j ++){
        if (j != i){
            res *= (tochka-x[begin+j]);
        }
    }
    return res/znam[i];
}

void block(double *randd, double *x, double znam[N], double mat[][N], int seg){
    int i, j, m, ind;
    double phi1, phi2;
    for (i = 0; i < N; i ++){
        for (j = 0; j < N; j ++){
            mat[i][j] = 0.0;
            for (m = 0; m < L; m ++){
                ind = seg*L + m; 
                phi1 = base(randd[ind], x, znam, i, seg);
                phi2 = base(randd[ind], x, znam, j, seg);
                mat[i][j] += phi1 * phi2;
            }
        }
    }        
}

void Gram(double **gram, double *randd, double *x, double znam[N]){
    int i, j, k;
    double mat[N][N];
    for (i = 0; i < M; i ++)
        for (j = 0; j < M; j ++)
            gram[i][j] = 0.0;

    for (k = 0; k < K; k ++){
        block(randd, x, znam, mat, k);
        for (i = 0; i < N; i ++){
            for (j = 0; j < N; j++){
                gram[i+(N-1)*k][j+(N-1)*k] += mat[i][j];
            }
        }
    }
}

void Vector_b (double *B, double *randd, double *x, double znam[N]){
    int i, m, k, ind;
    double phi, f;
    for (i = 0; i < M; i ++)
        B[i] = 0.0;
    for (k = 0; k < K; k ++){
        for (i = 0; i < N; i ++){
            ind = i+(N-1)*k;
            for (m = 0; m < L; m ++){
                phi = base(randd[k*L+m], x, znam, i, k);
                f = func(randd[k*L+m]);
                B[ind] += phi*f;
            }
        }
    }
}

void SLAU(double **gram, double *B, double *c){
    int i, j;
    Eigen::MatrixXd A(M, M);
    for (i = 0; i < M; i ++)
        for (j = 0; j < M; j ++)
            A(i, j) = gram[i][j];
    Eigen::VectorXd b(M);
    for (i = 0; i < M; i ++)
        b(i) = B[i];
    // решение слау -- находим коэффициенты для наилучшего приближения
    Eigen::VectorXd x = A.fullPivLu().solve(b);
    for (i = 0; i < M; i ++)
        c[i] = x(i);
}

double approximate(double *c, double *x, double znam[N], double tochka){
    double res = 0.0;
    int i, k, begin, end;

    for (k = 0; k < K; k ++){
        begin = (N-1)*k;
        end = begin + N-1;
        if (tochka >= x[begin] && tochka <= x[end]){
            for (i = 0; i < N; i ++)
                res += c[begin+i]*base(tochka, x, znam, i, k);
            break;
        }
    }
    return res;
}

int main(){
    int i, j, size = L*K;
    double h = (b-a)/(M-1);
    double *x = new double [M], *y = new double [size], *randd = new double [size];
    
    ofstream fout("gramm.txt");
    for (i = 0; i < M; i ++)
        x[i] = a + i*h;
    mass(randd, x, y);
    double znam[N] = {0.0};
    denominator(znam, x);
    
    double *B = new double [M];
    double *c = new double [M];
    double **gram = new double* [M];
    for (i = 0; i < M; i ++)
        gram[i] = new double [M];
    Gram(gram, randd, x, znam);
    Vector_b(B, randd, x, znam);
    SLAU(gram, B, c);

    ////////////////// Погрешности //////////////////
    ofstream fout2("deltas.txt");
    int n = M + 99*K;
    double tmp, step = h/100.0;
    double *in_random = new double [size];
    double *with_step = new double [n];
    double *y_1 = new double [n];
    for (i = 0; i < size; i ++)
        in_random[i] = approximate(c, x, znam, randd[i]);
    deltas(y, in_random, size, fout2, "In random points");
    for (i = 0; i < n; i ++){
        tmp = a + i*step;
        with_step[i] = approximate(c, x, znam, tmp);
        y_1[i] = func(tmp);
    }
    deltas(y_1, with_step, n, fout2, "With step h/100");
    
    ////////////////// Вывод //////////////////
    for (i = 0; i < M; i ++){
        for (j = 0; j < M; j ++)
            print(gram[i][j], 6, 15, fout);
        fout << "\n";
    }

    ofstream fout1("data.txt");
    step = (b-a)/(M_vis-1);
    fout1 << N << "  " << M << "  " << M_vis << endl;
    // Узлы
    for (i = 0; i < M; i ++){
        print(x[i], 6, 15, fout1);
        print(func(x[i]), 6, 15, fout1);
        fout1 << endl;
    }
    // Для графика
    for (i = 0; i < M_vis; i ++){
        tmp = a + i*step;
        print(tmp, 6, 15, fout1);
        print(func(tmp), 6, 15, fout1);
        print(approximate(c, x, znam, tmp), 6, 15, fout1);
        fout1 << endl;
    }
    
    delete [] x; delete [] y; delete [] randd;
    for (i = 0; i < M; i++)
        delete[] gram[i]; 
    delete[] gram; 
    delete [] B; delete [] c;
    delete [] in_random; delete [] with_step; delete [] y_1;
    fout.close(); fout1.close(); fout2.close();
    system("python 5.py");
    return 0;
}