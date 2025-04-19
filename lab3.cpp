#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include<iostream>
#include<vector>
#define N 10000

std::ostream& operator<<(std::ostream& out, const std::vector<double>& vec) {
    for (const auto & elem : vec) 
    {
        out << elem << " ";
    }
    out << std::endl;
    return out;
}

std::ostream& operator<<(std::ostream& file, const std::vector<std::vector<double>>& macierzatko) 
{
    for (const auto& row : macierzatko) 
    {
        for (const auto& elem : row) {
            file << elem << " "; 
        }
        file << "\n"; 
    }
    return file;
}

std::vector<double> matrixmulti(const std::vector<std::vector<double>>& A, std::vector<double> &r, double b) // macierz 2x2 razy wektor
{
    std::vector<double> result(A.size(), 0);
    for (int i = 0; i < A.size(); i++)
    {
        for (int j = 0; j < A[i].size(); j++)
        {
            result[i] += A[i][j] * r[j] *b;
        }
    }
    return result;
}

std::vector<std::vector<double>> R_matrix(double alfa) // macierz R
{
    std::vector<std::vector<double>> R(2, std::vector<double>(2));
    R[0][0] = cos(alfa);
    R[0][1] = -sin(alfa);
    R[1][0] = sin(alfa);
    R[1][1] = cos(alfa);

    return R;
}

std::vector<std::vector<double>> cov(std::vector<std::vector<double>> & wektor) // macierz kowariancji
{
    double x {0};
    double y{0};
    double xy{0};
    double xx{0};
    double yy{0};

    for(int i = 0; i < wektor.size(); i++)
    {
        x += wektor[i][0];
        y += wektor[i][1];
        xy += wektor[i][0] * wektor[i][1];
        xx += wektor[i][0] * wektor[i][0];
        yy += wektor[i][1] * wektor[i][1];
    }
    std::vector<std::vector<double>> kow(2, std::vector<double>(2));
    kow[0][0] = xx / wektor.size()- (x*x) / (wektor.size()*wektor.size());
    kow[0][1] = xy / wektor.size() - x / wektor.size() * y / wektor.size();
    kow[1][0] = xy / wektor.size() - x / wektor.size() * y / wektor.size();
    kow[1][1] = yy / wektor.size() - (y*y) / (wektor.size()*wektor.size());
    
    return kow;
}

double uniform()
{
    return rand() / double (RAND_MAX);;

}

std::vector<double> box_muller(double u1, double u2)
{   
    double X = std::sqrt(-2*log(1-u1))*cos(2*M_PI*u2);
    double Y = std::sqrt(-2*log(1-u1))*sin(2*M_PI*u2);
    double r = std::sqrt(X*X + Y*Y);
    std::vector<double> arr {X, Y, r};

    return arr;

}

std::vector<double> jednorodne_kolo(double X, double Y)
{
    double X_prim = X / std::sqrt(X*X + Y*Y);
    double Y_prim = Y / std::sqrt(X*X + Y*Y);
    std::vector<double> arr {X_prim, Y_prim, std::sqrt(X_prim*X_prim + Y_prim*Y_prim)};

    return arr;
}

std::vector<double> okrag(double X, double Y, double u1)
{
    double X_bis = std::sqrt(u1)*X;
    double Y_bis = std::sqrt(u1)*Y;
    std::vector<double> arr {X_bis, Y_bis, std::sqrt(X_bis*X_bis + Y_bis*Y_bis)};

    return arr;
}



int main()
{
    std::ofstream file("zad1.csv");
    std::ofstream file2("zad2.csv");
    std::ofstream file2b("zad2b.csv");
    std::ofstream file3("zad3.csv");
    std::ofstream file4("zad4.csv");

    std::vector<double> e1 = {1, 0};
    std::vector<double> e2 = {0, 1};

    std::vector<double> re1 = matrixmulti(R_matrix(M_PI/4), e1,1);
    std::vector<double> re2 = matrixmulti(R_matrix(M_PI/4), e2,0.2);
    std::vector<std::vector<double>> A_matrix(2, std::vector<double>(2));
    for (size_t i = 0; i < e1.size(); ++i) 
    {
        A_matrix[i][0] = re1[i];
        A_matrix[i][1] = re2[i];
    }

    std::vector<std::vector<double>> kow(2, std::vector<double>(2));
    std::vector<std::vector<double>> D_matrix;
    std::vector<std::vector<double>> E_matrix;;

    for(int i = 0; i <= N; i++)
    {   
        std::vector<double> A = box_muller(uniform(), uniform());
        std::vector<double>  B = jednorodne_kolo(A[0], A[1]);
        std::vector<double> C = okrag(B[0], B[1], uniform());
        std::vector<double> D = matrixmulti(A_matrix, C, 1);
        D_matrix.push_back(D);
        std::vector<double> E = matrixmulti(A_matrix, A, 1);
        E_matrix.push_back(E);

        file<<A<<"\n";
        file2<<B<<"\n";
        file2b<<C<<"\n";
        file3<<D<<"\n";
        file4<<E<<"\n";


    }
    //std::cout<<cov(D_matrix)<<"\n";
    std::cout<<cov(D_matrix)[0][1] / sqrt(cov(D_matrix)[0][0]* cov(D_matrix)[1][1])<<"\n";
    std::cout<<cov(E_matrix)[0][1] / sqrt(cov(E_matrix)[0][0]* cov(E_matrix)[1][1])<<"\n";
    std::cout<<A_matrix<<"\n";
    file.close();
    file2.close();
    file2b.close();
    file3.close();
    file4.close();

}