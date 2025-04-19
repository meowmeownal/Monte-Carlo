#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include<iostream>
#include<vector>
#include <iostream>
#include <tuple>
#include <algorithm>
#include <initializer_list>
#include <ctime>
#include<string>
#include <sstream>

#include <random>

std::mt19937 gen(time(NULL));  
std::uniform_real_distribution<double> dist(0.0, 1.0);

double uniform() 
{
    return dist(gen);
}


template<class T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& vec) {
    for (const auto & elem : vec) 
    {
        out << elem << " ";
    }
    out << std::endl;
    return out;
}

std::vector<double> kolo2d(double u1, double u2, double u3, double x_alfa, double y_alfa, double R_alfa)
{
    double X_1 = std::sqrt(-2*log(u1))*sin(2*M_PI*u2);
    double Y_1 = std::sqrt(-2*log(u1))*cos(2*M_PI*u2);

    double X = X_1 / std::sqrt(X_1*X_1 + Y_1*Y_1);
    double Y = Y_1 / std::sqrt(X_1*X_1 + Y_1*Y_1);

    double X_prim = X*std::sqrt(u3) * R_alfa + x_alfa;
    double Y_prim = Y*std::sqrt(u3) * R_alfa + y_alfa;

    std::vector<double> arr {X_prim, Y_prim, R_alfa, x_alfa, y_alfa};
    return arr;

}

bool punkt_w_kole(double x, double y, double x_alfa, double y_alfa, double R_alfa)
{
    double dx = x - x_alfa;  // x - x_alfa
    double dy = y - y_alfa;  // y - y_alfa
    return dx*dx + dy*dy <=  R_alfa* R_alfa;
}

std::vector<double> wobu(double R_alfa, double R_beta, double x_alfa, double x_beta, double y_alfa, double y_beta, int n)
{

    double sum  {0.0};
    double R2_alfa = M_PI * R_alfa*R_alfa;

    for (int i = 0; i < n; ++i) 
    {
        std::vector<double> punkt = kolo2d(uniform(), uniform(), uniform(), x_alfa, y_alfa, R_alfa);

        if (punkt_w_kole(punkt[0], punkt[1], x_beta, y_beta, R_beta)) sum += 1.0;
    }

    double mu1 = R2_alfa * (sum / n);
    double mu2 = R2_alfa * mu1;

    double sigma2 = (mu2 - mu1 * mu1) / n;
    double sigma = std::sqrt(sigma2);

    return {mu1, sigma};
}

//double C = 1/M_PI/std::sqrt(X*X + Y*Y)/std::sqrt(X*X + Y*Y)

int main()
{
    std::srand(time(NULL));
    std::vector<double> ra(2, 0);
    std::vector<double> rb (2, 0);
    double N {10000};
    double Ra = 1.0;
    double Rb = std::sqrt(2)*Ra;
    //double xA = Rb + Ra;
    double yA = 0;
    double xB = 0;
    double yB = 0;

    //----------------zad2-------------------------
    std::ofstream file1("./koloA.csv");
    std::ofstream file2("./koloB.csv");
    for(int i = 0; i < N; i++)
    {
        std::vector<double> koloA = kolo2d(uniform(), uniform(), uniform(), Rb+Ra, 0,Ra);
        file1<<koloA<<"\n";
        std::vector<double> koloB = kolo2d(uniform(), uniform(), uniform(), 0, 0,Rb);
        file2<<koloB<<"\n";
    }
    file1.close();
    file2.close();

    //----------------zad3----------------------- Generacja 1e6 punktow w kole A
    //std::vector<double> k{2,3,4,5,6};
    std::vector<std::vector<double>> A2_matrix;
    std::vector<std::vector<double>> A3_matrix;
    std::vector<std::vector<double>> B2_matrix;
    std::vector<std::vector<double>> B3_matrix;

    for(int i = 0; i < 1000000; i++)
    {
        std::vector<double> koloA2 = kolo2d(uniform(), uniform(), uniform(), Rb+ 0.5*Ra, 0,Ra);
        A2_matrix.push_back(koloA2);

        std::vector<double> koloB2 = kolo2d(uniform(), uniform(), uniform(), Rb+0.5*Ra, 0,Rb);
        B2_matrix.push_back(koloB2);

        std::vector<double> koloA3 = kolo2d(uniform(), uniform(), uniform(), 0, 0,Ra);
        A3_matrix.push_back(koloA3);

        std::vector<double> koloB3 = kolo2d(uniform(), uniform(), uniform(), 0, 0,Rb);
        B3_matrix.push_back(koloB3);
    }
    std::vector<int> n_values = {100, 1000, 10000, 100000, 1000000};
    std::vector<double> wynik(2,0);

    for (int t = 0; t < 4; ++t)
    {
        std::vector<double> kolo_beta(3);
        std::vector<double> kolo_alfa(3);

        if (t == 0)
        { 
            kolo_beta = {0, 0, Rb};
            kolo_alfa = {Rb + 1.0/2.0 * Ra, 0, Ra};
        }
        else if (t == 1) 
        { 
            kolo_beta = {0, 0, Rb};
            kolo_alfa = {0, 0, Ra};
        }
        else if (t == 2) 
        { 
            kolo_beta = {0, 0, Ra};
            kolo_alfa = {Rb + 1.0/2.0 * Ra, 0, Rb};
        }
        else if (t == 3) 
        { 
            kolo_beta = {0, 0, Ra};
            kolo_alfa = {0, 0, Rb};
        }

        std::ofstream file3("wyniki_podpunkt_" + std::to_string(t) + ".csv");
        for (int n : n_values)
        {
            wynik = wobu(kolo_alfa[2], kolo_beta[2], kolo_alfa[0], kolo_beta[0], kolo_alfa[1], kolo_beta[1], n);
            file3 << n << ", " << wynik[0] << ", " << wynik[1] << "\n";

        
        }
        file3.close();
    }
}