#include <cstdlib>
#include <cmath>
#include <fstream>
#include<iostream>
#include<string>
#include <algorithm>
#include <vector>
#define g1 0.8
#define g2 0.2
#define N 1000000

double f(double x)
{
    return 4.0/5.0* (1 + x - x*x*x);
}

double F(double x)
{
    return 4.0/5.0* (x + x*x/2 - x*x*x*x/4);
}

double uniform()
{
    return (double)rand()/(double)(RAND_MAX);

}

 double rozklad_zolozny(double U1, double U2)
 {
    if (U1 <= 4.0/5.0) return U2;
    else return std::sqrt(1.0 - std::sqrt(1.0 - U2));

 }

 double markow(double U1, double U2, double delta, double X_start) 
 {
    double X {};
    X = X_start + (2*U1 - 1)*delta;
    double p_acc = std::min( f(X) / f(X_start), 1.0);
   // double p_acc = std::min( (4.0/5.0* (1.0 + X - X*X*X)) / (4.0/5.0* (1.0 - X_start + X_start*X_start*X_start)), 1.0);   
    if(X <= 1 && X >= 0 && U2 <= p_acc) return X;
    else return X_start;
    
 }

double eliminacja() 
{
    while (true) 
    {
        double U1 = uniform();
        double G2 = 1.15 * uniform();
        if (G2 <= f(U1)) return U1;
    }
}

double chi_2(std::vector<double>& data) 
{
    int K = 10; //biny
    std::vector<int> probki(K, 0); //probki w kazdym binie
    std::vector<double> binowe_krance(K+1, 0);

    for (int i = 0; i <= K; i++) binowe_krance[i] = i / double(K);

    for (double x : data) 
    {
        for (int i = 0; i < K; i++) 
        {
            if (x >= binowe_krance[i] && x < binowe_krance[i + 1]) //gdzie nalezy x
            {
                probki[i]++;
                break;
            }
        }
    }
    probki[K - 1]++; 

    std::vector<double> teor(K); //przewidywane wartosci
    for (int i = 0; i < K; i++) {
        teor[i] = (F(binowe_krance[i + 1]) - F(binowe_krance[i])) * data.size();
    }

    double chi2 = 0.0;
    for (int i = 0; i < K; i++) {
        if (teor[i] > 0) {
            chi2 += (probki[i] - teor[i]) * (probki[i] - teor[i]) / teor[i];
        }
    }

    return chi2;
}

void test_chi(const std::string& plik)
{
    std::ifstream file(plik);

    std::vector<double> data;
    double value;

    while (file >> value) data.push_back(value);
    file.close();

    double chi2_stat = chi_2(data);
    
    double wartosc_krytyczna = 16.919; //dla alfa 0.05 i 9 ndf

    std::cout << "plik: "<< plik<<"\n";

    if (chi2_stat < wartosc_krytyczna) std::cout << "Nie ma podstaw aby odrzucic H0" << "chi^2: "<<chi2_stat<< "\n";
     else std::cout << "Sa podstawy do odrzucenia H0" << "chi^2: "<<chi2_stat<< "\n";
    
}

int main()
{
    srand(time(NULL));

    std::ofstream file1("rozkladzlozony.csv");
    for(int i = 0; i < N; i++)
    {
        double X = rozklad_zolozny(uniform(), uniform());

        file1<<X<<"\n"; 
        
    }
    file1.close();

    std::ofstream file2("eliminacja.csv");
    int k = 0;
    while (k < N) 
    {
        double X3 = eliminacja();
        file2 << X3 << "\n";
        k++;
    }
    file2.close();

    std::ofstream file3("markow1.csv");
    double X_start {0.5};
    for(int i = 0; i < N; i++)
    {
        double X2 = markow(uniform(), uniform(), 0.5, X_start);
        file3<<X2<<"\n";
        X_start = X2;
    }
    file3.close();

    std::ofstream file4("markow2.csv");
    for(int i = 0; i < N; i++)
    {
        double X2 = markow(uniform(), uniform(), 0.05, X_start);
        file4<<X2<<"\n";
        X_start = X2;
    }
    file4.close();



    test_chi("rozkladzlozony.csv");
    test_chi("eliminacja.csv");
    test_chi("markow1.csv");
    test_chi("markow2.csv");

    }