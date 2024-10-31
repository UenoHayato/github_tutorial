#include <iostream>
#include <cmath>
#include <complex>
#include <math.h>
#include <bits/stdc++.h>
#define _USE_MATH_DEFINES
using namespace std;

const double PI = 3.1415926535;

int main()
{
  double ω;
  double τn; // the electron recombination lifetime
  double f;  // frequency

  complex<long double> Q;
  complex<long double> Q2;

  double Dn; // the electron diffusion coeffcient
  double Ln; // the diffusion length
  double α;  // the absorption coefficient
  double d;  // the cell thickness

  double Rs = 10.0;                  // ressitance
  double Cg = 3.6 * pow(10.0, -9.0); // geometrical capacitance

  double Real_Q; // real part
  double Imag_Q; // imagnary part

  double Real_Q2; // real part
  double Imag_Q2; // imagnary part

  char filename[256]; // ファイル名

  cout << " the electron recombination lifetime (μ*s): ";
  cin >> τn;
  cout << endl;

  cout << " the electron diffusion coeffcient (cm^2/s): ";
  cin >> Dn;
  cout << endl;

  cout << " the absorption coefficient (nm) : ";
  cin >> α;
  cout << endl;

  cout << " the cell thickness (μ*m): ";
  cin >> d;
  cout << endl;

  cout << " diffusion length (μ*m): ";
  cin >> Ln;
  cout << endl;

  cout << "Please enter a file name:";
  cin >> filename;
  cout << endl;

  double n = 0.0;

   Ln = sqrt(Dn * τn);

  τn = τn / pow(10.0, 6.0);
  α = α / pow(10.0, 4.0);
  d = d / pow(10.0, 4.0);
 Ln = Ln / pow(10.0, 4.0);

  //α = pow(α, -1.0);

  FILE *fp;
  fopen_s(&fp, filename, "w");
  fprintf(fp, "f ,Real ,Imagnary ,Real_Q2,Imag_Q2 \n");

  while (n <= 10.0)
  {

    f = pow(10.0, n);

    ω = 2 * PI * f;

    complex< double> z(1.0, ω * τn);
    complex<long double> A(1.0, ω * Rs * Cg);

    z = sqrt(z);
    A = pow(A, -1);

    Q = (1.0 - exp(-α * d) * (exp((z * d) / Ln) + (z / (α * Ln) - 1.0) * sinh((z * d) / Ln))) / ((1.0 - pow((z / (Ln * α)), 2.0)) * cosh((z * d) / Ln));

    Q2 = A * Q;

    cout << Q << endl;

    double judgeQ = abs(Q);

    Real_Q = real(Q);
    Imag_Q = imag(Q);

    Real_Q2 = real(Q2);
    Imag_Q2 = imag(Q2);

    fprintf(fp, "%lf,%lf,%lf,%lf,%lf\n", f, Real_Q, Imag_Q, Real_Q2, Imag_Q2); // CSVファイルに上書き保存

    if (judgeQ == 0)
    {
      break; // 絶対値Qが0で終了
    }

    n = n + 0.01;
  }

  fclose(fp);
  return 0;
}
