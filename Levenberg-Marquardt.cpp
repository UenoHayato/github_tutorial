#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <bits/stdc++.h>
#include <math.h>
#include <complex>
#include <cmath>
#include <algorithm>
using namespace std;

const double PI = acos(-1);
char Readfilename[256]; // ファイル名

double f[200];
double realQ[200];
double imagQ[200];
double ω[200];
double ReQrecord[200];
double ImQrecord[200];

vector<string> split(string &input, char delimiter)
{
    istringstream stream(input);
    string field;
    vector<string> result;
    while (getline(stream, field, delimiter))
    {
        result.push_back(field);
    }
    return result;
}
int main()
{
    double table[3][1200];

    // データ読み込み
    double tau_n, L_n;

    cout << "Please write the name of the data file to be read: ";
    cin >> Readfilename;
    cout << endl;

    cout << "Initial value of τn : "; // 初期値の入力
    cin >> tau_n;
    cout << endl;

    cout << "Initial value of Ln: ";
    cin >> L_n;
    cout << endl;

    ifstream ifs(Readfilename);

    string line;
    int linenum = 0;

    while (getline(ifs, line))

    {
        int i = 0;
        vector<string> strvec = split(line, ',');

        while (i < strvec.size())
        {

            table[i][linenum] = stof(strvec.at(i));

            // printf("%lf", table[i][linenum]);
            i++;
        }

        cout << endl;
        linenum++;
    }
    ifs.close();

    // 各数列(f,realQ,imagQ)に格納

    for (int n = 0; n < linenum; n++)
    {

        f[n] = table[0][n];
        realQ[n] = table[1][n];
        imagQ[n] = table[2][n];

        ω[n] = 2 * PI * table[0][n];

        ReQrecord[n] = realQ[n];
        ImQrecord[n] = imagQ[n];

        //  printf("%lf,  %lf,  %lf \n", f[n], realQ[n], imagQ[n]);
        //  cout << ReQrecord[n]<<" "<<ω[n] << endl;
    }

    // Levenberg-Marquardt法スタート
    double d = 4 * pow(10, -5); // 0 * pow(10, -7); //cm

    double α = 1.52 * pow(10, 5); //(cm)^-1//  double RC = 7.23e-05;
    double Cg = 1.3 * pow(10, -8);
    double Rc = 40;
    double RC = Cg * Rc;

    // パラメータxの初期値の設定

    double a[2] = {tau_n, L_n}; // τn,Ln

    long double Rer[linenum] = {0}; // 残差6
    long double Imr[linenum] = {0}; // 残差
    long double r[2 * linenum];
    long double totalR[124000] = {0}; // 残差二乗和

    // 1.初期値残差の計算
    for (int p = 0; p < linenum; p++)
    {

        complex<double> z(1.0, ω[p] * a[0]);
        complex<double> A(1.0, ω[p] * RC);
        complex<long double> z1;
        z = sqrt(z);
        A = pow(A, -1);

        z1 = A * ((1.0 - exp(-α * d) * (exp((z * d) / a[1]) + (z / (α * a[1]) - 1.0) * sinh((z * d) / a[1]))) / ((1.0 - pow((z / (a[1] * α)), 2.0)) * cosh((z * d) / a[1])));
        Rer[p] = real(z1);
        Imr[p] = -imag(z1);

        //   cout<< z<<endl;
    }

    // 初期値の正規化
    long double max_value = Imr[0]; // 最初の要素を最大値として初期化
    for (int i = 1; i < linenum; i++)
    {
        if (Imr[i] > max_value)
        {
            max_value = Imr[i]; // 新しい最大値を見つけたら更新
        }
    }

    for (int p = 0; p < linenum; p++)
    {
        Rer[p] = ((Rer[p]) / max_value - ReQrecord[p]);
        Imr[p] = ((Imr[p]) / max_value - ImQrecord[p]);

        r[p] = Rer[p];
        r[p + linenum] = Imr[p];
    }

    // 重み関数
    long double weightingfactor[2 * linenum];

    for (int p1 = 0; p1 < linenum; p1++)
    {
        weightingfactor[p1] = 1 / (pow(ReQrecord[p1], 2) + pow(ImQrecord[p1], 2)); // 重み関数

        weightingfactor[p1 + linenum] = 1 / (pow(ReQrecord[p1], 2) + pow(ImQrecord[p1], 2));

        // cout << weightingfactor[linenum] << endl;
    }
    // cout << weightingfactor[linenum] << endl;

    //  初期値残差2乗和の計算
    for (int p = 0; p < 2 * linenum; p++)
    {
        totalR[0] = totalR[0] + (r[p] * r[p]) * weightingfactor[p];

        //  cout << totalR[0] << endl;
    }
    // cout << totalR[0] << endl;

    // 単位行列作成
    double I[2][2] = {0};
    for (int g = 0; g < 2; g++)
    {
        I[g][g] = 1.0;
    }
    long double λ = 1;

    long double update[2];
    long double jTj = 0;
    int LoopCount = 1;

    // 7.残差の二乗和f(x)がしきい値以下なら終了する。
    // 2.ヤコビアンj(x)の計算
    while (totalR[LoopCount - 1] > 0.0000008)
    {
        long double Rightpart[2] = {0};
        long double leftpartReverse[2][2] = {0};
        long double leftpart[2][2] = {0};
        long double Jacobian[2 * linenum][2] = {0};
        long double Jacobian22[2][2] = {0};

        for (int e = 0; e < linenum; e++)
        {
            complex<double> z(1.0, a[0] * ω[e]);
            complex<double> A(1.0, ω[e] * RC);
            z = sqrt(z);
            A = pow(A, -1.0);
            //  A = 1.0;
            complex<double> J1(0, 0);
            complex<double> J2(0, 0);

            J1 =
                A * (-((1.0 / cosh((d * z) / a[1]) * ((complex<double>(0, 0.5) * d * exp((d * z) / a[1]) * ω[e]) / (a[1] * z) + (complex<double>(0, 0.5) * d * ω[e] * (-1.0 + z / (a[1] * α)) * cosh((d * z) / a[1])) / (a[1] * z) + (complex<double>(0, 0.5) * ω[e] * sinh((d * z) / a[1])) / (a[1] * α * z))) / (exp(d * α) * (1.0 - complex<double>(1.0, a[0] * ω[e]) / (pow(a[1], 2.0) * pow(α, 2.0))))) +
                     (complex<double>(0, 1.0) * ω[e] / cosh((d * z) / a[1]) * (1.0 - (exp((d * z) / a[1]) + (-1.0 + z / (a[1] * α)) * sinh((d * z) / a[1])) / exp(d * α))) / (pow(a[1], 2) * pow(α, 2.0) * pow(1.0 - complex<double>(1.0, a[0] * ω[e]) / (pow(a[1], 2) * pow(α, 2)), 2)) -
                     (complex<double>(0, 0.5) * d * ω[e] / cosh((d * z) / a[1]) * (1.0 - (exp((d * z) / a[1]) + (-1.0 + z / (a[1] * α)) * sinh((d * z) / a[1])) / exp(d * α)) * sinh((d * z) / a[1])) / (a[1] * z * (1.0 - complex<double>(1.0, a[0] * ω[e]) / (pow(a[1], 2) * pow(α, 2)))) * cosh((d * z) / a[1]));

            J2 =
                A * (-((1.0 / cosh((d * z) / a[1]) * (-((d * exp((d * z) / a[1]) * z) / pow(a[1], 2.0)) - (d * z * (-1.0 + z / (a[1] * α)) * cosh((d * z) / a[1])) / pow(a[1], 2) - (z * sinh((d * z) / a[1])) / (pow(a[1], 2) * α))) /
                       (exp(d * α) * (1.0 - complex<double>(1.0, a[0] * ω[e]) / (pow(a[1], 2) * pow(α, 2.0))))) -
                     (2.0 * complex<double>(1.0, a[0] * ω[e]) / cosh((d * z) / a[1]) * (1.0 - (exp((d * z) / a[1]) + (-1.0 + z / (a[1] * α)) * sinh((d * z) / a[1])) / exp(d * α))) / (pow(a[1], 3) * pow(α, 2) * pow(1.0 - complex<double>(1.0, a[0] * ω[e]) / (pow(a[1], 2) * pow(α, 2)), 2)) +
                     (d * z / cosh((d * z) / a[1]) * (1.0 - (exp((d * z) / a[1]) + (-1.0 + z / (a[1] * α))) * sinh((d * z) / a[1]) / exp(d * α)) * sinh((d * z) / a[1])) / (pow(a[1], 2) * (1.0 - complex<double>(1.0, a[0] * ω[e]) / (pow(a[1], 2) * pow(α, 2)))) * cosh((d * z) / a[1]));

            Jacobian[e][0] = real(J1);

            Jacobian[e + linenum][0] = -imag(J1);

            Jacobian[e][1] = real(J2);

            Jacobian[e + linenum][1] = -imag(J2);

            // cout<< z <<endl;
            // cout << J1 << " " << J2 << endl;
            // cout << Jacobian[e][0] << " " <<Jacobian[e][1]  << endl;
        }

        //  cout << Jacobian[0][0] << " " << Jacobian[0][1] << " " << Jacobian[3][0] << " " << Jacobian[3][1] << endl;

        // 3.パラメータの更新

        for (int k3 = 0; k3 < 2; k3++)
        {
            for (int k1 = 0; k1 < 2; k1++)
            {
                jTj = 0.0;

                for (int k2 = 0; k2 < 2 * linenum; k2++)
                {
                    jTj = jTj + Jacobian[k2][k3] * Jacobian[k2][k1] * weightingfactor[k2];
                    // cout << setprecision(30) << jTj << endl;
                }

                // cout << jTj << endl;
                Jacobian22[k3][k1] = jTj;
            }
        }
        //  cout << Jacobian22[0][0] << endl;

        for (int k2 = 0; k2 < 2; k2++)
        {
            for (int k1 = 0; k1 < 2; k1++)
            {

                leftpart[k2][k1] = Jacobian22[k2][k1] + λ * I[k2][k1];

                // cout <<λ<< endl;
                // cout <<Jacobian22[k2][k1] << endl;
            }
        }
        // cout << setprecision(30)<< leftpart[0][0] <<" "<< leftpart[0][1] <<" "<< leftpart[1][0] <<" "<< leftpart[1][1] << endl;

        // 逆行列にする

        leftpartReverse[0][0] = leftpart[1][1] / (leftpart[0][0] * leftpart[1][1] - leftpart[0][1] * leftpart[1][0]);
        leftpartReverse[1][1] = leftpart[0][0] / (leftpart[0][0] * leftpart[1][1] - leftpart[0][1] * leftpart[1][0]);

        leftpartReverse[0][1] = -leftpart[0][1] / (leftpart[0][0] * leftpart[1][1] - leftpart[0][1] * leftpart[1][0]);
        leftpartReverse[1][0] = -leftpart[1][0] / (leftpart[0][0] * leftpart[1][1] - leftpart[0][1] * leftpart[1][0]);

        // 逆行列の完了

        //   cout << (leftpart[0][0] * pow(10, -80)) * (leftpart[1][1] * pow(10, -80))  << " " << leftpart[0][0] <<" "<<leftpart[1][1] <<setprecision(30)<< endl;
        // cout << leftpartReverse[0][0] << " " << leftpartReverse[1][1] << " " << leftpartReverse[0][1] << " " << leftpartReverse[1][0] << endl;
        // cout<< (leftpart[0][0] * leftpart[1][1] - leftpart[0][1] * leftpart[1][0])<<endl;

        // 右辺掛け算

        for (int k1 = 0; k1 < 2; k1++)
        {
            Rightpart[k1] = 0.0;

            for (int k2 = 0; k2 < 2 * linenum; k2++)
            {
                Rightpart[k1] = Rightpart[k1] + (r[k2] * Jacobian[k2][k1]) * weightingfactor[k2];
            }
        }
        //  cout <<setprecision(20)<< Rightpart[1] << " " << Rightpart[0] << endl;

        for (int k1 = 0; k1 < 2; k1++)
        {
            update[k1] = 0;

            for (int k2 = 0; k2 < 2; k2++)
            {
                update[k1] = update[k1] - Rightpart[k2] * leftpartReverse[k1][k2];

                // cout << Rightpart[k2] << " " << leftpartReverse[k2][k1] << endl;
            }
            // cout << update[k1] << endl;
        }
        // cout << setprecision(20) << update[1] << " " << update[0] << endl;

        // パラメータの更新量の計算終了

        // 4.パラメータｘの更新

        a[0] = a[0] + update[0];
        a[1] = a[1] + update[1];

        //  cout<<setprecision(30) <<a[0]<<" "<< a[1]<< endl;

        // 5.更新したパラメータで残差二乗和f(x)を求める。

        for (int p = 0; p < linenum; p++)
        {
            complex<double> z(1.0, ω[p] * a[0]);
            complex<double> A(1.0, RC * ω[p]);
            z = sqrt(z);
            A = pow(A, -1);
            // A = 1.0;
            complex<long double> z2(0, 0);
            z2 = A * (((1.0 - exp(-α * d) * (exp((z * d) / a[1]) + (z / (α * a[1]) - 1.0) * sinh((z * d) / a[1]))) / ((1.0 - pow((z / (a[1] * α)), 2.0)) * cosh((z * d) / a[1]))));

            Rer[p] = real(z2);
            Imr[p] = -imag(z2);

            //  cout << z/a[1] << endl;
        }

        long double max_value = Imr[0];
        for (int e = 1; e < linenum; e++)
        {
            if (Imr[e] > max_value)
            {
                max_value = Imr[e]; // 新しい最大値を見つけたら更新
            }
        }
        // cout<<max_value  <<endl;
        for (int p = 0; p < linenum; p++)
        {
            Rer[p] = (Rer[p] / max_value - ReQrecord[p]);
            Imr[p] = (Imr[p] / max_value - ImQrecord[p]);

            r[p] = Rer[p];
            r[p + linenum] = Imr[p];
        }

       
        // 更新後の残差2乗和の計算
        for (int p = 0; p < 2 * linenum; p++)
        {
            totalR[LoopCount] = totalR[LoopCount] + (r[p] * r[p])* weightingfactor[p] ;

            // cout <<  Rer[p] * weightingfactor[p] * Rer[p] + Imr[p] * weightingfactor[p] * Imr[p] << endl;
            // cout << setprecision(30) << totalR[LoopCount] << endl;/ (linenum - 2 - 1)
        }

        //   cout << setprecision(30) << totalR[LoopCount] - totalR[LoopCount - 1] << endl;
        // cout << update[1] << endl;

        // 6.減衰係数λ決定方法にしたがって、λの更新

        // cout << λ << endl;
        /*         if (totalR[LoopCount - 1] >= totalR[LoopCount])
                  {
                      λ = λ * 10;
                  }

                  else
                  {
                      λ = λ / 5;
                  };*/
       // λ = 0;
        for (int p = 0; p < 2 * linenum; p++)
        {
            λ = λ + ((r[p] * r[p]) + 0.001);
        }

        LoopCount = 1 + LoopCount;
        cout << setprecision(20) << a[0] << " " << a[1] << " " << LoopCount << " " << λ << " " << totalR[LoopCount - 1] << " " << update[1] << " " << update[0] << endl;
    }
    //  cout << a[0] << " " << a[1] << endl;

    return 0;
}
