#include <iostream>
#include <string>
#include <stdio.h>
#include <ctime>
#include <math.h>
#include <random>
#include <cassert>
#include <chrono>
#include <algorithm>
#include <thread>
#include <mutex>
#include <iostream>
#include <fstream>

double N = ui->Grid_Size_X->text().toDouble(); //количество узлов
double M = ui->Grid_Size_Y->text().toDouble();

std::vector<double> a(N), b(N), c(N),f(N), alpha(N), beta(N) ; //для прогонки
std::vector<std::vector<double>> P(N, std::vector<double>(M)); //матрица давлений
std::vector<std::vector<double>> P_old(N, std::vector<double>(M));
std::vector<std::vector<double>> Kxy(N, std::vector<double>(M)); //матрица коэффициента фильтрации
std::vector<std::vector<double>> Kz(N, std::vector<double>(M)); //матрица коэффициента фильтрации
std::vector<std::vector<double>> Q(N, std::vector<double>(M)); // матрица количества проппанта
std::vector<std::vector<double>> V(N, std::vector<double>(M)); // матрица модулей скорости
std::vector<std::vector<double>> HZ(N, std::vector<double>(M)); // матрица раскрытия
std::vector<std::vector<double>> A(N, std::vector<double>(M));
std::vector<std::vector<double>> B(N, std::vector<double>(M));
double dQ= ui->DQ->text().toDouble(); //скорость удаления проппанта


double dt = ui->dt->text().toDouble();
double T = ui->Time->text().toDouble();
double T0 = ui->Time_0->text().toDouble(); //время нагнетания давления
double P0 = ui->Enter_Pressure->text().toDouble();
double Pz = ui->Rock_Pressure->text().toDouble();

double Lx0 = ui->L0_x->text().toDouble(); //начало пустоты
double LxN = ui->L0_xm->text().toDouble(); // конец пустоты
double Ly0 = ui->L0_y->text().toDouble(); //начало пустоты
double LyN = ui->L0_ym->text().toDouble(); // конец пустоты
double Propp_pore = ui->Number_Proppant_Pore->text().toDouble();

double FF = ui->Friction_Force->text().toDouble();
double Propp = ui->Number_Proppant->text().toDouble();
double Propp_max = 1 ;
double gammaz=ui->K_filtration_z->text().toDouble();   //коэффициент фильтрации через породы
double gammaxy=ui->K_filtration_0->text().toDouble();  //начальный коэффициент фильтрации по ОХ
double HZ0 = ui->HZ->text().toDouble();
double l = ui->L_size->text().toDouble();
double mult_P = ui->Pressure_mult->text().toDouble();
double mult_gamma = ui->Filtration_mult->text().toDouble();
double Viscosity = ui->Viscosity->text().toDouble();

int time_steps=round(T/dt) ;
int time_steps0=round(T0/dt) ;

double trans=0 ; // переменная переноса проппанта


double h=l/N ; //размер шага, м
double atm=101325 ; //атмосфера в паскали
double ro = 850 ; //плотность жидкости кг/м^3
double Bh = 0.0001/atm ; //объемная упругость жидкости
double nu=7000000 ; //вязкость Па*с
nu=nu*Viscosity ;


P0=P0*atm*mult_P ; //давление на входе, Па
Pz=Pz*atm*mult_P ; //давление пластовое, Па
FF=FF*atm ; //предельное давление, сдерживающее проппант
int timecount=1 ;

QVector<double> x(N), y0(N), y1(N), y2(N), x3(time_steps), y3(time_steps), y4(time_steps) ;
std::vector<double> P_enter(time_steps)  ;
std::vector<double> Q_delete(time_steps)  ;
int M05=M/2 ;

P_enter[0]=Pz ;
for (int j=1; j<time_steps; j++)  //распределение давления на входе
    {
    Q_delete[j]=0 ;
    if (j<time_steps0)
        {
        P_enter[j] = P_enter[j-1]-(Pz-P0)/time_steps0 ;
        }
    else
        {
        P_enter[j]= P0 ;
        }
    }

for (int j=0; j<M; j++)  //начальное распределение проппанта
    {
    for (int i=0; i<N; i++)
        {
        if ((i>Lx0)&(i<LxN)&(j>Ly0)&(j<LyN))
            {
            Q[i][j]=Propp_pore ;
            }
        else
            {
            Q[i][j]=Propp ;
            }
        }
    }

for (int j=0; j<M; j++)  //задание коэффициентов
    {
    for (int i=0; i<N; i++)
        {
        HZ[i][j]=HZ0 ;
        Kxy[i][j]=gammaxy*mult_gamma ;
        Kz[i][j]=gammaz*mult_gamma ;
        A[i][j]=Kxy[i][j]/(nu*Bh*(Propp+(1-Propp)*Q[i][j]));
        B[i][j]=Kz[i][j]/(ro*Bh);
        P[i][j]=Pz;
        }
    }


//конец задания начальных условий

for (double k=dt; k<=T; k+=dt)  //цикл по времени
    {
    P[0][M05]=P_enter[timecount] ;
    for (int j=0; j<M; j++)  //обновление давления на предыдущем шаге
        {
        for (int i=0; i<N; i++)
            {
            P_old[i][j]=P[i][j];
            }
        }


    for (int j=0; j<M-1; j++) //Прогон по ОХ
        {
        if (j==M05)
            {
            a[0]=0 ;
            b[0]=1 ;
            c[0]=0 ;
            f[0]=P_old[0][j] ;// P_enter[timecount]
            }
        else
            {
            a[0]=0 ;
            b[0]=1 ;
            c[0]=0 ;
            f[0]=P_old[0][j] ; // b=1,c=1,f=0 - за коэффициентами еще надо последить
            }

        for (int i=1; i<N-1; i++)
            {
            a[i]=(A[i-1][j]+A[i][j])/(2*h*h) ;
            b[i]=(A[i-1][j]+A[i][j])/(2*h*h)+(A[i+1][j]+A[i][j])/(2*h*h)+2*B[i][j]+1/dt ;
            c[i]=(A[i+1][j]+A[i][j])/(2*h*h) ;
            f[i]=2*B[i][j]*Pz+P_old[i][j]/dt ;
            }

        a[N-1]=0 ;
        b[N-1]=1 ;
        c[N-1]=0 ;
        f[N-1]=Pz ;

        alpha[0]=c[0]/b[0] ;
        beta[0]=f[0]/b[0] ;

        for (int i=1; i<N; i++)
            {
            alpha[i]=c[i]/(b[i]-a[i]*alpha[i-1]) ;
            beta[i]=(f[i]+a[i]*beta[i-1])/(b[i]-a[i]*alpha[i-1]) ;
            }
        P[N-1][j]=beta[N-1] ;

        for (int i=N-2; i>=0; i--)
            {
            P[i][j]=alpha[i]*P[i+1][j]+beta[i] ;
            }
        }   //конец прогонки по ОХ



    for (int j=0; j<M; j++)  //обновление давления на предыдущем шаге
        {
        for (int i=0; i<N; i++)
            {
            P_old[i][j]=P[i][j];
            P[i][0]=Pz;
            P[i][M-1]=Pz;
            }
        }


    for (int j=0; j<M-1; j++)  //перераспределение проппанта
        {
        for (int i=0; i<N-1; i++)
            {
            if (((P[i][j]-P[i+1][j])>FF)&(Q[i][j]<Propp_max))
                {
                trans=dQ*Q[i+1][j] ;
                Q[i+1][j]-=trans ;
                Q_delete[timecount]+=trans ;
                }
            if (((P[i][j]-P[i+1][j])<-FF)&(Q[i+1][j]<Propp_max))
                {
                trans=dQ*Q[i][j] ;
                Q[i][j]-=trans ;
                Q_delete[timecount]+=trans ;
                }
            }
        }

        for (int j=0; j<M; j++)  //задание коэффициентов
            {
            for (int i=0; i<N; i++)
                {
                A[i][j]=Kxy[i][j]/(nu*Bh*(Propp+(1-Propp)*Q[i][j]));
                B[i][j]=Kz[i][j]/(ro*Bh);
                }
            }


/*
        for (int j=1; j<N-1; j++) //Прогон по ОУ
            {
            // от 0 до M/2
            a[0]=0 ;
            b[0]=1 ;
            c[0]=0 ;
            f[0]=Pz ;

            for (int i=1; i<M05-1; i++)
                {
                a[i]=(A[j][i-1]+A[j][i])/(2*h*h) ;
                b[i]=(A[j][i-1]+A[j][i])/(2*h*h)+(A[j][i+1]+A[j][i])/(2*h*h)+2*B[j][i] ;
                c[i]=(A[j][i+1]+A[j][i])/(2*h*h) ;
                f[i]=2*B[j][i]*P_old[j][M05] ;
                }

            a[M05-1]=0 ;
            b[M05-1]=1 ;
            c[M05-1]=0 ;
            f[M05-1]=P_old[j][M05] ;

            alpha[0]=c[0]/b[0] ;
            beta[0]=f[0]/b[0] ;

            for (int i=1; i<M05; i++)
                {
                alpha[i]=c[i]/(b[i]-a[i]*alpha[i-1]) ;
                beta[i]=(f[i]+a[i]*beta[i-1])/(b[i]-a[i]*alpha[i-1]) ;
                }

            P[j][M05-1]=beta[M05-1] ;

            for (int i=M05-2; i>=0; i--)
                {
                P[j][i]=alpha[i]*P[j][i+1]+beta[i] ;
                }
            }//конец прогонки по ОУ


        P[0][M05]=P_enter[timecount] ;
        for (int j=0; j<N-1; j++) //Прогон по ОУ
            {
            // От M/2 до M
            a[M05]=0 ;
            b[M05]=1 ;
            c[M05]=0 ;
            f[M05]=P_old[j][M05] ; // b=1,c=1,f=0 - за коэффициентами еще надо последить

            for (int i=M05+1; i<M-1; i++)
                {
                a[i]=(A[j][i-1]+A[j][i])/(2*h*h) ;
                b[i]=(A[j][i-1]+A[j][i])/(2*h*h)+(A[j][i+1]+A[j][i])/(2*h*h)+2*B[j][i] ;
                c[i]=(A[j][i+1]+A[j][i])/(2*h*h)
