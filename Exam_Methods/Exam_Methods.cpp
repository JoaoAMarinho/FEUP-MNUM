#define _USE_MATH_DEFINES

#include <iostream>
#include<vector>
#include <utility>
#include <math.h>
#include <cmath>
using namespace std;

double erro = 0.000001;
int cnt;

//Functions:
double f(double x) {
    return pow(x,7)+0.5*x-0.5;
}
double fderived(double x) {
    return 4 * pow(x, 3) + 6 * pow(x, 2) - 1;
}

//função que isola o x=g(x)
double g(double x) {
    return x;
}

//Função 1 e 2 e suas derivadas parciais
double f1(double x, double y) {
    return 2 * pow(x,2) - x * y - 5 * x + 1;
}
double dxf1(double x, double y) {
    return 4 * x - y - 5;
}
double dyf1(double x, double y) {
    return -x;
}
double f2(double x, double y) {
    return x + 3 * log(x) - pow(y, 2);
}
double dxf2(double x, double y) {
    return 3 / x + 1;
}
double dyf2(double x, double y) {
    return -2 * y;
}
//Função x=g1 e y=g2 tal que as suas derivadas parciais são <=1
double g1(double x, double y) {
    return sqrt((x * y + 5 * x - 1) / 2);
}
double g2(double x, double y) {
    return sqrt(x + 3 * log(x));
}

//Fórmulas de recorrência
double x_recor(double y, double z) {
    return (7 - y - z) / 3;
}
double y_recor(double x, double z) {
    return (4 - x - 2 * z) / 4;
}
double z_recor(double x, double y) {
    return (5 - 2 * y) / 5;
}

//Function Quadratura e cubatura:
double func(double x) {
    return sqrt(1+pow(1.5*pow(M_E,1.5*x),2));
}
double func(double x, double y) {
    return pow(M_E, y-x);
}

//Funções diferenciais
//Uma equação
double diff(double x, double y) {
    return pow(x, 2) + pow(y, 2);
}
//Sistema de equações
double diff_y(double x, double y, double z) {
    return z;
}
double diff_z(double x, double y, double z) {
    return x-3*z-2*y;
}

//Função para otimização uni
double f_uni(double x) {
    return 4*pow(x,3)+2;
}

//Otimização multi
double grad(double x, double y, char op) {
    if (op == 'x'){
        return 2*x+2;
    }
    else if (op == 'y') {
        return 2*y-8;
    }
}
double f_multi(double x, double y) {
    return pow(x+1, 2) + pow(y-4, 2);
}
double Hessiana_inverted_times_grad(double x, double y, char op) {
    if (op == 'x') {
        return x+1;
    }
    else if (op == 'y') {
        return y-4;
    }
}


//---------------------------------------------------------------------------------------------------------------------------------------------
//2- Zeros de uma função real

//2.1 Isolamento das raízes
/*
Num intervalo [a,b] se f(a)*f(b)<0, então existe um número impar de zeros.
Caso f(a)*f(b)>0, então existe um número par de zeros (incluindo nenhuma).

Usar máxima e ver no intervalo (possíveis valores)
*/

//Métodos intervalares
//2.2 Método da Bissecção
/*
Critério de precisão absoluta:
|x1 −x2| ≤ ε
Consiste em parar quando o intervalo que contêm a raiz for menor que um dado valor (pequeno)
*/
void Bissection(double a, double b) {
    double m;
    do {
        m = (a + b) / 2;
        if (f(a) * f(m) < 0) {
            b = m;
        }
        else {
            a = m;
        }
    } while (abs(a - b) > erro);
    cout << "Zero value: " << m << endl;
}
/*
Apesar da sua extrema robustez, o método da bissecção tem um ponto fraco que é frequentemente ignorado: no caso de uma função descontínua
pode ocorrer um polo de primeira ordem e, então, o método da bissecção reduzirá o intervalo inicial a um pequeno intervalo na vizinhança desse polo,
sem denunciar o facto, por isso, após a terminação da iteração, devemos sempre prever um teste dos valores de f(an)
e de f(bn) antes de aceitarmos o intervalo como enquadrante de um zero.
*/

//2.3 Método da Corda
/*
Usar critério com o membro seginte - o anterior
*/
void Corda(double a, double b) {
    double m,m1;
    do {
        m = (a*f(b) - b*f(a)) / (f(b)-f(a));
        if (f(a) * f(m) < 0) {
            b = m;
        }
        else {
            a = m;
        }
        m1= (a * f(b) - b * f(a)) / (f(b) - f(a));
    } while (abs(m1-m) > erro);
    cout << "Zero value: " << m << endl;
}

//Métodos não intervalares
//2.4 Método da tangente
//Ou de Newton

//Passar o valor seguinte para xn+1=xn-f(xn)/f'(xn);
void Newton(double guess) {
    double xn, xn1=guess;
    do {
        xn = xn1;
        xn1 = xn - f(xn) / fderived(xn);
    } while (abs(xn1-xn)>erro);
    cout << "Zero value: " << xn1 << endl;
}
/*
Quando funciona bem, o método da tangente é excelente, até porque, nas vizinhanças da raiz, tende, a
cada iteração, a dobrar o número de algarismos exactos da solução. Porém, tem muitas limitações.
A existência dessas limitações faz com que o método da tangente deva ser evitado.
*/

//2.5 Método de iteração de Picard-Peano
/*
Neste método:
1- Descobrir G(x) tal que x=G(x)
2- Verificar se |G'(guess)| < 1, para convergir
3- Usar a guess se 2.
*/
void Picar_Peano(double guess) {
    double xn, xn1 = guess;
    do {
        xn = xn1;
        xn1 = g(xn);
    } while (abs(xn1 - xn) > erro);
    cout << "Zero value: " << xn1 << endl;
}

//2.6 Resolução Iterativa de Sistemas de Equações Não Lineares

/*
Newton:
1- Achar derivadas parciais
2- Usar fórmula de recorrência
3- Paragem para ambos os valores
*/
void Newton_nonLinear(double x_guess, double y_guess) {
    double xn=x_guess, yn=y_guess;
    do {
        x_guess=xn;
        y_guess=yn;

        xn = x_guess - (f1(x_guess, y_guess) * dyf2(x_guess, y_guess) - f2(x_guess, y_guess) * dyf1(x_guess, y_guess)) / (dxf1(x_guess, y_guess) * dyf2(x_guess, y_guess) - dxf2(x_guess, y_guess) * dyf1(x_guess, y_guess));
        yn = y_guess - (f2(x_guess, y_guess) * dxf1(x_guess, y_guess) - f1(x_guess, y_guess) * dxf2(x_guess, y_guess)) / (dxf1(x_guess, y_guess) * dyf2(x_guess, y_guess) - dxf2(x_guess, y_guess) * dyf1(x_guess, y_guess));


    } while (abs(xn-x_guess)>erro || abs(yn - y_guess) > erro);
    cout << "Zero value- x:" << xn << " y:" << yn <<endl;
}

/*
Picard Peano:
1- Achar g1(x,y) e g2(x,y) tal que: |(dx e dy)g1(x_guess,y_guess)| e |(dx e dy)g2(x_guess,y_guess)| são <=1
2- Substituir xn por G1(x,y) e yn por G2(x,y);
3- Paragem para ambos os valores
*/
void Picard_Peano_nonLinear(double x_guess, double y_guess) {
    double xn = x_guess, yn=y_guess;
    do {
        x_guess = xn;
        y_guess = yn;
        xn = g1(x_guess,y_guess);
        yn = g2(x_guess, y_guess);
    } while (abs(xn - x_guess) > erro || abs(yn - y_guess) > erro);
    cout << "Zero value- x:" << xn << " y:" << yn << endl;
}


//---------------------------------------------------------------------------------------------------------------------------------------------
//3- Sistemas de Equações Lineares

//3.1 Eliminação Gaussiana 
/*
Notemos por A a matriz quadrada dos coeficientes e por B e X as matrizes-colunas dos termos independentes e das incógnitas.
Fica portanto A.X=B; Se |A|!=0 então o sistema tem solução única.

Maxima: Invert(matrix(A)).B=X;
*/

//3.2 O Erro no Método de Gauss
/*
Estailidade Externa: Erro nos dados, erros nos coeficientes e nos termos independentes;
Resolver a equação: A * erroX = erroB - erroA * X (erroB=[0.5 0.5 0.5] por exemplo, erroA é matrix 3x3 com erro 0.5)
Fica então A * erroX = K
Maxima: Invert(matrix(A)).K = erroX 
Conforme os resultados ver as incógnitas mais afetadas
*/

/*
Estabilidade Interna: Erro cometido por arredondamento e truncatura.
Maxima: B-A.X = delta(Erro da máquina);
void Estabilidade_Int()
{
    deltaX = B[0] - (x * A[0][0] + y * A[0][1] + z * A[0][2]);
    deltaY = B[1] - (x * A[1][0] + y * A[1][1] + z * A[1][2]);
    deltaZ = B[2] - (x * A[2][0] + y * A[2][1] + z * A[2][2]);
}
*/

//3.5 Método de Cholesky
/*
Tendo: A*X=B
Procuremos representar A na forma do produto de uma matriz triangular inferior, L, por uma matriz triangular superior de diagonal unitária, U: A=L.U
1- Colocar L e U com tudo 0 (tamanho de A)
2- Colocar a diagonal de U com 1's (if row==collumn, então colocar 1)
3- Ajeitar o resto da matrix L e U
4- Fazer eliminação de Gauss para L*Y=B, achando Y
5- Fazer eliminação de Gauss para U*X=Y, achando o X
*/
void Khaletsky(vector<vector<double>> A, vector<double> B) {
    //Criar matrix L do tamanho de A e com tudo a 0
    vector<vector<double>> L(A.size(), vector<double>(B.size(), 0));
    //Criar matrix U do tamanho de A e com tudo a 0
    vector<vector<double>> U(A.size(), vector<double>(B.size(), 0));

    //Colocar a diagonal de U a 1
    for (int i = 0; i < U.size(); i++)
    {
        for (int j = 0; j < U[i].size(); j++)
        {
            if (i == j)
            {
                U[i][j] = 1;
            }
        }
    }

    //Ajeitar o L e o U
    for (int i = 0; i < A.size(); i++)
    {
        for (int j = 0; j < A.size(); j++)
        {
            if (i >= j)
            {
                double x = 0;
                for (int k = 0; k < j; k++)
                {
                    x += L[i][k] * U[k][j];
                }
                L[i][j] = A[i][j] - x;
            }
            else if (j > i)
            {
                double x = 0;
                for (int k = 0; k < i; k++)
                {
                    x += L[i][k] * U[k][j];
                }
                U[i][j] = (A[i][j] - x) / L[i][i];

            }

        }

    }
    //Achar y da forma: L*Y=B
    vector< double> y(B.size(), 0);
    for (int i = 0; i < L.size(); i++)
    {
        double n = 0;
        for (int j = 0; j < L[i].size(); j++)
        {
            if (i != j)
            {
                n += L[i][j] * y[j];
            }
        }
        y[i] = (B[i] - n) / L[i][i];
    }
    //Achar x da forma: U*X=Y
    vector<double> x(B.size(), 0);
    for (int i = L.size() - 1; i > -1; i--)
    {
        double n = 0;
        for (int j = L[i].size() - 1; j > -1; j--)
        {
            if (i != j)
            {
                n += U[i][j] * x[j];
            }
        }
        x[i] = (y[i] - n) / U[i][i];
    }

    cout << "x: " << x[0] << "\ny: " << x[1] << "\nz: " << x[2] << endl;
}

//3.6 Métodos iterativos
/*
Um método iterativo de resolução de sistemas de equações lineares pode ser vantajoso em relação a um método directo, quando, por exemplo, o número de equações é muito grande,
e então os erros de truncatura e arredondamento propagados ao longo da aplicação de um método directo podem não compensar os erros afectos a uma solução aproximada resultante
da aplicação de um método iterativo.
*/
//3.6.1 Método de Gauss - Jacobi
/*
Este método é uma transposição do método de Picard_Peano
1- Achar fórmulas de recorrência(Livro): exemplo de x=(d - b*y - c*z)/a  
2- Verificar se o critério de convergência é menor ou igual a 1: exemplo deltax = (|b|+|c|)/|a|, atenção aos módulos
3- Criar algoritmo substituindo: xn+1 = x_recor(xn) e os restantes
*/
void GaussJacobi(double x, double y, double z) {
    double xn=x, yn=y, zn=z;
    do {
        x = xn;
        y = yn;
        z = zn;
        xn = x_recor(y, z);
        yn = y_recor(x, z);
        zn = z_recor(x, y);
    } while (abs(xn - x) > 0.0001 | abs(yn - y) > 0.0001 | abs(zn - z) > 0.0001);
    cout << "x: " << xn << "\ny: " << yn << "\nz: " << zn << endl;
}

//3.6.2 Método de Gauss - Seidel
/*
Igual ao Gauss Jacobi mas com valores mais recentes
*/
void GaussSeidel(double x, double y, double z) {
    double xn = x, yn = y, zn = z;
    do {
        x = xn;
        y = yn;
        z = zn;
        xn = x_recor(y, z);
        yn = y_recor(xn, z);
        zn = z_recor(xn, yn);
    } while (abs(xn - x) > 0.00001 | abs(yn - y) > 0.00001 | abs(zn - z) > 0.00001);
    cout << "x: " << xn << "\ny: " << yn << "\nz: " << zn << endl;
}

//---------------------------------------------------------------------------------------------------------------------------------------------
//4 Quadratura e Cubatura

//Integração Numérica simples ou quadratura
//4.2 Regra dos Trapézios
/*
Método de 2ª Ordem

Neste método, a ideia central consiste em substituir, em cada intervalo, o arco da curva pela sua corda,
calculando em seguida a área sob a poligonal assim definida.

1- Achar o h, h=(b-a)/n e n=(b-a)/h
2- Achar o resultado, S = h/2 * (f(a)+f(b)+ somatório(2*F(a+i*h)); (somatório começa com i=1 até n)
*/
double Trapezios(double n, double a, double b) {
    //Achar h (passo)
    double h = (b - a) / n;
    //Somatório
    double s = 0;
    for (int i = 1; i < n; i++)
    {
        s += 2 * func(a + i * h);
    }
    double S = h / 2 * (func(a) + s + func(b));
    cout << "St value: " << S << endl;
    return S;
}

//4.3 Regra de Simpson
/*
Método de 4ª Ordem

Só pode ser utilizado se o n for par

Neste método, a ideia central consiste em substituir, em cada intervalo, o arco da curva pelas parábolas definidas
por cada trio de pontos consecutivos.

1- Achar o h, h=(b-a)/n e n=(b-a)/h
2- Achar o resultado, S = h/3 * (f(a)+f(b) + somatório(4*F(a+i*h) + somatório(2*F(a+i*h)); (somatório começa com i=1 até 2*n-1)/(somatório começa com i=2 até 2*n-2)
*/
double Simpson(double n, double a, double b) {
    //Achar o h (passo)
    double h = (b - a) / n;
    double s = 0;
    for (int i = 1; i < n; i++)
    {
        if (i % 2 == 0)
        {
            s += 2 * func(a + i * h);
        }
        else
        {
            s += 4 * func(a + i * h);
        }
    }
    double S = h / 3 * (func(a) + s + func(b));
    cout << "Ss value: " << S << endl;
    return S;
}

//Controlo do erro

//Quociente de Convergência (Trapezios, Simpson)
/*
Fórmula Geral: Qc= (S'-S)/(S''-S') tem de ser igual a 2^n (n==ordem do método)
Caso se verifique o valor de Qc adequado, então o h (passo) utilizado também o é.
De forma a aproximar o Qc fazer n*2 (ou seja, h/2);
*/
void Quocient_Conv(double n, double a, double b, bool trap_or_simps) {
    double q;
    if (trap_or_simps)
    {
        q = (Trapezios(2 * n, a, b) - Trapezios(n, a, b)) / (Trapezios(4 * n, a, b) - Trapezios(2 * n, a, b));
        cout << "QC trapezio: " << q << endl;
    }
    else
    {
        q = (Simpson(2 * n, a, b) - Simpson(n, a, b)) / (Simpson(4 * n, a, b) - Simpson(2 * n, a, b));
        cout << "QC simpson: " << q << endl;
    }
}

//Error (Trapezios, Simpson)
/*
Fórmula Geral: Error= (S''-S')/(2^n-1) (n==ordem do método)
Só pode ser utilizado para um Qc correto
*/
void Error(double n, double a, double b, bool trap_or_simps) {
    double e;
    if (trap_or_simps)
    {
        e = (Trapezios(4 * n, a, b) - Trapezios(2 * n, a, b)) / 3;
        cout << "Error trapezio: " << e << endl;
    }
    else
    {
        e = (Simpson(4 * n, a, b) - Simpson(2 * n, a, b)) / 15;
        cout << "Error simpson: " << e << endl;
    }
}

//Integração Numérica dupla ou cubatura
//4.6 Cubatura
/*
Trapézios
1- St= hx*hy/4* [Somatório dos valores do vértice + 2*Somatório dos valores nos pontos médios + 4*Somatório do valor do ponto médio]
Neste caso nx=2 e ny=2, de forma a dar 4 parcelas;
Se n dividir em quadrados da forma 2^n, fazer a soma das diversas malhas.
*/
double Trapezios(double x0, double x1, double y0, double y1, double nx, double ny) {
    double sum_v = 0;
    double sum_pontos_i = 0;
    //Versão 1
    vector<pair<double, double>> vertices = { make_pair(x0,y0),make_pair(x0,y1), make_pair(x1,y0), make_pair(x1,y1) };
    vector<pair<double, double>> pontos_int = { make_pair(x0,(y0 + y1) / 2),make_pair((x0 + x1) / 2,y0), make_pair((x0 + x1) / 2,y1), make_pair(x1,(y0 + y1) / 2) };
    pair<double, double> central = make_pair((x0 + x1) / 2, (y0 + y1) / 2);
    for (int i = 0; i < vertices.size(); i++)
    {
        sum_v += func(vertices[i].first, vertices[i].second);
        sum_pontos_i += func(pontos_int[i].first, pontos_int[i].second);
    }
    //Versão 2
    /*
    sum_v=func(x0,y0)+func(x0,y1)+func(x1,y0)+func(x1,y1);
    sum_pontos_i=func(x0,(y0+y1)/2)+func((x0+x1)/ 2,y0)+func((x0+x1)/ 2,y1)+func(x1,(y0+y1)/ 2);
    */
    double sum_central = func(central.first, central.second);
    double hx = (x1 - x0) / nx;
    double hy = (y1 - y0) / ny;
    return (hx * hy / 4) * (sum_v + 2 * sum_pontos_i + 4 * sum_central);
}

/*
Simpsons
1- St= hx*hy/9* [Somatório dos valores do vértice + 4*Somatório dos valores nos pontos médios + 16*Somatório do valor do ponto médio]
Neste caso nx=2 e ny=2, de forma a dar 4 parcelas;
Se n dividir em quadrados da forma 2^n, fazer a soma das diversas malhas.
*/
double Simpson(double x0, double x1, double y0, double y1, double nx, double ny) {
    double sum_v = 0;
    double sum_pontos_i = 0;
    //Versão 1
    vector<pair<double, double>> vertices = { make_pair(x0,y0),make_pair(x0,y1), make_pair(x1,y0), make_pair(x1,y1) };
    vector<pair<double, double>> pontos_int = { make_pair(x0,(y0+y1) / 2),make_pair( (x0+x1)/ 2,y0), make_pair( (x0+x1)/ 2,y1), make_pair(x1,(y0+y1)/ 2) };
    pair<double, double> central = make_pair( (x0+x1) / 2, (y0+y1) / 2);
    for (int i = 0; i < vertices.size(); i++)
    {
        sum_v += func(vertices[i].first, vertices[i].second);
        sum_pontos_i += func(pontos_int[i].first, pontos_int[i].second);
    }
    //Versão 2
    /*
    sum_v=func(x0,y0)+func(x0,y1)+func(x1,y0)+func(x1,y1);
    sum_pontos_i=func(x0,(y0+y1)/2)+func((x0+x1)/ 2,y0)+func((x0+x1)/ 2,y1)+func(x1,(y0+y1)/ 2);
    */
    double sum_central = func(central.first, central.second);
    double hx = (x1 - x0) / nx;
    double hy = (y1 - y0) / ny;
    return ((hx * hy / 9) * (sum_v + 4 * sum_pontos_i + 16 * sum_central));
}

//Quociente de Convergência (Trapezios integral duplo, Simpson integral duplo)
//Igual ao quadratura mas aceita 2 n's
void Quocient_Conv(double x0, double x1, double y0, double y1, double nx, double ny, bool trap_or_simps) {
    double q;
    if (trap_or_simps)
    {
        q = (Trapezios(x0, x1, y0, y1, 2 * nx, 2 * ny) - Trapezios(x0, x1, y0, y1, nx, ny)) / (Trapezios(x0, x1, y0, y1, 4 * nx, 4 * ny) - Trapezios(x0, x1, y0, y1, 2 * nx, 2 * ny));
        cout << "QC trapezio: " << q << endl;
    }
    else
    {
        q = (Simpson(x0, x1, y0, y1, 2 * nx, 2 * ny) - Simpson(x0, x1, y0, y1, nx, ny)) / (Simpson(x0, x1, y0, y1, 4 * nx, 4 * ny) - Simpson(x0, x1, y0, y1, 2 * nx, 2 * ny));
        cout << "QC simpson: " << q << endl;
    }
}
//Erro
//Igual ao quadratura mas aceita 2 n's (Trapezios integral duplo, Simpson integral duplo)
void Error(double x0, double x1, double y0, double y1, double nx, double ny, bool trap_or_simps) {
    double e;
    if (trap_or_simps)
    {
        e = (Trapezios(x0, x1, y0, y1, 4 * nx, 4 * ny) - Trapezios(x0, x1, y0, y1, 2 * nx, 2 * ny)) / 3;
        cout << "Error trapezio: " << e << endl;
    }
    else
    {
        e = (Simpson(x0, x1, y0, y1, 4 * nx, 4 * ny) - Simpson(x0, x1, y0, y1, 2 * nx, 2 * ny)) / 15;
        cout << "Error simpson: " << e << endl;
    }
}

//---------------------------------------------------------------------------------------------------------------------------------------------
//5 Integração de equações diferenciais ordinárias

//5.2 Método de Euler
/*
Método de 1ª Ordem (Qc=2 e E = S´´- S')

Este método corresponde em termos analíticos, a utilizar a fórmula dos acréscimos finitos, isto é, um
desenvolvimento em série de Taylor limitado à primeira ordem.

1- Achar f(x,y), corresponde à derivada de y
2- Algoritmo:
{xn+1 = xn+h;
{yn+1 = yn+h*f(x,y);
*/
double Euler(double x0, double y0, double xf, double h) {
    double derived = 0;
    while (abs(xf - x0) > 0.000001) {
        derived = diff(x0, y0);
        x0 += h;
        y0 = y0 + h * derived;
    }
    return y0;
}

//5.4 Métodos de Runge - Kutta

//5.4.1 Método de Runge - Kuta de Segunda Ordem
/*
No caso mais simples, um método de Runge-Kutta de segunda ordem, muitas vezes designado por RK2,
pode visualiza-se do seguinte modo: (ver livro)

1-Achar y' ou f(x,y), é o mesmo
2-Achar f(x+h/2,y+y'*h/2)
3- Algoritmo:
{xn+1 = xn+h;
{yn+1 = yn+h*f(x+h/2,y+f(x,y)*h/2);
*/
double RK2(double x0, double y0, double xf, double h) {
    double incremento = 0;
    while (abs(xf - x0) > 0.000001) {
        incremento = diff(x0 + h / 2, y0 + diff(x0, y0) * h / 2);
        x0 += h;
        y0 += h * incremento;
    }
    return y0;
}

//5.4.2 Método de Runge - Kuta de Quarta Ordem
/*
É possível construir, mediante um raciocínio semelhante, um método de Runge-Kutta de terceira ordem,
mas o interesse é pequeno porque o trabalho de cálculo envolvido é pouco menor que o de um método
de Runge-Kutta de quarta ordem, conhecido por RK4.

1-Achar y' ou f(x,y), é o mesmo
2-Calcular o delta1=h*f(x,y);
3-A partir de delta1, achar delta2=h*f(x+h/2 , y+delta1/2);
4-Achar delta3=h*f(x+h/2 , y+delta2/2);
5-Achar delta4=h*f(x+h , y+delta3);
6- Achar incremento de y = delta1/6 + delta2/3 + delta3/3 + delta4/6;
7-Algoritmo
{xn+1= xn+h
{yn+1 = yn+incremento
*/
double RK4(double x0, double y0, double xf, double h) {
    double delta1, delta2, delta3, delta4;
    while (abs(xf-x0)>0.000001) {
        delta1 = h * diff(x0, y0);
        delta2 = h * diff(x0 + h / 2, y0 + delta1 / 2);
        delta3 = h * diff(x0 + h / 2, y0 + delta2 / 2);
        delta4 = h * diff(x0 + h, y0 + delta3);
        x0 += h;
        y0 += delta1 / 6 + delta2 / 3 + delta3 / 3 + delta4 / 6;
    }
    return y0;
}

//Quociente de Convergência (Euler, RK2, RK4)
void Quocient_Conv(double x0, double y0, double xf, double h, string op) {
    double q;
    if (op == "e")
    {
        q = (Euler(x0, y0, xf, h / 2) - Euler(x0, y0, xf, h)) / (Euler(x0, y0, xf, h / 4) - Euler(x0, y0, xf, h / 2));
        cout << "QC euler: " << q << endl;
    }
    else if (op == "rk2")
    {
        q = (RK2(x0, y0, xf, h / 2) - RK2(x0, y0, xf, h)) / (RK2(x0, y0, xf, h / 4) - RK2(x0, y0, xf, h / 2));
        cout << "QC rk2: " << q << endl;
    }
    else if (op == "rk4")
    {
        q = (RK4(x0, y0, xf, h / 2) - RK4(x0, y0, xf, h)) / (RK4(x0, y0, xf, h / 4) - RK4(x0, y0, xf, h / 2));
        cout << "QC rk4: " << q << endl;
    }
}
//Erro(Euler, RK2, RK4)
void Error(double x0, double y0, double xf, double h, string op) {
    double e;
    if (op == "e")
    {
        e = Euler(x0, y0, xf, h / 4) - Euler(x0, y0, xf, h / 2);
        cout << "Error euler: " << e << endl;
    }
    else if (op == "rk2")
    {
        e = (RK2(x0, y0, xf, h / 4) - RK2(x0, y0, xf, h / 2)) / 3;
        cout << "Error rk2: " << e << endl;
    }
    else if (op == "rk4")
    {
        e = (RK4(x0, y0, xf, h / 4) - RK4(x0, y0, xf, h / 2)) / 15;
        cout << "Error rk4: " << e << endl;
    }
}

//5.6 Sistemas de Equações e Equações de Ordem Superior

//Sistemas de equações diferenciais:
/*
Euler : Igual ao Euler normal mas com incremento no z
*/
double Euler(double x0, double y0, double z0, double xf, double h, char op) {
    double increment_y, increment_z;
    while (abs(xf - x0) > 0.000001) {
        increment_y = diff_y(x0, y0, z0);
        increment_z = diff_z(x0, y0, z0);
        x0 += h;
        y0 += h * increment_y;
        z0 += h * increment_z;
    }
    cout << "y: " << y0 << " z: " << z0 << endl;
    if (op == 'y') return y0;
    else if (op == 'z') return z0;
    else return -1;
}
/*
RK2 : Igual ao RK2 normal mas com incremento no z
*/
double RK2(double x0, double y0, double z0, double xf, double h, char op)//Runge-Kutta 2ª ordem 
{
    double increment_y, increment_z;
    while (abs(xf - x0) > 0.000001) {
        increment_y = diff_y(x0 + h / 2, y0 + diff_y(x0, y0,z0) * h / 2, z0 + diff_z(x0, y0, z0) * h / 2);
        increment_z = diff_z(x0 + h / 2, y0 + diff_y(x0, y0, z0) * h / 2, z0 + diff_z(x0, y0, z0) * h / 2);
        x0 += h;
        y0 += h * increment_y;
        z0 += h * increment_z;
    }
    cout << "y: " << y0 << " z: " << z0 << endl;
    if (op == 'y') return y0;
    else if (op == 'z') return z0;
    else return -1;
}
/*
RK4 : Igual ao RK4 normal mas com incremento no z, os deltas vêm por ordem:
deltay_1 >> deltaz_1 >> deltay_2 ...
*/
double RK4(double x0, double y0, double z0, double xf, double h, char op)//Runge-Kutta 4ª ordem 
{
    double delta_y1, delta_y2, delta_y3, delta_y4;
    double delta_z1, delta_z2, delta_z3, delta_z4;
    while (abs(xf - x0) > 0.000001) {
        delta_y1 = h * diff_y(x0, y0, z0);
        delta_z1 = h * diff_z(x0, y0, z0);
        delta_y2 = h * diff_y(x0 + h / 2, y0 + delta_y1 / 2, z0 + delta_z1 / 2);
        delta_z2 = h * diff_z(x0 + h / 2, y0 + delta_y1 / 2, z0 + delta_z1 / 2);
        delta_y3 = h * diff_y(x0 + h / 2, y0 + delta_y2 / 2, z0 + delta_z2 / 2);
        delta_z3 = h * diff_z(x0 + h / 2, y0 + delta_y2 / 2, z0 + delta_z2 / 2);
        delta_y4 = h * diff_y(x0 + h, y0 + delta_y3, z0 + delta_z3);
        delta_z4 = h * diff_z(x0 + h, y0 + delta_y3, z0 + delta_z3);
        x0 += h;
        y0 += delta_y1 / 6 + delta_y2 / 3 + delta_y3 / 3 + delta_y4 / 6;
        z0 += delta_z1 / 6 + delta_z2 / 3 + delta_z3 / 3 + delta_z4 / 6;
    }
    cout << "y: " << y0 << " z: " << z0 << endl;
    if (op == 'y') return y0;
    else if (op == 'z') return z0;
    else return -1;
}

//Quociente de Convergência Sistema de Equações (Euler, RK2, RK4)
void Quocient_Conv(double x0, double y0, double z0, double xf, double h, string op) {
    double q;
    if (op == "e")
    {
        q = (Euler(x0, y0, z0, xf, h / 2, 'y') - Euler(x0, y0, z0, xf, h, 'y')) / (Euler(x0, y0, z0, xf, h / 4, 'y') - Euler(x0, y0, z0, xf, h / 2, 'y'));
        cout << "QC euler y: " << q << endl;
        q = (Euler(x0, y0, z0, xf, h / 2, 'z') - Euler(x0, y0, z0, xf, h, 'z')) / (Euler(x0, y0, z0, xf, h / 4, 'z') - Euler(x0, y0, z0, xf, h / 2, 'z'));
        cout << "QC euler z: " << q << endl;
    }
    else if (op == "rk2")
    {
        q = (RK2(x0, y0, z0, xf, h / 2, 'y') - RK2(x0, y0, z0, xf, h, 'y')) / (RK2(x0, y0, z0, xf, h / 4, 'y') - RK2(x0, y0, z0, xf, h / 2, 'y'));
        cout << "QC rk2 y: " << q << endl;
        q = (RK2(x0, y0, z0, xf, h / 2, 'z') - RK2(x0, y0, z0, xf, h, 'z')) / (RK2(x0, y0, z0, xf, h / 4, 'z') - RK2(x0, y0, z0, xf, h / 2, 'z'));
        cout << "QC rk2 z: " << q << endl;
    }
    else if (op == "rk4")
    {
        q = (RK4(x0, y0, z0, xf, h / 2, 'y') - RK4(x0, y0, z0, xf, h, 'y')) / (RK4(x0, y0, z0, xf, h / 4, 'y') - RK4(x0, y0, z0, xf, h / 2, 'y'));
        cout << "QC rk4 y: " << q << endl;
        q = (RK4(x0, y0, z0, xf, h / 2, 'z') - RK4(x0, y0, z0, xf, h, 'z')) / (RK4(x0, y0, z0, xf, h / 4, 'z') - RK4(x0, y0, z0, xf, h / 2, 'z'));
        cout << "QC rk4 z: " << q << endl;
    }
}
//Erro Sistema de Equações(Euler, RK2, RK4)
void Error(double x0, double y0, double z0, double xf, double h, string op) {
    double e;
    if (op == "e")
    {
        e = Euler(x0, y0, z0, xf, h / 4, 'y') - Euler(x0, y0, z0, xf, h / 2, 'y');
        cout << "Error euler y: " << e << endl;
        e = Euler(x0, y0, z0, xf, h / 4, 'z') - Euler(x0, y0, z0, xf, h / 2, 'z');
        cout << "Error euler z: " << e << endl;
    }
    else if (op == "rk2")
    {
        e = (RK2(x0, y0, z0, xf, h / 4, 'y') - RK2(x0, y0, z0, xf, h / 2, 'y')) / 3;
        cout << "Error rk4 y: " << e << endl;
        e = (RK2(x0, y0, z0, xf, h / 4, 'z') - RK2(x0, y0, z0, xf, h / 2, 'z')) / 3;
        cout << "Error rk4 z: " << e << endl;
    }
    else if (op == "rk4")
    {
        e = (RK4(x0, y0, z0, xf, h / 4, 'y') - RK4(x0, y0, z0, xf, h / 2, 'y')) / 15;
        cout << "Error rk4 y: " << e << endl;
        e = (RK4(x0, y0, z0, xf, h / 4, 'z') - RK4(x0, y0, z0, xf, h / 2, 'z')) / 15;
        cout << "Error rk4 z: " << e << endl;
    }
}

//Equações de Ordem Superior:
/*
Não se pode aplicar os métodos diretamente
Exemplo:
No caso de uma derivada de 2ª ordem: y''+ 3*y'+ 2*y= x
1- Temos de fazer uma mudança de variável - y'= z (ou seja y'' = z') de forma a obter uma equação de 1ª ordem
2- Substituir na equação diferencial dada: z' + 3*z + 2*y = x
3- Compor o sistema
{y' = z
{z' = x - 3*z - 2*y
4- Resolver o sistema normalmente
*/

//---------------------------------------------------------------------------------------------------------------------------------------------
//6 Optimização

//6.3 As técnicas concretas

//6.3.1 Pesquisa unidimensional
/*
O mais elementar problema de optimização consiste na pesquisa do extremo de uma função de uma só
variável. Os critérios analíticos são bem conhecidos, mas são numerosos os casos em que, por falta de
uma definição analítica conveniente da função objectivo, ou de algum ou alguns dos seus componentes,
há que recorrer a métodos numéricos.
*/

//6.3.2 Métodos Intervalares

/*
Regra Áurea
O método da secção áurea utiliza esta condição. Com efeito, começando com o intervalo [x1, x2] em que
se sabe estar o mínimo, escolhemos:
1- x3 é tal que x3 - x1 = A*(x2-x1)
2- x4 é tal que x4 - x1 = B*(x2-x1)
3- Se for mínimo -> se f(x3) < f(x4) então o mínimo está entre [x1,x4]
-> se f(x3) > f(x4) então o mínimo está entre [x3,x2]
4- Se for máximo trocar os sinais de comparação das funções em x3 e x4
*/
double B = (sqrt(5) - 1) / 2; //Número de ouro
double A = pow(B, 2);         //A=B^2
double Aurea(double x1, double x2, char op) { //op == 'm'(mínimo), "M" (máximo)
    double x3, x4;
    while (abs(x1 - x2) > 0.000001) {
        x3 = (x1 + A * (x2 - x1));
        x4 = (x1 + B * (x2 - x1));
        if (op == 'm') {
            if (f_uni(x3) < f_uni(x4))
            {
                x1 = x1;
                x2 = x4;
            }
            else
            {
                x1 = x3;
                x2 = x2;
            }
        }
        else if (op == 'M') {
            if (f_uni(x3) >= f_uni(x4))
            {
                x1 = x1;
                x2 = x4;
            }
            else
            {
                x1 = x3;
                x2 = x2;
            }
        }
    }
    return x1;
}

//6.3.3 Pesquisa multidimensional

/*
Método do Gradiente
Ao contrário dos métodos analíticos, que tentam ir directamente ao valor desejado, os métodos numéricos, iterativos, baseiam-se no princípio de dar sucessivos passos
descendentes até encontrar o ponto mais baixo possível.

1- xn+1 = xn - h*nabla, (h multiplicador decidido pelo operador e o nabla é o gradiente)
2- Achar o gradiente
3- Algoritmo:
Começar com um h razoável (1 por exemplo)
Se f(xn+1) < f(xn), aumentar o passo, h*2
se não voltar atrás, h/2
*/
void Gradiente(double x, double y) {
    double xn = x;
    double yn = y;
    double h = 1;
    while (true) {
        xn = x - h * grad(x, y, 'x');
        yn = y - h * grad(x, y, 'y');
        if (abs(xn - x) <= 0.01 || abs(yn - y) <= 0.01)
        {
            break;
        }
        if (f_multi(xn, yn) < f_multi(x, y))
        {
            h = h * 2;
            x = xn;
            y = yn;
        }
        else {
            h = h / 2;
        }
    };
    cout << "x: " << xn << "\ny: " << yn << endl;
}

//6.3.4 Método da quádrica
/*
Não usa h diretamente
1- Achar o gradiente
2- Achar a inversa da hessiana
3- xn+1 = xn - H^-1*nabla (H é a hessiana de f(x,y) e nabla o gradiente) (H=[d^2f/d^2x, ...])
4- Usar o maxima: hessian(f(x,y),[x,y]) obtemos a hessiana
5- invert(%).matrix(grad(x,y))
6- ratsimp(%) e obtemos então, H^-1*nabla
7- Algoritmo: Caso a função cresça (f(xn,yn) > f(x,y)) alterar o ponto inicial
*/
void Quadrica(double x, double y) {
    double xn=x;
    double yn=y;
    do {
        x = xn;
        y = yn;
        xn = x - Hessiana_inverted_times_grad(x, y, 'x');
        yn = y - Hessiana_inverted_times_grad(x, y, 'y');
        //Parte possivelmente desnecesária
        if (abs(xn - x) <= pow(10, -4) && abs(yn - y) <= pow(10, -4))
        {
            break;
        }
        if (f_multi(xn, yn) - f_multi(x, y) > 0)
        {
            cout << "Error!!" << endl;
            break;
        }
    } while (abs(xn - x) > pow(10, -4) || abs(yn - y) > pow(10, -4));

    cout << "x: " << xn << "\ny: " << yn << endl;
}

//6.3.5 Método de Levenberg - Marquardt
/*
Método que não volta atrás ao contrário do gradiente
A ideia notável que ocorreu separadamente a Kenneth Levenberg e a Donald Marquardt foi a de combinar
os dois métodos anteriores no mesmo passo, fazendo:
h = hquad(H^-1*nabla) +λ.nabla
1- Começar com λ grande e ir diminuindo de forma a obter um resultado
2- xn+1= xn - hLM
3- Passos da quadrica
4- Se a função decrescer diminuir o lambda, lambda/2
else aumentar lambda, lambda*2
*/
void Levenberg_Marquardt(double x, double y, double lambda) {
    double xn = x;
    double yn = y;
    do {
        x = xn;
        y = yn;
        double hLM = Hessiana_inverted_times_grad(xn, yn, 'x') + lambda * grad(xn, yn, 'x');
        xn = x - hLM;
        hLM = Hessiana_inverted_times_grad(xn, yn, 'y') + lambda * grad(xn, yn, 'y');
        yn = y - hLM;
        if (f_multi(xn, yn) - f_multi(x, y) < 0)
        {
            lambda /= 2;
        }
        else
        {
            lambda *= 2;
        }
    } while (abs(xn - x) > pow(10, -4) || abs(yn - y) > pow(10, -4));

    cout << "x: " << xn << "\ny: " << yn << endl;
}


int main()
{
    cout << "Hope you enjoy!" << endl;
}
