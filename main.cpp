#include <iostream>
#include <algorithm>
#include <iomanip>
#include <string>
#include <iostream>
#include <fstream>

#include "parser.h"

std::string ReplaceAll(std::string str, const std::string& from, const std::string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    return str;
}

// equation = -((1,804 * z^2) / (1000)) - (490 / 1000) 

char * substitute_values_to_eq(std::string equation, 
                                    double x,
                                    double y,
                                    double z) {
    std::string res;
    res = ReplaceAll(equation, "x", std::to_string(x));
    res = ReplaceAll(res,      "y", std::to_string(y));
    res = ReplaceAll(res,      "z", std::to_string(z));

    const int len = res.length();
    char* char_array = new char[len + 1];
    strcpy(char_array, res.c_str()); 

    return char_array;
}

void print_table(double **table, int ii, std::ofstream& file) {
    const int xi_pos = 0, yi_pos = 1, zi_pos = 2;

    // writing to console output
    std::cout << std::setw(15) << "i" 
              << std::setw(15) << "xi" 
              << std::setw(15) << "yi" 
              << std::setw(15) << "zi" << std::endl;

    // writing to file
    file << 'x' << ' ' 
         << 'y' << ' ' 
         << 'z' << std::endl;

    for (size_t i = 0; i < ii; ++i) {
        // writing to console output
        std::cout << std::setw(15) << i 
                  << std::setw(15) << table[i][xi_pos]
                  << std::setw(15) << table[i][yi_pos]
                  << std::setw(15) << table[i][zi_pos] << std::endl;

        // writing to file
        file << table[i][xi_pos] << ' ' 
             << table[i][yi_pos] << ' ' 
             << table[i][zi_pos] << std::endl;
    }
}

void runge_kutta(std::ofstream& file,
                 std::string equation1,
                 std::string equation2,
                 int    n,
                 double h, 
                 double x0, 
                 double y0, 
                 double z0) {

    // bool x_in_eq = equation.find('x') == std::string::npos;
    // bool y_in_eq = equation.find('y') == std::string::npos;
    // bool z_in_eq = equation.find('z') == std::string::npos;

    int i = 1;
    
    // initialize 2d array with xi, yi, zi values
    double** values = new double*[n+1];
    for (size_t j = 0; j < n+1; j++){
        values[j] = new double[3](); // zeros ?
    }
    const int xi_pos = 0, yi_pos = 1, zi_pos = 2;

    // set initial values
    values[0][xi_pos] = x0;
    values[0][yi_pos] = y0;
    values[0][zi_pos] = z0;

    double zi, yi, k1, k2, k3, k4, q1, q2, q3, q4;
    while (i < n + 1) {
        // y_i = y_(i-1) + h/6 * (k_1 + 2*k_2 + 2*k_3 + k_4)
        // k_1 = f(x_(i-1),       y_(i-1))
        // k_2 = f(x_(i-1) + h/2, y_(i-1) + h/2 * k_1)
        // k_3 = f(x_(i-1) + h/2, y_(i-1) + h/2 * k_2)
        // k_4 = f(x_(i-1) + h,   y_(i-1) + h   * k_3)
        parser obj;
        q1 = obj.eval_exp(substitute_values_to_eq(equation1, values[i-1][xi_pos], 
                                                             values[i-1][yi_pos], 
                                                             values[i-1][zi_pos]));

        k1 = obj.eval_exp(substitute_values_to_eq(equation2, values[i-1][xi_pos], 
                                                             values[i-1][yi_pos], 
                                                             values[i-1][zi_pos]));

        q2 = obj.eval_exp(substitute_values_to_eq(equation1, values[i-1][xi_pos] + h/2, 
                                                             values[i-1][yi_pos] + k1*h/2, 
                                                             values[i-1][zi_pos] + q1*h/2));
                            
        k2 = obj.eval_exp(substitute_values_to_eq(equation2, values[i-1][xi_pos] + h/2, 
                                                             values[i-1][yi_pos] + k1*h/2, 
                                                             values[i-1][zi_pos] + q1*h/2));

        q3 = obj.eval_exp(substitute_values_to_eq(equation1, values[i-1][xi_pos] + h/2, 
                                                             values[i-1][yi_pos] + k2*h/2, 
                                                             values[i-1][zi_pos] + q2*h/2));
                            
        k3 = obj.eval_exp(substitute_values_to_eq(equation2, values[i-1][xi_pos] + h/2, 
                                                             values[i-1][yi_pos] + k2*h/2, 
                                                             values[i-1][zi_pos] + q2*h/2));   

        q4 = obj.eval_exp(substitute_values_to_eq(equation1, values[i-1][xi_pos] + h, 
                                                             values[i-1][yi_pos] + k3*h, 
                                                             values[i-1][zi_pos] + q3*h));
                            
        k4 = obj.eval_exp(substitute_values_to_eq(equation2, values[i-1][xi_pos] + h, 
                                                             values[i-1][yi_pos] + k3*h, 
                                                             values[i-1][zi_pos] + q3*h));      

        zi = values[i-1][zi_pos] + (h/6) * (q1 + 2*q2 + 2*q3 + q4);
        yi = values[i-1][yi_pos] + (h/6) * (k1 + 2*k2 + 2*k3 + k4);

        values[i][xi_pos] = values[i-1][xi_pos] + h; 
        values[i][yi_pos] = yi;
        values[i][zi_pos] = zi;

        i++;
    }

    print_table(values, n + 1, file);
}

int main(){
    std::ofstream myfile;
    myfile.open("example.txt");

    // std::cout << "Hello World!" << std::endl;
    std::string eq1 = "-0.001804 * z^2 - 0.49";
    std::string eq2 = "z";
    // ans: y = e^(2 * x) * (-0.00341797*x^4 + 0.0625*x^3 - 0.25*x^2 + 0.75*x - 9.75)
    runge_kutta(myfile, eq1, eq2, 1000, 0.1, 0, 0, 50);

    myfile.close();
    return 0;
}
