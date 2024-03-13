#include <iostream>
#include <algorithm>
#include <iomanip>
#include <string>
#include <iostream>
#include <fstream>

#include "parser.h"


// replace all ocurrences of a string with a string
std::string ReplaceAll(std::string str, const std::string& from, const std::string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    return str;
}


// replace x, y and z with values
char* substitute_values_to_eq(std::string& equation, 
                              double x,
                              double y,
                              double z) {
    std::string res;
    res = ReplaceAll(equation, "x", std::to_string(x));
    res = ReplaceAll(res,      "y", std::to_string(y));
    res = ReplaceAll(res,      "z", std::to_string(z));

    // turning std::string to a c-style string for parser
    const int len = res.length();
    char* char_array = new char[len + 1];
    strcpy(char_array, res.c_str()); 

    return char_array;
}


// write the result values to console and file
void print_table(double **table, int ii, std::ofstream& file) {
    const int xi_pos = 0, yi_pos = 1, zi_pos = 2; // declare positions

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


// runge kutta (4) method implementation
void runge_kutta(std::ofstream& file,
                 std::string&   equation1,
                 std::string&   equation2,
                 int    n,
                 double h, 
                 double x0, 
                 double y0, 
                 double z0) {
    int i = 1; // initial values are given (i = 0)
    
    // initialize 2d array with xi, yi, zi values
    double** values = new double*[n+1];
    for (size_t j = 0; j < n+1; j++){
        values[j] = new double[3](); // fill with zeros
    }
    const int xi_pos = 0, yi_pos = 1, zi_pos = 2; // declare positions

    // set initial values
    values[0][xi_pos] = x0;
    values[0][yi_pos] = y0;
    values[0][zi_pos] = z0;

    double zi, yi, k1, k2, k3, k4, q1, q2, q3, q4;
    while (i < n + 1) {
        // z' = U(x, y, z)
        // y' = V(x, y, z)

        // z_i = z_(i-1) + h/6 * (q_1 + 2*q_2 + 2*q_3 + q_4)

        // y_i = y_(i-1) + h/6 * (k_1 + 2*k_2 + 2*k_3 + k_4)

        // q_1 = U(x_(i-1),       y_(i-1),             z_(i-1))
        // q_2 = U(x_(i-1) + h/2, y_(i-1) + h/2 * k_1, z_(i-1) + h/2 * q_1)
        // q_3 = U(x_(i-1) + h/2, y_(i-1) + h/2 * k_2, z_(i-1) + h/2 * q_2)
        // q_4 = U(x_(i-1) + h,   y_(i-1) + h   * k_3, z_(i-1) + h   * q_3)

        // k_1 = V(x_(i-1),       y_(i-1),             z_(i-1))
        // k_2 = V(x_(i-1) + h/2, y_(i-1) + h/2 * k_1, z_(i-1) + h/2 * q_1)
        // k_3 = V(x_(i-1) + h/2, y_(i-1) + h/2 * k_2, z_(i-1) + h/2 * q_2)
        // k_4 = V(x_(i-1) + h,   y_(i-1) + h   * k_3, z_(i-1) + h   * q_3)

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

        // add values to the table
        values[i][xi_pos] = values[i-1][xi_pos] + h; 
        values[i][yi_pos] = yi;
        values[i][zi_pos] = zi;

        i++;
    }

    print_table(values, n + 1, file);
}

int main(){
    std::ofstream myfile;
    myfile.open("output.txt");

    // 1000 * x'' + 1.804 * (x')^2 + 490 = 0
    // z = x'
    //
    // => z' = -0.001804 * z^2 - 0.49 
    //    x' = z
    //
    // x(0) = 0
    // z(0) = 50

    std::string eq1 = "-0.001804 * z^2 - 0.49";
    std::string eq2 = "z";

    runge_kutta(myfile, eq1, eq2, 1000, 0.1, 0, 0, 50);

    myfile.close();
    return 0;
}
