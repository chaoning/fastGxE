#include "string_utils.hpp"
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <algorithm>
#include <cmath>

using std::string;
using std::vector;
using std::istringstream;
using std::ostringstream;
using std::cout;
using std::endl;

// Process a line by removing comments, trimming spaces, and replacing certain delimiters with spaces
// line: input string to be processed
// comment_str: comment prefix (default is "#", can be "//" or "%")
void process_line(string &line, const string comment_str) {
    // Replace tab, comma, semicolon, carriage return, and newline with space
    for (char &c : line) {
        if (c == '\t' || c == ',' || c == ';' || c == '\r' || c == '\n') {
            c = ' ';
        }
    }

    // Trim leading and trailing spaces
    line.erase(0, line.find_first_not_of(" "));
    line.erase(line.find_last_not_of(" ") + 1);

    // Remove comments
    size_t comment_pos = line.find_first_of(comment_str);
    if (comment_pos != string::npos) {
        line.erase(comment_pos);
    }
}

// Split a string by spaces into a vector of strings
vector<string> split_string(const std::string &str) {
    vector<string> result;
    istringstream stream(str);
    string token;
    while (stream >> token) {
        result.push_back(token);
    }
    return result;
}

// Join a vector of strings using a given separator
string join_string(const std::vector<std::string> &str_vec, const std::string &split_str) {
    if (str_vec.empty()) return "";

    ostringstream result;
    for (size_t i = 0; i < str_vec.size(); ++i) {
        if (i > 0) result << split_str;
        result << str_vec[i];
    }
    return result.str();
}

std::string join_string(const Eigen::VectorXd &Vec, const std::string &split_str){
    std::ostringstream oss;
    for (size_t i = 0; i < Vec.size(); ++i) {
        if (i > 0) oss << split_str;  // Add space between elements
        oss << Vec[i];
    }
    return oss.str();
}

// Convert string to double with error handling
double string_to_double(const std::string &str, char *endptr) {
    double result = strtod(str.c_str(), &endptr);
    if (*endptr != '\0') {
        cout << "ERROR: '" << str << "' contains invalid characters: '" << endptr 
             << "', cannot be converted to double." << endl;
        exit(EXIT_FAILURE);
    }
    if (!std::isfinite(result)) {
        cout << "ERROR: The number is not finite: " << result << endl;
        exit(EXIT_FAILURE);
    }
    return result;
}


// Check if a string represents a NaN value
bool is_nan(const std::string &str, const std::vector<std::string> &nan_vec) {
    return std::find(nan_vec.begin(), nan_vec.end(), str) != nan_vec.end();
}



// Convert string to long long integer with error handling
long long string_to_longlong(const std::string &str, char *endptr) {
    long long result = strtoll(str.c_str(), &endptr, 10);
    if (*endptr != '\0') {
        cout << "ERROR: '" << str << "' contains invalid characters: '" << endptr 
             << "', cannot be converted to long long." << endl;
        exit(EXIT_FAILURE);
    }
    return result;
}


// Convert double to string with specified decimal places
string double_to_string(double value, int decimal) {
    ostringstream stream;
    stream << std::fixed << std::setprecision(decimal) << value;
    return stream.str();
}
