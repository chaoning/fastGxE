
#include <ctime>
#include <iostream>
#include <string>
#include <set>
#include "random_num.hpp"
#include "string_utils.hpp"
using std::cout;
using std::endl;
using std::vector;
using std::string;


vector<long long> GenerateDiffNumber(long long min, long long max, long long num)
{
    if(num > max - min){
        cout << num << " exceeds the upper limit" << endl;
        exit(1);
    }
    long long rnd;
    vector<long long> diff;
    vector<long long> tmp;//存储剩余的数
    //初始化
    for(long long i = min; i < max; i++ )
    {
        tmp.push_back(i);
    }
    srand((unsigned)time(0)); //初始化随机数种子
    random_shuffle(tmp.begin(), tmp.end());
    for(long long i = 0 ; i < num ; i++)
    {
        diff.push_back(tmp[i]);
    }
    return diff;
}


void GenerateDiffPairNumber(long long min, long long max, long long num, vector<long long>& res_vec0, vector<long long>& res_vec1)
{
    if(num > (max - min) * (max - min - 1) / 2){
        cout << num << " exceeds the upper limit" << endl;
        exit(1);
    }

    vector<long long> tmp_vec;
    for(long long i = min; i < max; i++ )
    {
        tmp_vec.push_back(i);
    }


    srand((unsigned)time(0)); 
    std::set<std::string> random_set;
    for(long long i = 0; i < 100; i++){
        random_shuffle(tmp_vec.begin(), tmp_vec.end());
        vector<long long> tmp_vec1 = tmp_vec;
        random_shuffle(tmp_vec.begin(), tmp_vec.end());
        vector<long long> tmp_vec2 = tmp_vec;

        for(long long k = 0; k < tmp_vec1.size(); k++){
            long long a = tmp_vec1[k];
            long long b = tmp_vec2[k];
            if(a < b){
                random_set.insert(std::to_string(a) + " " + std::to_string(b));
            }else if (a > b){
                random_set.insert(std::to_string(b) + " " + std::to_string(a));
            }
        }
        if(random_set.size() >= num) break;
    }

    if(random_set.size() < num){
        cout << "Cannot produce enough random number pairs!" << endl;
        exit(1);
    }

    vector<string> random_vec(random_set.begin(), random_set.end());
    random_shuffle(random_vec.begin(), random_vec.end());
    res_vec0.clear();
    res_vec1.clear();
    for(long long i = 0; i < num; i++){
        vector<string> tmp_vec = split_string(random_vec[i]);
        res_vec0.push_back(atoll(tmp_vec[0].c_str()));
        res_vec1.push_back(atoll(tmp_vec[1].c_str()));
    }
}


void GenerateDiffPairNumber2(long long min, long long max, long long num, vector<long long>& res_vec0, vector<long long>& res_vec1)
{
    if(num > (max - min) * (max - min - 1)){
        cout << num << " exceeds the upper limit" << endl;
        exit(1);
    }

    vector<long long> tmp_vec;
    for(long long i = min; i < max; i++ )
    {
        tmp_vec.push_back(i);
    }


    srand((unsigned)time(0)); 
    std::set<std::string> random_set;
    for(long long i = 0; i < 100; i++){
        random_shuffle(tmp_vec.begin(), tmp_vec.end());
        vector<long long> tmp_vec1 = tmp_vec;
        random_shuffle(tmp_vec.begin(), tmp_vec.end());
        vector<long long> tmp_vec2 = tmp_vec;

        for(long long k = 0; k < tmp_vec1.size(); k++){
            long long a = tmp_vec1[k];
            long long b = tmp_vec2[k];
            if(a != b){
                random_set.insert(std::to_string(a) + " " + std::to_string(b));
            }
        }
        if(random_set.size() >= num) break;
    }

    if(random_set.size() < num){
        cout << "Cannot produce enough random number pairs!" << endl;
        exit(1);
    }

    vector<string> random_vec(random_set.begin(), random_set.end());
    random_shuffle(random_vec.begin(), random_vec.end());
    res_vec0.clear();
    res_vec1.clear();
    for(long long i = 0; i < num; i++){
        vector<string> tmp_vec = split_string(random_vec[i]);
        res_vec0.push_back(atoll(tmp_vec[0].c_str()));
        res_vec1.push_back(atoll(tmp_vec[1].c_str()));
    }
}

