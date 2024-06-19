#include <iostream>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <random>
#include <chrono>
#include <ctime>
#include "MurmurHash3.h"

// using namespace std;

std::chrono::high_resolution_clock::time_point lm_begin_update, lm_end_update; //更新
std::chrono::high_resolution_clock::time_point lm_begin_estimate, lm_end_estimate; //估计

class LM
{
    
    public:

        LM(int sketch_size);
        void Update(double *data, int data_num);
        void EstimateCard(); 
        
        double hash(int data_id, int counter_id);

        int sketch_size;
        double *sketch;
        double estimated_card;
        double update_time;
        double estimation_time;
        uint32_t *seed;

};

LM::LM(int sketch_size)
{
    std::random_device rd;
    this->estimated_card = 0.0;
    this->sketch_size = sketch_size;
    this->sketch = new double[sketch_size];
    this->seed = new uint32_t[sketch_size];
    for(int i = 0; i < sketch_size; i++){
        this->sketch[i] = 100000.0;
        this->seed[i] = rd();
    }
    this->update_time = 0.0;
    this->estimation_time = 0.0;
}

inline double LM::hash(int data_id, int counter_id)
{
    double hash_value = 0.0;
    uint32_t hash_result;
    std::string key = std::to_string(data_id) + "|" + std::to_string(counter_id);
    MurmurHash3_x86_32(&key, sizeof(key), this->seed[counter_id], &hash_result);

    hash_value = (double)hash_result / (double)UINT32_MAX;
    return hash_value;
}


void LM::Update(double *data, int data_num){
    
    lm_begin_update = std::chrono::high_resolution_clock::now();

    double u = 0.0, exp = 0.0;
    for (int i = 0; i < data_num; i++){
        for (int j = 0; j < this->sketch_size; j++) {
            u = LM::hash(i, j);
            exp = -log(u) / data[i];
            this->sketch[j] = std::min(this->sketch[j], exp);
        }
    } 

    lm_end_update = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time1 = lm_end_update - lm_begin_update;
    this->update_time = time1.count();
}

void LM::EstimateCard(){

    lm_begin_estimate = std::chrono::high_resolution_clock::now();

    double lambda = 0.0;
    for (int i = 0; i < this->sketch_size; i++) {
        lambda += this->sketch[i];
    }
    this->estimated_card = double((this->sketch_size - 1) / lambda);

    lm_end_estimate = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time2 = lm_end_estimate - lm_begin_estimate;
    this->estimation_time = time2.count();

}
