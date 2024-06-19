
#include <iostream>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <random>
#include <vector>
#include <string>
#include <chrono>
#include <ctime>
#include "MurmurHash3.h"
using namespace std;

std::uniform_real_distribution<double> dis2(0.0, 1.0);
std::uniform_int_distribution<int> uni2(0, RAND_MAX);
std::chrono::high_resolution_clock::time_point fg_begin_update, fg_end_update;
std::chrono::high_resolution_clock::time_point fg_begin_estimate, fg_end_estimate;

inline int argmax(double* array, int m)
{
    int idx = 0, maxv = array[0];
    for (int i = 1; i < m; i++) {
        if (array[i] > maxv) {
            maxv = array[i];
            idx = i;
        }
    }
    return idx;
}

inline double max_value(double* array, int m)
{
    int idx = 0; 
    double maxv = array[0];
    for (int i = 1; i < m; i++) {
        if (array[i] > maxv) {
            maxv = array[i];
            idx = i;
        }
    }
    return maxv;
}

class FastExp
{
    public:

        FastExp(int sketch_size);
        void Update(double *data, int data_num);
        void EstimateCard();

        double hash(int data_id, int counter_id);

        int sketch_size;
        double *sketch;
        uint32_t *pi;
        uint32_t *pii;
        double estimated_card;
        double update_time;
        double estimation_time;
        uint32_t *seed;
        uint32_t global_seed;

};

FastExp::FastExp(int sketch_size)
{
    std::random_device rd;
    this->sketch_size = sketch_size;
    this->sketch = new double[sketch_size];
    this->pii = new uint32_t[sketch_size];
    this->seed = new uint32_t[sketch_size];

    for (int i = 0; i < sketch_size; i++) {
        this->sketch[i] = 100000;
        this->pii[i] = i;
        this->seed[i] = rd();
    }
    this->update_time = 0.0;
    this->estimation_time = 0.0;
    this->global_seed = rd();
}

inline double FastExp::hash(int data_id, int counter_id)
{
    double hash_value = 0.0;
    uint32_t hash_result;
    std::string key = std::to_string(data_id) + "|" + std::to_string(counter_id);
    MurmurHash3_x86_32(&key, sizeof(key), this->seed[counter_id], &hash_result);

    hash_value = (double)hash_result / (double)UINT32_MAX;
    return hash_value;
}

void FastExp::Update(double *data, int n)
{
    fg_begin_update = std::chrono::high_resolution_clock::now();

    int j_max = 0, jj = 0;
    double r = 0.0, u = 0.0;
    double MAX = 100000;
    int z=0;

    for (int t = 0; t < n; t++) {

        this->pi = this->pii;
        bool updateMAX = false;
        r = 0.0;
        for (int i = 0; i < this->sketch_size; i++) {

            u = FastExp::hash(t, i);
            r = r - log(u) / (data[t] * (this->sketch_size - i));
            if (r > MAX) {
                break; 
            }

            jj = rand() % (this->sketch_size - i) + i;
            std::swap(this->pi[i], this->pi[jj]);
            int j=this->pi[i];

            if ( MAX == this->sketch[j] ) {updateMAX = true;}
            
            this->sketch[j] = std::min(r, this->sketch[j]);
        }

        if (updateMAX) { MAX = max_value(this->sketch, this->sketch_size); 
        }
    }

    fg_end_update = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time1 = fg_end_update - fg_begin_update;
    this->update_time = time1.count();
    
}

void FastExp::EstimateCard(){
    
    fg_begin_estimate = chrono::high_resolution_clock::now();
    double lambda = 0.0;

    for (int i = 0; i < this->sketch_size; i++){
        lambda += this->sketch[i];
    }
    this->estimated_card = (double)(this->sketch_size - 1) / lambda;

    fg_end_estimate = chrono::high_resolution_clock::now();
    std::chrono::duration<double> time2 = fg_end_estimate - fg_begin_estimate;
    this->estimation_time = time2.count();
}
