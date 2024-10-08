#include <ctime>
#include <cmath>
#include <algorithm>
#include <random>
#include <string>
#include <chrono>
#include "MurmurHash3.h"
#include "PackedVector.hpp"

std::uniform_int_distribution<int> TempUni(0, RAND_MAX);
std::chrono::high_resolution_clock::time_point qdynneg_begin_update, qdynneg_end_update;
std::random_device temp_rd;

class QDynNeg{

    public:

        QDynNeg(int m, int b, int seed);
        void Update(double *data, int data_num);

        uint32_t hash(int data_id, int counter_id);
        uint32_t hash2(int data_id, int counter_id);

        int r_max;
        int r_min;
        int offset;
        int sketch_size;
        int register_size;
        int* T;
        double estimated_card;
        int range;
        double update_time;
        double estimation_time;
        uint32_t *seed;
        int new_seed;

        PackedVector qdynneg;
};

QDynNeg::QDynNeg(int m, int b, int seed) : qdynneg(b, m)
{
    this->range = pow(2, b);
    this->r_max = this->range - 1;
    this->r_min = 0;
    this->offset = pow(2, b-1) - 1;
    this->estimated_card = 0.0;
    this->register_size = b;
    this->sketch_size = m;
    this->seed = new uint32_t[m];
    for (int i = 0; i < m; i++){
        this->seed[i] = temp_rd();
    }
    this->new_seed = seed;

    int table_size = this->range;
    this->T = new int[table_size];
    for (int i = 0; i < table_size; i++){
        this->T[i] = 0;
    }
    this->T[0] = this->sketch_size;

    this->update_time = 0.0;
    this->estimation_time = 0.0;

}

inline uint32_t QDynNeg::hash(int data_id, int counter_id)
{
    uint32_t hash_result;
    std::string key = std::to_string(data_id) + "|" + std::to_string(counter_id);
    MurmurHash3_x86_32(&key, sizeof(key), this->new_seed, &hash_result);

    return hash_result;
}

inline uint32_t QDynNeg::hash2(int data_id, int counter_id)
{
    uint32_t hash_result;
    std::string key = std::to_string(data_id) + "|" + std::to_string(counter_id);
    MurmurHash3_x86_32(&key, sizeof(key), this->seed[0], &hash_result);

    return hash_result;
}

void QDynNeg::Update(double *data, int data_num)
{
    qdynneg_begin_update = std::chrono::high_resolution_clock::now();

    uint32_t counter_id = 0;
    double u = 0.0, r = 0.0;
    int32_t y = 0, origin_val = 0;
    double qR = 0.0, tmp = 1.0;

    for (int i = 0; i < data_num; i++) {

        if (data[i] == 0.0) {
            continue;
        }

        uint32_t hv = QDynNeg::hash(i, counter_id);

        counter_id = hv % this->sketch_size;
        u = (double)hv / (double)UINT32_MAX;

        r = -log(u) / abs(data[i]); //ABS data[i]
        //r = 1- pow(1-u,(double)1.0/data[i]);
        y = floor(-log2(r)) + this->offset;
        y = std::min(std::max(y, this->r_min), this->r_max);

        origin_val = qdynneg.get(counter_id);
        if (y > origin_val){

            // update qR
            tmp = 0.0;
            for (int j = 0; j < this->range; j++){
                if(T[j]) tmp += this->T[j] * exp(-abs(data[i]) * (pow(2, this->offset-j-1))); 
            }
            qR = 1 - tmp / this->sketch_size;

            // update estimation value
            this->estimated_card += (data[i]) / qR;

            // update T
            if (this->T[origin_val] > 0){ // otherwise 
                this->T[origin_val] -= 1;
                this->T[y] += 1;
            }
            else{ // initialize
                this->T[y] += 1;
            }
            
            if (this->r_min < y && y < this->r_max) {
                qdynneg.set(counter_id, y);
            } else if (y >= this->r_max) {
                qdynneg.set(counter_id, this->r_max);
            } else {
                continue;
            }
        }
    }

    qdynneg_end_update = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time1 = qdynneg_end_update - qdynneg_begin_update;
    this->update_time = time1.count();
    this->estimation_time = this->update_time;
}