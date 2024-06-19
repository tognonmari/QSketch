#include <iostream>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <random>
#include <vector>
#include "MurmurHash3.h"

using namespace std;

chrono::high_resolution_clock::time_point fg2_begin_update, fg2_end_update;
chrono::high_resolution_clock::time_point fg2_begin_estimate, fg2_end_estimate;

struct NextCUST{
    double x_i;
    int c;
};

inline int argmax(vector<double> GBM_sketch){
    int j_star=0;
    double max = GBM_sketch[0];
    for (int i = 0; i < GBM_sketch.size(); i++) {
        if (max < GBM_sketch[i]){
            max = GBM_sketch[i];
            j_star = i;
        }
    }   
    return j_star;
}

class FastGM
{

    public:
        FastGM(int n, int k);
        void AscEXP_STM(int i, int k, int z);
        void Update(double* data, int n);
        void EstimateCard();

        double hash(int data_id, int counter_id);

        vector<NextCUST> CUSTVec;
        vector<int> PermutationVec;
        vector<double> GBM_sketch;
        double estimated_card;
        int sketch_size;
        double update_time;
        double estimation_time;
        vector<uint32_t> seed;

};

FastGM::FastGM(int n, int k)
{
    NextCUST CUST;
    random_device rd;

    this->CUSTVec.resize(n);
    this->PermutationVec.resize(n*k);
    this->GBM_sketch.resize(k);
    this->seed.resize(k);
    this->sketch_size = k;
    this->estimated_card = 0.0;

    for(int i=0;i<n;i++) {
        CUST.x_i = 0;
        this->CUSTVec[i] = CUST;

        for(int j=0;j<k;j++) {
            this->PermutationVec[i*k+j] = j;
        }
    }

    for(int i=0;i<k;i++) {
        this->GBM_sketch[i] = -1;
        this->seed[i] = rd();
    }
    this->update_time = 0.0;
    this->estimation_time = 0.0;
}

inline double FastGM::hash(int data_id, int counter_id)
{
    double hash_value = 0.0;
    uint32_t hash_result;
    // long id = ((long)counter_id << 32) + data_id;
    // const long *key = & id;
    std::string key = std::to_string(data_id) + "|" + std::to_string(counter_id);
    // MurmurHash3_x86_32(&key, sizeof(key), this->seed[0], &hash_result);
    MurmurHash3_x86_32(&key, sizeof(key), this->seed[counter_id], &hash_result);

    hash_value = (double)hash_result / (double)UINT32_MAX;
    return hash_value;
}

void FastGM::AscEXP_STM(int i, int k, int z){

    double u=FastGM::hash(i, z-1);
    double x=0;
    int j=0;

    j = (rand() % (k-z+1))+z-1;
    x = -(log(u))/(k+1-z);

    swap(this->PermutationVec[i * k + z-1], this->PermutationVec[i * k + j]);

    this->CUSTVec[i].x_i += x;
    this->CUSTVec[i].c = this->PermutationVec[i * k + z-1];
}

void FastGM::Update(double* data, int n){

    fg2_begin_update = chrono::high_resolution_clock::now();

    int k = this->sketch_size;
    uint32_t delta = 0, k_star = k;   
    int j_star = 0, Prune_Module_Flag = false;
    for (int i=0;i<n;i++){
        double v_i = data[i];
        if(v_i == 0) continue;
        for (int l = 1; l < k+1; l++){
            AscEXP_STM(i, k, l);
            double b_i = this->CUSTVec[i].x_i /(v_i);
            int c = this->CUSTVec[i].c;
            if (Prune_Module_Flag == false){
                if (this->GBM_sketch[c] < 0 ){
                        this->GBM_sketch[c] = b_i;
                        k_star = k_star - 1;
                        if (k_star == 0){
                            Prune_Module_Flag = true;
                            j_star = argmax(this->GBM_sketch);
                        }
                    }
                    else if (b_i < this->GBM_sketch[c]){
                        this->GBM_sketch[c] = b_i;
                    }
            }
            else{
                if (this->GBM_sketch[j_star] < b_i) break;
                else if (b_i < this->GBM_sketch[c]){
                    this->GBM_sketch[c] = b_i;
                    if(c == j_star){
                        j_star = argmax(this->GBM_sketch);
                    }
                }
            }
        }
    }

    fg2_end_update = chrono::high_resolution_clock::now();
    chrono::duration<double> time1 = fg2_end_update - fg2_begin_update;
    this->update_time = time1.count();
}

void FastGM::EstimateCard(){

    fg2_begin_estimate = chrono::high_resolution_clock::now();

    double lambda = 0.0;
    for (int i=0; i<this->GBM_sketch.size(); i++){
        lambda += this->GBM_sketch[i];
    }
    this->estimated_card = (double)(this->GBM_sketch.size() - 1) / lambda;

    fg2_end_estimate = chrono::high_resolution_clock::now();
    std::chrono::duration<double> time2 = fg2_end_estimate - fg2_begin_estimate;
    this->estimation_time = time2.count();

}


