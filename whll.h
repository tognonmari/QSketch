#include <ctime>
#include <cmath>
#include <algorithm>
#include <random>
#include <string>
#include <chrono>
#include "MurmurHash3.h"
#include "PackedVector.hpp"
#include "windows.h"
std::chrono::high_resolution_clock::time_point whll_begin_update, whll_end_update;
std::chrono::high_resolution_clock::time_point whll_begin_estimate, whll_end_estimate;


class WHLLSketch
{
    public:

        
        WHLLSketch(int sketch_size, int register_size);
        
        void Update(double* data, int data_num);
        void EstimateCard();
        double getSumInversePow2();
        uint32_t hash(int data_id);

        uint32_t seed;
        int sketch_size;
        int register_size;
        int p;
        double estimated_card;
        double update_time;
        double estimation_time;
        double power2to32minusP;
        double alphaMM;
        PackedVector whll;
};
//naive implementation without constants
double WHLLSketch::getSumInversePow2(){

    double sum =0; 
    
    for (int t=0; t<this->sketch_size; t++){
        
        sum += pow(2, -(whll.get(t)));

    }
    
    return sum;

}

WHLLSketch::WHLLSketch(int sketch_size, int register_size): whll(register_size, sketch_size)
{
    std::random_device rd;
    this->sketch_size = sketch_size;
    this->seed = 1230;
    this->p= (int) log2(sketch_size);
    this->power2to32minusP = pow(2.0, 32-this->p);
    this->update_time = 0.0;
    this->estimation_time = 0.0;
    this-> alphaMM = alphaMM = 0.697; 
    
    
}

inline uint32_t WHLLSketch::hash(int data_id)
{
    uint32_t hash_result;
    std::string key = std::to_string(data_id);
    MurmurHash3_x86_32(&key, sizeof(key), this->seed, &hash_result);
    
    return hash_result;

}

void WHLLSketch::Update(double *data, int data_num)
{
    whll_begin_update = std::chrono::high_resolution_clock::now();

    for (int t=0; t<data_num; t++){

        int element_id = t;

        double weight = data[t];
        
        uint32_t hashcode = hash(t);
        
        int register_id = (int) (hashcode & (sketch_size - 1));
        //int register_id = (int) (hashcode%this->sketch_size);
        uint32_t w = hashcode >>this-> p;
        
        double unif01 = ((double)(w+1))/this->power2to32minusP;
        //double unif01 = ((double)(hashcode+1))/(double) UINT32_MAX;

        double h2Tilde = 1 - pow(1.0 - unif01, 1.0/weight);
        
        
        uint8_t to_be_inserted = ((uint8_t)floor(-log2(h2Tilde)))+1;
        
        

       if( to_be_inserted > this->whll.get(register_id)) {
            //std::cout << "here we have in register_id "<< this->whll.get(register_id)<<std::endl;
            this->whll.set(register_id, to_be_inserted);

       }
    }

    whll_end_update = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time1 = whll_end_update - whll_begin_update;
    this->update_time = time1.count();
}
 
void WHLLSketch::EstimateCard()
{
    whll_begin_estimate = std::chrono::high_resolution_clock::now();
    //Count the cardinality
    double sum = getSumInversePow2();
    this->estimated_card = (this->alphaMM * this->sketch_size * this->sketch_size)/sum;

    whll_end_estimate = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time2 = whll_end_estimate - whll_begin_estimate;
    this->estimation_time = time2.count();

}

