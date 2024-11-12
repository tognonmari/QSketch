#include <ctime>
#include <cmath>
#include <algorithm>
#include <string>
#include <chrono>
#include "MurmurHash3.h"
#include "PackedVector.hpp"

std::chrono::high_resolution_clock::time_point qs_begin_update, qs_end_update;
std::chrono::high_resolution_clock::time_point qs_begin_estimate, qs_end_estimate;
/**
 * Returns the index of @param array that stores the minimum value.
 * 
 * @param array Array of numerical values
 * @param m Length of @param array
 * 
 * @return index of minimum element stored in @param array
 */
inline int argmin(PackedVector array, int m) 
{
    int idx = 0, minv = array.get(0);
    for (int i = 1; i < m; i++) {
        if (array.get(i) < minv) {
            minv = array.get(i);
            idx = i;
        }
    }
    return idx;
}
/**
 * Returns initial value to be provided to Newton-Raphson method, as described in the original paper.
 * 
 * @param sketch Packedvector representing QSketch
 * @param m Number of registers of the sketch
 * 
 * @return Initial value for Newton-Raphson method.
 */
double InitialValue(PackedVector sketch, int m)
{
    double c0 = 0.0;
    double tmp_sum = 0.0;
 
    for(int i=0; i<m; i++) { 
        tmp_sum += pow(2, -sketch.get(i)); 
    }

    c0 = (double)(m-1) / tmp_sum;
    return c0;
}
/**
 * Evaluates the log-likelihood function of the sketch's registers distribution at cardinality @param w.
 * 
 * @param sketch Packedvector representing QSketch
 * @param k Sketch size
 * @param w Point upon which the log-likelihood is evaluated
 * 
 * @return log-likelihood value at point @param w.
 * 
 */
inline double ffunc(PackedVector sketch, uint32_t k, double w) 
{
    double res = 0;
    double e = 2.718282;
    for (int i = 0; i < k; ++i) {
        double x = pow(2.0, -sketch.get(i) - 1);
        double ex = pow(e, w * x);
        res += x * (2.0 - ex) / (ex - 1.0);
    }
    return res;
}

/**
 * Evaluates the derivative for the log-likelihood function of the sketch's registers distribution at cardinality @param w.
 * 
 * @param sketch Packedvector representing QSketch
 * @param k Sketch size
 * @param w Point upon which the log-likelihood derivative is evaluated
 * 
 * @return log-likelihood derivative value at point @param w.
 * 
 */
inline double dffunc(PackedVector sketch, uint32_t k, double w) 
{
    double res = 0;
    double e = 2.718282;
    for (int i = 0; i < k; ++i) {
        double x = pow(2.0, -sketch.get(i) - 1);
        double ex = pow(e, w * x);
        res += -x * x * ex * pow(ex - 1, -2);
    }
    return res;
}

/**
 * Implements Newton-Raphson method for estimating the weighted cardinality. 
 * Computation is stopped when the newly-derived value of the estimate differs from the last one for less than 10e-5.
 * 
 * @param sketch PackedVector representing QSketch
 * @param k Sketch size
 * @param c0 Initial value 
 * 
 * @return estimation for the weighted cardinality based on sketch @param sketch
 */

double Newton(PackedVector sketch, uint32_t k, double c0) 
{
    double err = 1e-5;
    double c1 = c0 - ffunc(sketch, k, c0) / dffunc(sketch, k, c0);

    while (abs (c1 - c0) > err) {
        c0 = c1;
        c1 = c0 - (double)ffunc(sketch, k, c0) / dffunc(sketch, k, c0);
    }
    return c1;
}
/**
 * Utility method that produces a string with the content of a sketch.
 * 
 * @param qs PackedVector repressenting the sketch
 * 
 * @return string representing the sketch
 */
inline std::string drawRegister(PackedVector qs){

    std::string output = "\n";
    for (int i=0; i<qs.size();i++){

        output += ("R["+ std::to_string(i) +"] = "+ std::to_string(qs.get(i))+ "\n");

    }

    return output;

}
/**
 * Class implementing QSketch algorithm.
 */
class QSketch
{
    public:

        QSketch(int sketch_size, int register_size);
        QSketch(int sketch_size, int register_size, int multiplier_for_random_seed);
        void Update(double *data, int data_num);
        void EstimateCard();
        double hash(int data_id, int counter_id);
        void Update(double *data, int data_num, int *keys);
        void UpdateRandOnly(auto* data, int data_num, auto* keys);
        void UpdateRandOnlyNoPerm(double* data, int data_num, int* keys);
        
        int r_max;              //Maximum value that can be stored in the registers: any value above this threshold is stored as r_max.
        int r_min;              //Minimum value that can be stores in the registers: any value below this threshold is stored as r_min.
        int offset;             //Initial value stored in the registers.
        int sketch_size;        //Sketch size.
        int register_size;      //Sketch's registers' size in bits.
        int range;              //

        uint32_t *pi;           //Array that stores the random permutation produced by Fisher-Yates shuffle.
        uint32_t *seed;         //Vectors of random seeds that defined hash functions h_j.
        double estimated_card;  //Estimated weighted cardinality according to registers of PackedVector qs.
        double update_time;     //Total time required for insertion, assuming the update method is invoked once.
        double estimation_time; //Total time required for the estimation, i.e. to run Newton-Raphson method.
        int multiplier_for_random_seed;
        PackedVector qs;        //Actual representation of the sketch

        
};

QSketch::QSketch(int sketch_size, int register_size) : qs(register_size, sketch_size)
{

    this->register_size = register_size;
    this->sketch_size = sketch_size;
    this->range = pow(2, register_size);
    this->r_max = this->range - 1;
    this->r_min = 0;
    this->seed = new uint32_t[sketch_size];

    for (int i = 0; i < sketch_size; i++) {
        
        this->seed[i] = i; 
    }
    this->update_time = 0.0;
    this->estimation_time = 0.0;

}

QSketch::QSketch(int sketch_size, int register_size, int multiplier_for_random_seed) : qs(register_size, sketch_size)
{
    
    this->register_size = register_size;
    this->sketch_size = sketch_size;
    this->range = pow(2, register_size);
    this->r_max = this->range - 1;
    this->r_min = 0;
    this->seed = new uint32_t[sketch_size];
    this->multiplier_for_random_seed = multiplier_for_random_seed;
    for (int i = 0; i < sketch_size; i++) {
        this->seed[i] = i; 
    }
    this->update_time = 0.0;
    this->estimation_time = 0.0;

    
}
/**
 * Given a key ( @param data_id ) for an element of the data stream, it hashes it to a uniform value in (0,1) according to the @param counter_id-th hash function.
 * The method makes use of MurmurHash3 as an hash function.
 * @param datat_id key to be hashed
 * @param counter_id identifier of the hash function that we want to compute, coinciding with the register where the resulting value may be stored. 
 * 
 * @return decimal hash value for key ( @param data_id) according to the hash function specified by @param counter_id
 */
inline double QSketch::hash(int data_id, int counter_id)
{
    double hash_value = 0.0;
    uint32_t hash_result;
    std::string key = std::to_string(data_id)+"||"+ std::to_string(data_id); 
    
    MurmurHash3_x86_32(&key, sizeof(key), this->seed[counter_id], &hash_result);
    
    hash_value = (double)hash_result / (double)UINT32_MAX;
    
    return hash_value;
}
/*
void QSketch::Update(double *data, int data_num)
{
    qs_begin_update = std::chrono::high_resolution_clock::now();

    int j_min = 0, jj = 0;
    double r = 0.0, u = 0.0;
    int32_t y = 0;
    for (int t = 0; t < data_num; t++) {
        if (data[t] == 0.0) {
            continue;
        }
        this->pi = this->pii;
        r = 0.0;
        for (int i = 0; i < this->sketch_size; i++) {

            u = QSketch::hash(t, i);
            r = r - log(u) / (data[t] * (this->sketch_size - i));
            y = floor(-log2(r));

            if (y <= qs.get(j_min)) { break; }
            
            jj = rand() % (this->sketch_size - i) + i; //altro problema
            
            std::swap(this->pi[i], this->pi[jj]);

            if (y > int32_t(qs.get(this->pi[i]))) {

                if (this->r_min < y && y < this->r_max) {
                    qs.set(this->pi[i], y);
                } else if (y >= this->r_max) {
                    qs.set(this->pi[i], this->r_max);
                } else {
                    continue;
                }

                if (this->pi[i] == j_min) { j_min = argmin(qs, this->sketch_size); }
            }
        }
    }

    qs_end_update = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time1 = qs_end_update - qs_begin_update;
    this->update_time = time1.count();
}*/
/**
 * Computes an estimate for the weighted cardinality through Newton-Raphson method and saves it in the instance's estimated_card variable.
 * It also updates the estimation_time parameter.
 */
void QSketch::EstimateCard()
{
    qs_begin_estimate = std::chrono::high_resolution_clock::now();
    double c0 = InitialValue(qs, this->sketch_size);
    this->estimated_card = Newton(qs, this->sketch_size, c0);

    qs_end_estimate = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time2 = qs_end_estimate - qs_begin_estimate;
    this->estimation_time = time2.count();

}
/*
void QSketch::Update(double* data, int data_num, int* keys){
    qs_begin_update = std::chrono::high_resolution_clock::now();

    int j_min = 0, jj = 0;
    double r = 0.0, u = 0.0;
    int32_t y = 0;
    for (int t = 0; t < data_num; t++) {

        std::cout << "I am adding key "<< keys[t]<<" with weight "<< data[t]<<std::endl;
        if (data[t] == 0.0) {
            continue;
        }
        this->pi = this->pii;
        std::cout <<"this->pii[0] "<<this->pii[0]<<std::endl;
        r = 0.0;
        srand((uint64_t)(keys[t]));
        for (int i = 0; i < this->sketch_size; i++) {

            u = QSketch::hash(keys[t], i);
            r = r - log(u) / (data[t] * (this->sketch_size - i));
            y = floor(-log2(r));

            if (y <= qs.get(j_min)) { break; }
            
            jj = rand() % (this->sketch_size - i) + i;//fy shuffle
            std::cout << "For key "<< keys[t]<<"and iteration "<<i<<" the jj value is "<< jj<<std::endl;
            std::swap(this->pi[i], this->pi[jj]);
            if (y > int32_t(qs.get(this->pi[i]))) {

                if (this->r_min < y && y < this->r_max) {
                    qs.set(this->pi[i], y); 
                    std::cout <<"I am writing number "<<y<< " to register "<< this->pi[i]<<std::endl;
                } else if (y >= this->r_max) {
                    qs.set(this->pi[i], this->r_max);
                } else {
                    continue;
                }

                if (this->pi[i] == j_min) { j_min = argmin(qs, this->sketch_size); }
            }

        }

        std::cout << "After addition the sketch is "<< drawRegister(this->qs)<<std::endl;
    }

    qs_end_update = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time1 = qs_end_update - qs_begin_update;
    this->update_time = time1.count();
}*/
/**
 * Updates the sketch by inserting a total of @param data_num entries with weightes stored in @param data and keys stored in @param keys.
 * Saves in the instance's parameter update_time the time this specific update took.
 * 
 *@param data A pointer to an array of non-negative weights.
 *@param data_num The length of @param data and @param keys or eventually the last element to be added.
 *@param keys A pointer to an array of keys.
 */
void QSketch::UpdateRandOnly(auto* data, int data_num, auto* keys){

    qs_begin_update = std::chrono::high_resolution_clock::now();

    int j_min = 0, jj = 0;
    double r = 0.0, u = 0.0;
    int32_t y = 0;

    //Allocate the array that stores the permutation.
    this->pi = (uint32_t*)calloc(this->sketch_size,sizeof(uint32_t));


    for (int t = 0; t < data_num; t++) {

        //If an element has weight 0 it is not added.
        if (data[t] < 10e-9) {
            continue;
        }
        //Initialize the permutation array.
        for(int cc = 0; cc<this->sketch_size; cc++){

            this->pi[cc] = cc;
        }
        
        r = 0.0;

        //Initialize the random generator with a function of the key, so as to obtain equal permutations for the same key.
        srand((uint64_t)(keys[t])); //Not safe if key is not numerical: to be improved.

        for (int i = 0; i < this->sketch_size; i++) {

            u = QSketch::hash(keys[t], i);
            r = r - log(u) / (data[t] * (this->sketch_size - i));
            y = floor(-log2(r));

            if (y <= qs.get(j_min)) { break; }
            
            jj = rand() % (this->sketch_size - i) + i;
            
            std::swap(this->pi[i], this->pi[jj]);

            if (y > int32_t(qs.get(this->pi[i]))) {

                if (this->r_min < y && y < this->r_max) {
                    qs.set(this->pi[i], y); 
                    
                } else if (y >= this->r_max) {
                    qs.set(this->pi[i], this->r_max);
                } else {
                    continue;
                }

                if (this->pi[i] == j_min) { j_min = argmin(qs, this->sketch_size); }
            }

        }

    }

    qs_end_update = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time1 = qs_end_update - qs_begin_update;
    this->update_time = time1.count();
    free(this->pi);

}