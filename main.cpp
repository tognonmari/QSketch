#include <fstream>

#include "lm.h"
#include "qsketch.h"
#include "fastgm.h"
#include "qdyn.h"
#include "fastexp.h"
#include "whll.h"
#include "qdynbeta.h"
#include "qdynneg.h"
using namespace std;

// int register_size = 8;

void qs_proc(double sketch_size, double* data, int data_num, int register_size)
{

    QSketch qs(sketch_size, register_size);
    qs.Update(data, data_num);
    qs.EstimateCard();

    cout << qs.estimated_card << endl;
    cout <<qs.update_time<<endl;
}

void lm_proc(double sketch_size, double *data, int data_num)
{
    LM lm(sketch_size);
    lm.Update(data, data_num);
    lm.EstimateCard();

    cout << "LM Results: " << lm.estimated_card << endl;
    cout << "LM TIme : "<< lm.update_time<<endl;
}

void fg_proc(double sketch_size, double* data, int data_num)
{
    FastGM fg(data_num, sketch_size);
    fg.Update(data, data_num);
    fg.EstimateCard();

    cout << "FastGM Results: " << fg.estimated_card << endl;
    cout << "FastGM TIme : "<< fg.update_time<<endl;
}

void fe_proc(double sketch_size, double* data, int data_num)
{
    FastExp fe(sketch_size);
    fe.Update(data, data_num);
    fe.EstimateCard();

    cout << "FastExp Results: " << fe.estimated_card << endl;
    cout << "FastExp TIme : "<< fe.update_time<<endl;
}

void qdyn_proc(double sketch_size, double *data, int data_num, int register_size)
{
    QDyn qdyn(sketch_size, register_size, 0);
    qdyn.Update(data, data_num);

    cout  << qdyn.estimated_card << endl;
    cout <<qdyn.update_time<<endl;
}

void whll_proc(double sketch_size, double *data, int data_num){

    WHLLSketch whll = WHLLSketch(sketch_size, 8); //hardcoded space for a single byte
    
    whll.Update(data,data_num);
    
    whll.EstimateCard();
    //cout << "Whll Results: " << whll.estimated_card <<endl;
    //cout << "Whll TIme : "<< whll.update_time<<endl;
    cout << whll.estimated_card<<endl;
    cout << whll.update_time<<endl;
}

void qdynbeta_proc(double sketch_size, double *data, int data_num, int register_size)
{
    QDynBeta qdynbeta(sketch_size, register_size, 0);
    qdynbeta.Update(data, data_num);

    cout << "QdynBeta Results: " << qdynbeta.estimated_card << endl;
    cout << "QdynBeta TIme : "<<qdynbeta.update_time<<endl;
}

void qdynneg_proc(double sketch_size, double *data, int data_num, int register_size){
    QDynNeg qdynbeta(sketch_size, register_size, 0);
    qdynbeta.Update(data, data_num);

    cout << qdynbeta.estimated_card << endl;
    cout <<qdynbeta.update_time<<endl;
}

int main(int argc, char *argv[]){
    
    int sketch_size = atoi(argv[1]);
    int data_num = atoi(argv[3]);
    int rep_num = atoi(argv[4]);
    string file_name = argv[2];
    int register_size = atoi(argv[5]);
    int index = 0;
    double actual_count;

    double* data = new double[data_num];
    ifstream file(file_name);
    if (file) {
        while(!file.eof()) {
            file>>data[index++];
            actual_count += data[index-1];
            if (index == data_num) break;
        }
        file.close();
    } else {
        cout << "No such files !!!" << endl;
    }

    cout <<actual_count<< endl;
    qs_proc(sketch_size, data, data_num, register_size);
    qdyn_proc(sketch_size, data, data_num, register_size);
    //qdynneg_proc(sketch_size, data, data_num, register_size);
    //fg_proc(sketch_size, data, data_num);
    //fe_proc(sketch_size, data, data_num);
    //lm_proc(sketch_size, data, data_num);
    whll_proc(sketch_size, data, data_num);
    //qdynbeta_proc(sketch_size, data, data_num, register_size);
    return 0;
}
