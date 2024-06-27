#include <fstream>

#include "lm.h"
#include "qsketch.h"
#include "fastgm.h"
#include "qdyn.h"
#include "fastexp.h"

using namespace std;

int register_size = 8;

void qs_proc(double sketch_size, double* data, int data_num)
{

    QSketch qs(sketch_size, register_size);
    qs.Update(data, data_num);
    qs.EstimateCard();

    cout << "Qskech Results: " << qs.estimated_card << endl;
}

void lm_proc(double sketch_size, double *data, int data_num)
{
    LM lm(sketch_size);
    lm.Update(data, data_num);
    lm.EstimateCard();

    cout << "LM Results: " << lm.estimated_card << endl;
}

void fg_proc(double sketch_size, double* data, int data_num)
{
    FastGM fg(data_num, sketch_size);
    fg.Update(data, data_num);
    fg.EstimateCard();

    cout << "FastGM Results: " << fg.estimated_card << endl;
}

void fe_proc(double sketch_size, double* data, int data_num)
{
    FastExp fe(sketch_size);
    fe.Update(data, data_num);
    fe.EstimateCard();

    cout << "FastExp Results: " << fe.estimated_card << endl;
}

void qdyn_proc(double sketch_size, double *data, int data_num)
{
    QDyn qdyn(sketch_size, register_size, 0);
    qdyn.Update(data, data_num);

    cout << "Qdyn Results: " << qdyn.estimated_card << endl;
}

int main(int argc, char *argv[]){
    
    int sketch_size = atoi(argv[1]);
    int data_num = atoi(argv[2]);
    int rep_num = atoi(argv[3]);
    string file_name = argv[4];

    int index = 0;
    double* data = new double[data_num];
    ifstream file(file_name);
    if (file) {
        while(!file.eof()) {
            file>>data[index++];
            if (index == data_num) break;
        }
        file.close();
    } else {
        cout << "No such files !!!" << endl;
    }

    qs_proc(sketch_size, data, data_num);
    qdyn_proc(sketch_size, data, data_num);
    fg_proc(sketch_size, data, data_num);
    fe_proc(sketch_size, data, data_num);
    lm_proc(sketch_size, data, data_num);

    return 0;
}
