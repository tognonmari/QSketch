import os
import subprocess
import pandas as pd
import threading

DATA_REPOSITORY_ROOT = "./data"

RESULTS_REPOSITORY_ROOT = "./results_diff_seed"

STREAM_LENGTHS = [10,100,1000,10000]#list(range(800000, 10 ** 6, 10 ** 5)) + list(range(10 ** 6, 10 ** 7, 10 ** 6)) #list(range(10,100, 10)) + list(range(100,10**3, 100)) + list(range(10**3,10 ** 4, 10 ** 3)) + list(range(10 ** 4, 10 ** 5, 10 ** 4)) + list(range(10 ** 5, 10 ** 6, 10 ** 5)) + list(range(10 ** 6, 10 ** 7, 10 ** 6))

WEIGHTS_DISTRIBUTIONS = ["uniform01", "uniform_integer_small_range", "uniform_integer_large_range"] #["unif_neg_range2", "unif_neg_range10","unif_neg_range100"]

EXECUTABLE = "./main.exe"

COL_NAMES =["ActualCount", "QS_est", "QS_time", "QDyn_est", "QDyn_time", "Whll_est", "Whll_time"] #["ActualCount", "QDynNeg_est", "QDynNeg_time"] 


def execute_tests(weight_distribution, res_repo_root, data_repo_root, str_lengths, exe_path, cols ):
    if not os.path.exists(f"{res_repo_root}/{weight_distribution}/"):
        os.makedirs(f"{res_repo_root}/{weight_distribution}/")

    for str_len in str_lengths:

        if not os.path.exists(f"{res_repo_root}/{weight_distribution}/{str_len}/"):
            os.makedirs(f"{res_repo_root}/{weight_distribution}/{str_len}/")
        
        instances_path = f"{data_repo_root}/{weight_distribution}/{str_len}/"
        for file in os.listdir(instances_path):

            instance_results = pd.DataFrame(columns = cols, index = range(4,15))
            
            for p in range(4,15):
                
                #Execute test on instance and collect results into a df then moved to .csv file named as the instance  
                str_exec = f"{exe_path} {int(2**p)} {instances_path + file} {str_len} 1 8"
                output = subprocess.run(str_exec, capture_output = True, text = True).stdout
                
                lines = output.split("\n")
                lines = lines [0:len(lines)-1]
                for i in range(0,len(lines)):
                    lines[i] = float(lines[i])
                
                instance_results.loc[p] = lines 
            
            #copy dataframe as csv to file
            instance_results.to_csv(f"{res_repo_root}/{weight_distribution}/{str_len}/"+ file)
        print(f"finished {str_len} for distribution {weight_distribution}")            
    



                    


if __name__ == "__main__":

    if not os.path.exists(RESULTS_REPOSITORY_ROOT):
        os.makedirs(RESULTS_REPOSITORY_ROOT)
    
    t1 = threading.Thread(target= execute_tests, args = (WEIGHTS_DISTRIBUTIONS[0], RESULTS_REPOSITORY_ROOT, DATA_REPOSITORY_ROOT, STREAM_LENGTHS, EXECUTABLE, COL_NAMES))
    t2 = threading.Thread(target= execute_tests, args = (WEIGHTS_DISTRIBUTIONS[1], RESULTS_REPOSITORY_ROOT, DATA_REPOSITORY_ROOT, STREAM_LENGTHS, EXECUTABLE, COL_NAMES))
    t3 = threading.Thread(target= execute_tests, args = (WEIGHTS_DISTRIBUTIONS[2], RESULTS_REPOSITORY_ROOT, DATA_REPOSITORY_ROOT, STREAM_LENGTHS, EXECUTABLE, COL_NAMES))

    t1.start()
    t2.start()
    t3.start()

    t1.join()
    t2.join()
    t3.join()


    """
    for wd in WEIGHTS_DISTRIBUTIONS:
        if not os.path.exists(f"{RESULTS_REPOSITORY_ROOT}/{wd}/"):
            os.makedirs(f"{RESULTS_REPOSITORY_ROOT}/{wd}/")

        for str_len in STREAM_LENGTHS:

            if not os.path.exists(f"{RESULTS_REPOSITORY_ROOT}/{wd}/{str_len}/"):
                os.makedirs(f"{RESULTS_REPOSITORY_ROOT}/{wd}/{str_len}/")
            
            instances_path = f"{DATA_REPOSITORY_ROOT}/{wd}/{str_len}/"
            for file in os.listdir(instances_path):

                instance_results = pd.DataFrame(columns = COL_NAMES)

                for p in range(4,15):
                    
                    #Execute test on instance and collect results into a df then moved to .csv file named as the instance  
                    str_exec = f"{EXECUTABLE} {int(2**p)} {instances_path + file} {str_len} 1 8"
    """





        



            
