# QSketch

This is the source code of QSketch in KDD 2024.

# Compilation
To compile the code, navigate to the directory containing the source files and run the following command in your terminal:

```shell
g++ -o main main.cpp MurmurHash3.cpp -std=c++17 -O3
```

# Running the Code
To run the executable, use the following command format:

```shell
./main <register_num> <file_name> <data_size> <repeated_times> <register_size>
```

register_num: Number of registers.

file_name: Input data file

data_size: Number of data elements.

repeated_times: Number of times the experiment is repeated for averaging results.

register_size: Size of each register in bits.
