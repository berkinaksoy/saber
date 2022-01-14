# SABER
SABER is a Mod-LWR based KEM submitted to NIST Post-Quantum Cryptography Process.
Poynomial multiplications are implemented using iterative two level of Karatsuba following up Toom-Cook algorithm with block recombination method. 

First, use "make clean" command to clean executable files.
Second, use "make all" command to compile to source codes.
Third, use "./test/test_kex" to run KEM operations in loop. 
