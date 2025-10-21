 Machine-learning-exchange-field for ab initio spin dynamics. 
 =============
 
The machine learning model is built to learn the magnetic potenital energy surface of noncolllinear spin systems. Based on this ML model, machine larning ab initio spin dynamics can be relaized efficiently. This project is built on the machine learning force field as implemented in the QUIP code.  
# How to install?
1. To run the program, the "descriptor for noncollinear spin.f95" should be incorporated into the source file "descriptor.f95" of QUIP code. The descriptor including the freedom of spin is written in the subroutine named "sosd_ex_calc". Descriptor for 3b term of exchange interaction is also provided （subroutine named as "sosd_3b_ex_calc" in "descriptor-3b version 2.f95").  Users should copy the files in the folder "src" and replace the cooresponding files in QUIP code. We have rewrote part of the QUIP code, including atom types and neighbour map list in the src folder to read the spin system training sets.
2. "make config" to choose what should be included in QUIP code.  Set "GAP" to yes.
3. make
4. use the "gap_fit" program to train the model and "quip" program to do prediction.

# Example
An example of the input and output files are included in the "training and test" folder.  The results for bcc Fe are shown below.


<img width="659" height="833" alt="image" src="https://github.com/user-attachments/assets/c3757412-c921-43e5-af69-0e9840db0e1b" />


