 Machine-learning-ab-initio-spin-dynamics 
 
The machine learning model is built to learn the magnetic potenital energy surface of noncolllinear spin systems. Based on this ML model, machine larning ab initio spin dynamics can be relaized efficiently. This project is built on the machine learning force field as implemented in the QUIP code.  
To run the program, the "descriptor for noncollinear spin.f95" should be incorporated into the source files of QUIP code. The descriptor including the freedom of spin is written in the subroutine named "sosd_ex". Descriptor for 3b term of exchange interaction is also provided （subroutine named as "sosd_3b_ex" in "descriptor-3b version 2.f95"). The correspoinding atom types and neighbour map list should be changed to read the spin system training sets, including the spin vector and atomic coordinates. 
An example of the input and output files are included in the "training and test" folder.  The results for bcc Fe are shown below.


![300-1000K (1)_page-0001](https://github.com/YuqiangGao/Machine-learning-ab-initio-spin-dynamics/assets/25586920/04ac578a-107d-4370-857f-8e28750038e3)

