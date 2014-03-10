#Error-Correcting Code Goods for MATLAB
**(ecc_goods)**
MATLAB scripts and functions for error-correcting codes


###berlekamp_welc.m
This function uses the Berlekamp-Welch Algorithm to recover the original, uncorrupted n-length message of a received and potentially corrupted message 

*r<sub>1</sub> r<sub>2</sub> r<sub>3</sub> ...r<sub>n+2k</sub>*,
                
which is expected to have at most k corruption errors.

The recovered message returned by this function is in the form

*m<sub>1</sub> m<sub>2</sub> ... m<sub>n</sub> *
 
and has length n.

Note that in this context, operations are done modulo p, where p is prime, and "messages" are finite-length vectors/sequences of integers (mod p).

###row_reduce_p.m
Row reduces a matrix to its reduced row echelon form, if possible, modulo some number p. If complete, reduced row echelon is not possible for the input matrix, the function row reduces as much as possible.

###egcd.m
A simple MATLAB implementation of Euclid's extended algorithm.  Useful for finding multiplicative inverses modulo some number p. Or for actually finding the GCD of two integers.
