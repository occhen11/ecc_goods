#Error-Correcting Code Goods (ecc_goods) for MATLAB

##berlekamp_welc.m
This function uses the Berlekamp-Welch Algorithm to recover the original, uncorrupted message of a received and potentially corrupted message 
                        r1 r2 r3 ...,
which is expected to have a known max number of corruption errors.

The recovered message returned by this function is in the form
                        m_1  m_2 ... m_n
for a message of length n.

Note that in this context, operations are done modulo p, where p is prime, and "messages" are finite-length vectors/sequences of integers (mod p).

##row_reduce_p.m
Row reduces a matrix to its reduced row echelon form, if possible, modulo some number p. If complete, reduced row echelon is not possible for the input matrix, the function row reduces as much as possible.

##egcd.m
A simple MATLAB implementation of Euclid's extended algorithm.  Useful for finding multiplicative inverses modulo some number p. Or for actually finding the GCD of two integers.
