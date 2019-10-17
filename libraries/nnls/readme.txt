The following two files are fast versions 
of the M-file NNLS supplied by Mathworks.

The files given here are modifications of 
the NNLS m-file to ease comparison. Our own
implementations were used in the paper 
referenced below.

The file FNNLS corresponds exactly to NNLS
except it requires crossproducts of X and y
as input instead of X and y. FNNLS corresponds
to the algorithm FNNLSa in Bro & de Jong, 
J. Chemom., 1997, 11, p. 393

The file FNNLSB corresponds to NNLS but 
requires additional input, namely the index
set defining which constraints are initially
active and which are passive. This is helpful
for speeding up in iterative algorithms where
similar problems are solved numerous times.
FNNLSB corresponds to the algorithm FNNLSb
in Bro & de Jong, J. Chemom., 1997, 11, p. 393

Feb. 1997
Rasmus bro
Chemometrics Group, Food Technology
Dept. Dairy and Food Science
Royal Vet. & Agricultural
DK-1958 Frederiksberg C
Denmark
rb@kvl.dk
http://newton.foodsci.kvl.dk/rasmus.html
