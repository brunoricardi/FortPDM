How to use PDM

PDM is a program that calculates the pair density matrix for a pair of atoms interacting through the Aziz potential or the Lennard-Jones potential.
The calculation is done by applying numerical convolutions to a partial wave expansion of the pair density matrix at high temperatures where a semiclassical approximation is used.

Input parameters (optimum experienced values)

nsq: number of squarings (8)
tau_max: highest temperature (10240 K to achieve 40 K)
nw: number of partial waves (40)
ng: number of grid points in radial coordinates (2000)
na: number of grid points in angular coordinate (40)
rn: renormalization factor to write results in file (10)
dr: radial grid spacing (0.004)
da: angular grid spacing (0.005)


Output is composed of two files.

d.txt: contains the elements of the pair density matrix that have cos(theta) = 1, ie theta = 0.
pdm.txt: contains the full pair density matrix

Output format in both files is as follows:

radial coordinate 1, radial coordinate 2, angular coordinate (cos), density matrix.



In order to analise the results, three more programs are needed.

ler-diag.f:
reads d.txt and returns the pair action of the diagonal elements of the pair density matrix, output is ep.txt.
number of lines must be corrected for different parameters
output format is
radial coordinate, end-point action


ler-qsz.f:
uses output of ler-diag.f and pdm.txt to compute the density matrix in the auxiliary variables q, s and z, output is qsz.txt.
number of lines must be corrected for different parameters
output format is
q variable, s variable, z variable, pair action (diagonal term discounted)


ler-qc.f:
separates in different files points with fixed q variable, this should be the last execution.
number of lines must be corrected for different parameters
output format is
q variable, s variable, z variable, pair action (diagonal term discounted)


To compute the u_{kj} functions, one adjusts the surfaces obtained in this last step to a polynomial expression in the variables s and z.


