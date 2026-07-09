### OPS_OpenSEM

OPS version of methods implemented in https://github.com/chhc-1/OpenSEM .

oSEM, MRSEM, ISEM have been implemented. Note the version of ISEM in https://github.com/chhc-1/OpenSEM is currently bugged - it should not produce the correct Reynold's Stress Tensor.

These scripts have been tested using the OPS-DSL library [1] and HDF5-1.14.6.

The methods in this repository have been tested using the OpenSBLI CFD code [2], and the code for those simulations are in the 'OpenSBLI_cases' folder.

SEM animations:
![oSEM:](./figures/oSEM_gif.gif)

![MRSEM:](./figures/MRSEM_gif.gif)


\*\*\*Potential discrepancies in the manipulation of random numbers may exist: during the project, `ops_fill_random_uniform` produced both +ve and -ve numbers, whilst a test script on Iridis 6 produces only +ve numbers; modify this as necesary. 

### References
[1] I. Z. Reguly, G. R. Mudalige and M. B. Giles, Loop Tiling in Large-Scale Stencil Codes at Run-Time with OPS, in IEEE Transactions on Parallel and Distributed Systems, vol. 29, no. 4, pp. 873-886, 1 April 2018, doi: 10.1109/TPDS.2017.2778161.

[2] David J. Lusher, Andrea Sansica, Neil D. Sandham, Jianping Meng, Bálint Siklósi, Atsushi Hashimoto, OpenSBLI v3.0: High-fidelity multi-block transonic aerofoil CFD simulations using domain specific languages on GPUs, Computer Physics Communications, Volume 307, 2025, 109406, ISSN 0010-4655, https://doi.org/10.1016/j.cpc.2024.109406.


