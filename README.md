<h1 align="center">
    <br>parallel_example
</h1>

This is a small parallel programming experience for testing what the author has learned in the parallel program designing.

There are mainly three topics which includes the calculation of pi、the test of SPECOMP2012 and the calculation of matrix.

# 0x1.For the calculation of pi
we try to get the specific value as accurate as possible. Besides, it's also our goal to speed up the
program with specific machine includeing AVX-float、AVX-double、SSE-float、SSE-double、serial and openmp versions.
Set the compling argvs as follows:
~~~~bash
gcc -mavx2 -msse2 -fopenmp -O3 pi.c -o result
~~~~
Get the normal experience according to the results:
~~~~text
AVX-float>AVX-double≈SSE-float>SSE-double>serial
~~~~
![Image text](https://github.com/JK-Star/images/blob/master/pi.png)

Details can be reached at [topic-pi.docx](https://github.com/JK-Star/parallel_examples/blob/master/doc/topic-pi.docx)

# 0x2.For the test of SPECOMP2012
As to the second, just test the SPECOMP2012 with different threads or setting the environment at dynamic ,static or guided.
Details can be reached at [topic-SPECOMP2012.docx](https://github.com/JK-Star/parallel_examples/blob/master/doc/topic-SPECOMP2012.docx)

# 0x3.For the calculation of matrix
The last, may be the most familiar to us.The versions includes tranpose、blocked、openmp+tranpose and MPI.
Details can be reached at [topic-matrix.docx](https://github.com/JK-Star/parallel_examples/blob/master/doc/topic-matrix.docx)
~~~~bash
gcc -mavx2 -msse2 -O3 -fopenmp matrix.c -o test
~~~~

Assuming the matrix $A * B = C$, all are square matrices, $A=1_{n*n}, B=2_{n*n}$, $n=1024$, the result about time consumption is as follows:

|    Time: ms     |     1     |     2     |     3     |     4     |  average  |
| :-------------: | :-------: | :-------: | :-------: | :-------: | :-------: |
|     origin      | 6602757.0 | 6284082.0 | 6060442.0 | 6298664.0 | 6311486.3 |
|    tranpose     | 914533.0  | 873473.0  | 931000.0  | 1010096.0 | 932275.5  |
|     blocked     | 6349623.0 | 5452277.0 | 6361220.0 | 6153948.0 | 6079267.0 |
| openmp+tranpose | 724606.5  | 796924.2  | 995741.2  | 809043.0  | 831578.7  |

![Image text](https://github.com/JK-Star/images/blob/master/mpi.png)



# Acknowledgement
If anything wrong you find, welcome to correct me. Or you can reach me at Minglikesworking@163.com.




