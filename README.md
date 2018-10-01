# Table of Content
- [Table of Content](#table-of-content)
- [1. My Changes](#1-my-changes)
- [2. How to Use](#2-how-to-use)
- [3. Contact](#3-contact)
- [4. Information about original repository](#4-information-about-original-repository)

# 1. My Changes
* this projects dependents on eigen3 for fitting Ax=b.
  
To avoid calculate $A^TA$ to increase condition number, here householder QR decomposition (from eigen3) was used.

$$\min_{x}\|Ax-b\|^2 $$

$$
\begin{equation}
\begin{split}
\|Ax-b\|^2 &= \|QRx-b\|^2\\
           &= \|Q(Rx-Q^Tb)\|^2 \text{, (Q is unitary)}\\
           &:= \| \begin{pmatrix}R_{diag}x\\0\end{pmatrix} - \begin{pmatrix}Q_1^Tb\\Q_2^Tb\end{pmatrix} \|^2 \text{, (thin QR)} \\
           &= \|R_{diag}x-Q_1^Tb\|^2 + \|Q_2^Tb\|^2\\
\end{split}\end{equation}
$$

$$\min_{x}\|Ax-b\|^2 \Rightarrow \min \|R_{diag}x-Q_1^Tb\|^2$$

finally we have $ x = R_{diag}^{-1}Q_1^Tb$, where $Q_1,R_{diag}$ are from thin QR decomposition.


<!-- 
$$
\begin{equation}
\begin{split}
(Ax-b)^T(Ax-b) &=x^T A^TAx - x^TA^Tb -b^T Ax + b^Tb\\
\frac{d(Ax-b)^T(Ax-b)}{dx} &=2A^TAx-2A^Tb =0\\
\end{split}\end{equation}
$$

Here QR decomposition was used to avoid calculate $A^TA$, which squares condition number.
$$
\begin{equation}
\begin{split}
    x  &= (A^TA)^{-1}A^Tb\\
       &= ((QR)^TQR)^{-1}A^Tb\\
       &= (R^TR)^{-1}A^Tb 
\end{split}\end{equation}
$$ -->

* add quadratic model (done). 
  (_green dots are inliner, red dots are curve definition points, empty dots are outliers, legends are the same for all figures_)

![RANSAC quadratic fitting example](./examples/quadratic_fitting.png)

* add polynomial model (done)

![RANSAC polynomial fitting example](./examples/polynomial_fitting.png)
 
* add linear least square model for over-deterministic problem (done). 
  _performance is the same with quadratic and polynomial implementation, since all used lookup table for point-to-curve distance calculation._

![lls fitting example](./examples/lls_fitting.png)

* regularized (L2) least square (working), [discussion](https://math.stackexchange.com/questions/2013160/qr-factorization-regularized-least-squares)

$$\min_{x}\|Ax-b\|_2^2 + \lambda \|x\|_2^2$$

$$
\begin{equation}
\begin{split}
(Ax-b)^T(Ax-b) + \lambda x^tx   &=(x^TA^T-b^T)(Ax-b) + \lambda x^tx \\
(\text{take derivative}) &=x^T A^TAx - x^TA^Tb -b^T Ax + b^Tb + 2\lambda x\\
                         &=2A^TAx-2A^Tb + 2\lambda x = 0\\
\end{split}\end{equation}
$$

$$
\begin{equation}
\begin{split}
    x  &= (A^TA+\lambda I)^{-1}A^Tb\\
       &= ((QR)^TQR+\lambda I)^{-1}A^Tb\\
       &= (R^TR+\lambda I)^{-1}A^Tb 
\end{split}\end{equation}
$$

* regularized (L1) least square (maybe), has problem with subgradient calculation, related [discussion](https://stsievert.com/blog/2015/12/09/inverse-part-2/)

$$\min_{x}\|Ax-b\|_2^2 + \lambda \|x\|_1$$

$$
\begin{equation}
\begin{split}
(Ax-b)^T(Ax-b) + \lambda \vert x \vert   &=(x^TA^T-b^T)(Ax-b) + \lambda \vert x \vert \\
(\text{take derivative}) &=x^T A^TAx - x^TA^Tb -b^T Ax + b^Tb + \lambda \text{sign}(x)\\
                         &=2A^TAx-2A^Tb + \lambda \text{sign}(x) = 0\\
\end{split}\end{equation}
$$


* polish codes (working)
* GPU support (maybe)
* add SVD for least squares (maybe)


# 2. How to Use
* openmp is required for multi-threads, see cmake file
* compile
```bash
$ mkdir build && cd build
$ cmake ..
$ make
```

* API, see examples

# 3. Contact
Autonomous Driving / Machine Learning Engineer

Author, Zhiliang Zhou, [Linkedin](https://www.linkedin.com/in/zhiliang-zhou/)


---
# 4. Information about original repository
GRANSAC: Multi-threaded generic RANSAC implemetation

Author, Srinath Sridhar (srinaths@umich.edu), Max Planck Institute for Informatics
