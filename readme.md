# Current issues
*比熱の計算が合わない

# Canonical Summation
$$\braket{\beta|\hat{A}|\beta} = {\rm e}^{-\beta N l} \sum_k \left(\frac{(N\beta)^{2k}}{(2k)!} \braket{k|\hat{A}|k} + \frac{(N\beta)^{2k + 1}}{(2k + 1)!} \braket{k|\hat{A}|k + 1}\right)$$

$$\braket{\beta|\hat{A}|\beta} = {\rm e}^{-\beta N l} \braket{0|0} \sum_k \left[\left(\prod_{j = 1}^k \frac{(N\beta)^2}{(2j - 1)2j} \frac{\braket{j|j}}{\braket{j-1|j-1}}\right) \left(\braket{\psi^{(k)}|\hat{A}|\psi^{(k)}} + \frac{N\beta}{2k + 1} \braket{\psi^{(k)}|\hat{A}(l - h)|\psi^{(k)}}\right)\right]$$
