# Current issues
- 比熱の計算が合わない
  - 横磁場イジング, 4サイトの場合絶対零度で比熱が発散suru
  - ハイゼンベルグは関数形は良いが2-3倍異なる
  - エネルギー密度はかなり良い一致、磁化率は多少のずれはあるものの悪くない
  - $\hat{h}^2$の演算子にエラーがあるのか？
- 遅い
  - HPhiのmTPQに比べるとだいぶ遅い($\chi^2$回の平均を取っていると見なしても)
  - ベンチマークをとってみたところボトルネックはCanonical FormにするときのSVD

# Canonical Summation
$$\braket{\beta|\hat{A}|\beta} = {\rm e}^{-\beta N l} \sum_k \left(\frac{(N\beta)^{2k}}{(2k)!} \braket{k|\hat{A}|k} + \frac{(N\beta)^{2k + 1}}{(2k + 1)!} \braket{k|\hat{A}|k + 1}\right)$$

$$\braket{\beta|\hat{A}|\beta} = {\rm e}^{-\beta N l} \braket{0|0} \sum_k \left[\left(\prod_{j = 1}^k \frac{(N\beta)^2}{(2j - 1)2j} \frac{\braket{j|j}}{\braket{j-1|j-1}}\right) \left(\braket{\psi^{(k)}|\hat{A}|\psi^{(k)}} + \frac{N\beta}{2k + 1} \braket{\psi^{(k)}|\hat{A}(l - h)|\psi^{(k)}}\right)\right]$$
