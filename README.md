# Current issues
- 比熱の計算が合わない
  - ~~横磁場イジング, 4サイトの場合絶対零度で比熱が発散する~~ $l - \hat{h}$の演算子に打ち間違いあり. 治った. ただし磁化率は発散
  - ハイゼンベルグは関数形は良いが2-3倍異なる
  - エネルギー密度はどちらもかなり良い一致
- 遅い
  - ベンチマークをとってみたところボトルネックはCanonical FormにするときのSVD

# Note
## Canonical Summation
$$\braket{\beta|\hat{A}|\beta} = {\rm e}^{-\beta N l} \braket{0|0} \sum_k \left[\left(\prod_{j = 1}^k \frac{(N\beta)^2}{(2j - 1)2j} \frac{\braket{j|j}}{\braket{j-1|j-1}}\right) \left(\braket{\psi^{(k)}|\hat{A}|\psi^{(k)}} + \frac{N\beta}{2k + 1} \braket{\psi^{(k)}|\hat{A}(l - h)|\psi^{(k)}}\right)\right]$$