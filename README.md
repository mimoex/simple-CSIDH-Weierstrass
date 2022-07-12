# simple-CSIDH-Weierstrass

###　説明
CSIDH原文では$y^2=x^3+Ax$を利用しているため，右辺の$A$のみを求めます．

この実装では同種写像暗号を直感的に理解するために，$A,B$両方のパラメータを求めます．
この場合，$Fp$の計算を行うだけです．

#### Curve parameter
$y^2=x^3+Ax+B$

$E_{[0]}: A=1, B=0$
