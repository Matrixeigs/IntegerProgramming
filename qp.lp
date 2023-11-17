\ Model qp
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  2 xyz + [ 2 xyz ^2 + 2 xyz * C1 + 2 C1 ^2 + 2 C1 * C2 + 2 C2 ^2 ] / 2 
Subject To
 c1: xyz + 2 C1 + 3 C2 >= 4
 c2: xyz + C1 >= 1
 c3_1: xyz >= 0
 c3_2: C1 >= 0
 c3_3: C2 >= 0
Bounds
 xyz free
 C1 free
 C2 free
End
