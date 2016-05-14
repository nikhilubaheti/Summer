function [B,dBdx] = fcn_expCBF_IntChainDeg4(x,params)

  B(1,1)=8*(x(1) - params(1))^2 + 8*(x(2) - params(2))^2 + 8*(x(3) - params(3))^2 - 8*params(4)^2 +...
          4*x(4)*(2*x(1) - 2*params(1)) + 2*x(7)*(2*x(1) - 2*params(1)) + x(10)*(2*x(1) - 2*params(1)) + 4*x(5)*...
         (2*x(2) - 2*params(2)) + 2*x(8)*(2*x(2) - 2*params(2)) + x(11)*(2*x(2) - 2*params(2)) + 4*x(6)*(2*...
         x(3) - 2*params(3)) + 2*x(9)*(2*x(3) - 2*params(3)) + x(12)*(2*x(3) - 2*params(3));

   dBdx(1,1)=8*x(4) + 4*x(7) + 2*x(10) + 16*x(1) - 16*params(1);
  dBdx(1,2)=8*x(5) + 4*x(8) + 2*x(11) + 16*x(2) - 16*params(2);
  dBdx(1,3)=8*x(6) + 4*x(9) + 2*x(12) + 16*x(3) - 16*params(3);
  dBdx(1,4)=8*x(1) - 8*params(1);
  dBdx(1,5)=8*x(2) - 8*params(2);
  dBdx(1,6)=8*x(3) - 8*params(3);
  dBdx(1,7)=4*x(1) - 4*params(1);
  dBdx(1,8)=4*x(2) - 4*params(2);
  dBdx(1,9)=4*x(3) - 4*params(3);
  dBdx(1,10)=2*x(1) - 2*params(1);
  dBdx(1,11)=2*x(2) - 2*params(2);
  dBdx(1,12)=2*x(3) - 2*params(3);

 