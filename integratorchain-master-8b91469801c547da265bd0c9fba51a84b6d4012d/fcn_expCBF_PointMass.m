function [B,dBdx] = fcn_expCBF_PointMass(x,params)

  B(1,1)=(x(1) - params(1))^2 + (x(2) - params(2))^2 + (x(3) - params(3))^2 - params(4)^2 + x(4)*(2*...
         x(1) - 2*params(1)) + x(5)*(2*x(2) - 2*params(2)) + x(6)*(2*x(3) - 2*params(3));

   dBdx(1,1)=2*x(4) + 2*x(1) - 2*params(1);
  dBdx(1,2)=2*x(5) + 2*x(2) - 2*params(2);
  dBdx(1,3)=2*x(6) + 2*x(3) - 2*params(3);
  dBdx(1,4)=2*x(1) - 2*params(1);
  dBdx(1,5)=2*x(2) - 2*params(2);
  dBdx(1,6)=2*x(3) - 2*params(3);

 