function [DTI,DTJ] = gradHOij(T)
% Andrew Wheeler 2016
      D16 = zeros(7,10);
      D16(1,1)=-4192853913297.768d0/1.0d12;
      D16(1,2)=5339040076397.708d0/1.0d12;
      D16(1,3)=-995327078300.349d0/1.0d12;
      D16(1,4)=-150859084799.5902d0/1.0d12;
      D16(2,1)=-403794985785.2408d0/1.0d12;
      D16(2,3)=431827299959.2449d0/1.0d12;
      D16(2,4)=-44851702678.4067d0/1.0d12;
      D16(2,5)=16819388504.4025d0/1.0d12;
      D16(3,1)=157062810031.9159d0/1.0d12;
      D16(3,2)=-900989848376.5991d0/1.0d12;
      D16(3,4)=959478115560.413d0/1.0d12;
      D16(3,5)=-250644037526.018d0/1.0d12;
      D16(3,6)=35092960310.2883d0/1.0d12;
      D16(4,1)=16968084745.7541d0/1.0d12;
      D16(4,2)=66702567274.9842d0/1.0d12;
      D16(4,3)=-683894144860.1648d0/1.0d12;
      D16(4,5)=753863390649.9535d0/1.0d12;
      D16(4,6)=-178653360538.646d0/1.0d12;
      D16(4,7)=25013462728.1191d0/1.0d12;
      D16(5,2)=-26519952061.4978d0/1.0d12;
      D16(5,3)=189413141579.3246d0/1.0d12;
      D16(5,4)=-799266426974.1559d0/1.0d12;
      D16(5,6)=799266426974.1559d0/1.0d12;
      D16(5,7)=-189413141579.3246d0/1.0d12;
      D16(5,8)=26519952061.4978d0/1.0d12;
%      
[ni,nj] = size(T);
DTI(1:ni,1:nj) = 0.0;
DTJ(1:ni,1:nj) = 0.0;


% body
I = 5:ni-4;
J = 1:nj;
DTI(I,J) = D16(5,2)*T(I-3,J) ...
    + D16(5,3)*T(I-2,J) ...
    + D16(5,4)*T(I-1,J) ...
    - D16(5,2)*T(I+3,J) ...
    - D16(5,3)*T(I+2,J) ...
    - D16(5,4)*T(I+1,J) ;
%
for II=1:4
DTI(II,J) = D16(II,1)*T(1,J) ...
    + D16(II,2)*T(2,J) ...
    + D16(II,3)*T(3,J) ...
    + D16(II,4)*T(4,J) ...
    + D16(II,5)*T(5,J) ...
    + D16(II,6)*T(6,J) ...
    + D16(II,7)*T(7,J);
%
end
%
for II=ni-3:ni
DTI(II,J) = -D16(ni-II+1,1)*T(ni,J) ...
    - D16(ni-II+1,2)*T(ni-1,J) ...
    - D16(ni-II+1,3)*T(ni-2,J) ...
    - D16(ni-II+1,4)*T(ni-3,J) ...
    - D16(ni-II+1,5)*T(ni-4,J) ...
    - D16(ni-II+1,6)*T(ni-5,J) ...
    - D16(ni-II+1,7)*T(ni-6,J);
end
%
%
I = 1:ni;
J = 5:nj-4;
DTJ(I,J) = D16(5,2)*T(I,J-3) ...
    + D16(5,3)*T(I,J-2) ...
    + D16(5,4)*T(I,J-1) ...
    - D16(5,2)*T(I,J+3) ...
    - D16(5,3)*T(I,J+2) ...
    - D16(5,4)*T(I,J+1) ;
%
for JJ=1:4
DTJ(I,JJ) = D16(JJ,1)*T(I,1) ...
    + D16(JJ,2)*T(I,2) ...
    + D16(JJ,3)*T(I,3) ...
    + D16(JJ,4)*T(I,4) ...
    + D16(JJ,5)*T(I,5) ...
    + D16(JJ,6)*T(I,6) ...
    + D16(JJ,7)*T(I,7);
end
%
for JJ=nj-3:nj
DTJ(I,JJ) = -D16(nj-JJ+1,1)*T(I,nj) ...
    - D16(nj-JJ+1,2)*T(I,nj-1) ...
    - D16(nj-JJ+1,3)*T(I,nj-2) ...
    - D16(nj-JJ+1,4)*T(I,nj-3) ...
    - D16(nj-JJ+1,5)*T(I,nj-4) ...
    - D16(nj-JJ+1,6)*T(I,nj-5) ...
    - D16(nj-JJ+1,7)*T(I,nj-6);
end

return