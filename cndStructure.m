%Note that the function cnd.m prepares the calculations per the documented formula for the other large program, whose output is shown above. The program stores in a structure all the sub-matrices of the main matrix and the associated norms, cond numbers etc.

clear all
d=100
M = rand(d);
cndStructureFunction(M)


function [dimM1,cnM1,normM1,normM1i,M1,M1i,M2,M2i,cnM2,normM2,normM2i,V1,V2,k0,b0,m0,VX,L1X,L1,F1,F2,p]=cndStructureFunction(M)

 disp '------------------'
 disp '--- Author: Raul E Garcia'
 disp '--- Email: siliconvalleypro1@yahoo.com' 
 disp '--- Program Name: cndStructureFuction'
 disp '------------------'

 p='fro';
 M1=M;
 M1i=inv(M1);
 [m,n] = size(M1)

 if abs(m-n)>0 
    disp 'Function can decompose condition numbers only for square matrices' 
     
 else
     dimM1=m
 end

 M2=M1(1:dimM1-1,1:dimM1-1); %-- An submatrix of An+1
 M2i=inv(M2);

 cnM1=cond(M1,p);
 cnM2=cond(M2,p);

 normM1=norm(M1,p);
 normM1i=norm(M1i,p);

 normM2=norm(M2,p);
 normM2i=norm(M2i,p);

 %% -- now the key - load in the nd x nd the previous matrix

 %-- resetted the 1st dim's to the above fixed random matrix

 %% -- create the edge vecors V1 (nx1) and V2 (1xn) on M2 last row/col
 V1a=M1(:,dimM1);
 V1=V1a(1:dimM1-1);%--skip las dia elem 1xn

 V2a=M1(dimM1,:);
 V2=V2a(1:dimM1-1);%--skip las dia elem nx1

 VX=V1*V2; %-- nxn matrix - nx1 * 1xn of edge vectors

 %% -- all values below are scalars
 k0=V2*M2i*V1; %-- a scalar - 1xn * nxn * nx1 then 1x1
 b0=M1(dimM1,dimM1); %-- last diag element of M1 - matrix with 1 more dim
 m0=1/(b0 - k0); %--main value, inverse matrix of scalar

 I=eye(dimM1-1);

 L1X=I + m0*VX*M2i;
 L1=M2i*L1X;

 %% -- now the matrices on the upper right, lower left corners
 F1= M2i*V1 ;
 F2= V2*M2i ;
 %% -- compute all the needed Norms
     
 end