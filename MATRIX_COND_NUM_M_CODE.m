function [dimM1,cnM1,normM1,normM1i,M1,M1i,M2,M2i,cnM2,normM2,normM2i,V1,V2,k0,b0,m0,VX,L1X,L1,F1,F2,p]=cnd(M)

 disp '------------------'
 disp '--- Author: Raul E Garcia'
 disp '--- Email: siliconvalleypro1@yahoo.com' 
 disp '--- Program Name: MATXCOND'
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
 
  

 clear all
 disp '------------------'
 disp '--- Author: Raul E Garcia'
 disp '--- Email: siliconvalleypro1@yahoo.com' 
 disp '--- Program Name: MATXCOND'
 disp '--- Version 7 - Date Updated: 1/13/13'
 disp '--- Last Upd: added plot of diagonals matrix'
 disp '--- This program decomposes a matrix condtion number into its causes'
 disp '--- We assume the matrix to decompose is in Matlab desktop with the name: condMatrix or a matrix file'
 disp '------------------'

% --- Uppercase variables are matrices or vectors,i=inverse, lowercase scalars

 % If the matrix is loaded in mem and the name is condMatrix then comment
 % the load command below, else we assume the matrix is saved in the curr
 % 
 % Matrices below generated and saved in C:\MATXCOND\ by prog: 
 % largeCondNumVer7.m

 fLoaded='C:\MATXCOND\condMatrix10001.mat';
 load (fLoaded);
 str1= ['--- Loaded Matrix from file: ',fLoaded];
 disp(str1)

 cm=condMatrix; %-- Matrix An+1
 m=cm;
 mi=inv(cm);

 %-- dcm = condition matrix dim
 dcm=length(cm);
 str2= ['--- The matrix dimension is: ',num2str(dcm)];
 disp(str2)
 %---
 p='fro';
 str3= ['--- Type of norm used: ',p];
 disp(str3)
 %---
 format long
 CondNum=round(cond(m,p));
 MatrixNorm1=round(norm(m,p));
 MatrixNorm2=round(norm(mi,p));
 str4= ['--- The Condition Number is: ',num2str(CondNum),', The norm: ',num2str(MatrixNorm1),', Inverse matrix Norm: ',num2str(MatrixNorm2)];
 disp(str4)
 %---

 format short
% fro=Frobenius norm

 dg1n=diag(m);
 figure
 plot(dg1n)
 title 'Diagonal Elements of Full Matrix (dg1n)'
 grid on

 dg1i=diag(mi);
 figure
 plot(dg1i)
 title 'Diagonal Elements of INVERSE of the Full Matrix (dg1i)'
 grid on

D=dcm:-1:2;

 %-- (decreasing) index domain for count of randome matrices for a fixed dimension input

sm=cm; %-- initial value given matrix, sm=submatrix 
 diffCond=zeros(1,dcm);

 disp '--- The matrix is: '
 m

for dm = D
               
         [dimM1,cnM1,normM1,normM1i,M1,M1i,M2,M2i,cnM2,normM2,normM2i,V1,V2,k0,b0,m0,VX,L1X,L1,F1,F2,p]=cnd(sm);

         %% NMD - Natural Matrix Decomposition into 2x2 matrices

         nmdstrucM1(dm)=struct('dimM1',dm,'M1',M1,'M2',M2,'V1',V1,'V2',V2,'b0',b0);

         nmdstrucM1i(dm)=struct('dimM1i',dm,'M1i',M1i,'M2i',M2i,'F1',F1,'F2',F2,'L1X',L1X,'L1',L1,'m0',m0);

         %% CND - Condition Number Descomposition
         cndstruc(dm) = struct ('dimM1',dimM1,'cnM1',cnM1,'normM1',normM1,'normM1i',normM1i,'M1',M1,'M1i',M1i,'M2',M2,'M2i',M2i,'cnM2',cnM2,'normM2',normM2, 'normM2i',normM2i,'V1',V1,'V2',V2,'k0',k0,'b0',b0,'m0',m0,'VX',VX,'L1X',L1X,'L1',L1,'F1',F1,'F2',F2,'p',p);
         
         normM1(dm)=norm(M1,p); %--given matrix An+1
         normM2(dm)=norm(M2,p); %-- summatrix An

         normM1Log(dm)=log10(normM1(dm));
         normM2Log(dm)=log10(normM2(dm));
         
         normV1(dm)=norm(V1,p);
         normV2(dm)=norm(V2,p);

         normV1Log(dm)=log10(norm(V1,p));
         normV2Log(dm)=log10(norm(V2,p));
         
         normB0(dm)=abs(b0); %--norm is abs value for scalar-last diag element
         
         normL1(dm)=norm(L1,p);
         normF1(dm)=norm(F1,p);
         normF2(dm)=norm(F2,p);
         
         normL1Log(dm)=log10(normL1(dm));
         normF1Log(dm)=log10(normF1(dm));
         normF2Log(dm)=log10(normF2(dm));

         %% next use the formula to compute the An+1 cond numner based on An
  
         n(1)=normM2(dm)^2;
         n(2)=normV1(dm)^2;
         n(3)=normV2(dm)^2;
         n(4)=b0^2;
         
         n(5)=sum(n(1:4));
         M1norm(1)=sqrt(n(5));
         
         ni(1)=normL1(dm)^2;
         ni(2)=(m0*normF1(dm))^2;
         ni(3)=(m0*normF2(dm))^2;
         ni(4)=m0^2;
         
         ni(5)=sum(ni(1:4));
         M1norm(2)=sqrt(ni(5));

         %% note that we have 4 terms since is the norm of a 2x2 (block) matrix
                
         condNumCalc(dm)=  M1norm(1) * M1norm(2)   ;

         %dm-1   -- submatrix dimension 

         V= norm(V1,p)^2 + norm(V2,p)^2;
         A1=norm(M2,p)^2; %-- summatrix An
         A2=norm(M2i,p)^2;

         mv=abs(m0*k0);

         factor1=cond(M2,p)^2 + A2*(V + b0^2);
         factor2=(dm-1) + 2*mv*sqrt((dm-1)*A2) + (mv^2)*A2 + (m0^2)*(V+1);

         upperBoundCN(dm)=factor1*factor2;
         
         condNum1(dm,1)=cond(M1,p); %% direct value - cond num
         condNum1(dm,2)=log10(condNum1(dm,1));

         condNum2(dm,1)=cond(M2,p); %% direct value - cond num
         condNum2(dm,2)=log10(condNum2(dm,1));
         
         cn(dm,1)=dm;
         cn(dm,2)=condNum1(dm);
         cn(dm,3)=condNum2(dm);
         cn(dm,4)=condNum1(dm)-condNum2(dm);
         cn(dm,5)=log10(condNum1(dm)/condNum2(dm));
         
         normsM1(dm,1)= M1norm(1);
         normsM1(dm,2)= M1norm(2);
         normsM1Log(dm,1)=log10(M1norm(1));
         normsM1Log(dm,2)=log10(M1norm(2));
         normsM1Log(dm,3)= normsM1Log(dm,1)* normsM1Log(dm,2);
    
         normB0(dm)=abs(b0); %--norm is abs value for scalar-last diag element
         
         normL1(dm)=norm(L1,p);
         normF1(dm)=norm(F1,p);
         normF2(dm)=norm(F2,p);
                 
         scalarsB0(dm,1)=b0;
         scalarsB0(dm,2)=log10(abs(b0));
         scalarsM0(dm,1)=m0;
         scalarsM0(dm,2)=log10(abs(m0));
         scalarsV0(dm,1)=k0;
         scalarsV0(dm,2)=log10(abs(k0));
         
         partSum(dm,1)= n(5);
         partSum(dm,2)= ni(5);
         
         for x=1:5
             sum1normsLog(dm,x)=log10(n(x));
             sum2normsLog(dm,x)=log10(ni(x));
         end    
                
       
         diffCond(dm)= abs(condNumCalc(dm) - condNum1(dm)) ;
                
         if diffCond(dm) > 1-eps
             disp '** Non zero difference**'
             diffCond(dm)
             break
         else
             M1;
             M2;
             disp 'ok...Matrix Dim: ',dm
         end
       
         % -- initilize to submatrix and repeat process
         sm=M2;
 end
  
 condNumberPlot1(sum1normsLog,false, CondNum, MatrixNorm1)
 condNumberPlot2(sum2normsLog,false, CondNum, MatrixNorm2)
 condNumSubPlot1(m, mi, sum1normsLog, sum2normsLog, CondNum, MatrixNorm1, MatrixNorm2)

 % specific plots

 m0plot(scalarsM0(:,1).^2,scalarsM0(:,2).*2)
 cnPlot( condNum1(:,1), condNum1(:,2),condNum2(:,1), condNum2(:,2))

 disp '** Done - List of Condition Numbers: '
 condNum1
 condNumCalc





 
function condNumberPlot1(y1,subplot, cn, mn)
 %CREATEFIGURE(Y1)
 %  Y1:  matrix of y data
 %  CN:  condition number
 %  MN:  matrix fro norm

  
 %  Auto-generated by MATLAB on 08-Oct-2010 17:03:02
  
 %% Create figure
 if ~subplot
    figure1 = figure;
 else
    figure1=gcf;  
 end

 %% Create axes
 axes1 = axes('Parent',figure1);
 title(axes1,['Decomposition of the Matrix Condition Number: ',num2str(cn),' Components - Norm of Matrix || A_n_+_1 || = ',num2str(mn)]);
 xlabel(axes1,'Dimension of Condition Number Natural Block Sub-Matrices');
 ylabel(axes1,'Log10 of Norm Component of Condition Number Squared');
 box(axes1,'on');
 hold(axes1,'all');
  
 %% Create mutliple lines using matrix input to plot
 plot1 = plot(y1);
 set(plot1(1),'Marker','o');
 set(plot1(2),'Marker','^');
 set(plot1(3),'Marker','v');
 set(plot1(4),'Marker','square');
 set(plot1(5),...
   'Color',[1 0 0],...
   'LineWidth',2,...
   'Marker','+');
  
 %% Create legend
 legend1 = legend(axes1,{'Log Norm M2 (submatrix An-1) Squared','Log Norm V1 Squared','Log Norm V2 Squared','Log Norm B0 Squared','Log Sum all Norms Squared'},'Location','SouthOutside');
 grid on



The next m-file finds a 100x100 matrix with a large cond number and saves it with the std naming format condMatrix + a 4 digit id

clear all

 dir = 'c:\MATXCOND\'
 fprefix = 'condMatrix' %%- used in MATXCOND program to read array
 fsufixNo='10003'

 d=100;
 large=10000000;
 trials=100000;

 p='fro'

 disp 'Finds a dim = ',num2str(d)
 disp 'Square 100 X 100 matrix with Large 2-norm condition number over: 10,000,000'

 ix=1:trials;

 for x=ix
     m=rand(d);
     if cond(m,p)>large
      break
     end
 end
 

 condMatrix=m;
 disp '---matrix:'
 condMatrix
 disp ' '
 disp '---Cond Number: '
 condNumber=cond(condMatrix,p)
 disp ' '
 disp '---Counter:'
 x
  
 % Use the, #### = serial num to identify matrix, start 10001
 f= [dir,fprefix,fsufixNo,'.mat']
 save(f ,fprefix)

end