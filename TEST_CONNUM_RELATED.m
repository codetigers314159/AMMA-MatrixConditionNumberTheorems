% m=rand(50);
% % Generate a single random number between -2 and 2
% r = unifrnd(-2, 2);
 % d=100;
 % N=d;
% 
% MU = rand(d, d);
% disp(MU);
% figure
% plot(MU)
% figure
% plot(inv(MU))

%complex mat from real mat and its inverse
d=100;
N=d;


% User prompt
defaultThreshold = 1e6;
prompt = 'Enter minimum Frobenius condition number threshold (default 1e6): ';
userInput = input(prompt);
if isempty(userInput)
    condThreshold = defaultThreshold;
else
    condThreshold = userInput;
end

% Matrix search
d = 100;
trials = 1000000;
found = false;

disp(['Searching for a random ', num2str(d), 'x', num2str(d), ' matrix with Frobenius cond > ', num2str(condThreshold)])

for x = 1:trials
    M = rand(d);
    if cond(M, 'fro') > condThreshold
        disp('Matrix found');
        found = true;
        break;
    end
end

if ~found
    disp(['No matrix found with cond > ', num2str(condThreshold), ' after ', num2str(trials), ' trials.']);
    return;
end
   

M1=M %rand(d,d);

invM1=inv(M1);

SumInvMat=zeros(d,d);
SumInvMatExp=zeros(d,d);

for ix1=1:N
    for ix2=1:N


        z = complex( M1(ix1,ix2),invM1(ix1,ix2) );
        
        SumInvMat(ix1,ix2) = z;
       
    
    end
end

cn=num2str( floor(cond(SumInvMat,"fro")))

figure

plot (SumInvMat )
xlabel("Re(z)")
ylabel("Im(z)")
title ("EXP of the Sum of Real Mat and its Inv as a Complex Matrix-Cond# " + cn)
grid on

for ix1=1:N
    for ix2=1:N


        z = exp(complex( M1(ix1,ix2),invM1(ix1,ix2)) );
        
        SumInvMatExp(ix1,ix2) = z;
       
    
    end
end



% To plot a complex matrix in MATLAB, you can visualize the real and imaginary parts separately or represent the complex numbers on the complex plane. Below are examples of how to achieve this:
% 
% 1. Plot Real vs Imaginary Parts (Complex Plane)
% You can treat the real part as the x-axis and the imaginary part as the y-axis.


% Flatten the matrix for plotting
Z_flat = SumInvMat(:);

% Plot on the complex plane
figure;
plot(real(Z_flat), imag(Z_flat), 'o');
xlabel('Real Part');
ylabel('Imaginary Part');
title("Complex Matrix on Complex Plane-Cond# " + cn);
grid on;


% Flatten the matrix for plotting
Z_flat = SumInvMatExp(:);

% Plot on the complex plane
figure;
plot(real(Z_flat), imag(Z_flat), 'o');
xlabel('Real Part');
ylabel('Imaginary Part');
title("Complex Matrix on Complex Plane-Cond# " + cn);
grid on;

% 2. Visualize Real and Imaginary Parts as Separate Surfaces
% You can use surf or imagesc to visualize the real and imaginary parts as 2D or 3D plots.
% Example complex matrix
Z = Z_flat;

% Plot real part
figure;
%subplot(1, 2, 1);
surf(real(SumInvMat));

title("Real Part-Cond# " + cn);
xlabel('X');
ylabel('Y');
zlabel('Value');

figure;
%subplot(1, 2, 1);
surf(imag(SumInvMat));

title("Imaginary Part-Cond# " + cn);
xlabel('X');
ylabel('Y');
zlabel('Value');

% % Plot imaginary part
% subplot(1, 2, 2);
% surf(imag(Z));
% title('Imaginary Part');
% xlabel('X');
% ylabel('Y');
% zlabel('Value');


% 3. Use imagesc for Heatmap Representation
% If you prefer a heatmap-style visualization:
% Matlab% Example complex matrix

% Plot real part

figure;
title ("Sum of Real Mat and its Inv as a Complex Matrix-Cond# " + cn)
subplot(1, 2, 1);
imagesc(real(Z));
colorbar;
title("Real Part - Cond# " + cn);

% Plot imaginary part
subplot(1, 2, 2);
imagesc(imag(Z));
colorbar;
title('Imaginary Part');


% These approaches allow you to explore the structure of your complex matrix effectively. Choose the one that best fits your needs!
% 

% figure
% for ix1=1:N
%     for ix2=1:N
% 
% 
%         z = complex( M1(ix1,ix2),invM1(ix1,ix2) );
%         hold on
%         polarplot(angle(z),abs(z),"o")
%     end
% end
% 
% figure;
% 
% for ix1=1:N
%     for ix2=1:N
% 
% 
%         x= M1(ix1,ix2);
%         y= invM1(ix1,ix2);
% 
%        z = complex( x,y );
% 
%         subplot(2,1,1); plot(x, abs(z));
%         subplot(2,1,2); plot(x, angle(z));
%     end
% end

% 
% for ix1=1:N
%     for ix2=1:N
% 
% 
%         x= M1(ix1,ix2);
%         y= invM1(ix1,ix2);
% 
%         subplot(2,1,1); plot(x, abs(y));
%         subplot(2,1,2); plot(x, angle(y));
%     end
% end



%-- other tests

d=N;
M3 = M1;
invM3=invM1;
M1Stats=zeros(d,d);

for ix1=1:N
   

        %z = complex( M1(ix1,ix2),invM1(ix1,ix2));
     
        M1Stats(ix1,1)=max(invM3(ix1,:))
        M1Stats(N,ix1)=min(invM3(ix1,:))
       
end

figure
grid on
title("Max and Min for Inverse Matrix M1Stats")
plot(M3Stats)