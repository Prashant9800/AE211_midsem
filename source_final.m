airfoil = 'naca0009';                                                           %airfoilname
alpha =0.0;                                                                    %angle of attack
vstr = 10;                                                                      % free stream velocity                    
Clen = 10;                                                                      %chord length
pthick = 100;                                                                   %percent thickness

%% Opening file
filetype = '.txt';
airfoil = strcat(airfoil,filetype);
Cxy = readmatrix('naca0009');
x= Cxy(:,1);
y= Cxy(:,2);
x = Clen*x;
y = (Clen*pthick*y)/100;
%% number of panels
numPAN = length(x)-1;
disp(numPAN);
num_point = length(x);
%%  Edge length 
edgelen = zeros(numPAN,1);
edge = zeros(numPAN,1);
for i = 1:1:numPAN
    edgelen(i) = (x(i+1)-x(i))^2 + (y(i+1)-y(i))^2;                         %edgelength square    
end
edgelen = sqrt(edgelen);
%% COMPUTE GEOMETRIC VARIABLES
XC   = zeros(numPAN,1);                                                     % Initialize control point X-coordinate array
YC   = zeros(numPAN,1);                                                     % Initialize control point Y-coordinate array
theta = zeros(numPAN,1);                                                    % Initialize panelangle array
for i = 1:1:numPAN                                                          
    XC(i)   = 0.5*(x(i)+x(i+1));                                          
    YC(i)   = 0.5*(y(i)+y(i+1));                                          
    dx      = x(i+1)-x(i);                                                
    dy      = y(i+1)-y(i);                                                
    theta(i) = atan2(dy,dx);
end
if(theta<0)
theta = theta + 2*pi;
end
%% vortex panel method
S = edgelen;
Cn1 = zeros(numPAN,numPAN);                                                 %integral result array 
Cn2 = zeros(numPAN,numPAN);
CT1 = zeros(numPAN,numPAN);
CT2 = zeros(numPAN,numPAN);
%% ARRAY calculation coffeicents intialized
for i = 1:1:numPAN 
    for j = 1:1:numPAN
            if (i == j)
            Cn1(i,j) = -1;
            Cn2(i,j) = 1 ;
            CT1(i,j) = 0.5*pi;
            CT2(i,j) = 0.5*pi;
            else
        midpoint_x = (x(i)+x(i+1))*0.5;
        midpoint_y = (y(i)+y(i+1))*0.5;
        DX = (midpoint_x-x(j));
        DY = (midpoint_y-y(j));
       A = -DX*cos(theta(i))-DY*sin(theta(j));
       B = DX^2 +DY^2;
       C = sin(theta(i)-theta(j));
       D = cos(theta(i)-theta(j));   
       E = DX*sin(theta(j)) - DY*cos(theta(j));
       S = edgelen(i);
       F = log(1 + ((S^2+2*A*S)/B)) ;
       G = atan2(E*S,B+A*S);
       P = DX*(sin(theta(i) -2*theta(j))) + DY*(cos(theta(i)-2*theta(j))) ;
       Q = DX*cos(theta(i)-2*theta(j)) + DY*sin(theta(i)-2*theta(j));
       Cn2(i) =real(( D + (0.5*Q*F/S) + ((A*C + D*E)*G/S)));
       Cn1(i,j) = real(((0.5*D*F) + (C*G) -Cn2(i,j)));
       CT2(i,j) =real( C + ((0.5*P*F)/S) + ((A*D - C*E)*(G/S)));
       CT1(i,j) = real(0.5*C*F - D*G - CT2(i,j));
    end
    end
end
%% defining system of linear equation
coff_matrix1 = zeros(num_point,num_point);
coff_matrix2 = zeros(num_point-1,num_point);
velocity = zeros(num_point,1);
for i = 1:1:num_point-1
    for j = 1:1:num_point
       if (j==1)
         coff_matrix1(i,j) =   Cn1(i,j);
       elseif  (j==num_point)
         coff_matrix1(i,j) =   Cn2(i,j-1); 
       else 
           coff_matrix1(i,j) =   Cn1(i,j) + Cn2(i,j-1);
       end
    end
    velocity(i) = sin(theta(i)-alpha);
end
for i = 1:1:num_point-1
    for j = 1:1:num_point
       if (j==1)
         coff_matrix2(i,j) =   CT1(i,j);
       elseif  (j==num_point)
         coff_matrix2(i,j) =   CT2(i,j-1); 
       else 
           coff_matrix2(i,j) =   CT1(i,j) + CT2(i,j-1);
       end
    end
    velocity(i) = sin(theta(i)-alpha);
end
%% defining kutta condition
coff_matrix1(num_point,1) = 1;
coff_matrix1(num_point,num_point) = 1;
velocity(num_point) = 0; 
%% solving 
gamma = linsolve(coff_matrix1,velocity);
%% Solve for Gamma and velocity/pressure                           
Cp = zeros(1,num_point);
for i = 1:1:num_point-1
    V(i) = cos(theta(i)-alpha);
    for j = 1:1:num_point-1
        V(i) = V(i) + coff_matrix2(i,j)*gamma(j);
        Cp(i) = 1 - (V(i))^2;
    end
end
figure(1)
plot(x,Cp);
set(gca,'Ydir','reverse');
xlabel('x/c');
ylabel('Coefficient of Pressure');
grid on;
grid minor;