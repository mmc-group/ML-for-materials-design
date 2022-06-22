clc; clear all; close all;

% Note: this code only allows visualizing ONE stiffness at a time.

% To-do: enter the 9-dimensional stiffness vector between single [] brackets
S = [0.3068, 0.0560, 0.1079, 0.2979, 0.1073, 0.6205, 0.1163, 0.1183, 0.0498]


%------------------------------------------------------------------------
% Code: Do not edit anything below this line

visualize_2d(S)

function [] = visualize_2d(S)

C = zeros(6,6);
C(1,1)=S(1);
C(1,2)=S(2);
C(1,3)=S(3);

C(2,1) = C(1,2);
C(2,2) = S(4);
C(2,3) = S(5);

C(3,1) = C(1,3);
C(3,2) = C(2,3);
C(3,3) = S(6);

C(4,4) = S(7);
C(5,5) = S(8);
C(6,6) = S(9);

Cvec = reshape(C,6,6)';
S = getCompliance(Cvec);

Emax = max([1/S(1,1,1,1),1/S(2,2,2,2),1/S(3,3,3,3)]);

figure
elastSurf3(S,0)
xlim(Emax*[-1,1])
ylim(Emax*[-1,1])
zlim(Emax*[-1,1])

set(gca,'FontSize',18)
xlabel('x')
ylabel('y')
zlabel('z')

end

function [] = elastSurf3(S,bounds)
n=200;
phi = linspace(0,pi,n);
theta = linspace(0,2*pi,n);
[x,y,z,E] = getE_3d(S,phi,theta);
surf(x,y,z,E,'EdgeColor','none');
minE = min(min(E));
maxE = max(max(E));
caxis([minE,maxE]);
camlight
material dull
hold on
axis equal
set(gca,'DataAspectRatio',[1 1 1])
colormap jet
end




function [x,y,z,E] = getE_3d(S,phi,theta)
E=zeros(length(phi),length(theta));
x=E;
y=E;
z=E;
for p = 1:length(phi)
   for t = 1:length(theta)
       [x(p,t),y(p,t),z(p,t),E(p,t)] = getE(S,phi(p),theta(t));
   end
end

end


function [x,y,z,E] = getE(S,phi,theta)
d=[cos(theta)*sin(phi),sin(theta)*sin(phi),cos(phi)];
E = 1/doubleContraction(S,d);
x = E*d(1);
y = E*d(2);
z = E*d(3);
end


function [x] = doubleContraction(S,d)
x=0;
for i=1:3
       for j=1:3
           for k=1:3
               for l=1:3
                   x = x+S(i,j,k,l)*d(i)*d(j)*d(k)*d(l);
               end
           end
       end
end
end
   

function [S] = getCompliance(vec)


map = [...
    1,6,5;...
    6,2,4;...
    5,4,3];

Cv = reshape(vec,[6,6]);
Cv = 0.5*(Cv+Cv');
Sv = inv(Cv);

S=zeros([3,3,3,3]);

for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                vr = map(i,j);
                vc = map(k,l);
                
                factor=1;
                if(vr>3)
                    factor = 2*factor;
                end
                if(vc>3)
                    factor = 2*factor;
                end
                
                S(i,j,k,l) = 1/factor * Sv(vr,vc);
                    
            end
        end
    end
end

end
