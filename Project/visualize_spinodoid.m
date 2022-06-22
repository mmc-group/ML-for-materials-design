clc; clear all; close all;

%res: resolution; bigger means more accurate geometry
%However, computations increase with cubic order of resolution
res = 51 %resolution

%density: density of the solid; between 0 and 1
%Geometry for density<=0.15 doesn't make sense because of disconnected domains
%For meaningful geometries, use density>=0.3
density = 0.5 % To-do

%beta: controls the pore size
%Higher values means smaller pore size
%A crude estimate of pore size = 1/beta
%Higher beta requires higher resolution in order to capture microstructures
beta = 5 % To-do

%thetaX, thetaY, thetaZ: anisotropy parameters
%each parameter is an angle between [0,90] degrees.
%For sensible geometries, use either 0 or between [15,90] degrees.
%In other words, don't use values between 0 and 15 degrees.
%e.g., thetaX = 0, thetaY = 15, thetaZ = 45 | good example
%e.g., thetaX = 10, thetaY = 15, thetaZ = 45 | bad example
%e.g., thetaX = 15, thetaY = 15, thetaZ = 45 | good example
% When no angle is zero: cubic topology
% When one angle is zero: columnar topology
% When two angles are zero: lamellar topology
% When either one of the angles is 90: isotropic topology
thetaX = 0 % To-do
thetaY = 0 % To-do
thetaZ = 90 % To-do

%name of output file
name = 'spinodoid_demo.stl'



%----------------------------------------------------------------
% Code: Do not edit below this line

disc=linspace(0,1,res);
x=zeros(res,res,res);
y=zeros(res,res,res);
z=zeros(res,res,res);
f=zeros(res,res,res);

for i=1:res
    for j=1:res
        for k=1:res
            x(i,j,k)=disc(i);
            y(i,j,k)=disc(j);
            z(i,j,k)=disc(k);
        end
    end
end

N=10000; 
g = 2*pi*rand(1,N);
n = zeros(N,3);
for w=1:N
    while(true)
        vec = 2*rand(1,3)-1;
            if(norm(vec)>1)
                continue
            end
        vec = vec/norm(vec);
        t1 = (abs(dot(vec,[1,0,0])) > cosd(thetaX));
        t2 = (abs(dot(vec,[0,1,0])) > cosd(thetaY));
        t3 = (abs(dot(vec,[0,0,1])) > cosd(thetaZ));
        if(t1 || (t2 || t3))
            n(w,:)=vec;
            break
        else
            continue
        end
    end
end

for w=1:N
    dp = n(w,1)*x + n(w,2)*y + n(w,3)*z;
    f = f + sqrt(2/N)*cos(2*pi*beta*dp+g(w));
end

levelset = sqrt(2)*erfinv(2*density-1);

isos = isosurface(x,y,z,f,levelset);
isoc = isocaps(x,y,z,f,levelset,'enclose','below');
out.faces=[isos.faces; isoc.faces + length(isos.vertices(:,1))];
out.vertices=[isos.vertices;isoc.vertices];

stl = triangulation(out.faces,out.vertices);
stlwrite(stl, name)
