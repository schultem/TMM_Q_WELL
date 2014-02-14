%Bound Energies in a Quantum Well Calculated Using TMM

%electron mass is taken as m = 0.067*m0 throughout the entire region
%(effective mass of Gallium Arsenide is 0.067)
m=.067*9.10938215*10^-31;
%An additional vector input could be defined for effective mass mi 
%throughout the homogenous potentials.

%reduced Plancks constant
hbar=1.054571628*10^-34;

%The homogenous potentials, the plots below help portray these graphically
V=[0.1 -.45];
for l = 1:28
    V=[V (-.45+(.45*l/28))];
end
N = length(V);

%The physical location in length for the 1D problem
z=-7*10^-9;
for l = 1:(N-1)
    z = [z (-7*10^-9 + (7*10^-9)*l/(N-1))];
end

%Develop the set of energies to use in finding the bound states
E=-.45 + .45/100;
for l=2:100*N
    E = [E (-.45 + .45*l/(100*N))];
end
T11=0;
for e=1:100*N
    %avoid states involving calcualtions of E=Vi, one of the estimated
    %potentials, which would result in ki=0 and the local probability
    %current density in the Vi segment would always equal zero.
    singular=0;
    for j = 1:N
        if E(e) == V(j)
            singular=1;
        end
    end
    if singular == 0,
        %k1 assumes E(1) < V(1) (bound)
        k=1i*sqrt(2*m*(V(1)-E(e))*1.60217646*10^-19)/hbar;
        %Developing the rest of the k vector k=[k1...ki...kn]
        for l = 2:N
            if E(e) > V(l)
                k = [k (sqrt(2*m*(E(e)-V(l))*1.60217646*10^-19)/hbar)];
            end
            if E(e)<V(l)
                k = [k (1i*sqrt(2*m*(V(l)-E(e))*1.60217646*10^-19)/hbar)];
            end
        end
        %Each Mi of the total transfer matrix T is cascaded. Effective electron mass m is
        %simplified out of the expressions because it is the same on the
        %numerator and the denominator in this example.  If m(i) is not the 
        %same as m(i+1) the 2 x 2 M matrix becomes(with the addition of an added 
        %m vector input):
        %Mi=[ ((1/2 + k(l+1)*m(i)/(2*k(l)*m(i+1))*exp(1i*(k(l+1)-k(l))*z(l+1)))  ...
        %     ((1/2 - k(l+1)*m(i)/(2*k(l)*m(i+1))*exp(-1i*(k(l+1)+k(l))*z(l+1)));...
        %     ((1/2 - k(l+1)*m(i)/(2*k(l)*m(i+1))*exp(1i*(k(l+1)+k(l))*z(l+1)))  ...
        %     ((1/2 + k(l+1)*m(i)/(2*k(l)*m(i+1))*exp(-1i*(k(l+1)-k(l))*z(l+1))) ];
        
        M1=[ ((1/2 + k(2)/(2*k(1)))*exp(1i*(k(2)-k(1))*z(2)))  ...
             ((1/2 - k(2)/(2*k(1)))*exp(-1i*(k(2)+k(1))*z(2)));...
             ((1/2 - k(2)/(2*k(1)))*exp(1i*(k(2)+k(1))*z(2)))  ...
             ((1/2 + k(2)/(2*k(1)))*exp(-1i*(k(2)-k(1))*z(2))) ];
        T=M1;
        for l=2:(N-1),
            Mi=[ ((1/2 + k(l+1)/(2*k(l)))*exp(1i*(k(l+1)-k(l))*z(l+1)))  ...
                 ((1/2 - k(l+1)/(2*k(l)))*exp(-1i*(k(l+1)+k(l))*z(l+1)));...
                 ((1/2 - k(l+1)/(2*k(l)))*exp(1i*(k(l+1)+k(l))*z(l+1)))  ...
                 ((1/2 + k(l+1)/(2*k(l)))*exp(-1i*(k(l+1)-k(l))*z(l+1))) ];
             T=T*Mi;
        end
        T11=[T11 T(1)];
    end
end

plot(1:2997,T11)
hold
plot(1:2997,0,'r')
xlabel('Energy Levels [0.45 eV/x]')
ylabel('T11 Eigen-energies')
title('T11=0 indicates a Bound Energy State in a Quantum Well')
%zero crossings for N=30: E(2966)=-0.0051  E(1576)=-0.2136

%interpolate potential vector using the z vector
Vhomogenous=V(1);
for l=2:N
    for m=1:10
        Vhomogenous = [Vhomogenous V(l)];
    end
end
Vinterpolated=[0.1 -0.45];
for l = 1:28
    Vinterpolated=[Vinterpolated (-.45+(.45*(l-1)/28))];
end

x=-7*10^-9:(7/290)*10^-9:0;
plot(x,Vhomogenous,'b')
hold
x=[-7*10^-9 -7*10^-9];
x=[x (-7*10^-9):(7/29)*10^-9:0];
plot(x,[Vinterpolated 0 0],'r')
x=-7*10^-9:(7/290)*10^-9:0;
plot(x,-0.2136,'b')%E1
plot(x,-0.0051,'b')%E2
xlabel('z position [nm]')
ylabel('Energy [eV]')
title('Bound Energies in a Quantum Well Calculated Using TMM')

%Exact answers
plot(x,-0.2099,'r')%E1
plot(x,-0.0064,'r')%E2




    
    
    
    
    
    
    
    
    
    
    