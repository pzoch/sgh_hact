% na podstawie kodu Bena Molla z https://benjaminmoll.com/codes/

clear all; 
clc; 
close all;

tic; % pomiar czasu
 
s = 2; % parameter funkcji uzytecznosci; u(c) = (c^(1-s))/(1-s) 

tau = 0.0;
T   = 0.0;
rho = 0.05;
yL = .1;
yH = .2;
y = [yL,yH];
lambdaL = 1.2;
lambdaH = 1.2;
lambda = [lambdaL,lambdaH];

r = 0.03;

amin = -0.15; 
amax = 5;

% tworzymy siatke aktywow
I= 1000;
a = linspace(amin,amax,I)';
da = (amax-amin)/(I-1);

aa = [a,a];
yy = ones(I,1)*y;


maxit= 100;
crit = 10^(-6);
Delta = 1000;

dVf = zeros(I,2);
dVb = zeros(I,2);
c = zeros(I,2);

Aswitch = [-speye(I)*lambda(1),speye(I)*lambda(1); speye(I)*lambda(2),-speye(I)*lambda(2)];


% zgadniecie v0

v0(:,1) = ((1-tau)*y(1)+T + r.*a).^(1-s)/(1-s)/rho;
v0(:,2) = ((1-tau)*y(2)+T + r.*a).^(1-s)/(1-s)/rho;

v = v0;
for n=1:maxit
    V = v;
    V_n(:,:,n)=V;
    
    %pochodne
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,:) = ((1-tau)*y+T + r.*amax).^(-s);
    
    dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVb(1,:) = ((1-tau)*y+T + r.*amin).^(-s); 
    
    I_concave = dVb > dVf;
    
    % konsumpcja i oszczednosci w kazdym wariancie
    cf = max(dVf,10^(-10)).^(-1/s);
    ssf = (1-tau)*yy+T + r.*aa - cf;
    
    cb = max(dVb,10^(-10)).^(-1/s);
    ssb = (1-tau)*yy+T + r.*aa - cb;
   
    c0 = (1-tau)*yy+T + r.*aa;
    
  
    If = ssf > 0;
    Ib = ssb < 0; 
    I0 = (1-If-Ib);
    
    c = cf.*If + cb.*Ib + c0.*I0;
    ss = ssf.*If + ssb.*Ib;
    u = c.^(1-s)/(1-s);
    
    % macierz przejscia
    X = -min(ssb,0)/da;
    Y = -max(ssf,0)/da + min(ssb,0)/da;
    Z = max(ssf,0)/da;
    
    A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
    A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
    
    if max(abs(sum(A,2)))>10^(-9)
       disp('niewlasciwa macierz')
       break
    end
    
    B = (1/Delta + rho)*speye(2*I) - A;

    u_stacked = [u(:,1);u(:,2)];
    V_stacked = [V(:,1);V(:,2)];
    
    b = u_stacked + V_stacked/Delta;
    V_stacked = B\b; %uklad rownan liniowych
    
    V = [V_stacked(1:I),V_stacked(I+1:2*I)];
    
    Vchange = V - v;
    v = V;

    dist(n) = max(max(abs(Vchange)));
    if dist(n)<crit
        disp('Zbieznosc funkcji wartosci, iteracja = ')
        disp(n)
        break
    end
end
toc;


% rozklad %

AT = A';
b = zeros(2*I,1);


i_fix = 1;
b(i_fix)=.1;
row = [zeros(1,i_fix-1),1,zeros(1,2*I-i_fix)];
AT(i_fix,:) = row;

%uklad rownan liniowych
gg = AT\b;
g_sum = gg'*ones(2*I,1)*da;
gg = gg./g_sum;

g = [gg(1:I),gg(I+1:2*I)];

check1 = g(:,1)'*ones(I,1)*da;
check2 = g(:,2)'*ones(I,1)*da;

c_stack = [c(:,1);c(:,2)];
g_stack = [c(:,1);c(:,2)];
a_stack = [a;a];
sd_c = std(c_stack,g_stack);
sd_a = std(a_stack,g_stack);

amax1 = 0.6;
amin1 = amin-0.01;

figure(1)
set(gca,'FontSize',16)
h1 = plot(a,c(:,1),'b',a,c(:,2),'r','LineWidth',2);
legend(h1,'c_L(a)','c_H(a)','Location','NorthEast')
text(-0.155,-.105,'$\underline{a}$','FontSize',16,'interpreter','latex')
xlabel('Aktywa, $a$','interpreter','latex')
ylabel('Konsumpcja, $c_i(a)$','interpreter','latex')
xlim([amin1 amax1])

figure(2)
set(gca,'FontSize',16)
h1 = plot(a,ss(:,1),'b',a,ss(:,2),'r','LineWidth',2);
legend(h1,'s_L(a)','s_H(a)','Location','NorthEast')
text(-0.155,-.105,'$\underline{a}$','FontSize',16,'interpreter','latex')
xlabel('Aktywa, $a$','interpreter','latex')
ylabel('Oszczednosci, $c_i(a)$','interpreter','latex')
xlim([amin1 amax1])

figure(3)
set(gca,'FontSize',16)
h1 = plot(a,g(:,1),'b',a,g(:,2),'r','LineWidth',2);
legend(h1,'g_L(a)','g_H(a)')
text(-0.155,-.12,'$\underline{a}$','FontSize',16,'interpreter','latex')
xlabel('Aktywa, $a$','interpreter','latex')
ylabel('Gestosc, $g_i(a)$','interpreter','latex')
xlim([amin1 amax1])
ylim([0 4])