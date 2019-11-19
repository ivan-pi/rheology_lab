
podatki = dlmread('1skupina_temp.csv',',',0,0)


omega = podatki(:,2)
G1 = podatki(:,3)
G2 = podatki(:,4)


%definicija funkcija G' in G''
G1func = @(g,lambda) g.*lambda.^2.*omega'.^2./(1+lambda.^2.*omega'.^2)
G2func = @(g,lambda ) g.*lambda.*omega'./(1+lambda.^2.*omega'.^2)



%vsota vrsti po g_i in lambda_i
G1func2 = @(g,lambda) sum(cell2mat(arrayfun(G1func,g,lambda,'UniformOutput',0)'),1)
G2func2 = @(g,lambda) sum(cell2mat(arrayfun(G2func,g,lambda,'UniformOutput',0)'),1)


%ta del je samo za risanje
omegaA = logspace(-2,3,124)';
G1funcA = @(g,lambda) g.*lambda.^2.*omegaA'.^2./(1+lambda.^2.*omegaA'.^2)
G2funcA = @(g,lambda) g.*lambda.*omegaA'./(1+lambda.^2.*omegaA'.^2)
G1func2A = @(g,lambda) sum(cell2mat(arrayfun(G1funcA,g,lambda,'UniformOutput',0)'),1)
G2func2A = @(g,lambda) sum(cell2mat(arrayfun(G2funcA,g,lambda,'UniformOutput',0)'),1)


%razlike med izracuna in izmerjeno vrednost
vectorG1 = @(x) (G1' - G1func2(x(1,:),x(2,:)));
vectorG2 = @(x) (G2' - G2func2(x(1,:),x(2,:)));


utez =  1 %utez v prid G' ali G''

%skupni vector razlik
vector = @(x) [vectorG1(x).^2,vectorG2(x).^2*utez]


%število Maxwellovih elementov
n = 6


%naredi resitev
resitev = lsqnonlin(vector,rand(2,n)*150,zeros(2,n),ones(2,n)*1000)


g = resitev(1,:);
lambda = resitev(2,:);

%malo posortiram rezultate
[lambda,red] = sort(lambda);

%nato jih se izpisemo
lambda = lambda*1
g = g(red)

subplot(1,2,1)
loglog(omega,G1,'r.',omegaA,G1func2A(g,lambda),'r-',omega,G2,'b.',omegaA,G2func2A(g,lambda),'b-')
grid on
xlabel('\omega [rad/s]')
legend('G'' meritev','G'' model','G'''' meritev','G'''' model','Location','Best')

subplot(1,2,2)
loglog(lambda,g,'b.-')
xlabel('\lambda')
ylabel('g')
grid on
