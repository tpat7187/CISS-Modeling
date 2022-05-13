clc;
clear all;
close all;

total_sites = 4
A = zeros(1,total_sites*2) ;
electron_sites = 2
A(1:2:2*electron_sites) = 1

lap = 3;
ions_per_lap = ceil(total_sites/lap);
a = 1; %radius
c = 1; %length
%laps*ionsperlap must equal total_sites

jj = sqrt(-1);

%unique permutations
n = size(A,2);
k = sum(A==1);
C = nchoosek(1:n,k);
m = size(C,1); 
perms = zeros(m,n);
perms(repmat((1-m:0)',1,k)+m*C) = 1;

for i = 1:length(perms)
    cat_index(i,1) = indexy(perms(i,:));
end

run("VSetup.m")

%electron-electron intractions 
U = 0;

%energy
t = 1;

%Spin Orbit Interaction Parameter
lambda = 0;



H = sparse(length(perms),length(perms));
TestU = sparse(length(perms),length(perms));
TestD = sparse(length(perms),length(perms));
TestS = sparse(length(perms),length(perms));


%density operators
for i = 1:2:total_sites*2
EE = sparse(diag(perms(:,i))*diag(perms(:,i+1)))*U;
H = H + EE;
end


%NN
for i = 1:total_sites-1
for ifw = 1:length(perms)
A = perms(ifw,:);
[a,H] = hopping(A,i,i+1,"up","up",cat_index,H,-t);
[a,H] = hopping(A,i,i+1,"down","down",cat_index,H,-t);
end
end


%NN boundary conditions
for ifw = 1:length(perms)
A = perms(ifw,:);
[a,H] = hopping(A,total_sites,1,"up","up",cat_index,H,-t);
[a,H] = hopping(A,total_sites,1,"down","down",cat_index,H,-t);
end



%next NN
for i = 1:total_sites-2
for ifw = 1:length(perms)
A = perms(ifw,:);
%V1
[a,H] = hopping(A,i,i+2,"up","down",cat_index,H,V1(i,i+2)*lambda*jj);
[a,H] = hopping(A,i,i+2,"down","up",cat_index,H,V1(i,i+2)*lambda*jj);
%V2
[a,H] = hopping(A,i,i+2,"up","down",cat_index,H,V2(i,i+2)*lambda*(-jj)*jj);
[a,H] = hopping(A,i,i+2,"down","up",cat_index,H,V2(i,i+2)*lambda*(-1)*(-jj)*jj);
%V3
[a,H] = hopping(A,i,i+2,"up","up",cat_index,H,V3(i,i+2)*lambda*jj);
[a,H] = hopping(A,i,i+2,"down","down",cat_index,H,V3(i,i+2)*lambda*(-1)*jj);
end
end


%Next NN Boundary Conditions
for ifw = 1:length(perms)
A = perms(ifw,:);
%V1
[a,H] = hopping(A,total_sites-1,1,"up","down",cat_index,H,V1(total_sites-1,1)*lambda*jj);
[a,H] = hopping(A,total_sites,2,"up","down",cat_index,H,V1(total_sites,2)*lambda*jj);
[a,H] = hopping(A,total_sites-1,1,"down","up",cat_index,H,V1(total_sites-1,1)*lambda*jj);
[a,H] = hopping(A,total_sites,2,"down","up",cat_index,H,V1(total_sites,2)*lambda*jj);
%V2
[a,H] = hopping(A,total_sites-1,1,"up","down",cat_index,H,V2(total_sites-1,1)*lambda*(-jj)*jj);
[a,H] = hopping(A,total_sites,2,"up","down",cat_index,H,V2(total_sites,2)*lambda*(-jj)*jj);
[a,H] = hopping(A,total_sites-1,1,"down","up",cat_index,H,V2(i,i+2)*lambda*(-1)*(-jj)*jj);
[a,H] = hopping(A,total_sites,2,"down","up",cat_index,H,V2(i,i+2)*lambda*(-1)*(-jj)*jj);
%V3
[a,H] = hopping(A,total_sites-1,1,"up","up",cat_index,H,V3(i,i+2)*lambda*jj);
[a,H] = hopping(A,total_sites,2,"up","up",cat_index,H,V3(i,i+2)*lambda*jj);
[a,H] = hopping(A,total_sites-1,1,"down","down",cat_index,H,V3(i,i+2)*lambda*(-1)*jj);
[a,H] = hopping(A,total_sites,2,"down","down",cat_index,H,V3(i,i+2)*lambda*(-1)*jj);
end

H_ctr = ctranspose(H);
H_final = H_ctr + H;


%Current Term 
for i = 1:total_sites-1
for ifw = 1:length(perms)
A = perms(ifw,:);
[a,TestU] = hopping(A,i,i+1,"up","up",cat_index,TestU,jj);
[a,TestD] = hopping(A,i,i+1,"down","down",cat_index,TestD,jj);

end
end

%Current Term Boundary Conditions
for i = 1:total_sites
for ifw = 1:length(perms)
A = perms(ifw,:);
[a,TestU] = hopping(A,total_sites,1,"up","up",cat_index,TestU,jj);
[a,TestD] = hopping(A,total_sites,1,"down","down",cat_index,TestD,jj);
end
end

%Spin Density Matrix
for i = 1:total_sites-1
for ifw = 1:length(perms)
A = perms(ifw,:);
[a,TestS] = hopping(A,i,i,"up","up",cat_index,TestS,1/2);
end
end


J_CtrU = ctranspose(TestU); 
J_finalU = TestU + J_CtrU;
J_CtrD = ctranspose(TestD);
J_finalD = TestD + J_CtrD;


J = linspace(-2,2,100);
Expectation_JU = zeros(1,length(J));
Expectation_JD = zeros(1,length(J));

for i = 1:length(Expectation_JU)
H_final_J = H_final + J(1,i)*(J_finalU);
[E_Vectors, E_Values] = eigs(H_final_J,2,'smallestreal');
Ground_State = E_Vectors(:,1);
Expectation_JU(1,i) = ctranspose(Ground_State)*J_finalU*Ground_State;
end

for i = 1:length(Expectation_JD)
H_final_J = H_final + J(1,i)*(J_finalD);
[E_Vectors, E_Values] = eig(full(H_final_J));
[E_values] = eig(full(H_final_J));
Ground_State = E_Vectors(:,1);
Expectation_JD(1,i) = ctranspose(Ground_State)*J_finalD*Ground_State;
end

%expectation value plot
figure(1)
hold on
plot(J,Expectation_JU)
plot(J,Expectation_JD)
title('Current Expectation Value vs Constant');
xlabel('J / t'); 
ylabel('Expectation value');
axis([-2 2 -5 5])
legend('JU','JD')
hold off



%function inputs: state,starting site, ending site, spin1, spin2, index,
%hamiltonian matrix, electron-electron interaction energy, constant,
%intermediate site (only valid for next NN hopping







