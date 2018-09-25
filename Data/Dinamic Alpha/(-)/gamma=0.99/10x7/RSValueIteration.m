%% Risk Value Iteration Algorithm implemented in River Domain
% Responsible for make value iteration based on Risk Escalar Function
% proposed by Mihatsch and Neuneier (2002)
% coded by Igor Oliveira Borges | advisory by Karina Valdivia Delgado and Valdinei Freire da Silva 
classdef RSValueIteration < handle
    properties
        epsilon
        gamma
        alpha
        Q
        V
        Pi
        M
        k
        A
        counter
        residua
        inicial
        alph
        epsil
        risk
        neighbour
    end
methods

%Constructor of agent           
function G = RSValueIteration(A,epsilon,gamma,alpha,M,k)
    G.gamma = gamma;
    G.A = A;
    if (G.A.isDinamicAlpha)
        G.alpha = A.ini_alpha;
    else
        G.alpha = alpha;
    end
    if (G.A.isDinamicRisk)
        G.k = A.ini_k;
        G.epsilon = A.ini_epsilon;
        G.alpha = A.ini_alpha;
    else
        G.k= k;
        G.epsilon = epsilon;
    end
    G.M = M;
    G.V = zeros(1,G.M.St);
    G.Pi = zeros(1,G.M.St);
    G.Q = zeros(G.M.St,G.M.A);
    G.residua = zeros(G.A.limit,1);
    G.inicial = zeros(G.A.limit,1);
    G.alph = zeros(G.A.limit,1);
    G.epsil = zeros(G.A.limit,1);
    G.risk = zeros(G.A.limit,1);
    G.neighbour = zeros(G.A.limit,1);
end

function [coord] = state2coord(G,s,Nx,Ny)
   s = s-1;
   coord(2) = mod(s,Ny);
   coord(1) = mod( floor(s/Ny) , Nx);
   coord = coord + 1;
end
        
function s = coord2state(G,coord,Nx,Ny)
    coord = coord - 1;
    s = coord(2) + coord(1)*Ny + 1;
end

function c = imprimir(G)
for s=1:G.M.St
maxV = max(G.Q(s,:)); %find max value in Q(s,:)
indice = find(G.Q(s,:)==maxV);
if size(indice)>0
    G.Pi(s) = indice(1); % put into policy pi(s) = max value of Q in s
end
if G.M.transict(s,:)==0
    if s == G.M.St-G.M.Ny+1
        G.Pi(s) = 0;
    else
        G.Pi(s) = -1; %if not mark with zero cause there is no valid action
    end
end
end  
vetor = zeros(G.M.Nx,G.M.Ny);
c = zeros(G.M.Nx,G.M.Ny);

for i=1:G.M.St
    coord = G.state2coord(i,G.M.Nx,G.M.Ny);
    vetor(coord(1),coord(2)) = G.Pi(i);
end

for i=1:G.M.St
coord = G.state2coord(i,G.M.Nx,G.M.Ny);

switch vetor(coord(1),coord(2))
    case -1 
        c(coord(1),coord(2))=43;
    case 0
        c(coord(1),coord(2))=111;  
    case 1
        c(coord(1),coord(2))=8593;
    case 2
        c(coord(1),coord(2))=8595;
    case 3
        c(coord(1),coord(2))=8594;
	case 4	
	    c(coord(1),coord(2))=8592;
    otherwise
        c(coord(1),coord(2))=0;
end
end  
i=G.M.Ny;
fileID = fopen(G.A.print,'w');
encoded='';
while (i>=1)
    linha = (c(:,i))';
    encoded = unicode2native(char(linha), 'UTF-8');
    fwrite(fileID,encoded,'uint8');
    fprintf(fileID,'\n');
    i=i-1; 
end
fprintf(fileID,'\n');
fwrite(fileID,num2str(G.counter),'uint8');
fclose(fileID);
end

function y = escalar(G,x,k)
if x > 0
    y = (1 - k) *x;
else
    y = (1 + k) *x;
end
end

%Method for update Q values using Value Iteration
function RSValuation(G)
%inicialization
V0 = zeros(1,G.M.St); V1 = zeros(1,G.M.St);
G.Q(:,:) = zeros(G.M.St,G.M.A); old_Q = zeros(G.M.St,G.M.A);
G.counter = 0;
for i=1:G.M.St
 for a=1:G.M.A
    if sum(G.M.T(i,a,:))==0 && i~=(G.M.Nx-1)*G.M.Ny+1
        G.Q(i,a)=-Inf;
    end
 end
end
G.Q((G.M.Nx-1)*G.M.Ny+1,:) = ones; residual = Inf(1,G.M.St); old_residual = zeros(1, G.M.St);
while max(residual(:)) > G.epsilon %iterates VI while max(V1-V0)>epsilon
    old_residual(:) = residual(:);
    old_Q(:) = G.Q(:);
    for s=1:G.M.St           %states
        for a=1:G.M.A        %actions
            j = (nonzeros(find(G.M.T(s,a,:))))';
            if s~=(G.M.Nx-1)*G.M.Ny+1 && numel(j)>0 
                soma = 0;        %clean acumulated value
                for p=1:numel(j)
                    delta = G.M.recomp(s,a) + (G.gamma * (max(old_Q(j(p),:))) - old_Q(s,a));
                    rho = G.escalar(delta, G.k);
                    soma = soma + G.M.T(s,a,j(p)) * rho;
                end
                G.Q(s,a) = G.Q(s,a) + G.alpha * soma;    
            end
        end
        V1(s) = max(G.Q(s,:)); %find VI
        residual(s) = abs((V1(s)-V0(s))/V0(s)); %calculates residual = |V1-V0|
    end
    G.counter = G.counter + 1;
    if (G.A.isDinamicAlpha)
        if (max(residual(:)) > max(old_residual(:)))
           G.alpha = G.A.alpha; 
        end
    end
    
    G.residua(G.counter) = max(residual(:));
    G.inicial(G.counter) = V1(1);
    G.neighbour(G.counter) = V1((G.M.Nx-1)*G.M.Ny+2);
    G.epsil(G.counter) = G.epsilon;
    
    
    G.alph(G.counter) = G.alpha;
    G.risk(G.counter) = G.k;
    
    if (G.counter == G.A.limit)
        disp('EXIT BY LIMIT EXCEPTION');
        break;
    end
    
    V0(:) = V1(:); %copy V1 to V0
    V1(:) = zeros; %clean V1
    if (G.A.isDinamicRisk)
        if (max(residual(:)) <= G.epsilon && G.k~=G.A.k && G.epsilon~=G.A.epsilon)
            G.k = G.A.k;
            G.epsilon = G.A.epsilon;
            G.alpha = G.A.alpha;
        end
    end
end
G.imprimir();  %print policy    
end
end
end

