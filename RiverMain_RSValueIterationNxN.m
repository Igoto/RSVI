%% Main file to invoke River Domain
%% Using RS Value Iteration     
function RiverMain_RSValueIterationNxN
ks = [0.99 0.8 0.5 0 -0.5 -0.8 -0.99]; %sets risk values (-1,1) which k<0 is expected a prone policy and k>0 is averse in k=0 is equal VA

for w=1:length(ks)

    limit = 250000;
    
    isDinamicRisk = true; %set dinamic value
    
    isDinamicAlpha = false; %set dinamic alpha
    
    k = ks(w); %factor of risk in (-1,1) range -1 is prone attitude +1 is averse 0 is neutral
    
    alpha = (1 +abs(k))^-1; %stepsize

    epsilon = 0.00001; %precision to stop planning algorithm

    
    if isDinamicAlpha == true
        ini_alpha = 1; 
    else
        ini_alpha = alpha;  
    end
    
    
    %if Dinamic is enable initial value is used after next value is used
    if isDinamicRisk == true
        ini_k = 0;
        ini_epsilon = 0.1;
        ini_alpha = (1 +abs(ini_k))^-1;
    else
        ini_k = k;
        ini_epsilon = epsilon;
    end
        
    scenario = false; %true (+) scenario [1 in goal and 0 otherwise] or false (-) scenario [0 in goal and -1 otherwise]

    p = 0.2; %probability to move in river ProbMoveInsideRiver

    nA = 4; %number of actions (N,S,W,E)

    %size of grid instance
    Nx = 10; Ny = 7;

    gamma = 0.99; %discount parameter for convergence

    if (scenario)
       scene = '(+)';
    else
       scene = '(-)';
    end
    
    %directory for save files
    mkdir_if_not_exist(strcat('Data/',scene,'/','gamma=',num2str(gamma),'/',num2str(Nx),'x',num2str(Ny)));

    agentStorage =  strcat('Data/',scene,'/','gamma=',num2str(gamma),'/',num2str(Nx),'x',num2str(Ny),'/','Agent','_','RS','_','VI','_',num2str(Nx),'x',num2str(Ny),'_','k=',num2str(k),'_','p=',num2str(p),'_e=0.00001.mat'); %local for storage instance for agents simulated (could be reused in a new execution)

    print = strcat('Data/',scene,'/','gamma=',num2str(gamma),'/',num2str(Nx),'x',num2str(Ny),'/','Print','_','RS','_','VI','_',num2str(Nx),'x',num2str(Ny),'_','k=',num2str(k),'_','p=',num2str(p),'_e=0.00001.txt'); %local for storage instance for agents simulated (could be reused in a new execution)

    %vinculates set of information for agent
    A = setAgent(agentStorage, print, epsilon, alpha, gamma, k,  ini_k, ini_alpha, ini_epsilon, isDinamicRisk, isDinamicAlpha, limit);

    %vinculates all probabilities and transiction of river matrix
    M = setRiver(Nx, Ny, p, nA, scenario); 

    %instance
    agent = RSValueIteration(A, A.epsilon, A.gamma, A.alpha, M, A.k); %vinculates instance for RSValueIteration

    %time counting
    t=0; e=0;

    t = cputime;
    agent.RSValuation(); %call for RSValueIteration

    e = cputime - t;

    %save agent after planning
    save(A.agentStorage,'agent','e','M','A');
    end
end

function mkdir_if_not_exist(dirpath)
    if dirpath(end) ~= '/', dirpath = [dirpath '/']; end
    if (exist(dirpath, 'dir') == 0), mkdir(dirpath); end
end

function A = setAgent(agentStorage, print, epsilon, alpha, gamma, k, ini_k, ini_alpha, ini_epsilon, isDinamicRisk, isDinamicAlpha, limit)
    A.agentStorage = agentStorage;
    A.print = print;
    A.epsilon = epsilon;
    A.alpha = alpha;
    A.gamma = gamma;
    A.k = k;
    A.ini_k = ini_k;
    A.ini_alpha = ini_alpha;
    A.ini_epsilon = ini_epsilon;
    A.isDinamicRisk = isDinamicRisk;
    A.isDinamicAlpha = isDinamicAlpha;
    A.limit = limit;
end