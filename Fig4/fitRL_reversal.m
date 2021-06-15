%%%Minimizes negative llh 10 times using different initial values. The results with
%%%lowest llh are output as the best parameters.

function [param,AIC]= fitRL_reversal(data,DiscLength,RecallLength)

global D paramcount QStor QTest PTest QAnimal Model DL RL sheet BetaStor ; D=data; DL=DiscLength; RL=RecallLength;
QTest=cell(10,1);
PTest=cell(10,1);
BetaStor=cell(10,1);
for iter=1:10
    %%%Initialize number of parameters to minimize over
    pmin=zeros(1,paramcount);
    pmax=ones(1,paramcount);
    %pmax = [];
    %for z = 1:paramcount
    %    pmax = [pmax 0.75];
    %end
    
   %%%Chooses a random initial value to run fmincon
   init=pmin+rand(1,length(pmin)).*(pmax-pmin);
   %init(init==0) = 0.00001;
   %%%Uses Matlab builtin optimization function fmincon to genereate a
   %%%vector of fit parameters. Localmin outputs the corresponding llh for
   %%%these parameters.
   QStor=iter;
   [bestparam,localmin]=fmincon(@llhRL,init,[],[],[],[],pmin,pmax);
   
   %%%Stores outputs for all 10 fmincon tests as a matrix
   results(iter,:)=[bestparam,localmin]
   
end
%%% Finds the set of parameters for which llh was lowest and outputs these
%%% outputs these parameters as best values.
[~,index]=min(results(:,end));
param= results(index(1),:);
Model.Qs{QAnimal}=QTest{index(1)};
Model.Ps{QAnimal}=PTest{index(1)};
Model.QsDisc{QAnimal}=Model.Qs{QAnimal}(1:DiscLength,:);
Model.PsDisc{QAnimal}=Model.Ps{QAnimal}(1:DiscLength,:);
Model.QsRec{QAnimal}=Model.Qs{QAnimal}(DiscLength+1:DiscLength+RecallLength,:);
Model.PsRec{QAnimal}=Model.Ps{QAnimal}(DiscLength+1:DiscLength+RecallLength,:);
Model.QsRev{QAnimal}=Model.Qs{QAnimal}(DiscLength+RecallLength+1:end,:);
Model.PsRev{QAnimal}=Model.Ps{QAnimal}(DiscLength+RecallLength+1:end,:);
AIC=2*paramcount+2*param(1,paramcount+1);
Model.AICc{QAnimal}=AIC +(2*(paramcount)^2+2*paramcount)/(length(Model.Animals)-paramcount-1);
Model.BIC{QAnimal}=2*log(length(Model.Animals))*paramcount+2*param(1,paramcount+1);
end


function llh=llhRL(par)
%%%Uses Model Dependent Likelihood Function to Determine llh to pass into
%%%FitRL.

global Model modelStrxr

switch Model.name
 
    case modelStrxr.RLTwoParam.name
        llh=llhRLTwoParam(par);
        
    case modelStrxr.RLSplitAlphaRev.name
        llh=llhRLSplitAlphaRev(par);
        
    case modelStrxr.RLSplitAlphaBetaRev.name
        llh=llhRLSplitAlphaBetaRev(par);
     
    case modelStrxr.RLSplitBetaRev.name
        llh=llhRLSplitBetaRev(par);
        
end
end

function llh= llhRLTwoParam(par)
global D DL RL Model QStor QTest PTest sheet;
beta= par(1);
alpha= par(2);

DiscLength = DL;
RecallLength = RL;

data=D;
if char(sheet) == 'D2_mcherry'
    Q0 = Model.initialvalue_D2_mcherry;
elseif char(sheet) == 'D2___hm3dq'
    Q0 = Model.initialvalue_D2___hm3dq;
elseif char(sheet) == 'D2___hm4di'
    Q0 = Model.initialvalue_D2___hm4di;
end


   
Q = Q0;
T = length(data);
llh = 0;

epsilon = .00001;


for t =1:T
    
    if t == DiscLength + RecallLength + 1
        Q(4) = (sum(Q)/100) * 0.6757;
    end
    
    softmax = exp(beta*Q)/sum(exp(beta*Q));
    esoftmax = (1-epsilon)*softmax + epsilon/4;
    
    a=data(t);
    llh = llh + log(esoftmax(a));
    
    
    if t <= DiscLength + RecallLength
        r = 100*(a==1);
    else
        r = 100*(a==2);
    end
    Q(a) = Q(a) + alpha*(r-Q(a));
    data(t)=a;
    Qs(t,:)=Q;
    Ps(t,:)=esoftmax;
end
QTest{QStor}=Qs;
PTest{QStor}=Ps;
%%%Function Output (negative so we can minimize instead of maximize)
llh=-llh;
end

function llh= llhRLSplitAlphaRev(par)
global D DL RL Model QStor QTest PTest sheet;
beta=par(1);
discalpha=par(2);
revalpha=par(3);

DiscLength = DL;
RecallLength = RL;

data=D;
if char(sheet) == 'D2_mcherry'
    Q0 = Model.initialvalue_D2_mcherry;
elseif char(sheet) == 'D2___hm3dq'
    Q0 = Model.initialvalue_D2___hm3dq;
elseif char(sheet) == 'D2___hm4di'
    Q0 = Model.initialvalue_D2___hm4di;
end

Q=Q0;
T=length(data);
llh=0;
epsilon=.00001;

for t=1:length(data)
    
    if t == DiscLength + RecallLength + 1
        Q(4) = (sum(Q)/100) * 0.6757;
    end
  
    softmax = exp(beta*Q)/sum(exp(beta*Q));
    esoftmax = (1-epsilon)*softmax + epsilon/4;
   
    a=data(t);
    llh = llh + log(esoftmax(a));
    
    
    if t <= DiscLength + RecallLength
        r = 100*(a==1); 
        Q(a) = Q(a) + discalpha*(r-Q(a));
    else
        r = 100*(a==2);
        Q(a) = Q(a) + revalpha*(r-Q(a));
    end
    data(t)=a;
    Qs(t,:)=Q;
    Ps(t,:)=esoftmax;  
end
QTest{QStor}=Qs;
PTest{QStor}=Ps;
llh=-llh;
end

function llh= llhRLSplitBetaRev(par)
global D DL RL Model QStor QTest PTest sheet;
discbeta=par(1);
alpha=par(2);
revbeta=par(3);

DiscLength = DL;
RecallLength = RL;

data=D;
if char(sheet) == 'D2_mcherry'
    Q0 = Model.initialvalue_D2_mcherry;
elseif char(sheet) == 'D2___hm3dq'
    Q0 = Model.initialvalue_D2___hm3dq;
elseif char(sheet) == 'D2___hm4di'
    Q0 = Model.initialvalue_D2___hm4di;
end


Q=Q0;
T=length(data);
llh=0;
epsilon=.00001;

for t=1:length(data)
    
    if t == DiscLength + RecallLength + 1
        Q(4) = (sum(Q)/100) * 0.6757;
    end
    
    if t <= DiscLength + RecallLength
        softmax = exp(discbeta*Q)/sum(exp(discbeta*Q));
        esoftmax = (1-epsilon)*softmax + epsilon/4;
    else
        softmax = exp(revbeta*Q)/sum(exp(revbeta*Q));
        esoftmax = (1-epsilon)*softmax + epsilon/4;
    end
    a=data(t);
    llh = llh + log(esoftmax(a));
    
    
    if t <= DiscLength + RecallLength
        r = 100*(a==1); 
    else
        r = 100*(a==2);
    end
    Q(a) = Q(a) + alpha*(r-Q(a));
    data(t)=a;
    Qs(t,:)=Q;
    Ps(t,:)=esoftmax;  
end
QTest{QStor}=Qs;
PTest{QStor}=Ps;
llh=-llh;
end

function llh= llhRLSplitAlphaBetaRev(par)
global D DL RL Model QStor QTest PTest sheet OSwitch ASwitch;
discbeta=par(1);
discalpha=par(2);
revbeta=par(3);
revalpha=par(4);

DiscLength = DL;
RecallLength = RL;

data=D;
if char(sheet) == 'D2_mcherry'
    Q0 = Model.initialvalue_D2_mcherry;
elseif char(sheet) == 'D2___hm3dq'
    Q0 = Model.initialvalue_D2___hm3dq;
elseif char(sheet) == 'D2___hm4di'
    Q0 = Model.initialvalue_D2___hm4di;
end

Q=Q0;
T=length(data);
llh=0;
epsilon=.00001;
ASwitch = [];
OSwitch = [];
rewards1 = [];
rewards2 = [];
rewards3 = [];
rewards4 = [];
for t=1:length(data)
    
    if t == DiscLength + RecallLength + 1
        Q(4) = (sum(Q)/100) * 0.6757;
    end
    
    if t <= DiscLength + RecallLength
        softmax = exp(discbeta*Q)/sum(exp(discbeta*Q));
    else
        softmax = exp(revbeta*Q)/sum(exp(revbeta*Q));
    end
    esoftmax = (1-epsilon)*softmax + epsilon/4;
   
    a=data(t);
    llh = llh + log(esoftmax(a));
    
    if t <= DiscLength + RecallLength
        r = 100*(a==1); 
        Q(a) = Q(a) + discalpha*(r-Q(a));
    else
        r = 100*(a==2);
        Q(a) = Q(a) + revalpha*(r-Q(a));
    end
    s = r/100;
    
    if t <= DiscLength + RecallLength
        if a == 1 
            rewards1 = [rewards1 s];
        elseif a == 2
            rewards2 = [rewards2 s];
        elseif a == 3
            rewards3 = [rewards3 s];
        elseif a == 4
            rewards4 = [rewards4 s];
        end
    end
    
    if t > DiscLength + RecallLength 
        if isempty(ASwitch) && a ~= 1
            ASwitch = t;
        end
        
        allrewards = [rewards1 rewards2 rewards3 rewards4];
        if mean(rewards1) <= mean(allrewards)
            OSwitch = t;
        end
    end
    
    data(t)=a;
    Qs(t,:)=Q;
    Ps(t,:)=esoftmax;  
end



QTest{QStor}=Qs;
PTest{QStor}=Ps;
llh=-llh;
end


