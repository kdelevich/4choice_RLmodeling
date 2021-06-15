function [data,param,AIC,dTTC,recTTC,revTTC]= generate_and_recover
global   modelStrxr Model D sheet QAnimal

%%% Selects Simulator to use. Model passed in by Four_Choice_Simulator. See
%%% generate_and_recoverTwoParamDisc and llhTwoParamDisc for detailed
%%% comments.

switch Model.parameters.name
    
    case modelStrxr.RLTwoParam.name
        [data,param,AIC,dTTC,recTTC,revTTC]=generate_and_recoverRLTwoParam;
        
    case modelStrxr.RLSplitBetaRev.name
        [data,param,AIC,dTTC,recTTC,revTTC]=generate_and_recoverRLSplitBetaRev;
        
    case modelStrxr.RLSplitAlphaRev.name
        [data,param,AIC,dTTC,recTTC,revTTC]=generate_and_recoverRLSplitAlphaRev;
        
    case modelStrxr.RLSplitAlphaBetaRev.name
        [data,param,AIC,dTTC,recTTC,revTTC]=generate_and_recoverRLSplitAlphaBetaRev;
        
end

end


function [data,param,AIC,dTTC,recTTC,revTTC]= generate_and_recoverRLTwoParam
global Q0 alpha beta Model D Discdata paramcount sheet QAnimal


paramcount=2;
T=300;
data=[];
Discdata=[];
recalldata=[];
revdata=[];
Q = Q0;
for t =1:T
    
    softmax = exp(beta*Q);
    softmax=softmax/sum(softmax);
    cdf = [0 cumsum(softmax)];
    
    a= find(cdf<rand); a=a(end);
    
    r = 100*(a==1);                    
    Q(a) = Q(a) + alpha*(r-Q(a));
    Discdata(t)=a;
    %%%This is the Checking Condition to see if TTC has been met%%%
    if t>=10                                   
        checkmatrix=Discdata((t-9):t);               
        checkmatrixcounts=sum(checkmatrix==1);
        if checkmatrixcounts>=8             
            break                                
        end
    end  
end
for t =1:T
    
    softmax = exp(beta*Q);
    softmax=softmax/sum(softmax);
    cdf = [0 cumsum(softmax)];
    
    a= find(cdf<rand); a=a(end);
    
   
    r = 100*(a==1);                    
    Q(a) = Q(a) + alpha*(r-Q(a));
    recalldata(t)=a;
    %%%This is the Checking Condition to see if TTC has been met%%%
    if t>=10                                   
        checkmatrix=recalldata((t-9):t);               
        checkmatrixcounts=sum(checkmatrix==1);
        if checkmatrixcounts>=8             
            break                                
        end
    end  
end
Q(4) = (sum(Q)/100) * 0.6757;
for t =1:T
    
    softmax = exp(beta*Q);
    softmax=softmax/sum(softmax);
    cdf = [0 cumsum(softmax)];
    
    a= find(cdf<rand); a=a(end);
    
 
    r = 100*(a==2);                    
    Q(a) = Q(a) + alpha*(r-Q(a));
    revdata(t)=a;
    %%%This is the Checking Condition to see if TTC has been met%%%
    if t>=10                                   
        checkmatrix=revdata((t-9):t);               
        checkmatrixcounts=sum(checkmatrix==2);
        if checkmatrixcounts>=8             
            break                                
        end
    end  
end

data=[Discdata, recalldata, revdata];
%%%Function Output%%%
dTTC=length(Discdata);
recTTC=length(recalldata);
revTTC=length(revdata);
[param,AIC] = fitRL_recover(data,dTTC,recTTC);
end

function [data,param,AIC,dTTC,recTTC,revTTC]= generate_and_recoverRLSplitAlphaRev
global Q0 discalpha beta revalpha Model D Discdata paramcount sheet QAnimal


paramcount=3;
T=300;
data=[];
Discdata=[];
recalldata=[];
revdata=[];
Q = Q0;
for t =1:T
    
    softmax = exp(beta*Q);
    softmax=softmax/sum(softmax);
    cdf = [0 cumsum(softmax)];
    
    a= find(cdf<rand); a=a(end);
    
    r = 100*(a==1);                    
    Q(a) = Q(a) + discalpha*(r-Q(a));
    Discdata(t)=a;
    %%%This is the Checking Condition to see if TTC has been met%%%
    if t>=10                                   
        checkmatrix=Discdata((t-9):t);               
        checkmatrixcounts=sum(checkmatrix==1);
        if checkmatrixcounts>=8             
            break                                
        end
    end  
end
for t =1:T
    
    softmax = exp(beta*Q);
    softmax=softmax/sum(softmax);
    cdf = [0 cumsum(softmax)];
    
    a= find(cdf<rand); a=a(end);
    
   
    
    r = 100*(a==1);                    
    Q(a) = Q(a) + discalpha*(r-Q(a));
    recalldata(t)=a;
    %%%This is the Checking Condition to see if TTC has been met%%%
    if t>=10                                   
        checkmatrix=recalldata((t-9):t);               
        checkmatrixcounts=sum(checkmatrix==1);
        if checkmatrixcounts>=8             
            break                                
        end
    end  
end
Q(4) = (sum(Q)/100) * 0.6757;
for t =1:T
    
    softmax = exp(beta*Q);
    softmax=softmax/sum(softmax);
    cdf = [0 cumsum(softmax)];
    
    a= find(cdf<rand); a=a(end);
    
    r = 100*(a==2);                    
    Q(a) = Q(a) + revalpha*(r-Q(a));
    revdata(t)=a;
    %%%This is the Checking Condition to see if TTC has been met%%%
    if t>=10                                   
        checkmatrix=revdata((t-9):t);               
        checkmatrixcounts=sum(checkmatrix==2);
        if checkmatrixcounts>=8             
            break                                
        end
    end  
end

data=[Discdata, recalldata, revdata];
%%%Function Output%%%
dTTC=length(Discdata);
recTTC=length(recalldata);
revTTC=length(revdata);
[param,AIC] = fitRL_recover(data,dTTC,recTTC);
end

function [data,param,AIC,dTTC,recTTC,revTTC]= generate_and_recoverRLSplitBetaRev
global Q0 alpha discbeta revbeta Model D Discdata paramcount sheet QAnimal


paramcount=3;
T=300;
data=[];
Discdata=[];
recalldata=[];
revdata=[];
Q = Q0;
for t =1:T
    
    softmax = exp(discbeta*Q);
    softmax=softmax/sum(softmax);
    cdf = [0 cumsum(softmax)];
    
    a= find(cdf<rand); a=a(end);
    
   
    
    r = 100*(a==1);                    
    Q(a) = Q(a) + alpha*(r-Q(a));
    Discdata(t)=a;
    %%%This is the Checking Condition to see if TTC has been met%%%
    if t>=10                                   
        checkmatrix=Discdata((t-9):t);               
        checkmatrixcounts=sum(checkmatrix==1);
        if checkmatrixcounts>=8             
            break                                
        end
    end  
end
for t =1:T
    
    softmax = exp(discbeta*Q);
    softmax=softmax/sum(softmax);
    cdf = [0 cumsum(softmax)];
    
    a= find(cdf<rand); a=a(end);
    
   
    
    r = 100*(a==1);                    
    Q(a) = Q(a) + alpha*(r-Q(a));
    recalldata(t)=a;
    %%%This is the Checking Condition to see if TTC has been met%%%
    if t>=10                                   
        checkmatrix=recalldata((t-9):t);               
        checkmatrixcounts=sum(checkmatrix==1);
        if checkmatrixcounts>=8             
            break                                
        end
    end  
end
for t =1:T
    
    softmax = exp(revbeta*Q);
    softmax=softmax/sum(softmax);
    cdf = [0 cumsum(softmax)];
    
    a= find(cdf<rand); a=a(end);
    
   
    
    r = 100*(a==2);                    
    Q(a) = Q(a) + alpha*(r-Q(a));
    revdata(t)=a;
    %%%This is the Checking Condition to see if TTC has been met%%%
    if t>=10                                   
        checkmatrix=revdata((t-9):t);               
        checkmatrixcounts=sum(checkmatrix==2);
        if checkmatrixcounts>=8             
            break                                
        end
    end  
end

data=[Discdata, recalldata, revdata];
%%%Function Output%%%
dTTC=length(Discdata);
recTTC=length(recalldata);
revTTC=length(revdata);
[param,AIC] = fitRL_recover(data,dTTC,recTTC);
end

function [data,param,AIC,dTTC,recTTC,revTTC]= generate_and_recoverRLSplitAlphaBetaRev
global Q0 discalpha discbeta revbeta revalpha Model D Discdata paramcount sheet QAnimal


paramcount=4;
T=300;
data=[];
Discdata=[];
recalldata=[];
revdata=[];
Q = Q0;
for t =1:T
    
    softmax = exp(discbeta*Q);
    softmax=softmax/sum(softmax);
    cdf = [0 cumsum(softmax)];
    
    a= find(cdf<rand); a=a(end); 
    
    r = 100*(a==1);                    
    Q(a) = Q(a) + discalpha*(r-Q(a));
    Discdata(t)=a;
    %%%This is the Checking Condition to see if TTC has been met%%%
    if t>=10                                   
        checkmatrix=Discdata((t-9):t);               
        checkmatrixcounts=sum(checkmatrix==1);
        if checkmatrixcounts>=8             
            break                                
        end
    end  
end
for t =1:T
    
    softmax = exp(discbeta*Q);
    softmax=softmax/sum(softmax);
    cdf = [0 cumsum(softmax)];
    
    a= find(cdf<rand); a=a(end);
    
    r = 100*(a==1);                    
    Q(a) = Q(a) + discalpha*(r-Q(a));
    recalldata(t)=a;
    %%%This is the Checking Condition to see if TTC has been met%%%
    if t>=10                                   
        checkmatrix=recalldata((t-9):t);               
        checkmatrixcounts=sum(checkmatrix==1);
        if checkmatrixcounts>=8             
            break                                
        end
    end  
end
Q(4) = (sum(Q)/100) * 0.6757;
for t =1:T
    
    softmax = exp(revbeta*Q);
    softmax=softmax/sum(softmax);
    cdf = [0 cumsum(softmax)];
    
    a= find(cdf<rand); a=a(end);
    
   
    
    r = 100*(a==2);                    
    Q(a) = Q(a) + revalpha*(r-Q(a));
    revdata(t)=a;
    %%%This is the Checking Condition to see if TTC has been met%%%
    if t>=10                                   
        checkmatrix=revdata((t-9):t);               
        checkmatrixcounts=sum(checkmatrix==2);
        if checkmatrixcounts>=8             
            break                                
        end
    end  
end

data=[Discdata, recalldata, revdata];
%%%Function Output%%%
dTTC=length(Discdata);
recTTC=length(recalldata);
revTTC=length(revdata);
[param,AIC] = fitRL_recover(data,dTTC,recTTC);
end

%%% Minimizes negative llh 10 times using different initial values. The results with
%%% lowest llh are output as the best parameters.
function [param,AIC]= fitRL(data)
global paramcount D 

D=data;

for iter=1:10
    %%%Initialize number of parameters to minimize over
    pmin=zeros(1,paramcount);
    pmax=ones(1,paramcount);

   %%%Chooses a random initial value to run fmincon
   init=pmin+rand(1,length(pmin)).*(pmax-pmin);
   
   %%%Uses Matlab builtin optimization function fmincon to genereate a
   %%%vector of fit parameters. Localmin outputs the corresponding llh for
   %%%these parameters.
   [bestparam,localmin]=fmincon(@llhRL,init,[],[],[],[],pmin,pmax);
   
   %%%Stores outputs for all 10 fmincon tests as a matrix
   results(iter,:)=[bestparam,localmin]
end
%%% Finds the set of parameters for which llh was lowest and outputs these
%%% outputs these parameters as best values.
[~,index]=min(results(:,end));
param= results(index(1),:);

%%% Calculates AIC Score. Note that AIC penalizes based on number of
%%% parameters, so paramcount must be taken into account.
AIC=2*paramcount+2*param(1,paramcount+1);
end

function llh=llhRL(par)
%%% Uses Model Dependent Likelihood Function to Determine llh to pass into
%%% FitRL.

global Model modelStrxr

switch Model.parameters.name
    
    case modelStrxr.RLTwoParam.name
        llh=llhRLTwoParam(par);
        
    case modelStrxr.RLSplitBeta.name
        llh=llhRLSplitBeta(par);
  
    case modelStrxr.RLSplitAlphaForget.name
        llh=llhRLSplitAlphaForget(par);

    case modelStrxr.RLSplitImpulsive.name
        llh=llhRLSplitImpulsive(par);
  
    case modelStrxr.RLForget.name
        llh=llhRLForget(par);

    case modelStrxr.RLSplitBetaForget.name
        llh=llhRLSplitBetaForget(par);
      
end
end
