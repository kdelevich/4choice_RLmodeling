%%% Runs generate_and_recover to simulate 100 trials per mouse. Simulated
%%% Data is re-fit using llh and fitRL subfuncitons in
%%% generate_and_recover. All simulation information is stored in .mat
%%% files (one per cohort).

%%% Begin Function
function Four_Choice_Simulator

%%% Clear Workspace and Load Models
clear all
LoadModels_reversal;
global Model modelStrxr sheet 


%%% User input. Load one of the .mat files corresponding to the model you 
%%% want to run. It does not matter which cohort you choose for the .mat
%%% file; all cohorts will be simulated.
Model=load('RLTwoParam_modelingparameters061421_D2_mcherry.mat');


%%% Decides which Four_Choice to run based on model name.
switch Model.parameters.name
    
    case modelStrxr.RLTwoParam.name
        Four_Choice_SimulatorRLTwoParam;
        
    case modelStrxr.RLSplitAlphaRev.name
        Four_Choice_SimulatorRLSplitAlphaRev;
        
    case modelStrxr.RLSplitBetaRev.name
        Four_Choice_SimulatorRLSplitBetaRev;
        
    case modelStrxr.RLSplitAlphaBetaRev.name
        Four_Choice_SimulatorRLSplitAlphaBetaRev;
end

    
end


function Four_Choice_SimulatorRLTwoParam

global Model alpha beta Q0 sheet QAnimal

cohorts_list={'D2_mcherry','D2___hm3dq','D2___hm4di'};
progressbar('Cohort Progress','Animal Progress','Simulation Progress')
for sheet = cohorts_list
    Model=load(sprintf('RLTwoParam_modelingparameters061421_%s.mat',char(sheet)));
    datasize=size(Model.parameters.Animals);
    
for j=1:datasize(2)
    progressbar([],[],[])
    beta=Model.parameters.beta(j);
    alpha=Model.parameters.alpha(j);
    QAnimal=j;
    if char(sheet) == 'D2_mcherry'
        Q0=Model.parameters.initialvalue_D2_mcherry;
        u=1;
    elseif char(sheet) == 'D2___hm3dq'
        Q0=Model.parameters.initialvalue_D2___hm3dq;
        u=2;
    elseif char(sheet) == 'D2___hm4di'
        Q0=Model.parameters.initialvalue_D2___hm4di;
    end
    for i=1:100
        
        [datastrx{i},paramstrx{i},AICstrx(i),dTTC(i),recTTC(i),revTTC(i)]=generate_and_recover;
        simTTC(i)=(length(datastrx{i}));
        progressbar([],[],i/100)
    end
    
   
    Model.parameters.recovereddata{j}=datastrx;
    Model.parameters.simalpha(j)= mean(cellfun(@(x) (x(2)), paramstrx));
    Model.parameters.simbeta(j)=  mean(cellfun(@(x) (x(1)), paramstrx));
    Model.parameters.simTTC(:,j)=simTTC;
    Model.parameters.simdTTC(:,j)=dTTC;
    Model.parameters.simrecTTC(:,j)=recTTC;
    Model.parameters.simrevTTC(:,j)=revTTC;
    Model.parameters.simTTCavg(j)=mean(simTTC);
    Model.parameters.simDiscTTCavg(j)=mean(dTTC);
    Model.parameters.simRecallTTCavg(j)=mean(recTTC);
    Model.parameters.simDiscRecTTCavg(j)=mean(dTTC)+mean(recTTC);
    Model.parameters.simRevTTCavg(j)=mean(revTTC);
   
    progressbar([],((j-.00001)/datasize(2)),[])
end
parameters=Model.parameters;
    save([sprintf('%s_fit_and_recover_parameters_%s.mat',Model.parameters.name,char(sheet))],'parameters');
    progressbar(.33*u,0,0)
end
progressbar(1,1,1)
end

function Four_Choice_SimulatorRLSplitAlphaBetaRev

global Model discalpha discbeta revalpha revbeta Q0 sheet QAnimal

cohorts_list={'D2_mcherry','D2___hm3dq','D2___hm4di'};
progressbar('Cohort Progress','Animal Progress','Simulation Progress')
for sheet = cohorts_list
    Model=load(sprintf('RLSplitAlphaBetaRev_modelingparameters061421_%s.mat',char(sheet)));
    datasize=size(Model.parameters.Animals);
for j=1:datasize(2)
    progressbar([],[],[])
    discbeta=refit_Model.parameters.discbeta(j);
    discalpha=refit_Model.parameters.discalpha(j);
    revalpha=refit_Model.parameters.revalpha(j);
    revbeta=refit_Model.parameters.revbeta(j);
    QAnimal=j;
    if char(sheet) == 'D2_mcherry'
        Q0=Model.parameters.initialvalue_D2_mcherry;
        u=1;
    elseif char(sheet) == 'D2___hm3dq'
        Q0=Model.parameters.initialvalue_D2___hm3dq;
        u=2;
    elseif char(sheet) == 'D2___hm4di'
        Q0=Model.parameters.initialvalue_D2___hm4di;
    end
    for i=1:100
        
        [datastrx{i},paramstrx{i},AICstrx(i),dTTC(i),recTTC(i),revTTC(i)]=generate_and_recover;
        simTTC(i)=(length(datastrx{i}));
        progressbar([],[],i/100)
    end
    
   
    Model.parameters.recovereddata{j}=datastrx;
    Model.parameters.simdiscalpha(j)= mean(cellfun(@(x) (x(2)), paramstrx));
    Model.parameters.simdiscbeta(j)=  mean(cellfun(@(x) (x(1)), paramstrx));
    Model.parameters.simrevbeta(j)=  mean(cellfun(@(x) (x(3)), paramstrx));
    Model.parameters.simrevalpha(j)=  mean(cellfun(@(x) (x(4)), paramstrx));
    Model.parameters.simTTC(:,j)=simTTC;
    Model.parameters.simdTTC(:,j)=dTTC;
    Model.parameters.simrecTTC(:,j)=recTTC;
    Model.parameters.simrevTTC(:,j)=revTTC;
    Model.parameters.simTTCavg(j)=mean(simTTC);
    Model.parameters.simDiscTTCavg(j)=mean(dTTC);
    Model.parameters.simRecallTTCavg(j)=mean(recTTC);
    Model.parameters.simDiscRecTTCavg(j)=mean(dTTC)+mean(recTTC);
    Model.parameters.simRevTTCavg(j)=mean(revTTC);
   
    progressbar([],((j-.00001)/datasize(2)),[])
end
parameters=Model.parameters;
    save([sprintf('%s_fit_and_recover_parameters_%s.mat',Model.parameters.name,char(sheet))],'parameters');
    progressbar(.33*u,0,0)
end
progressbar(1,1,1)
end

function Four_Choice_SimulatorRLSplitAlphaRev

global Model discalpha beta revalpha Q0 sheet QAnimal

cohorts_list={'D2_mcherry','D2___hm3dq','D2___hm4di'};
progressbar('Cohort Progress','Animal Progress','Simulation Progress')
for sheet = cohorts_list
    Model=load(sprintf('RLSplitAlphaRev_modelingparameters061421_%s.mat',char(sheet)));
    datasize=size(Model.parameters.Animals);
    
for j=1:datasize(2)
    progressbar([],[],[])
    beta=Model.parameters.beta(j);
    discalpha=Model.parameters.discalpha(j);
    revalpha=Model.parameters.revalpha(j);
    QAnimal=j;
    if char(sheet) == 'D2_mcherry'
        Q0=Model.parameters.initialvalue_D2_mcherry;
        u=1;
    elseif char(sheet) == 'D2___hm3dq'
        Q0=Model.parameters.initialvalue_D2___hm3dq;
        u=2;
    elseif char(sheet) == 'D2___hm4di'
        Q0=Model.parameters.initialvalue_D2___hm4di;
    end
    for i=1:100
        
        [datastrx{i},paramstrx{i},AICstrx(i),dTTC(i),recTTC(i),revTTC(i)]=generate_and_recover;
        simTTC(i)=(length(datastrx{i}));
        progressbar([],[],i/100)
    end
    
   
    Model.parameters.recovereddata{j}=datastrx;
    Model.parameters.simdiscalpha(j)= mean(cellfun(@(x) (x(2)), paramstrx));
    Model.parameters.simbeta(j)=  mean(cellfun(@(x) (x(1)), paramstrx));
    Model.parameters.simrevalpha(j)=  mean(cellfun(@(x) (x(3)), paramstrx));
    Model.parameters.simTTC(:,j)=simTTC;
    Model.parameters.simdTTC(:,j)=dTTC;
    Model.parameters.simrecTTC(:,j)=recTTC;
    Model.parameters.simrevTTC(:,j)=revTTC;
    Model.parameters.simTTCavg(j)=mean(simTTC);
    Model.parameters.simDiscTTCavg(j)=mean(dTTC);
    Model.parameters.simRecallTTCavg(j)=mean(recTTC);
    Model.parameters.simDiscRecTTCavg(j)=mean(dTTC)+mean(recTTC);
    Model.parameters.simRevTTCavg(j)=mean(revTTC);
   
    progressbar([],((j-.00001)/datasize(2)),[])
end
parameters=Model.parameters;
    save([sprintf('%s_fit_and_recover_parameters_%s.mat',Model.parameters.name,char(sheet))],'parameters');
    progressbar(.33*u,0,0)
end
progressbar(1,1,1)
end

function Four_Choice_SimulatorRLSplitBetaRev

global Model alpha revbeta discbeta Q0 sheet QAnimal

cohorts_list={'D2_mcherry','D2___hm3dq','D2___hm4di'};
progressbar('Cohort Progress','Animal Progress','Simulation Progress')
for sheet = cohorts_list
    Model=load(sprintf('RLSplitBetaRev_modelingparameters061321_%s.mat',char(sheet)));
    datasize=size(Model.parameters.Animals);
    
for j=1:datasize(2)
    progressbar([],[],[])
    alpha=Model.parameters.alpha(j);
    discbeta=Model.parameters.discbeta(j);
    revbeta=Model.parameters.revbeta(j);
    QAnimal=j;
    if char(sheet) == 'D2_mcherry'
        Q0=Model.parameters.initialvalue_D2_mcherry;
        u=1;
    elseif char(sheet) == 'D2___hm3dq'
        Q0=Model.parameters.initialvalue_D2___hm3dq;
        u=2;
    elseif char(sheet) == 'D2___hm4di'
        Q0=Model.parameters.initialvalue_D2___hm4di;
    end
    for i=1:100
        
        [datastrx{i},paramstrx{i},AICstrx(i),dTTC(i),recTTC(i),revTTC(i)]=generate_and_recover;
        simTTC(i)=(length(datastrx{i}));
        progressbar([],[],i/100)
    end
    
   
    Model.parameters.recovereddata{j}=datastrx;
    Model.parameters.simalpha(j)= mean(cellfun(@(x) (x(2)), paramstrx));
    Model.parameters.simdiscbeta(j)=  mean(cellfun(@(x) (x(1)), paramstrx));
    Model.parameters.simrevbeta(j)=  mean(cellfun(@(x) (x(3)), paramstrx));
    Model.parameters.simTTC(:,j)=simTTC;
    Model.parameters.simdTTC(:,j)=dTTC;
    Model.parameters.simrecTTC(:,j)=recTTC;
    Model.parameters.simrevTTC(:,j)=revTTC;
    Model.parameters.simTTCavg(j)=mean(simTTC);
    Model.parameters.simDiscTTCavg(j)=mean(dTTC);
    Model.parameters.simRecallTTCavg(j)=mean(recTTC);
    Model.parameters.simDiscRecTTCavg(j)=mean(dTTC)+mean(recTTC);
    Model.parameters.simRevTTCavg(j)=mean(revTTC);
   
    progressbar([],((j-.00001)/datasize(2)),[])
end
parameters=Model.parameters;
    save([sprintf('%s_fit_and_recover_parameters_%s.mat',Model.parameters.name,char(sheet))],'parameters');
    progressbar(.33*u,0,0)
end
progressbar(1,1,1)
end

