%%%Fits Parameters Using FitRL Function
function FitRL_Test
clear all
LoadModels_reversal;
global modelStrxr Model paramcount DiscLength RecallLength RevLength QAnimal AnimalList sheet

%%%User Input Here--Change to desired model to run. Check names of Models
%%%by looking in LoadModels function.
Model=modelStrxr.RLTwoParam;


%%%Finds number of parameters used by model. Output as global to fitRL.
paramcount=Model.paramcount;


%%% Names Excel File AllCohortTrialHistories that holds FourChoice Data for
%%% all cohorts
file='/Users/kdelevic/Dropbox/D2 DMS OVX project/eLife submission/transparent reporting info/modeling for github/Fig4/Fig4_TrialHistories.xlsx';

%%% Names Cohorts for fitting.
cohorts_list={'D2_mcherry','D2___hm3dq','D2___hm4di'}; %Odd OVX name to have same length name as other cohorts
AnimalIter=0;
%%%Run FitRL over the cohorts named above
for sheet = cohorts_list
    AnimalIter=AnimalIter+1;
    %%% Clear all relevant variables from previous cohort
    clear  param AIC QEndDisc QEndRecall QEndRev PEndDisc PEndRecall PEndRev TTC DiscTTC RecallTTC RevTTC TrialHistoryDisc TrialHistoryRecall TrialHistoryRev data DiscLength RecallLength 
    Model.Qs={};
    Model.Ps={};
    Model.QsDisc={};
    Model.PsDisc={};
    Model.QsRec={};
    Model.PsRec={};
    Model.QsRev={};
    Model.PsRev={};
    Model.data={};
    Model.DiscErrors=[];
    Model.RecallErrors=[];
    Model.RevErrors=[];
    Model.Disc4error=[];
    Model.Disc2error=[];
    Model.Disc3error=[];
    Model.Recall4error=[];
    Model.Recall2error=[];
    Model.Recall3error=[];
    Model.PerseverativeErrors=[];
    Model.RegressiveErrors=[];
    Model.NovelErrors=[];
    Model.IrrelevantErrors=[];
    
    Model.AICc={};
    Model.BIC={};
    
    Model.DiscCumulativeSumReward = [];
    Model.RecallCumulativeSumReward = [];
    Model.RevCumulativeSumReward = [];
    Model.DiscRewardPlacement = [];
    Model.RecallRewardPlacement = [];
    Model.RevRewardPlacement = [];
    
    Model.Animals=AnimalList{AnimalIter};
    %%% Pulls Cohort data from AllCohortTrialHistories excel sheet
    TrialHistoryDisc=xlsread(file,[char(sheet),'_Disc']);
    TrialHistoryRecall=xlsread(file,[char(sheet),'_Recall']);
    TrialHistoryRev=xlsread(file,[char(sheet),'_Rev']);
    
    %%% Gives dimensions of matrix containing cohort data. 
    SizeDisc = size(TrialHistoryDisc);
    SizeRev = size(TrialHistoryRev);
    if SizeDisc(1) > SizeRev(1)
        datasize=SizeDisc;
    else
        datasize=[SizeRev(1),SizeDisc(2)];
    end

%%% Runs fitRL for every animal in the cohort
for i=1:datasize(2)

    
    
    %%%Concatenates Discrimination and Recall Trial Histories and Stores
    %%%Raw Data
    data=[TrialHistoryDisc(:,i);TrialHistoryRecall(:,i);TrialHistoryRev(:,i)]; 
    
    
    %%%Removes All NaNs and omissions from trial history 
    
    data(isnan(data))= [];
    Model.rawdata{i}=data;
    data(data==0) = [];
    
    
    %%% Determine DiscLength by clearing NaNs and omissions from
    %%% discrimination only. DiscLength used to switch from discrimination
    %%% to recall in fitRL and to store TTCs.
    DiscriminationTrials=TrialHistoryDisc(:,i);
    DiscriminationTrials(isnan(DiscriminationTrials))=[];
    RawDiscLength=length(DiscriminationTrials);
    Model.DiscErrors{i}=sum(DiscriminationTrials>1);
    Model.Disc4error{i}=sum(DiscriminationTrials==4);
    Model.Disc2error{i}=sum(DiscriminationTrials==2);
    Model.Disc3error{i}=sum(DiscriminationTrials==3);
    DiscriminationTrials(DiscriminationTrials==0)=[];
    RecallTrials=TrialHistoryRecall(:,i);
    RecallTrials(isnan(RecallTrials))=[];
    RecallTrials(RecallTrials==0)=[];
    Model.RecallErrors{i}=sum(RecallTrials>1);
    Model.Recall4error{i}=sum(RecallTrials==4);
    Model.Recall2error{i}=sum(RecallTrials==2);
    Model.Recall3error{i}=sum(RecallTrials==3);
    RecallLength=length(RecallTrials);
    DiscLength=length(DiscriminationTrials);
    
    ReversalTrials=TrialHistoryRev(:,i);
    ReversalTrials(isnan(ReversalTrials))=[];
    RawRevLength=length(ReversalTrials);
    Model.RevErrors{i}=sum(ReversalTrials==1) + sum(ReversalTrials>2);
    n = 0;
    RevT = ReversalTrials.';
    for k = RevT
        n = n + 1;
        if k == 2
            break
        end
    end
    RegressiveTrials=ReversalTrials(n:end);
    PerseverativeTrials=ReversalTrials(1:n-1);
    Model.PerseverativeErrors{i}=sum(PerseverativeTrials==1);
    Model.RegressiveErrors{i}=sum(RegressiveTrials==1);
    Model.NovelErrors{i}=sum(ReversalTrials==4);
    Model.IrrelevantErrors{i}=sum(ReversalTrials==3);
    ReversalTrials(ReversalTrials==0)=[];
    RevLength=length(ReversalTrials);
    
    %%%follow accumulation of reward across behavior
    rewardsum = 0;
    discrewardsum = [];
    discrewardplacement = [];
    for k = 1:DiscLength
        if DiscriminationTrials(k) == 1
            rewardsum = rewardsum + 1;
            discrewardplacement = [discrewardplacement (k/DiscLength)*100];
        end 
        discrewardsum(k) = rewardsum;
    end
    rewardsum = 0;
    recallrewardsum = [];
    recallrewardplacement = [];
    for k = 1:RecallLength
        if RecallTrials(k) == 1
            rewardsum = rewardsum + 1;
            recallrewardplacement = [recallrewardplacement (k/RecallLength)*100];
        end 
        recallrewardsum(k) = rewardsum;
    end
    rewardsum = 0;
    revrewardsum = [];
    revrewardplacement = [];
    for k = 1:RevLength
        if ReversalTrials(k) == 2
            rewardsum = rewardsum + 1;
            revrewardplacement = [revrewardplacement (k/RevLength)*100];
        end 
        revrewardsum(k) = rewardsum;
    end
    discrewardsum = discrewardsum.';
    discrewardplacement = discrewardplacement.';
    recallrewardsum = recallrewardsum.';
    recallrewardplacement = recallrewardplacement.';
    revrewardsum = revrewardsum.';
    revrewardplacement = revrewardplacement.';
    
    Model.DiscCumulativeSumReward{i} = discrewardsum;
    Model.DiscRewardPlacement{i} = discrewardplacement;
    Model.RecallCumulativeSumReward{i} = recallrewardsum;
    Model.RecallRewardPlacement{i} = recallrewardplacement;
    Model.RevCumulativeSumReward{i} = revrewardsum;
    Model.RevRewardPlacement{i} = revrewardplacement;
    
    %%%QAnimal Set to Store Q History for each animal
    QAnimal=i;
    
    %%%boolean for fitRL to know whether this is initial fit or refitting
    %%%with priors
    refit_bool = 0;
    %%% Run fitRL and store parameters and AIC for each mouse
    [param(i,:),AIC(i)]=fitRL_reversal(data,DiscLength,RecallLength);
   
    
    %%% Stores TTCs and Ending Discrimination Odor Values
    Model.data{i}=data;

    TTC(i)=length(data);
    DiscTTC(i)=DiscLength;
    RecallTTC(i)=TTC(i)-DiscLength-RevLength;
    DiscRecTTC(i)=DiscTTC(i)+RecallTTC(i);
    RevTTC(i)=RevLength;
    QEndDisc(i,:)=Model.Qs{QAnimal}(DiscLength,:);
    PEndDisc(i,:)=Model.Ps{QAnimal}(DiscLength,:);
    QEndRecall(i,:)=Model.Qs{QAnimal}(DiscLength+RecallLength,:);
    PEndRecall(i,:)=Model.Ps{QAnimal}(DiscLength+RecallLength,:);
    QEndRev(i,:)=Model.Qs{QAnimal}(end,:);
    PEndRev(i,:)=Model.Ps{QAnimal}(end,:);
    
     
    %%%Replace values with NaN if mouse disqualified in any previous days.
    %%%This only works if data is input as a single 0 for any day following
    %%%disqualification.
    if sum(DiscriminationTrials==1) == 0
        Model.DiscErrors{i}=nan;
        Model.Disc4error{i}=nan;
        Model.Disc2error{i}=nan;
        Model.Disc3error{i}=nan;
        
        
        DiscTTC(i)=nan;
        QEndDisc(i,:)=nan;
        PEndDisc(i,:)=nan;
    end
    
    if sum(RecallTrials==1) == 0
        Model.RecallErrors{i}=nan;
        Model.Recall4error{i}=nan;
        Model.Recall2error{i}=nan;
        Model.Recall3error{i}=nan;
        RecallTTC(i)=nan;
    end
    
    if sum(ReversalTrials==2)==0
        Model.RevErrors{i}=nan;
        Model.PerseverativeErrors{i}=nan;
        Model.RegressiveErrors{i}=nan;
        Model.NovelErrors{i}=nan;
        Model.IrrelevantErrors{i}=nan;
        RevTTC(i)=nan;
    end

end


%%% Stores individual parameters into cell structure. Parameters are
%%% model dependent, so helper functions used to store parameters.
switch Model.name
    
        
    case modelStrxr.RLTwoParam.name
        FitRL_RLTwoParam(param);
        
    case modelStrxr.RLSplitAlphaRev.name
        FitRL_RLSplitAlphaRev(param);
 
    case modelStrxr.RLSplitBetaRev.name
        FitRL_RLSplitBetaRev(param);

    case modelStrxr.RLSplitAlphaBetaRev.name
        FitRL_RLSplitAlphaBetaRev(param);
    
    
end
%%% Stores model-independent variables into cell structure.
Model.TTC=TTC;
Model.DiscTTC=DiscTTC;
Model.RecallTTC=RecallTTC;
Model.DiscRecTTC=DiscRecTTC;
Model.RevTTC=RevTTC;
Model.llh=param(:,end);
Model.AIC=AIC;
Model.QEndDisc=QEndDisc;
Model.PEndDisc=PEndDisc;
Model.QEndRecall=QEndRecall;
Model.PEndRecall=PEndRecall;
Model.QEndRev=QEndRev;
Model.PEndRev=PEndRev;

%%% Saves to .mat files
parameters=Model;
save([sprintf('%s_modelingparameters061421_%s.mat',Model.name,char(sheet))],'parameters');
end
end


function FitRL_RLTwoParam(param)
global Model 

Model.beta=param(:,1);
Model.alpha=param(:,2);
end


function FitRL_RLSplitAlphaRev(param)
global Model 

Model.beta=param(:,1);
Model.discalpha=param(:,2);
Model.revalpha=param(:,3);
end


function FitRL_RLSplitBetaRev(param)
global Model

Model.discbeta=param(:,1);
Model.alpha=param(:,2);
Model.revbeta=param(:,3);

end


function FitRL_RLSplitAlphaBetaRev(param)
global Model 

Model.discbeta=param(:,1);
Model.discalpha=param(:,2);
Model.revbeta=param(:,3);
Model.revalpha=param(:,4);
end

