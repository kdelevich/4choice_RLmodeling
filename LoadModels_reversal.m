%%%Places Models into structure array for access in generate_and_recover
function LoadModels_reversal
global modelStrxr AnimalList

%%%AnimalList%%%
AnimalList{1}={'M1'	'M2' 'M5' 'M6' 'M7'	'M10' 'M11' 'M13' 'M16' 'M17' 'M20' 'M21' 'M24' 'M25' 'M26'};
AnimalList{2}={'M3'	'M4' 'M8' 'M9' 'M12' 'M14' 'M15' 'M18' 'M19' 'M22' 'M23' 'M27' 'M28' 'M29' 'M30' 'M31'};

%%% FINAL USED MODELS
%%% initial odor values [O1 O2 O3 O4]

modelStrxr.RLTwoParam.name='RLTwoParam';
modelStrxr.RLTwoParam.initialvalue_F_Sham=[18.55 11.29 2.42 67.74];
modelStrxr.RLTwoParam.initialvalue_F__OVX=[18.55 11.29 2.42 67.74];
modelStrxr.RLTwoParam.paramcount=2;

modelStrxr.RLSplitBetaRev.name='RLSplitBetaRev';
modelStrxr.RLSplitBetaRev.initialvalue_F_Sham=[18.55 11.29 2.42 67.74];
modelStrxr.RLSplitBetaRev.initialvalue_F__OVX=[18.55 11.29 2.42 67.74];
modelStrxr.RLSplitBetaRev.paramcount=3;

modelStrxr.RLSplitAlphaRev.name='RLSplitAlphaRev';
modelStrxr.RLSplitAlphaRev.initialvalue_F_Sham=[18.55 11.29 2.42 67.74];
modelStrxr.RLSplitAlphaRev.initialvalue_F__OVX=[18.55 11.29 2.42 67.74];
modelStrxr.RLSplitAlphaRev.paramcount=3;

modelStrxr.RLSplitAlphaBetaRev.name='RLSplitAlphaBetaRev';
modelStrxr.RLSplitAlphaBetaRev.initialvalue_F_Sham=[18.55 11.29 2.42 67.74];
modelStrxr.RLSplitAlphaBetaRev.initialvalue_F__OVX=[18.55 11.29 2.42 67.74];
modelStrxr.RLSplitAlphaBetaRev.paramcount=4;

%%%END Final Models

end