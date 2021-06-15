%%%Places Models into structure array for access in generate_and_recover
function LoadModels_reversal
global modelStrxr AnimalList

%%%AnimalList%%%
AnimalList{1}={'M1'	'M2' 'M3' 'M4' 'M5' 'M6' 'M7' 'M8'  'M9'};
AnimalList{2}={'M15' 'M16' 'M17' 'M18' 'M19' 'M20'};
AnimalList{3}={'M10' 'M11' 'M12' 'M13' 'M14'}

%%% FINAL USED MODELS
%%% initial odor values [O1 O2 O3 O4]

modelStrxr.RLTwoParam.name='RLTwoParam';
modelStrxr.RLTwoParam.initialvalue_D2_mcherry=[26.25	5	1.25	67.5];
modelStrxr.RLTwoParam.initialvalue_D2___hm3dq=[26.25	5	1.25	67.5];
modelStrxr.RLTwoParam.initialvalue_D2___hm4di=[26.25	5	1.25	67.5];
modelStrxr.RLTwoParam.paramcount=2;

modelStrxr.RLSplitBetaRev.name='RLSplitBetaRev';
modelStrxr.RLSplitBetaRev.initialvalue_D2_mcherry=[26.25	5	1.25	67.5];
modelStrxr.RLSplitBetaRev.initialvalue_D2___hm3dq=[26.25	5	1.25	67.5];
modelStrxr.RLSplitBetaRev.initialvalue_D2___hm4di=[26.25	5	1.25	67.5];
modelStrxr.RLSplitBetaRev.paramcount=3;

modelStrxr.RLSplitAlphaRev.name='RLSplitAlphaRev';
modelStrxr.RLSplitAlphaRev.initialvalue_D2_mcherry=[26.25	5	1.25	67.5];
modelStrxr.RLSplitAlphaRev.initialvalue_D2___hm3dq=[26.25	5	1.25	67.5];
modelStrxr.RLSplitAlphaRev.initialvalue_D2___hm4di=[26.25	5	1.25	67.5];

modelStrxr.RLSplitAlphaRev.paramcount=3;

modelStrxr.RLSplitAlphaBetaRev.name='RLSplitAlphaBetaRev';
modelStrxr.RLSplitAlphaBetaRev.initialvalue_D2_mcherry=[26.25	5	1.25	67.5];
modelStrxr.RLSplitAlphaBetaRev.initialvalue_D2___hm3dq=[26.25	5	1.25	67.5];
modelStrxr.RLSplitAlphaBetaRev.initialvalue_D2___hm4di=[26.25	5	1.25	67.5];

modelStrxr.RLSplitAlphaBetaRev.paramcount=4;

%%%END Final Models

end