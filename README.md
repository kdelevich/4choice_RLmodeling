# 4choice_RLmodeling

Step-by-step guide to fit trial-by-trial multiple choice reversal task behavior:
1) Run Find_Initial_Values.m to calculate the initial odor values (1—4); calculated from the first 4 choices during the Discrimination phase
2) Open LoadModels_reversal.m; if necessary, input animals names and initial odor values calculated from step 1.
3) Open FitRL_Test.m and input model selection. Set directory for xlsx file that contains trial history data. Input .mat file name on line 277 if desired.
4) Run FitRL_Test.m

Step-by-step guide to generate trial data using best-fit parameters:
1) Open Four_Choice_Simulator.m and input .mat file corresponding to parameters for desired model. It does not matter which cohort you choose for the .mat
file; all cohorts will be simulated.
2) 
