discrimination_choice_FSham=xlsread('/Users/kdelevic/Dropbox/D2 DMS OVX project/eLife submission/transparent reporting info/modeling for github/TrialHistories','F_Sham_Disc');
discrimination_choice_FOVX=xlsread('/Users/kdelevic/Dropbox/D2 DMS OVX project/eLife submission/transparent reporting info/modeling for github/TrialHistories','F__OVX_Disc');
datasize1=size(discrimination_choice_FSham)
datasize2=size(discrimination_choice_FOVX)

for i=1:datasize1(2)

data1=discrimination_choice_FSham(:,i);
data1(isnan(data1)) = [];
data1(data1==0)=[];

data1stxr{i}=data1

end

for i=1:datasize2(2)

data2=discrimination_choice_FOVX(:,i);
data2(isnan(data2)) = [];
data2(data2==0)=[];

data2stxr{i}=data2

end

data11= cellfun(@(x) (x(1)), data1stxr);
data12= cellfun(@(x) (x(2)), data1stxr);
data13= cellfun(@(x) (x(3)), data1stxr);
data14= cellfun(@(x) (x(4)), data1stxr);
fulldata1=[data11 data12 data13 data14];

data21= cellfun(@(x) (x(1)), data2stxr);
data22= cellfun(@(x) (x(2)), data2stxr);
data23= cellfun(@(x) (x(3)), data2stxr);
data24= cellfun(@(x) (x(4)), data2stxr);
fulldata2=[data21 data22 data23 data24];


trialcounts=tabulate([fulldata1 fulldata2]);