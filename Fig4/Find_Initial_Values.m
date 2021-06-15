discrimination_choice_D2mcherry=xlsread('/Users/kdelevic/Dropbox/D2 DMS OVX project/eLife submission/transparent reporting info/modeling for github/Fig4/Fig4_TrialHistories','D2_mcherry_Disc');
discrimination_choice_D2hm3dq=xlsread('/Users/kdelevic/Dropbox/D2 DMS OVX project/eLife submission/transparent reporting info/modeling for github/Fig4/Fig4_TrialHistories','D2___hm3dq_Disc');
discrimination_choice_D2hm4di=xlsread('/Users/kdelevic/Dropbox/D2 DMS OVX project/eLife submission/transparent reporting info/modeling for github/Fig4/Fig4_TrialHistories','D2___hm4di_Disc');

datasize1=size(discrimination_choice_D2mcherry)
datasize2=size(discrimination_choice_D2hm3dq)
datasize3=size(discrimination_choice_D2hm4di)


for i=1:datasize1(2)

data1=discrimination_choice_D2mcherry(:,i);
data1(isnan(data1)) = [];
data1(data1==0)=[];

data1stxr{i}=data1

end

for i=1:datasize2(2)

data2=discrimination_choice_D2hm3dq(:,i);
data2(isnan(data2)) = [];
data2(data2==0)=[];

data2stxr{i}=data2

end


for i=1:datasize3(2)

data3=discrimination_choice_D2hm4di(:,i);
data3(isnan(data3)) = [];
data3(data3==0)=[];

data3stxr{i}=data3

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

data31= cellfun(@(x) (x(1)), data3stxr);
data32= cellfun(@(x) (x(2)), data3stxr);
data33= cellfun(@(x) (x(3)), data3stxr);
data34= cellfun(@(x) (x(4)), data3stxr);
fulldata3=[data31 data32 data33 data34];


trialcounts=tabulate([fulldata1 fulldata2 fulldata3]);