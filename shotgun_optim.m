clc;
clear;
format long
tic
chromosomes=generate_chromosomes(100);
dt=[];
dtBig=[]

for i=1:10
    dt=[]
    for i=1:100
        [dV1, dV2, dV3, dV4, dV5,  dV6]=evaluate_chromosome(chromosomes(i,:));
        total=norm(dV2-dV1)+norm(dV3-dV2)+norm(dV4-dV3)+norm(dV5-dV4)+norm(dV6-dV5);
        dt=[dt;total];
    end
    dtBig=[dtBig dt];
end
%% Plots
dtBig=dtBig';
figure(1)
hold on;
mins=[];
dims=[];
for i=1:size(dtBig,2)
    dt=round([dtBig(:,i) ones(length(dtBig(:,i)),1)*i]);
    dims=[dims; round([dtBig(:,i) ones(length(dtBig(:,i)),1)*i])];
%    scatter(dt(:,2),dt(:,1),'s','b','flat');
    d=min(dt(:,1));
    mins=[mins;d i];
end
plot(mins(:,2),mins(:,1),'b','LineWidth',1);
