%Author: Asheesh Sharma
%License: UIA (Use It Anywhere)
%Disclaimer: The UIA licesne applies to any part of this code except for
%any where it is specifically mentioned. You are given this software for
%free so dont try to sell it and do not bother me if some thing is broken.
%Finally, to all the students out there. The challenge is not to complete a
%task. It is to not plagiarize. :D 
%Usual stuff
clc;
clear;
format long
number_chromosomes=500;
ga_iterations=100;
elite_percentage=0.1;
mutation_percentage=0.4;
best=[10^50,0,0,0,0,0,0,0,0]; %This is a dummy chromosome which is replaced
%by the best chromosome every time. This can be replace by an initial
%solution.The first number is the fitness, set really high so it can be
%replace by something sensible.

%Generate initial population. Look inside the function to see the ranges.
%These were decided emperically. A better way would have been something
%like the interval approach. But this was faster since the time window for
%the mission is already small (Something withing 2 year since epoch).
chromosomes=generate_chromosomes(number_chromosomes);
best_population=[];
display_chromosomes=[];
display_ies=[];
display_C3s=[];
display_vectors=[];
display_phies=[];
min_dv=[];
ss=[];
for i=1:ga_iterations
    %Evalues the fitness of current popoulation
    fitness=[];
    C3=[];
    vectors=[];
    for j=1:size(chromosomes,1)
        [total,C3_total,vector]=evaluate_chromosome(chromosomes(j,:));
        C3=[C3;C3_total];
        fitness=[fitness;total];
        vectors=[vectors;vector];
    end
    S=(max(fitness)-min(fitness))/min(fitness);%Scatter function
    W2=(i/ga_iterations)*100; %Scatter weight
    W1=100-W2; %delta V weight
    fitness_real=fitness; %Actual delta V
    ss=[ss;(S*W2)+90];
    fitness=(W1*fitness)+(S*W2);
    elites=[];
    %Elite are kept as they are and used to populate next generation.
    % For this we need the sorted indces
    [~,I]=sort(fitness,'ascend');
    display_chromosomes=[display_chromosomes; fitness_real];
    if best(1)>fitness(I(1));
         best=[fitness_real(I(1)) chromosomes(I(1),:)];
         best_population=[best_population;best];
         display_C3s=[display_C3s; min(C3)];
         display_ies=[display_ies;ones(size(fitness_real,1),1)*i];
         display_phies=[display_phies;ones(size(fitness_real,1),1)*i];
         display_vectors=[display_vectors;vectors];
         min_dv=[min_dv;vectors(I(1),:)];
    else
         display_ies=[display_ies;ones(size(fitness_real,1),1)*i];
    end
    % next we get the "elite_percentage" worth elite chromosome from the
    % population
    elite=chromosomes(I(1:round(size(I,1)*elite_percentage)),:);
    while ~isempty(elite)
        eI=randi(size(elite,1),2,1); %generate two random integers
        %Cross over them parents
        [c1,c2]=PerformCrossover(elite(eI(1),:),elite(eI(2),:));
        %Evaluate the fitness of new chromosomes and parents.
        [t1,~]=evaluate_chromosome(c1);
        [t2,~]=evaluate_chromosome(c2);
        [e1,~]=evaluate_chromosome(elite(eI(1),:));
        [e2,~]=evaluate_chromosome(elite(eI(2),:));
        chromos=[c1;c2;elite(eI(1),:);elite(eI(2),:)];
        chromos_fit=[t1;t2;e1;e2];
        [~,io]=sort(chromos_fit,'ascend');
        child1=chromos(io(1),:);
        child2=chromos(io(2),:);
        elite(eI',:)=[];
        %Populate the new population
        elites=[elites;child1;child2];
    end
    chromosomes(I(1:round(size(I,1)*elite_percentage)),:)=[];
    %Mutation 
    mutated=mutation(chromosomes,mutation_percentage);
    chromosomes=[elites;mutated];
end
figure(1)
hold on
ylim([80 100]);
txt1 =num2str(['Minimum |\DeltaV|=',num2str(round(min(best_population(:,1)))),' Km/s']);
txt2=num2str(['Minimum C_3=',num2str(round(min(display_C3s)*10^-3)),'\times 10^{-3} Km^2/s^2'])
title({'GA \DeltaV minimisation results';txt1;txt2})
xlabel('iterations')
ylabel('|\DeltaV| (Km/s)')
h1=scatter(display_ies,display_chromosomes,5,'s','flat')
d=unique(display_ies);
dp=[];
for i=1:size(d,1)
    dp=[dp; d(i) min(display_chromosomes((find(display_ies == d(i)))))];
end
h2=plot(dp(:,1),dp(:,2),'r','LineWidth',1);
legend([h1,h2],{'Individual |\DeltaV|','[Minimum |\DeltaV|','spread'});
hold off;
figure(2)
display_optim(best(2:9));