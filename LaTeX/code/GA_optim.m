%Usual stuff
clc;
clear;
format long
number_chromosomes=100;
ga_iterations=500;
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
for i=1:ga_iterations
    %Evalues the fitness of current popoulation
    fitness=[];
    for j=1:size(chromosomes,1)
        [total,C3_total]=evaluate_chromosome(chromosomes(j,:));
        fitness=[fitness;total];
    end
    elites=[];
    %Elite are kept as they are and used to populate next generation.
    % For this we need the sorted indces
    [~,I]=sort(fitness,'ascend');
    if best(1)>fitness(I(1));
         best=[fitness(I(1)) chromosomes(I(1),:)];
         best_population=[best_population;best];
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
plot(best_population(:,1));
