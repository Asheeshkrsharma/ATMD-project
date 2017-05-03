function chromosomes=mutation(chromosomes,mutation_rate)
    for i=1:size(chromosomes)
        chromoeo=generate_chromosomes(1);%Get a dummy chromosome, from which an element will be choosen
        chromosome=chromosomes(i,:);
        for j=1:size(chromosome,2)
            if rand(1,1) > mutation_rate
                chromosome(j)=chromoeo(j);
            end
        end
        chromosomes(i,:)=chromosome;
    end
end