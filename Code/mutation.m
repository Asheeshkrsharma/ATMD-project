%Author: Asheesh Sharma
%License: UIA (Use It Anywhere)
%Disclaimer: The UIA licesne applies to any part of this code except for
%any where it is specifically mentioned. You are given this software for
%free so dont try to sell it and do not bother me if some thing is broken.
%Finally, to all the students out there. The challenge is not to complete a
%task. It is to not plagiarize. :D 
%PS: This is based on a book called Orbital dynamics. You can find it for
%free online.
function chromosomes=mutation(chromosomes,mutation_rate)
    for i=1:size(chromosomes)
        chromoeo=generate_chromosomes(1);%Get a dummy chromosome, from which an element will be choosen
        chromosome=chromosomes(i,:);
        for j=1:size(chromosome,2)
            if rand(1,1) < mutation_rate
                chromosome(j)=chromoeo(j);
            end
        end
        chromosomes(i,:)=chromosome;
    end
end