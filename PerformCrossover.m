function [ newChromosome1 , newChromosome2 ] = PerformCrossover(c1,c2)
    pos1=randi(length(c1)-1) ;
    newChromosome1 = [ c1(1:pos1) c2(pos1+1: length(c1))] ;
    newChromosome2 = [ c2(1:pos1) c1(pos1+1: length(c1))] ;
end