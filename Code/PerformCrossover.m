%Author: Asheesh Sharma
%License: UIA (Use It Anywhere)
%Disclaimer: The UIA licesne applies to any part of this code except for
%any where it is specifically mentioned. You are given this software for
%free so dont try to sell it and do not bother me if some thing is broken.
%Finally, to all the students out there. The challenge is not to complete a
%task. It is to not plagiarize. :D 
%PS: This is based on a book called Orbital dynamics. You can find it for
%free online.
function [ newChromosome1 , newChromosome2 ] = PerformCrossover(c1,c2)
    evaluate_chromosome(c1);
    evaluate_chromosome(c2);
    pos1=randi(length(c1)-1) ;
    newChromosome1 = [ c1(1:pos1) c2(pos1+1: length(c1))] ;
    newChromosome2 = [ c2(1:pos1) c1(pos1+1: length(c1))] ;
end