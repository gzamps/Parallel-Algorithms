clear all;
close all;
%easy periptwsh 1000 Nodes

size =1000000;



%-
%dianysma Po
Po = rand(size,1);
Po = Po/5;


%-
%pinakas sundesewn


for i=1:size
    pe=randperm(size,15);
    
    A(i,1:15)=pe(1:15);

end

   
%-
%dianusma E
E=rand(size,1);

fid = fopen('P1000000.bin', 'w', 'l');
fwrite(fid, Po, 'double');
fclose(fid);

fid = fopen('E1000000.bin', 'w', 'l');
fwrite(fid, E, 'double');
fclose(fid);


fid = fopen('G1000000.bin', 'w', 'l');
fwrite(fid, A,'integer*4');
fclose(fid);



