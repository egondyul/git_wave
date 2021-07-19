function FracCaseMatExport
freq=10000;%really start frequency
dfreq=20000;

for i=1:1 %session
    for j=1:2 %cases
        for k=1:2 %frequency
            out(k,j)=ResExport(freq,j,[num2str(i) '/']);
            freq=freq+dfreq;
        end  
        freq=10000;
    end
end

save([num2str(i) '.mat'],'out')
clear out;
end