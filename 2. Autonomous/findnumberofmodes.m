function NT=findnumberofmodes(sizeT,typeT)
% Number of modes NT in each function dimension of a tensor of size sizeT and type typeT.

sizeT=[sizeT,ones(1,length(typeT)-length(sizeT))];
notpha=find(~strcmp(typeT,"Pha"));
tt=typeT(notpha);
NT=sizeT(notpha)-1;
findfourier=find(strcmp("Fou",tt));
NT(findfourier)=NT(findfourier)/2;
findsine=find(strcmp("Sin",tt));
NT(findsine)=NT(findsine)+1;

end
