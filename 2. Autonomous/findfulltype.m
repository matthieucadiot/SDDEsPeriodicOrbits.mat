function [fulltype,newsizeT]=findfulltype(typeT,sizeT)
% Full type of a tensor of type typeT ("Pha" dimensions added if necessary) and its size in this
% full type (sizeT with ones inserted if necessary).

typeT=reshape(typeT,[1,length(typeT)]);
Npha0=length(find(strcmp(typeT,"Pha")));
dimT=length(sizeT);
if dimT==2 && sizeT(2)==1
    dimT=1;
end
Ntt=length(typeT);
if dimT<Ntt
    % trailing singleton dimensions (e.g. N = 1 for "Sin")
    sizeT=[sizeT,ones(1,Ntt-dimT)];
    dimT=Ntt;
elseif dimT>Ntt && dimT<2*Ntt
    % trailing singleton dimensions (e.g. N = 1 for "Sin")
    sizeT=[sizeT,ones(1,2*Ntt-dimT)];
    dimT=2*Ntt;
end
fulltype=typeT;
newsizeT=sizeT;
if Npha0==0
    if dimT==2*Ntt
        fulltype=[typeT,"Pha",typeT,"Pha"];
        newsizeT=[sizeT(1:length(sizeT)/2),1,sizeT(length(sizeT)/2+1:length(sizeT)),1];
    elseif dimT==Ntt
        fulltype=[typeT,"Pha"];
        newsizeT=[sizeT,1];
    end
elseif Npha0==1 && dimT==2*Ntt
    fulltype=[typeT,typeT];
end

end
