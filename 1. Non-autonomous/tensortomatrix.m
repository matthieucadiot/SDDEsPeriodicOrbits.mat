function [MT,phaseT] = tensortomatrix(T,typeT)
% Reshapes a tensor T of type typeT, representing a vector or a matrix, into the matrix MT;
% phaseT is obtained by transposing only the "Pha" dimensions (useful when MT is a vector).
% The first "Pha" dimension separates the output (row) and input (column) dimensions. If typeT
% contains function types only, T is a column vector if it has at most length(typeT) dimensions
% (an n x 1 array has one), and a matrix of type [typeT,typeT] if it has more (at most
% 2*length(typeT)).

[fulltype,sizeT]=findfulltype(typeT,size(T));
T=reshape(T,sizeT);

Npha=length(find(strcmp(fulltype,"Pha")));
if Npha==1
    % a vector is a matrix with one column
    fulltype=[fulltype,"Pha"];
    sizeT=[sizeT,1];
end

% output and input dimensions
findphase=find(strcmp(fulltype,"Pha"));
tout=fulltype(1:findphase(1)-1);
tin=fulltype(findphase(1)+1:length(fulltype)-1);
MMout=length(tout);
MMin=length(tin);
sizeT=[sizeT,ones(1,MMout+MMin+2-length(sizeT))]; % ones added for the case MMin = 0

% reshape into a matrix
MTcol=prod(sizeT(1:MMout+1));
MTrow=prod(sizeT(MMout+2:MMout+MMin+2));
MT=reshape(T,[MTcol,MTrow]);

% reshape into the phase-transposed matrix
phaseT=permute(T,[1:MMout,MMout+MMin+2,MMout+2:MMout+MMin+1,MMout+1]);
szphaseT=[size(phaseT),ones(1,MMout+MMin+2-length(size(phaseT)))];
phaseTcol=prod(szphaseT(1:MMout+1));
phaseTrow=prod(szphaseT(MMout+2:MMout+MMin+2));
phaseT=reshape(phaseT,[phaseTcol,phaseTrow]);

end
