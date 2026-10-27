function [Ma,Ta,typeTa,ttout,L] = multiplicationmatrix(a,typea,Nmat,typein)
% [MATRIX,TENSOR,typeT,TYPEOUT,L] = multiplicationmatrix(V,TYPE,N,TYPEIN)
% Matrix MATRIX and tensor TENSOR (of type typeT) of the (left) multiplication by V, a tensor of
% Fourier ("Fou"), Chebyshev ("Che"), Taylor ("Tay"), cosine ("Cos") or sine ("Sin") coefficients
% encoded by TYPE; TYPEOUT is the output type.
% Only TYPE = "Fou" is used here: multiplicationmatrix(V,"Fou",N) is the (2N+1)x(2N+1) Toeplitz
% matrix (V_{n-k})_{|n|,|k|<=N} of the multiplication by sum_n V_n e^{int} (e.g. L_2 of (62)).
%   N      : orders of the output (first row) and input (second row) in each dimension; a single
%            row is used for both (default: size(V)-1)
%   TYPEIN : input type; by default "Cos" acts on "Cos", and "Sin" acts on "Sin" (output "Cos")
% The "Pha" dimensions of V, if any, store the L(1) x L(2) entries of a Jacobian (otherwise
% L = [1;1]).

% full type, with the two "Pha" dimensions
Npha=length(find(strcmp(typea,"Pha")));
typea=reshape(typea,[1,length(typea)]);
if Npha==0
    fulltype=[typea,"Pha","Pha"];
elseif Npha==1
    fulltype=[typea,"Pha"];
elseif Npha==2
    fulltype=typea;
else
    error('the type is not compatible with the tensor')
end
tt=fulltype(~strcmp(fulltype,"Pha"));
MM=length(tt);

if ~exist('typein') || isempty(typein)
    ttin=tt;
else
    ttin=typein(~strcmp(typein,"Pha"));
    ttin=reshape(ttin,[1,length(ttin)]);
end

% sizes
lsizea=length(size(a));
sizea=[size(a),ones(1,MM+2-lsizea)];
sza=size(a,1:MM);
L=size(a,MM+1:MM+2);
a=reshape(a,[prod(sza),L]);
prodsza=[cumprod(sza(2:MM),"reverse"),1];

Na = findnumberofmodes(sza,tt);
findfourier=find(strcmp("Fou",tt));
findsine=find(strcmp("Sin",tt));

if ~exist('Nmat','var') || isempty(Nmat)
    Nmat=sza-1;
end
if size(Nmat,1)==1
    % square case
    Nmat(2,:)=Nmat(1,:);
end

ttout=tt;
cossintosin = ( strcmp("Cos",tt) & strcmp("Sin",ttin) );
sinsintocos = ( strcmp("Sin",tt) & strcmp("Sin",ttin) );
sincostosin = ( strcmp("Sin",tt) & strcmp("Cos",ttin) );
ttout(cossintosin)="Sin";
ttout(sinsintocos)="Cos";

sizeM = Nmat+1;
sizeM(:,findfourier)=2*Nmat(:,findfourier)+1;
sizeM(:,findsine)=Nmat(:,findsine);
sizeM(:,cossintosin)=sizeM(:,cossintosin)-1;
sizeM(1,sinsintocos)=sizeM(1,sinsintocos)+1;
sizeM(2,sincostosin)=sizeM(2,sincostosin)+1;

% preallocate
useintervals = isa(a(1),'intval');
Ta = altzeros([prod(sizeM(1,:)),prod(sizeM(2,:)),L],useintervals);

for d=1:L(1)
    for dd=1:L(2)
        % each entry of the Jacobian separately
        Mad=a(:,d,dd);
        sizeMad=[1;1];
        for m=1:MM
            % recursively over the dimensions
            Nam=Na(m);
            K=Nmat(:,m);
            switch tt(m)
            % Fourier
                case {"Fou"}
                    Kmax=max(K);
                    NK=min(Nam,2*Kmax);
                    % toeplitz matrix of indices
                    T_K = toeplitz([(Nam+1:Nam+1+NK) zeros(1,2*Kmax-NK)]',[(Nam+1:-1:Nam+1-NK) zeros(1,2*Kmax-NK)]);
                    T_K = T_K(Kmax+1+(-K(1):K(1)),Kmax+1+(-K(2):K(2)));

                    % reshape with next dimension in first dimension of the tensor
                    Mad = reshape(Mad,[sizeMad(1:2)',sizea(m),prodsza(m)]);
                    % concatenate zero matrix
                    Madzero = altzeros([size(Mad,1:2),1,size(Mad,4)],useintervals);
                    % (cat(3,.,.) does not always work for intervals with singleton dimensions)
                    Mad = permute([permute(Madzero,[3 1 2 4]);permute(Mad,[3 1 2 4])],[2 3 1 4]);
                    if ~useintervals
                        % turn into cell array of matrices
                        Mad = mat2cell(Mad, sizeMad(1), sizeMad(2), ones(2*Nam+2,1), size(Mad,4));
                        % Toeplitz arrangement of the blocks, back to a matrix
                        Mad = cell2mat(Mad(T_K+1));
                    else
                        % for intervals: midpoints and radii separately
                        Mad_mid = mat2cell(mid(Mad), sizeMad(1), sizeMad(2), ones(2*Nam+2,1), size(Mad,4));
                        Mad_mid = cell2mat(Mad_mid(T_K+1));
                        Mad_rad = mat2cell(rad(Mad), sizeMad(1), sizeMad(2), ones(2*Nam+2,1), size(Mad,4));
                        Mad_rad = cell2mat(Mad_rad(T_K+1));
                        Mad = midrad(Mad_mid,Mad_rad);
                    end
                case {"Che","Cos"}
                % Chebyshev and Cos
                    NK1=min(Nam,K(1));
                    NK2=min(Nam,K(2));
                    T_K = toeplitz([1:NK1+1 zeros(1,K(1)-NK1)]',[1:NK2+1 zeros(1,K(2)-NK2)]);
                    % also the hankel part
                    Kmax=max([Nam,K(1),K(2)]);
                    H_K = hankel([1:Nam+1 zeros(1,Kmax-Nam)]');
                    H_K = H_K(1:K(1)+1,1:K(2)+1);
                    H_K(:,1)=0;
                    sH=1;
                    if strcmp(ttin(m),"Sin")
                        % if working on "Sin"
                        T_K=T_K(2:K(1)+1,2:K(2)+1);
                        H_K=H_K(2:K(1)+1,2:K(2)+1);
                        sH=-1;
                    end
                    Mad = reshape(Mad,[sizeMad(1:2)',sizea(m),prodsza(m)]);
                    Madzero = altzeros([size(Mad,1:2),1,size(Mad,4)],useintervals);
                    Mad = permute([permute(Madzero,[3 1 2 4]);permute(Mad,[3 1 2 4])],[2 3 1 4]);
                    if ~useintervals
                        Mad = mat2cell(Mad, sizeMad(1), sizeMad(2), ones(Nam+2,1), size(Mad,4));
                        Mad = cell2mat(Mad(T_K+1))+sH*cell2mat(Mad(H_K+1));
                    else
                        Mad_mid = mat2cell(mid(Mad), sizeMad(1), sizeMad(2), ones(Nam+2,1), size(Mad,4));
                        Mad_mid_T = cell2mat(Mad_mid(T_K+1));
                        Mad_mid_H = cell2mat(Mad_mid(H_K+1));
                        Mad_rad = mat2cell(rad(Mad), sizeMad(1), sizeMad(2), ones(Nam+2,1), size(Mad,4));
                        Mad_rad_T = cell2mat(Mad_rad(T_K+1));
                        Mad_rad_H = cell2mat(Mad_rad(H_K+1));
                        % add rigorously
                        Mad = midrad(Mad_mid_T,Mad_rad_T)+sH*midrad(Mad_mid_H,Mad_rad_H);
                    end
                case {"Sin"}
                % Sine
                    NK1=min(Nam,K(1));
                    NK2=min(Nam,K(2));
                    % lower triangular
                    T1_K = toeplitz([0:NK1 zeros(1,K(1)-NK1)]',[zeros(1,K(2)+1)]);
                    % upper triangular
                    T2_K = toeplitz([zeros(K(1)+1,1)],[0:NK2 zeros(1,K(2)-NK2)]);
                    % also the hankel part
                    Kmax=max([Nam,K(1)+1,K(2)+1]);
                    H_K = hankel([0:Nam zeros(1,Kmax-Nam)]');
                    if strcmp(ttin(m),"Sin")
                        % if working on "Sin"
                        T1_K = T1_K(1:K(1)+1,2:K(2)+1);
                        T2_K = T2_K(1:K(1)+1,2:K(2)+1);
                        H_K = H_K(1:K(1)+1,2:K(2)+1);
                        sH=-1;
                    else
                        % if working on "Cos"
                        T1_K = T1_K(2:K(1)+1,1:K(2)+1);
                        T2_K = T2_K(2:K(1)+1,1:K(2)+1);
                        H_K = H_K(2:K(1)+1,1:K(2)+1);
                        H_K(:,1)=0;
                        sH=1;
                    end
                    Mad = reshape(Mad,[sizeMad(1:2)',sizea(m),prodsza(m)]);
                    Madzero = altzeros([size(Mad,1:2),1,size(Mad,4)],useintervals);
                    Mad = permute([permute(Madzero,[3 1 2 4]);permute(Mad,[3 1 2 4])],[2 3 1 4]);
                    % cell2mat does not work nicely if T_K is a vector
                    if size(T1_K,2)==1
                        flagup=true;
                        T1_K=repmat(T1_K,1,2);
                        T2_K=repmat(T2_K,1,2);
                        H_K=repmat(H_K,1,2);
                    else
                        flagup=false;
                    end
                    if ~useintervals
                        Mad = mat2cell(Mad, sizeMad(1), sizeMad(2), ones(Nam+1,1), size(Mad,4));
                        Mad = cell2mat(Mad(T1_K+1))-cell2mat(Mad(T2_K+1))+sH*cell2mat(Mad(H_K+1));
                    else
                        Mad_mid = mat2cell(mid(Mad), sizeMad(1), sizeMad(2), ones(Nam+1,1), size(Mad,4));
                        Mad_mid_T1 = cell2mat(Mad_mid(T1_K+1));
                        Mad_mid_T2 = cell2mat(Mad_mid(T2_K+1));
                        Mad_mid_H = cell2mat(Mad_mid(H_K+1));
                        Mad_rad = mat2cell(rad(Mad), sizeMad(1), sizeMad(2), ones(Nam+1,1), size(Mad,4));
                        Mad_rad_T1 = cell2mat(Mad_rad(T1_K+1));
                        Mad_rad_T2 = cell2mat(Mad_rad(T2_K+1));
                        Mad_rad_H = cell2mat(Mad_rad(H_K+1));
                        % subtract rigorously
                        Mad = midrad(Mad_mid_T1,Mad_rad_T1)-midrad(Mad_mid_T2,Mad_rad_T2);
                        Mad = Mad+sH*midrad(Mad_mid_H,Mad_rad_H);
                    end
                    if flagup
                        Mad=Mad(:,1:prod(sizeM(2,:)));
                    end
                case {"Tay"}
                % Taylor
                    NK1=min(Nam,K(1));
                    % lower triangular
                    T_K = toeplitz([1:NK1+1 zeros(1,K(1)-NK1)]',[1 zeros(1,K(2))]);
                    Mad = reshape(Mad,[sizeMad(1:2)',sizea(m),prodsza(m)]);
                    Madzero = altzeros([size(Mad,1:2),1,size(Mad,4)],useintervals);
                    Mad = permute([permute(Madzero,[3 1 2 4]);permute(Mad,[3 1 2 4])],[2 3 1 4]);
                    if ~useintervals
                        Mad = mat2cell(Mad, sizeMad(1), sizeMad(2), ones(Nam+2,1), size(Mad,4));
                        Mad = cell2mat(Mad(T_K+1));
                    else
                        Mad_mid = mat2cell(mid(Mad), sizeMad(1), sizeMad(2), ones(Nam+2,1), size(Mad,4));
                        Mad_mid = cell2mat(Mad_mid(T_K+1));
                        Mad_rad = mat2cell(rad(Mad), sizeMad(1), sizeMad(2), ones(Nam+2,1), size(Mad,4));
                        Mad_rad = cell2mat(Mad_rad(T_K+1));
                        Mad = midrad(Mad_mid,Mad_rad);
                    end
            end
            % size of the matrix at the current depth
            sizeMad = sizeMad.*sizeM(:,m);
        end
        Ta(:,:,d,dd)=Mad;
    end
end

% final reshaping
Ta2=Ta;
Ta=reshape(Ta,[sizeM(1,:),sizeM(2,:),sizea(MM+1:MM+2)]);
Ta=permute(Ta,[1:MM,2*MM+1,MM+1:2*MM,2*MM+2]);
if Npha==0 && all(strcmp(ttin,ttout))
    Ta=reshape(Ta,[sizeM(1,:),sizeM(2,:)]);
    typeTa=tt;
else
    typeTa=[ttout,"Pha",ttin,"Pha"];
    ttout=[ttout,"Pha"];
end
Ma=tensortomatrix(Ta,typeTa);
L=reshape(L,[length(L),1]);

end

function z = altzeros(n,useintervals)
% zeros(n), as an intval array if useintervals

if length(n)==1
  disp('please check code; possible problem in altzeros; you have to bracket dimensions')
end

z = zeros(n);
if useintervals
  z=intval(z);
end

end
