function [cormacellcoll, reg, dirsa, corcells] = getkernelcormat2(kernelMats, wi, nb, numcomp, numperm, reggrid, gsdval)


% calculate apects of the window, like cent (center index), offset (offset
% to the left), and offsetright.  This depends parity of window length.
if ~mod(nb, 2)
    %even
    cent=floor(mean(1:nb));
    offset=(wi/2) - 1;
    offsetright=wi/2;
else
    %odd
    cent=mean(1:nb);
    offset=floor(wi/2);
    offsetright=offset;
end
lasty=size(kernelMats, 2)-nb+1;
lastcent=lasty+cent-1;

firstnum=cent-offset;
lastnum=lastcent-offset;

dircell=cell(2, length(firstnum:wi:lastnum));


%Determine regularization parameter
reg=reggrid;
regnorm=zeros(length(reg),1);
pacorma=@processcorma;
randpermsize=size(kernelMats{1,1}, 2);
parfor inda=1:length(regnorm)
    reguse=reg(inda);
    cormacellcoll=cell(2, 1);
    cormacellcoll{1}=pacorma(1:randpermsize, reguse, 1);
    cormacellcoll{2}=pacorma(randperm(randpermsize), reguse, 1);
    regnorm(inda)=norm(cormacellcoll{1}-cormacellcoll{2}, 2);
end
[~,tmp]=max(regnorm);
reg=reg(tmp);





cormacellcoll=cell(numperm+1, 1);
cormacellcoll{1}=processcorma(1:size(kernelMats{1,1}, 2), reg);
dirsa=dircell;
randpermsize=size(kernelMats{1,1}, 2);
if numperm>10
    
    pacorma=@processcorma;
    parfor ind=1:numperm
        randoperm=randperm(randpermsize);
        cormacellcoll{ind+1}=pacorma(randoperm, reg);
        %            sum((randoperm-(1:randpermsize))==0)
    end
else
    for ind=1:numperm
        cormacellcoll{ind+1}=processcorma(randperm(randpermsize), reg);
    end
end



%%
    function [cormacell]=processcorma(v, rega, getrflag)
        switch nargin
            case 2
                getrflag=0;
        end
        k=0;
        
        numr=zeros(size(kernelMats{1,1},1), length(1:wi:lasty));
        for lind=1:wi:lasty
            
            k=k+1;
            finindi=lind+nb-1;
            indi=lind:finindi;
            kA=kernelMats{1,lind};
            kB=kernelMats{2,lind};
            for mind=2:length(indi)
                kA=kA+kernelMats{1,indi(mind)};
                kB=kB+kernelMats{2,indi(mind)};
            end
            kA=kA(v,v);
            [nalph, nbet, r]=kcanonca_reg_ver2(kA, kB, gsdval, rega, 0); %.999995
            
            numr(1:length(r),k)=r(end:-1:1);
            
            dircell{1,k}=nalph(:,end:-1:(end-numcomp+1));
            dircell{2,k}=nbet(:,end:-1:(end-numcomp+1));
            
        end
        
        if getrflag==1
            cormacell=numr;
        else
            corX=zeros(size(kernelMats{1,1}, 2), size(kernelMats, 2), numcomp);
            Rx=zeros(numcomp, numcomp, size(kernelMats, 2));
            varx=zeros(numcomp, size(kernelMats, 2));
            vary=varx;
            Ry=Rx;
            corY=corX;
            finnum=cent+offsetright;
            
            for lind=1:finnum
                [corX(:,lind,:), Rx(:,:,lind), varx(:,lind)]=calculateArrayQuantities(1, lind, 1, v);
                [corY(:,lind,:), Ry(:,:,lind), vary(:,lind)]=calculateArrayQuantities(1, lind, 2, v);
            end
            if finnum ~= size(kernelMats, 2)
                for lind=(finnum+1):(size(kernelMats, 2)-(finnum+1))
                    use=floor((lind-firstnum)/wi)+1;
                    [corX(:,lind,:), Rx(:,:,lind), varx(:,lind)]=calculateArrayQuantities(use, lind, 1, v);
                    [corY(:,lind,:), Ry(:,:,lind), vary(:,lind)]=calculateArrayQuantities(use, lind, 2, v);
                end
                for lind=(size(kernelMats, 2)-finnum):size(kernelMats,2)
                    [corX(:,lind,:), Rx(:,:,lind), varx(:,lind)]=calculateArrayQuantities(size(dircell, 2), lind, 1, v);
                    [corY(:,lind,:), Ry(:,:,lind), vary(:,lind)]=calculateArrayQuantities(size(dircell, 2), lind, 2, v);
                end
            end

            cormacell=corr(corX(:,:,1), corY(:,:,1));
            corcells=cell(2,1);
            corcells{1} = corX;
            corcells{2} = corY;
            if numcomp>1
                
                for jind=2:numcomp
                    tx=zeros(size(kernelMats{1,1}, 2),size(kernelMats,2));
                    ty=tx;
                    corXjind=tx;
                    corYjind=ty;
                    for knd=1:size(kernelMats,2)
                        tx(:,knd)=squeeze(corX(:,knd,1:(jind-1)))*Rx(1:(jind-1), jind, knd);
                        ty(:,knd)=squeeze(corY(:,knd,1:(jind-1)))*Rx(1:(jind-1), jind, knd);
                        corXjind(:,knd)=squeeze(corX(:,knd,jind))*Rx(jind, jind, knd);
                        corYjind(:,knd)=squeeze(corY(:,knd,jind))*Ry(jind, jind, knd);
                    end
                    papa=(getcov(corXjind, corYjind)+getcov(corXjind, ty)+ getcov(tx, corYjind))./sqrt(varx(jind,:)*vary(jind,:)');
                    cormacell=cormacell+papa;
                end
                
            end
            
            
        end
    end
    function [covmat]=getcov(matA, matB)
        matA=bsxfun(@minus,matA,mean(matA,1)); %%% zero-mean
        matB=bsxfun(@minus,matB,mean(matB,1)); %%% zero-mean
        covmat=abs(matA'*matB/size(matA, 1));
    end

    function [A, B, varset]=calculateArrayQuantities(i, j, k, v)
        tmpX=(dircell{k,i}'*kernelMats{k,j}(v,v))';
        varset=var(tmpX);
        [A, B]=qr(tmpX);
        A=A(:,1:size(B,2));
        B=B(1:size(B, 2), :);
    end
end

