
function [S,atomInfo]=H_sparse_gen(t)



atomInfo=[];
iInds =[];
jInds =[];
atomInd = 0;

for nnzCol=1:t


      tpoint = 1;
    
    for i=1:t-nnzCol+1
        atomInd = atomInd+1;
        

        
        AtomStart =  tpoint;
        AtomEnd  =  tpoint+nnzCol-1;
        tpoint = tpoint+1;
    
        atomInfo(atomInd,1) =  atomInd; 
        atomInfo(atomInd,2) =  AtomStart; 
        atomInfo(atomInd,3) =  AtomEnd;
        % concat new col inds
        jInds = [jInds,AtomStart:AtomEnd ];
        
        % row ind are not changeing per atom 
        rowInds = atomInd.*ones(1,nnzCol);
        iInds = [iInds,rowInds];
        
    
    end
end

    oneVec= ones(1,length(   iInds));
    
    S = sparse(iInds,jInds,oneVec);


end