
function [badAtoms] = badAtomFinder(atomInfo,atomEndInfo,atomStartInfo,SVfirstOc,SVs,EVfirstOc,EVs,CS,CE,numAtoms)

   
      
       %% find atom start between CE and CS; 
     startOfstarts = SVfirstOc(SVs==CS);

     ind =find(SVs==CE);
     endOfstarts = SVfirstOc( ind);

     % this can occur if CE is greater than any start 
     if  isempty(endOfstarts)
 
          % use all possible starts after startOfstarts 
          endOfstarts = length(atomStartInfo);
     end
    
     %% Atoms with start value equal to CE
         % get when number of atoms with this end 
      if ind==length(SVs)
           atomsSVEqCE = atomStartInfo(SVfirstOc(ind):end,1);
      else
            atomsSVEqCE = atomStartInfo(SVfirstOc(ind):SVfirstOc(ind+1),1);
      end



        
     atomsWithBadStarts = [atomStartInfo( startOfstarts:endOfstarts,1); atomsSVEqCE];





       %% find atom end between CS and CE; 
       startOfends = EVfirstOc(EVs==CE);

       inds =find(EVs==CS);
     

         % this can occur if CE is greater than any start 
      if  isempty( inds)
 
          % use all possible starts after startOfstarts 
          endOfends = length(atomEndInfo);
      else
            endOfends = EVfirstOc(inds);

      end

           %% Atoms with end value equal to CS
         % get when number of atoms with this end 
      if ind==1
           atomsEVEqCS = atomEndInfo(EVfirstOc(ind):end,1);
      else
            atomsEVEqCS = atomEndInfo(EVfirstOc(ind):EVfirstOc(ind+1),1);
      end
         
        
        
      atomsWithBadEnds = [atomEndInfo( startOfends: endOfends,1);atomsEVEqCS];

        %% atoms who start before CS and end after CE 
      atomsWithStartsBfCS = atomStartInfo( 1:startOfstarts-1,1);
      atomsWithEndsAfterCE = atomEndInfo( 1:  startOfends-1,1);
      atomsSpanningWholeBand = intersect( atomsWithEndsAfterCE  ,atomsWithStartsBfCS);

      badAtoms =unique( [atomsWithBadStarts; atomsWithBadEnds;  atomsSpanningWholeBand ]);


end

