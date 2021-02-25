function [TxSecondProjWallj, TxSecondReflecWallj] = ImageTx(tempRef,wall)
    for i = 1:size(wall.xyz1,1) % only works for first Tx 
        TxSecondProjWallj.xyz(:,:,i) = repmat((dot((wall.xyz1 - repmat(tempRef.xyz(i,:,1),size(wall.xyz1,1),1)),...
            wall.unitNormal.xyz,2)./dot(wall.unitNormal.xyz,wall.unitNormal.xyz,2)),1,3).* wall.unitNormal.xyz + ...
            repmat(tempRef.xyz(i,:,1),size(wall.unitNormal.xyz,1),1); % only works for 1 TX so (:,:,1)
        TxSecondReflecWallj.xyz(:,:,i) = repmat(tempRef.xyz(i,:,1),size(wall.unitNormal.xyz,1),1) + 2.*...
            (TxSecondProjWallj.xyz(:,:,i) -repmat(tempRef.xyz(i,:,1),size(wall.unitNormal.xyz,1),1));
    end
end