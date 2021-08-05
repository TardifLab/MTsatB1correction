%% Lookup table approach to MTsat with arbitrary readouts. 
function MTsat = calcMTsatThruLookupTable(inputMTw_img, b1, T1_1, mask_1, M0_1, echoSpacing, numExcitation, TD, flip)

% i = 1;
% inputMTw_img = Dual_comb(:,:,:,i);
% b1 = B1_gauss(:,:,:,i);
% T1_1 = T1(:,:,:,i);
% mask_1 = mask(:,:,:,i);
% M0_1 = M0(:,:,:,i);
% % echoSpacing, numExcitation, TD, flip
% 
% i = 60;j = 45; k = 40; 


MTsat = zeros(size(mask_1));
MTdrop_table = 0:0.005:0.40;
simSig = zeros(size(MTdrop_table));

tic
for i = 1:size(mask_1,1)
    for j = 1:size(mask_1,2)
        for k = 1:size(mask_1,3)
            if mask_1(i,j,k) > 0 && T1_1(i,j,k) > 500 && T1_1(i,j,k) < 5000 && b1(i,j,k) > 0.4% use for masked data
                
                for z = 1:size(MTdrop_table,2)
                    simSig(z) = MTrage_sig_eqn_v5(echoSpacing, flip, T1_1(i,j,k), TD, numExcitation, M0_1(i,j,k), MTdrop_table(z), b1(i,j,k), 1);  
                                                
                end
                
                MTsat(i,j,k) = interp1( simSig, MTdrop_table, inputMTw_img(i,j,k));

            end
        end
    end
end

toc
