%% global settings
fig_prefix = 'lifetime';
N_S = 110;
tau = 0.777660157519;
end_t = 0.9;

% species renaming
% "npropyloxy", "nRO"
% "npropylooh", "nROOH"
% "npropyloo", "nROO"
% "npropyl", "nR"
% "ipropyloxy", "iRO"
% "ipropylooh", "iROOH"
% "ipropyloo", "iROO"
% "ipropyl", "iR"
% "well_1", "O2QOOH1"
% "prod_1", "OQ" + '\prime' + "OOH1"
% "frag_1", "OQ" + '\prime' + "O1"

spe_name_latex = {'HE', 'AR', 'N_2', 'H', 'H_2', 'HEG', 'HCFG', 'Halt', 'O', 'O_2', 'OH', 'H_2O', 'HO_2', 'H_2O_2', 'CO', 'CO_2', 'HCO', 'CH_2O', 'formyloxy', 'formylperoxy', 'formylooh', 'C', 'CH', 'CH_2(S)', 'CH_2', 'CH_3', 'CH_4', 'CH_3OO', 'CH_3OOH', 'CH_2OH', 'CH_3O', 'CH_3OH', 'OCH_2OOH', 'HOCH_2OO', 'HOCH_2OOH', 'HOCH_2O', 'C_2H_2', 'C_2H_3', 'C_2H_4', 'C_2H_5', 'C_2H_6', 'ketene', 'HCCO', 'ethynol', 'acetaldehyde', 'acetyl', 'vinoxy', 'acetylperoxy', 'acetyloxy', 'ethoxy', 'CH_3CH_2OO', 'CH_3CH_2OOH', 'CH_2CH_2OOH', 'ethanol', 'CH_2CH_2OH', 'CH_3CHOH', 'oxirane', 'oxiranyl', 'ethenol', 'C_3H_6', 'nR', 'iR', 'C_3H_8', 'allyl', 'propen1yl', 'propen2yl', 'acetone', 'CH_3OCH_3', 'CH_3OCH_2', 'acrolein', 'CH_2CHCO', 'CH_3CHCO', 'allyloxy', 'allyl-alcohol', 'propanal', 'propionyl', 'nRO', 'iRO', 'nROO', 'nROOH', 'iROO', 'iROOH', 'propen1ol', 'propen2ol', 'CH_3CO_3H', 'O_2C_2H_4OH', 'propoxide', 'QOOH_1', 'QOOH_2', 'QOOH_3', 'O_2QOOH_1', 'well_2', 'well_3', 'well_5', 'OQ^\primeOOH_1', 'prod_2', 'prod_3', 'prod_4', 'prod_5', 'prod_6', 'prod_7', 'OQ^\primeO_1', 'frag_3', 'frag_4', 'frag_5', 'propen1oxy', 'propen2oxy', 'glyoxal', 'vinoxylmethyl', 'formylethyl'};
reaction_name_latex = {'H+O_2 \rightarrow O+OH', 'O+OH \rightarrow H+O_2', 'O+H_2 \rightarrow H+OH', 'H+OH \rightarrow O+H_2', 'H_2+OH \rightarrow H_2O+H', 'H_2O+H \rightarrow H_2+OH', 'O+H_2O \rightarrow OH+OH', 'OH+OH \rightarrow O+H_2O', 'H_2+M \rightarrow H+H+M', 'H+H+M \rightarrow H_2+M', 'H_2+AR \rightarrow H+H+AR', 'H+H+AR \rightarrow H_2+AR', 'H_2+HE \rightarrow H+H+HE', 'H+H+HE \rightarrow H_2+HE', 'O+O+M \rightarrow O_2+M', 'O_2+M \rightarrow O+O+M', 'O+O+AR \rightarrow O_2+AR', 'O_2+AR \rightarrow O+O+AR', 'O+O+HE \rightarrow O_2+HE', 'O_2+HE \rightarrow O+O+HE', 'O+H+M \rightarrow OH+M', 'OH+M \rightarrow O+H+M', 'H+OH+M \rightarrow H_2O+M', 'H_2O+M \rightarrow H+OH+M', 'H+O_2(+M) \rightarrow HO_2(+M)', 'HO_2(+M) \rightarrow H+O_2(+M)', 'HO_2+H \rightarrow H_2+O_2', 'H_2+O_2 \rightarrow HO_2+H', 'HO_2+H \rightarrow OH+OH', 'OH+OH \rightarrow HO_2+H', 'HO_2+O \rightarrow O_2+OH', 'O_2+OH \rightarrow HO_2+O', 'HO_2+OH \rightarrow H_2O+O_2', 'H_2O+O_2 \rightarrow HO_2+OH', 'HO_2+HO_2 \rightarrow H_2O_2+O_2', 'H_2O_2+O_2 \rightarrow HO_2+HO_2', 'H_2O_2(+M) \rightarrow OH+OH(+M)', 'OH+OH(+M) \rightarrow H_2O_2(+M)', 'H_2O_2+H \rightarrow H_2O+OH', 'H_2O+OH \rightarrow H_2O_2+H', 'H_2O_2+H \rightarrow HO_2+H_2', 'HO_2+H_2 \rightarrow H_2O_2+H', 'H_2O_2+O \rightarrow OH+HO_2', 'OH+HO_2 \rightarrow H_2O_2+O', 'H_2O_2+OH \rightarrow HO_2+H_2O', 'HO_2+H_2O \rightarrow H_2O_2+OH', 'CO+O(+M) \rightarrow CO_2(+M)', 'CO_2(+M) \rightarrow CO+O(+M)', 'CO+O_2 \rightarrow CO_2+O', 'CO_2+O \rightarrow CO+O_2', 'CO+HO_2 \rightarrow CO_2+OH', 'CO_2+OH \rightarrow CO+HO_2', 'CO+OH \rightarrow CO_2+H', 'CO_2+H \rightarrow CO+OH', 'HCO+M \rightarrow H+CO+M', 'H+CO+M \rightarrow HCO+M', 'HCO+O_2 \rightarrow CO+HO_2', 'CO+HO_2 \rightarrow HCO+O_2', 'HCO+H \rightarrow CO+H_2', 'CO+H_2 \rightarrow HCO+H', 'HCO+O \rightarrow CO+OH', 'CO+OH \rightarrow HCO+O', 'HCO+OH \rightarrow CO+H_2O', 'CO+H_2O \rightarrow HCO+OH', 'HCO+O \rightarrow CO_2+H', 'CO_2+H \rightarrow HCO+O', 'HCO+HO_2 \rightarrow CO_2+OH+H', 'CO_2+OH+H \rightarrow HCO+HO_2', 'HCO+CH_3 \rightarrow CO+CH_4', 'CO+CH_4 \rightarrow HCO+CH_3', 'HCO+HCO \rightarrow H_2+CO+CO', 'H_2+CO+CO \rightarrow HCO+HCO', 'HCO+HCO \rightarrow CH_2O+CO', 'CH_2O+CO \rightarrow HCO+HCO', 'HCO+O_2 \rightarrow formylperoxy', 'formylperoxy \rightarrow HCO+O_2', 'CH_2O+formylperoxy \rightarrow HCO+formylooh', 'HCO+formylooh \rightarrow CH_2O+formylperoxy', 'formylooh \rightarrow formyloxy+OH', 'formyloxy+OH \rightarrow formylooh', 'H+CO_2+M \rightarrow formyloxy+M', 'formyloxy+M \rightarrow H+CO_2+M', 'CH_2O+M \rightarrow HCO+H+M', 'HCO+H+M \rightarrow CH_2O+M', 'CH_2O+M \rightarrow CO+H_2+M', 'CO+H_2+M \rightarrow CH_2O+M', 'CH_2O+H \rightarrow HCO+H_2', 'HCO+H_2 \rightarrow CH_2O+H', 'CH_2O+O \rightarrow HCO+OH', 'HCO+OH \rightarrow CH_2O+O', 'CH_2O+OH \rightarrow HCO+H_2O', 'HCO+H_2O \rightarrow CH_2O+OH', 'CH_2O+O_2 \rightarrow HCO+HO_2', 'HCO+HO_2 \rightarrow CH_2O+O_2', 'CH_2O+HO_2 \rightarrow HCO+H_2O_2', 'HCO+H_2O_2 \rightarrow CH_2O+HO_2', 'CH_2O+CH_3 \rightarrow HCO+CH_4', 'HCO+CH_4 \rightarrow CH_2O+CH_3', 'CH_2O+HO_2 \rightarrow OCH_2OOH', 'OCH_2OOH \rightarrow CH_2O+HO_2', 'OCH_2OOH \rightarrow HOCH_2OO', 'HOCH_2OO \rightarrow OCH_2OOH', 'HOCH_2OO+HO_2 \rightarrow HOCH_2OOH+O_2', 'HOCH_2OOH+O_2 \rightarrow HOCH_2OO+HO_2', 'HOCH_2O+OH \rightarrow HOCH_2OOH', 'HOCH_2OOH \rightarrow HOCH_2O+OH', 'CH_3+O \rightarrow CH_2O+H', 'CH_2O+H \rightarrow CH_3+O', 'CH_3+O_2 \rightarrow CH_3O+O', 'CH_3O+O \rightarrow CH_3+O_2', 'CH_3+HO_2 \rightarrow CH_3O+OH', 'CH_3O+OH \rightarrow CH_3+HO_2', 'CH_3+HO_2 \rightarrow CH_4+O_2', 'CH_4+O_2 \rightarrow CH_3+HO_2', 'CH_3+CH_3(+M) \rightarrow C_2H_6(+M)', 'C_2H_6(+M) \rightarrow CH_3+CH_3(+M)', 'CH_3+H(+M) \rightarrow CH_4(+M)', 'CH_4(+M) \rightarrow CH_3+H(+M)', 'CH_4+H \rightarrow CH_3+H_2', 'CH_3+H_2 \rightarrow CH_4+H', 'CH_4+O \rightarrow CH_3+OH', 'CH_3+OH \rightarrow CH_4+O', 'CH_4+OH \rightarrow CH_3+H_2O', 'CH_3+H_2O \rightarrow CH_4+OH', 'CH_4+HO_2 \rightarrow CH_3+H_2O_2', 'CH_3+H_2O_2 \rightarrow CH_4+HO_2', 'CH_3+CH_3OH \rightarrow CH_4+CH_3O', 'CH_4+CH_3O \rightarrow CH_3+CH_3OH', 'CH_3O+CH_3 \rightarrow CH_2O+CH_4', 'CH_2O+CH_4 \rightarrow CH_3O+CH_3', 'CH_3O+H \rightarrow CH_2O+H_2', 'CH_2O+H_2 \rightarrow CH_3O+H', 'CH_3+O_2(+M) \rightarrow CH_3OO(+M)', 'CH_3OO(+M) \rightarrow CH_3+O_2(+M)', 'CH_3OO+CH_2O \rightarrow CH_3OOH+HCO', 'CH_3OOH+HCO \rightarrow CH_3OO+CH_2O', 'CH_4+CH_3OO \rightarrow CH_3+CH_3OOH', 'CH_3+CH_3OOH \rightarrow CH_4+CH_3OO', 'CH_3OH+CH_3OO \rightarrow CH_2OH+CH_3OOH', 'CH_2OH+CH_3OOH \rightarrow CH_3OH+CH_3OO', 'CH_3OO+CH_3 \rightarrow CH_3O+CH_3O', 'CH_3O+CH_3O \rightarrow CH_3OO+CH_3', 'CH_3OO+HO_2 \rightarrow CH_3OOH+O_2', 'CH_3OOH+O_2 \rightarrow CH_3OO+HO_2', 'CH_3OO+CH_3OO \rightarrow CH_2O+CH_3OH+O_2', 'CH_2O+CH_3OH+O_2 \rightarrow CH_3OO+CH_3OO', 'CH_3OO+CH_3OO \rightarrow O_2+CH_3O+CH_3O', 'O_2+CH_3O+CH_3O \rightarrow CH_3OO+CH_3OO', 'CH_3OO+H \rightarrow CH_3O+OH', 'CH_3O+OH \rightarrow CH_3OO+H', 'CH_3OO+O \rightarrow CH_3O+O_2', 'CH_3O+O_2 \rightarrow CH_3OO+O', 'CH_3OO+OH \rightarrow CH_3OH+O_2', 'CH_3OH+O_2 \rightarrow CH_3OO+OH', 'CH_3OOH \rightarrow CH_3O+OH', 'CH_3O+OH \rightarrow CH_3OOH', 'CH_2OH+M \rightarrow CH_2O+H+M', 'CH_2O+H+M \rightarrow CH_2OH+M', 'CH_2OH+H \rightarrow CH_2O+H_2', 'CH_2O+H_2 \rightarrow CH_2OH+H', 'CH_2OH+H \rightarrow CH_3+OH', 'CH_3+OH \rightarrow CH_2OH+H', 'CH_2OH+O \rightarrow CH_2O+OH', 'CH_2O+OH \rightarrow CH_2OH+O', 'CH_2OH+OH \rightarrow CH_2O+H_2O', 'CH_2O+H_2O \rightarrow CH_2OH+OH', 'CH_2OH+O_2 \rightarrow CH_2O+HO_2', 'CH_2O+HO_2 \rightarrow CH_2OH+O_2', 'CH_2OH+HO_2 \rightarrow CH_2O+H_2O_2', 'CH_2O+H_2O_2 \rightarrow CH_2OH+HO_2', 'CH_2OH+HCO \rightarrow CH_3OH+CO', 'CH_3OH+CO \rightarrow CH_2OH+HCO', 'CH_2OH+HCO \rightarrow CH_2O+CH_2O', 'CH_2O+CH_2O \rightarrow CH_2OH+HCO', '2CH_2OH \rightarrow CH_3OH+CH_2O', 'CH_3OH+CH_2O \rightarrow 2CH_2OH', 'CH_2OH+CH_3O \rightarrow CH_3OH+CH_2O', 'CH_3OH+CH_2O \rightarrow CH_2OH+CH_3O', 'CH_2OH+HO_2 \rightarrow HOCH_2O+OH', 'HOCH_2O+OH \rightarrow CH_2OH+HO_2', 'CH_3O+M \rightarrow CH_2O+H+M', 'CH_2O+H+M \rightarrow CH_3O+M', 'CH_3O+H \rightarrow CH_3+OH', 'CH_3+OH \rightarrow CH_3O+H', 'CH_3O+O \rightarrow CH_2O+OH', 'CH_2O+OH \rightarrow CH_3O+O', 'CH_3O+OH \rightarrow CH_2O+H_2O', 'CH_2O+H_2O \rightarrow CH_3O+OH', 'CH_3O+O_2 \rightarrow CH_2O+HO_2', 'CH_2O+HO_2 \rightarrow CH_3O+O_2', 'CH_3O+HO_2 \rightarrow CH_2O+H_2O_2', 'CH_2O+H_2O_2 \rightarrow CH_3O+HO_2', 'CH_3O+CO \rightarrow CH_3+CO_2', 'CH_3+CO_2 \rightarrow CH_3O+CO', 'CH_3O+HCO \rightarrow CH_3OH+CO', 'CH_3OH+CO \rightarrow CH_3O+HCO', '2CH_3O \rightarrow CH_3OH+CH_2O', 'CH_3OH+CH_2O \rightarrow 2CH_3O', 'OH+CH_3(+M) \rightarrow CH_3OH(+M)', 'CH_3OH(+M) \rightarrow OH+CH_3(+M)', 'H+CH_2OH(+M) \rightarrow CH_3OH(+M)', 'CH_3OH(+M) \rightarrow H+CH_2OH(+M)', 'H+CH_3O(+M) \rightarrow CH_3OH(+M)', 'CH_3OH(+M) \rightarrow H+CH_3O(+M)', 'CH_3OH+H \rightarrow CH_2OH+H_2', 'CH_2OH+H_2 \rightarrow CH_3OH+H', 'CH_3OH+H \rightarrow CH_3O+H_2', 'CH_3O+H_2 \rightarrow CH_3OH+H', 'CH_3OH+O \rightarrow CH_2OH+OH', 'CH_2OH+OH \rightarrow CH_3OH+O', 'CH_3OH+OH \rightarrow CH_3O+H_2O', 'CH_3O+H_2O \rightarrow CH_3OH+OH', 'CH_3OH+OH \rightarrow CH_2OH+H_2O', 'CH_2OH+H_2O \rightarrow CH_3OH+OH', 'CH_3OH+O_2 \rightarrow CH_2OH+HO_2', 'CH_2OH+HO_2 \rightarrow CH_3OH+O_2', 'CH_3OH+HCO \rightarrow CH_2OH+CH_2O', 'CH_2OH+CH_2O \rightarrow CH_3OH+HCO', 'CH_3OH+HO_2 \rightarrow CH_2OH+H_2O_2', 'CH_2OH+H_2O_2 \rightarrow CH_3OH+HO_2', 'CH_3OH+CH_3 \rightarrow CH_2OH+CH_4', 'CH_2OH+CH_4 \rightarrow CH_3OH+CH_3', 'CH_3O+CH_3OH \rightarrow CH_3OH+CH_2OH', 'CH_3OH+CH_2OH \rightarrow CH_3O+CH_3OH', 'CH_3+CH_3 \rightarrow H+C_2H_5', 'H+C_2H_5 \rightarrow CH_3+CH_3', 'CH_4+CH_2 \rightarrow CH_3+CH_3', 'CH_3+CH_3 \rightarrow CH_4+CH_2', 'CH_3+OH \rightarrow CH_2+H_2O', 'CH_2+H_2O \rightarrow CH_3+OH', 'CH_3+CH_2 \rightarrow C_2H_4+H', 'C_2H_4+H \rightarrow CH_3+CH_2', 'CH_2+H(+M) \rightarrow CH_3(+M)', 'CH_3(+M) \rightarrow CH_2+H(+M)', 'CH_2+O \rightarrow HCO+H', 'HCO+H \rightarrow CH_2+O', 'CH_2+OH \rightarrow CH_2O+H', 'CH_2O+H \rightarrow CH_2+OH', 'CH_2+H_2 \rightarrow H+CH_3', 'H+CH_3 \rightarrow CH_2+H_2', 'CH_2+O_2 \rightarrow HCO+OH', 'HCO+OH \rightarrow CH_2+O_2', 'CH_2+HO_2 \rightarrow CH_2O+OH', 'CH_2O+OH \rightarrow CH_2+HO_2', 'CH_2+CO(+M) \rightarrow ketene(+M)', 'ketene(+M) \rightarrow CH_2+CO(+M)', 'CH_2+CH_2 \rightarrow C_2H_2+H_2', 'C_2H_2+H_2 \rightarrow CH_2+CH_2', 'CH_2(S)+M \rightarrow CH_2+M', 'CH_2+M \rightarrow CH_2(S)+M', 'CH_2(S)+H_2O \rightarrow CH_2+H_2O', 'CH_2+H_2O \rightarrow CH_2(S)+H_2O', 'CH_2(S)+CO \rightarrow CH_2+CO', 'CH_2+CO \rightarrow CH_2(S)+CO', 'CH_2(S)+CO_2 \rightarrow CH_2+CO_2', 'CH_2+CO_2 \rightarrow CH_2(S)+CO_2', 'CH_2(S)+AR \rightarrow CH_2+AR', 'CH_2+AR \rightarrow CH_2(S)+AR', 'CH_2(S)+O \rightarrow CO+H_2', 'CO+H_2 \rightarrow CH_2(S)+O', 'CH_2(S)+O \rightarrow HCO+H', 'HCO+H \rightarrow CH_2(S)+O', 'CH_2(S)+OH \rightarrow CH_2O+H', 'CH_2O+H \rightarrow CH_2(S)+OH', 'CH_2(S)+H_2 \rightarrow CH_3+H', 'CH_3+H \rightarrow CH_2(S)+H_2', 'CH_2(S)+O_2 \rightarrow H+OH+CO', 'H+OH+CO \rightarrow CH_2(S)+O_2', 'CH_2(S)+O_2 \rightarrow CO+H_2O', 'CO+H_2O \rightarrow CH_2(S)+O_2', 'CH_2(S)+CO_2 \rightarrow CH_2O+CO', 'CH_2O+CO \rightarrow CH_2(S)+CO_2', 'CH_2(S)+C_2H_6 \rightarrow CH_3+C_2H_5', 'CH_3+C_2H_5 \rightarrow CH_2(S)+C_2H_6', 'CH_2(S)+ketene \rightarrow C_2H_4+CO', 'C_2H_4+CO \rightarrow CH_2(S)+ketene', 'H+HCCO \rightarrow CH_2(S)+CO', 'CH_2(S)+CO \rightarrow H+HCCO', 'CH_2(S)+H \rightarrow CH+H_2', 'CH+H_2 \rightarrow CH_2(S)+H', 'CH_2+H \rightarrow CH+H_2', 'CH+H_2 \rightarrow CH_2+H', 'CH_2+OH \rightarrow CH+H_2O', 'CH+H_2O \rightarrow CH_2+OH', 'CH+O_2 \rightarrow HCO+O', 'HCO+O \rightarrow CH+O_2', 'CH+H \rightarrow C+H_2', 'C+H_2 \rightarrow CH+H', 'CH+O \rightarrow CO+H', 'CO+H \rightarrow CH+O', 'CH+OH \rightarrow HCO+H', 'HCO+H \rightarrow CH+OH', 'CH+H_2O \rightarrow H+CH_2O', 'H+CH_2O \rightarrow CH+H_2O', 'CH+CO_2 \rightarrow HCO+CO', 'HCO+CO \rightarrow CH+CO_2', 'HCCO+M \rightarrow CH+CO+M', 'CH+CO+M \rightarrow HCCO+M', 'CH+CH_2O \rightarrow H+ketene', 'H+ketene \rightarrow CH+CH_2O', 'CH+HCCO \rightarrow CO+C_2H_2', 'CO+C_2H_2 \rightarrow CH+HCCO', 'CH+CH_4 \rightarrow C_2H_4+H', 'C_2H_4+H \rightarrow CH+CH_4', 'C_2H_6+CH \rightarrow C_2H_5+CH_2', 'C_2H_5+CH_2 \rightarrow C_2H_6+CH', 'C_2H_5+H(+M) \rightarrow C_2H_6(+M)', 'C_2H_6(+M) \rightarrow C_2H_5+H(+M)', 'C_2H_6+H \rightarrow C_2H_5+H_2', 'C_2H_5+H_2 \rightarrow C_2H_6+H', 'C_2H_6+O \rightarrow C_2H_5+OH', 'C_2H_5+OH \rightarrow C_2H_6+O', 'C_2H_6+OH \rightarrow C_2H_5+H_2O', 'C_2H_5+H_2O \rightarrow C_2H_6+OH', 'C_2H_6+O_2 \rightarrow C_2H_5+HO_2', 'C_2H_5+HO_2 \rightarrow C_2H_6+O_2', 'C_2H_6+CH_3 \rightarrow C_2H_5+CH_4', 'C_2H_5+CH_4 \rightarrow C_2H_6+CH_3', 'C_2H_6+HO_2 \rightarrow C_2H_5+H_2O_2', 'C_2H_5+H_2O_2 \rightarrow C_2H_6+HO_2', 'C_2H_6+CH_3OO \rightarrow C_2H_5+CH_3OOH', 'C_2H_5+CH_3OOH \rightarrow C_2H_6+CH_3OO', 'C_2H_6+CH_3O \rightarrow C_2H_5+CH_3OH', 'C_2H_5+CH_3OH \rightarrow C_2H_6+CH_3O', 'H+C_2H_4(+M) \rightarrow C_2H_5(+M)', 'C_2H_5(+M) \rightarrow H+C_2H_4(+M)', 'H_2+CH_3OO \rightarrow H+CH_3OOH', 'H+CH_3OOH \rightarrow H_2+CH_3OO', 'H_2+CH_3CH_2OO \rightarrow H+CH_3CH_2OOH', 'H+CH_3CH_2OOH \rightarrow H_2+CH_3CH_2OO', 'C_2H_4+C_2H_4 \rightarrow C_2H_5+C_2H_3', 'C_2H_5+C_2H_3 \rightarrow C_2H_4+C_2H_4', 'CH_3+C_2H_5 \rightarrow CH_4+C_2H_4', 'CH_4+C_2H_4 \rightarrow CH_3+C_2H_5', 'C_2H_5+H \rightarrow C_2H_4+H_2', 'C_2H_4+H_2 \rightarrow C_2H_5+H', 'C_2H_5+O \rightarrow acetaldehyde+H', 'acetaldehyde+H \rightarrow C_2H_5+O', 'C_2H_5+HO_2 \rightarrow ethoxy+OH', 'ethoxy+OH \rightarrow C_2H_5+HO_2', 'CH_3OO+C_2H_5 \rightarrow CH_3O+ethoxy', 'CH_3O+ethoxy \rightarrow CH_3OO+C_2H_5', 'ethoxy+O_2 \rightarrow acetaldehyde+HO_2', 'acetaldehyde+HO_2 \rightarrow ethoxy+O_2', 'CH_3+CH_2O \rightarrow ethoxy', 'ethoxy \rightarrow CH_3+CH_2O', 'acetaldehyde+H \rightarrow ethoxy', 'ethoxy \rightarrow acetaldehyde+H', 'C_2H_5+O_2 \rightarrow CH_3CH_2OO', 'CH_3CH_2OO \rightarrow C_2H_5+O_2', 'CH_3CH_2OO+CH_2O \rightarrow CH_3CH_2OOH+HCO', 'CH_3CH_2OOH+HCO \rightarrow CH_3CH_2OO+CH_2O', 'CH_4+CH_3CH_2OO \rightarrow CH_3+CH_3CH_2OOH', 'CH_3+CH_3CH_2OOH \rightarrow CH_4+CH_3CH_2OO', 'CH_3OH+CH_3CH_2OO \rightarrow CH_2OH+CH_3CH_2OOH', 'CH_2OH+CH_3CH_2OOH \rightarrow CH_3OH+CH_3CH_2OO', 'CH_3CH_2OO+HO_2 \rightarrow CH_3CH_2OOH+O_2', 'CH_3CH_2OOH+O_2 \rightarrow CH_3CH_2OO+HO_2', 'C_2H_6+CH_3CH_2OO \rightarrow C_2H_5+CH_3CH_2OOH', 'C_2H_5+CH_3CH_2OOH \rightarrow C_2H_6+CH_3CH_2OO', 'CH_3CH_2OOH \rightarrow ethoxy+OH', 'ethoxy+OH \rightarrow CH_3CH_2OOH', 'C_2H_5+O_2 \rightarrow CH_2CH_2OOH', 'CH_2CH_2OOH \rightarrow C_2H_5+O_2', 'C_2H_5+O_2 \rightarrow C_2H_4+HO_2', 'C_2H_4+HO_2 \rightarrow C_2H_5+O_2', 'C_2H_5+O_2 \rightarrow oxirane+OH', 'oxirane+OH \rightarrow C_2H_5+O_2', 'C_2H_5+O_2 \rightarrow acetaldehyde+OH', 'acetaldehyde+OH \rightarrow C_2H_5+O_2', 'CH_2CH_2OOH \rightarrow CH_3CH_2OO', 'CH_3CH_2OO \rightarrow CH_2CH_2OOH', 'CH_3CH_2OO \rightarrow acetaldehyde+OH', 'acetaldehyde+OH \rightarrow CH_3CH_2OO', 'CH_3CH_2OO \rightarrow C_2H_4+HO_2', 'C_2H_4+HO_2 \rightarrow CH_3CH_2OO', 'CH_3CH_2OO \rightarrow oxirane+OH', 'oxirane+OH \rightarrow CH_3CH_2OO', 'CH_2CH_2OOH \rightarrow oxirane+OH', 'oxirane+OH \rightarrow CH_2CH_2OOH', 'CH_2CH_2OOH \rightarrow C_2H_4+HO_2', 'C_2H_4+HO_2 \rightarrow CH_2CH_2OOH', 'CH_2CH_2OOH \rightarrow acetaldehyde+OH', 'acetaldehyde+OH \rightarrow CH_2CH_2OOH', 'oxirane \rightarrow CH_3+HCO', 'CH_3+HCO \rightarrow oxirane', 'oxirane \rightarrow acetaldehyde', 'acetaldehyde \rightarrow oxirane', 'oxirane+OH \rightarrow oxiranyl+H_2O', 'oxiranyl+H_2O \rightarrow oxirane+OH', 'oxirane+H \rightarrow oxiranyl+H_2', 'oxiranyl+H_2 \rightarrow oxirane+H', 'oxirane+HO_2 \rightarrow oxiranyl+H_2O_2', 'oxiranyl+H_2O_2 \rightarrow oxirane+HO_2', 'oxirane+CH_3OO \rightarrow oxiranyl+CH_3OOH', 'oxiranyl+CH_3OOH \rightarrow oxirane+CH_3OO', 'oxirane+CH_3CH_2OO \rightarrow oxiranyl+CH_3CH_2OOH', 'oxiranyl+CH_3CH_2OOH \rightarrow oxirane+CH_3CH_2OO', 'oxirane+CH_3 \rightarrow oxiranyl+CH_4', 'oxiranyl+CH_4 \rightarrow oxirane+CH_3', 'oxirane+CH_3O \rightarrow oxiranyl+CH_3OH', 'oxiranyl+CH_3OH \rightarrow oxirane+CH_3O', 'oxiranyl \rightarrow acetyl', 'acetyl \rightarrow oxiranyl', 'oxiranyl \rightarrow vinoxy', 'vinoxy \rightarrow oxiranyl', 'CH_3+HCO \rightarrow acetaldehyde', 'acetaldehyde \rightarrow CH_3+HCO', 'acetaldehyde+H \rightarrow acetyl+H_2', 'acetyl+H_2 \rightarrow acetaldehyde+H', 'acetaldehyde+O \rightarrow acetyl+OH', 'acetyl+OH \rightarrow acetaldehyde+O', 'acetaldehyde+OH \rightarrow acetyl+H_2O', 'acetyl+H_2O \rightarrow acetaldehyde+OH', 'acetaldehyde+O_2 \rightarrow acetyl+HO_2', 'acetyl+HO_2 \rightarrow acetaldehyde+O_2', 'acetaldehyde+CH_3 \rightarrow acetyl+CH_4', 'acetyl+CH_4 \rightarrow acetaldehyde+CH_3', 'acetaldehyde+HO_2 \rightarrow acetyl+H_2O_2', 'acetyl+H_2O_2 \rightarrow acetaldehyde+HO_2', 'CH_3OO+acetaldehyde \rightarrow CH_3OOH+acetyl', 'CH_3OOH+acetyl \rightarrow CH_3OO+acetaldehyde', 'acetaldehyde+acetylperoxy \rightarrow acetyl+CH_3CO_3H', 'acetyl+CH_3CO_3H \rightarrow acetaldehyde+acetylperoxy', 'acetaldehyde+OH \rightarrow vinoxy+H_2O', 'vinoxy+H_2O \rightarrow acetaldehyde+OH', 'acetyl(+M) \rightarrow CH_3+CO(+M)', 'CH_3+CO(+M) \rightarrow acetyl(+M)', 'acetyl+H \rightarrow ketene+H_2', 'ketene+H_2 \rightarrow acetyl+H', 'acetyl+O \rightarrow ketene+OH', 'ketene+OH \rightarrow acetyl+O', 'acetyl+CH_3 \rightarrow ketene+CH_4', 'ketene+CH_4 \rightarrow acetyl+CH_3', 'acetyl+O_2 \rightarrow acetylperoxy', 'acetylperoxy \rightarrow acetyl+O_2', 'acetylperoxy+HO_2 \rightarrow CH_3CO_3H+O_2', 'CH_3CO_3H+O_2 \rightarrow acetylperoxy+HO_2', 'H_2O_2+acetylperoxy \rightarrow HO_2+CH_3CO_3H', 'HO_2+CH_3CO_3H \rightarrow H_2O_2+acetylperoxy', 'CH_4+acetylperoxy \rightarrow CH_3+CH_3CO_3H', 'CH_3+CH_3CO_3H \rightarrow CH_4+acetylperoxy', 'CH_2O+acetylperoxy \rightarrow HCO+CH_3CO_3H', 'HCO+CH_3CO_3H \rightarrow CH_2O+acetylperoxy', 'C_2H_6+acetylperoxy \rightarrow C_2H_5+CH_3CO_3H', 'C_2H_5+CH_3CO_3H \rightarrow C_2H_6+acetylperoxy', 'CH_3CO_3H \rightarrow acetyloxy+OH', 'acetyloxy+OH \rightarrow CH_3CO_3H', 'acetyloxy+M \rightarrow CH_3+CO_2+M', 'CH_3+CO_2+M \rightarrow acetyloxy+M', 'ketene+H \rightarrow vinoxy', 'vinoxy \rightarrow ketene+H', 'vinoxy+O_2 \rightarrow CH_2O+CO+OH', 'CH_2O+CO+OH \rightarrow vinoxy+O_2', 'ketene+H \rightarrow CH_3+CO', 'CH_3+CO \rightarrow ketene+H', 'ketene+H \rightarrow HCCO+H_2', 'HCCO+H_2 \rightarrow ketene+H', 'ketene+O \rightarrow CH_2+CO_2', 'CH_2+CO_2 \rightarrow ketene+O', 'ketene+O \rightarrow HCCO+OH', 'HCCO+OH \rightarrow ketene+O', 'ketene+OH \rightarrow HCCO+H_2O', 'HCCO+H_2O \rightarrow ketene+OH', 'ketene+OH \rightarrow CH_2OH+CO', 'CH_2OH+CO \rightarrow ketene+OH', 'HCCO+OH \rightarrow H_2+CO+CO', 'H_2+CO+CO \rightarrow HCCO+OH', 'HCCO+O \rightarrow H+CO+CO', 'H+CO+CO \rightarrow HCCO+O', 'HCCO+O_2 \rightarrow OH+CO+CO', 'OH+CO+CO \rightarrow HCCO+O_2', 'C_2H_3+H(+M) \rightarrow C_2H_4(+M)', 'C_2H_4(+M) \rightarrow C_2H_3+H(+M)', 'C_2H_4(+M) \rightarrow C_2H_2+H_2(+M)', 'C_2H_2+H_2(+M) \rightarrow C_2H_4(+M)', 'C_2H_4+H \rightarrow C_2H_3+H_2', 'C_2H_3+H_2 \rightarrow C_2H_4+H', 'C_2H_4+O \rightarrow CH_3+HCO', 'CH_3+HCO \rightarrow C_2H_4+O', 'C_2H_4+O \rightarrow vinoxy+H', 'vinoxy+H \rightarrow C_2H_4+O', 'C_2H_4+OH \rightarrow C_2H_3+H_2O', 'C_2H_3+H_2O \rightarrow C_2H_4+OH', 'C_2H_4+CH_3 \rightarrow C_2H_3+CH_4', 'C_2H_3+CH_4 \rightarrow C_2H_4+CH_3', 'C_2H_4+O_2 \rightarrow C_2H_3+HO_2', 'C_2H_3+HO_2 \rightarrow C_2H_4+O_2', 'C_2H_4+CH_3O \rightarrow C_2H_3+CH_3OH', 'C_2H_3+CH_3OH \rightarrow C_2H_4+CH_3O', 'C_2H_4+CH_3OO \rightarrow C_2H_3+CH_3OOH', 'C_2H_3+CH_3OOH \rightarrow C_2H_4+CH_3OO', 'C_2H_4+CH_3CH_2OO \rightarrow C_2H_3+CH_3CH_2OOH', 'C_2H_3+CH_3CH_2OOH \rightarrow C_2H_4+CH_3CH_2OO', 'C_2H_4+acetylperoxy \rightarrow C_2H_3+CH_3CO_3H', 'C_2H_3+CH_3CO_3H \rightarrow C_2H_4+acetylperoxy', 'C_2H_4+CH_3OO \rightarrow oxirane+CH_3O', 'oxirane+CH_3O \rightarrow C_2H_4+CH_3OO', 'C_2H_4+CH_3CH_2OO \rightarrow oxirane+ethoxy', 'oxirane+ethoxy \rightarrow C_2H_4+CH_3CH_2OO', 'C_2H_4+HO_2 \rightarrow oxirane+OH', 'oxirane+OH \rightarrow C_2H_4+HO_2', 'C_2H_2+H(+M) \rightarrow C_2H_3(+M)', 'C_2H_3(+M) \rightarrow C_2H_2+H(+M)', 'C_2H_3+O_2 \rightarrow HCO+CH_2O', 'HCO+CH_2O \rightarrow C_2H_3+O_2', 'C_2H_3+O_2 \rightarrow HO_2+C_2H_2', 'HO_2+C_2H_2 \rightarrow C_2H_3+O_2', 'C_2H_3+O_2 \rightarrow O+vinoxy', 'O+vinoxy \rightarrow C_2H_3+O_2', 'CH_3+C_2H_3 \rightarrow CH_4+C_2H_2', 'CH_4+C_2H_2 \rightarrow CH_3+C_2H_3', 'C_2H_3+H \rightarrow C_2H_2+H_2', 'C_2H_2+H_2 \rightarrow C_2H_3+H', 'C_2H_3+OH \rightarrow C_2H_2+H_2O', 'C_2H_2+H_2O \rightarrow C_2H_3+OH', 'C_2H_2+O_2 \rightarrow HCCO+OH', 'HCCO+OH \rightarrow C_2H_2+O_2', 'C_2H_2+O \rightarrow CH_2+CO', 'CH_2+CO \rightarrow C_2H_2+O', 'C_2H_2+O \rightarrow HCCO+H', 'HCCO+H \rightarrow C_2H_2+O', 'C_2H_2+OH \rightarrow ketene+H', 'ketene+H \rightarrow C_2H_2+OH', 'C_2H_2+OH \rightarrow CH_3+CO', 'CH_3+CO \rightarrow C_2H_2+OH', 'OH+C_2H_2 \rightarrow H+ethynol', 'H+ethynol \rightarrow OH+C_2H_2', 'H+ethynol \rightarrow H+ketene', 'H+ketene \rightarrow H+ethynol', 'ethanol(+M) \rightarrow CH_2OH+CH_3(+M)', 'CH_2OH+CH_3(+M) \rightarrow ethanol(+M)', 'ethanol(+M) \rightarrow C_2H_5+OH(+M)', 'C_2H_5+OH(+M) \rightarrow ethanol(+M)', 'ethanol(+M) \rightarrow C_2H_4+H_2O(+M)', 'C_2H_4+H_2O(+M) \rightarrow ethanol(+M)', 'ethanol(+M) \rightarrow acetaldehyde+H_2(+M)', 'acetaldehyde+H_2(+M) \rightarrow ethanol(+M)', 'ethanol+O_2 \rightarrow CH_2CH_2OH+HO_2', 'CH_2CH_2OH+HO_2 \rightarrow ethanol+O_2', 'ethanol+O_2 \rightarrow CH_3CHOH+HO_2', 'CH_3CHOH+HO_2 \rightarrow ethanol+O_2', 'ethanol+OH \rightarrow CH_2CH_2OH+H_2O', 'CH_2CH_2OH+H_2O \rightarrow ethanol+OH', 'ethanol+OH \rightarrow CH_3CHOH+H_2O', 'CH_3CHOH+H_2O \rightarrow ethanol+OH', 'ethanol+OH \rightarrow ethoxy+H_2O', 'ethoxy+H_2O \rightarrow ethanol+OH', 'ethanol+H \rightarrow CH_2CH_2OH+H_2', 'CH_2CH_2OH+H_2 \rightarrow ethanol+H', 'ethanol+H \rightarrow CH_3CHOH+H_2', 'CH_3CHOH+H_2 \rightarrow ethanol+H', 'ethanol+H \rightarrow ethoxy+H_2', 'ethoxy+H_2 \rightarrow ethanol+H', 'ethanol+HO_2 \rightarrow CH_2CH_2OH+H_2O_2', 'CH_2CH_2OH+H_2O_2 \rightarrow ethanol+HO_2', 'ethanol+HO_2 \rightarrow CH_3CHOH+H_2O_2', 'CH_3CHOH+H_2O_2 \rightarrow ethanol+HO_2', 'ethanol+HO_2 \rightarrow ethoxy+H_2O_2', 'ethoxy+H_2O_2 \rightarrow ethanol+HO_2', 'ethanol+CH_3OO \rightarrow CH_2CH_2OH+CH_3OOH', 'CH_2CH_2OH+CH_3OOH \rightarrow ethanol+CH_3OO', 'ethanol+CH_3OO \rightarrow CH_3CHOH+CH_3OOH', 'CH_3CHOH+CH_3OOH \rightarrow ethanol+CH_3OO', 'ethanol+CH_3OO \rightarrow ethoxy+CH_3OOH', 'ethoxy+CH_3OOH \rightarrow ethanol+CH_3OO', 'ethanol+O \rightarrow CH_2CH_2OH+OH', 'CH_2CH_2OH+OH \rightarrow ethanol+O', 'ethanol+O \rightarrow CH_3CHOH+OH', 'CH_3CHOH+OH \rightarrow ethanol+O', 'ethanol+O \rightarrow ethoxy+OH', 'ethoxy+OH \rightarrow ethanol+O', 'ethanol+CH_3 \rightarrow CH_2CH_2OH+CH_4', 'CH_2CH_2OH+CH_4 \rightarrow ethanol+CH_3', 'ethanol+CH_3 \rightarrow CH_3CHOH+CH_4', 'CH_3CHOH+CH_4 \rightarrow ethanol+CH_3', 'ethanol+CH_3 \rightarrow ethoxy+CH_4', 'ethoxy+CH_4 \rightarrow ethanol+CH_3', 'ethanol+C_2H_5 \rightarrow CH_2CH_2OH+C_2H_6', 'CH_2CH_2OH+C_2H_6 \rightarrow ethanol+C_2H_5', 'ethanol+C_2H_5 \rightarrow CH_3CHOH+C_2H_6', 'CH_3CHOH+C_2H_6 \rightarrow ethanol+C_2H_5', 'C_2H_4+OH \rightarrow CH_2CH_2OH', 'CH_2CH_2OH \rightarrow C_2H_4+OH', 'CH_3CHOH+M \rightarrow acetaldehyde+H+M', 'acetaldehyde+H+M \rightarrow CH_3CHOH+M', 'O_2C_2H_4OH \rightarrow CH_2CH_2OH+O_2', 'CH_2CH_2OH+O_2 \rightarrow O_2C_2H_4OH', 'O_2C_2H_4OH \rightarrow OH+CH_2O+CH_2O', 'OH+CH_2O+CH_2O \rightarrow O_2C_2H_4OH', 'CH_3CHOH+O_2 \rightarrow acetaldehyde+HO_2', 'acetaldehyde+HO_2 \rightarrow CH_3CHOH+O_2', 'acetone(+M) \rightarrow acetyl+CH_3(+M)', 'acetyl+CH_3(+M) \rightarrow acetone(+M)', 'acetone+OH \rightarrow propen2oxy+H_2O', 'propen2oxy+H_2O \rightarrow acetone+OH', 'acetone+H \rightarrow propen2oxy+H_2', 'propen2oxy+H_2 \rightarrow acetone+H', 'acetone+O \rightarrow propen2oxy+OH', 'propen2oxy+OH \rightarrow acetone+O', 'acetone+CH_3 \rightarrow propen2oxy+CH_4', 'propen2oxy+CH_4 \rightarrow acetone+CH_3', 'acetone+CH_3O \rightarrow propen2oxy+CH_3OH', 'propen2oxy+CH_3OH \rightarrow acetone+CH_3O', 'acetone+O_2 \rightarrow propen2oxy+HO_2', 'propen2oxy+HO_2 \rightarrow acetone+O_2', 'acetone+HO_2 \rightarrow propen2oxy+H_2O_2', 'propen2oxy+H_2O_2 \rightarrow acetone+HO_2', 'acetone+CH_3OO \rightarrow propen2oxy+CH_3OOH', 'propen2oxy+CH_3OOH \rightarrow acetone+CH_3OO', 'propen2oxy \rightarrow ketene+CH_3', 'ketene+CH_3 \rightarrow propen2oxy', 'iR+O \rightarrow acetone+H', 'acetone+H \rightarrow iR+O', 'acetone+H \rightarrow iRO', 'iRO \rightarrow acetone+H', 'iRO+O_2 \rightarrow acetone+HO_2', 'acetone+HO_2 \rightarrow iRO+O_2', 'C_2H_3+HCO \rightarrow acrolein', 'acrolein \rightarrow C_2H_3+HCO', 'acrolein+H \rightarrow CH_2CHCO+H_2', 'CH_2CHCO+H_2 \rightarrow acrolein+H', 'acrolein+O \rightarrow CH_2CHCO+OH', 'CH_2CHCO+OH \rightarrow acrolein+O', 'acrolein+H \rightarrow C_2H_4+HCO', 'C_2H_4+HCO \rightarrow acrolein+H', 'acrolein+O \rightarrow ketene+HCO+H', 'ketene+HCO+H \rightarrow acrolein+O', 'acrolein+OH \rightarrow CH_2CHCO+H_2O', 'CH_2CHCO+H_2O \rightarrow acrolein+OH', 'acrolein+O_2 \rightarrow CH_2CHCO+HO_2', 'CH_2CHCO+HO_2 \rightarrow acrolein+O_2', 'acrolein+HO_2 \rightarrow CH_2CHCO+H_2O_2', 'CH_2CHCO+H_2O_2 \rightarrow acrolein+HO_2', 'acrolein+CH_3 \rightarrow CH_2CHCO+CH_4', 'CH_2CHCO+CH_4 \rightarrow acrolein+CH_3', 'acrolein+C_2H_3 \rightarrow CH_2CHCO+C_2H_4', 'CH_2CHCO+C_2H_4 \rightarrow acrolein+C_2H_3', 'acrolein+CH_3O \rightarrow CH_2CHCO+CH_3OH', 'CH_2CHCO+CH_3OH \rightarrow acrolein+CH_3O', 'acrolein+CH_3OO \rightarrow CH_2CHCO+CH_3OOH', 'CH_2CHCO+CH_3OOH \rightarrow acrolein+CH_3OO', 'C_2H_3+CO \rightarrow CH_2CHCO', 'CH_2CHCO \rightarrow C_2H_3+CO', 'CH_2CHCO+O_2 \rightarrow vinoxy+CO_2', 'vinoxy+CO_2 \rightarrow CH_2CHCO+O_2', 'CH_2CHCO+O \rightarrow C_2H_3+CO_2', 'C_2H_3+CO_2 \rightarrow CH_2CHCO+O', 'C_2H_5+HCO \rightarrow propanal', 'propanal \rightarrow C_2H_5+HCO', 'propanal+H \rightarrow propionyl+H_2', 'propionyl+H_2 \rightarrow propanal+H', 'propanal+O \rightarrow propionyl+OH', 'propionyl+OH \rightarrow propanal+O', 'propanal+OH \rightarrow propionyl+H_2O', 'propionyl+H_2O \rightarrow propanal+OH', 'propanal+CH_3 \rightarrow propionyl+CH_4', 'propionyl+CH_4 \rightarrow propanal+CH_3', 'propanal+HO_2 \rightarrow propionyl+H_2O_2', 'propionyl+H_2O_2 \rightarrow propanal+HO_2', 'propanal+CH_3O \rightarrow propionyl+CH_3OH', 'propionyl+CH_3OH \rightarrow propanal+CH_3O', 'propanal+CH_3OO \rightarrow propionyl+CH_3OOH', 'propionyl+CH_3OOH \rightarrow propanal+CH_3OO', 'propanal+C_2H_5 \rightarrow propionyl+C_2H_6', 'propionyl+C_2H_6 \rightarrow propanal+C_2H_5', 'propanal+ethoxy \rightarrow propionyl+ethanol', 'propionyl+ethanol \rightarrow propanal+ethoxy', 'propanal+CH_3CH_2OO \rightarrow propionyl+CH_3CH_2OOH', 'propionyl+CH_3CH_2OOH \rightarrow propanal+CH_3CH_2OO', 'propanal+O_2 \rightarrow propionyl+HO_2', 'propionyl+HO_2 \rightarrow propanal+O_2', 'propanal+acetylperoxy \rightarrow propionyl+CH_3CO_3H', 'propionyl+CH_3CO_3H \rightarrow propanal+acetylperoxy', 'propanal+C_2H_3 \rightarrow propionyl+C_2H_4', 'propionyl+C_2H_4 \rightarrow propanal+C_2H_3', 'C_2H_5+CO \rightarrow propionyl', 'propionyl \rightarrow C_2H_5+CO', 'CH_3OCH_3(+M) \rightarrow CH_3+CH_3O(+M)', 'CH_3+CH_3O(+M) \rightarrow CH_3OCH_3(+M)', 'CH_3OCH_3+OH \rightarrow CH_3OCH_2+H_2O', 'CH_3OCH_2+H_2O \rightarrow CH_3OCH_3+OH', 'CH_3OCH_3+H \rightarrow CH_3OCH_2+H_2', 'CH_3OCH_2+H_2 \rightarrow CH_3OCH_3+H', 'CH_3OCH_3+O \rightarrow CH_3OCH_2+OH', 'CH_3OCH_2+OH \rightarrow CH_3OCH_3+O', 'CH_3OCH_3+HO_2 \rightarrow CH_3OCH_2+H_2O_2', 'CH_3OCH_2+H_2O_2 \rightarrow CH_3OCH_3+HO_2', 'CH_3OCH_3+CH_3OO \rightarrow CH_3OCH_2+CH_3OOH', 'CH_3OCH_2+CH_3OOH \rightarrow CH_3OCH_3+CH_3OO', 'CH_3OCH_3+CH_3 \rightarrow CH_3OCH_2+CH_4', 'CH_3OCH_2+CH_4 \rightarrow CH_3OCH_3+CH_3', 'CH_3OCH_3+O_2 \rightarrow CH_3OCH_2+HO_2', 'CH_3OCH_2+HO_2 \rightarrow CH_3OCH_3+O_2', 'CH_3OCH_3+CH_3O \rightarrow CH_3OCH_2+CH_3OH', 'CH_3OCH_2+CH_3OH \rightarrow CH_3OCH_3+CH_3O', 'CH_3OCH_3+formylperoxy \rightarrow CH_3OCH_2+formylooh', 'CH_3OCH_2+formylooh \rightarrow CH_3OCH_3+formylperoxy', 'CH_3OCH_2 \rightarrow CH_2O+CH_3', 'CH_2O+CH_3 \rightarrow CH_3OCH_2', 'CH_3OCH_2+CH_3O \rightarrow CH_3OCH_3+CH_2O', 'CH_3OCH_3+CH_2O \rightarrow CH_3OCH_2+CH_3O', 'CH_3OCH_2+CH_2O \rightarrow CH_3OCH_3+HCO', 'CH_3OCH_3+HCO \rightarrow CH_3OCH_2+CH_2O', 'CH_3OCH_2+acetaldehyde \rightarrow CH_3OCH_3+acetyl', 'CH_3OCH_3+acetyl \rightarrow CH_3OCH_2+acetaldehyde', 'nROO+C_3H_8 \rightarrow nROOH+nR', 'nROOH+nR \rightarrow nROO+C_3H_8', 'iROO+C_3H_8 \rightarrow iROOH+nR', 'iROOH+nR \rightarrow iROO+C_3H_8', 'nROO+C_3H_8 \rightarrow nROOH+iR', 'nROOH+iR \rightarrow nROO+C_3H_8', 'iROO+C_3H_8 \rightarrow iROOH+iR', 'iROOH+iR \rightarrow iROO+C_3H_8', 'nROOH \rightarrow nRO+OH', 'nRO+OH \rightarrow nROOH', 'iROOH \rightarrow iRO+OH', 'iRO+OH \rightarrow iROOH', 'C_3H_8(+M) \rightarrow CH_3+C_2H_5(+M)', 'CH_3+C_2H_5(+M) \rightarrow C_3H_8(+M)', 'nR+H \rightarrow C_3H_8', 'C_3H_8 \rightarrow nR+H', 'iR+H \rightarrow C_3H_8', 'C_3H_8 \rightarrow iR+H', 'C_3H_8+O_2 \rightarrow iR+HO_2', 'iR+HO_2 \rightarrow C_3H_8+O_2', 'C_3H_8+O_2 \rightarrow nR+HO_2', 'nR+HO_2 \rightarrow C_3H_8+O_2', 'H+C_3H_8 \rightarrow H_2+iR', 'H_2+iR \rightarrow H+C_3H_8', 'H+C_3H_8 \rightarrow H_2+nR', 'H_2+nR \rightarrow H+C_3H_8', 'C_3H_8+O \rightarrow iR+OH', 'iR+OH \rightarrow C_3H_8+O', 'C_3H_8+O \rightarrow nR+OH', 'nR+OH \rightarrow C_3H_8+O', 'C_3H_8+OH \rightarrow nR+H_2O', 'nR+H_2O \rightarrow C_3H_8+OH', 'C_3H_8+OH \rightarrow iR+H_2O', 'iR+H_2O \rightarrow C_3H_8+OH', 'C_3H_8+HO_2 \rightarrow iR+H_2O_2', 'iR+H_2O_2 \rightarrow C_3H_8+HO_2', 'C_3H_8+HO_2 \rightarrow nR+H_2O_2', 'nR+H_2O_2 \rightarrow C_3H_8+HO_2', 'CH_3+C_3H_8 \rightarrow CH_4+iR', 'CH_4+iR \rightarrow CH_3+C_3H_8', 'CH_3+C_3H_8 \rightarrow CH_4+nR', 'CH_4+nR \rightarrow CH_3+C_3H_8', 'iR+C_3H_8 \rightarrow nR+C_3H_8', 'nR+C_3H_8 \rightarrow iR+C_3H_8', 'C_2H_3+C_3H_8 \rightarrow C_2H_4+iR', 'C_2H_4+iR \rightarrow C_2H_3+C_3H_8', 'C_2H_3+C_3H_8 \rightarrow C_2H_4+nR', 'C_2H_4+nR \rightarrow C_2H_3+C_3H_8', 'C_2H_5+C_3H_8 \rightarrow C_2H_6+iR', 'C_2H_6+iR \rightarrow C_2H_5+C_3H_8', 'C_2H_5+C_3H_8 \rightarrow C_2H_6+nR', 'C_2H_6+nR \rightarrow C_2H_5+C_3H_8', 'C_3H_8+allyl \rightarrow nR+C_3H_6', 'nR+C_3H_6 \rightarrow C_3H_8+allyl', 'C_3H_8+allyl \rightarrow iR+C_3H_6', 'iR+C_3H_6 \rightarrow C_3H_8+allyl', 'C_3H_8+CH_3O \rightarrow nR+CH_3OH', 'nR+CH_3OH \rightarrow C_3H_8+CH_3O', 'C_3H_8+CH_3O \rightarrow iR+CH_3OH', 'iR+CH_3OH \rightarrow C_3H_8+CH_3O', 'CH_3OO+C_3H_8 \rightarrow CH_3OOH+nR', 'CH_3OOH+nR \rightarrow CH_3OO+C_3H_8', 'CH_3OO+C_3H_8 \rightarrow CH_3OOH+iR', 'CH_3OOH+iR \rightarrow CH_3OO+C_3H_8', 'CH_3CH_2OO+C_3H_8 \rightarrow CH_3CH_2OOH+nR', 'CH_3CH_2OOH+nR \rightarrow CH_3CH_2OO+C_3H_8', 'CH_3CH_2OO+C_3H_8 \rightarrow CH_3CH_2OOH+iR', 'CH_3CH_2OOH+iR \rightarrow CH_3CH_2OO+C_3H_8', 'C_3H_8+acetylperoxy \rightarrow iR+CH_3CO_3H', 'iR+CH_3CO_3H \rightarrow C_3H_8+acetylperoxy', 'C_3H_8+acetylperoxy \rightarrow nR+CH_3CO_3H', 'nR+CH_3CO_3H \rightarrow C_3H_8+acetylperoxy', 'C_3H_8+formylperoxy \rightarrow nR+formylooh', 'nR+formylooh \rightarrow C_3H_8+formylperoxy', 'C_3H_8+formylperoxy \rightarrow iR+formylooh', 'iR+formylooh \rightarrow C_3H_8+formylperoxy', 'H+C_3H_6 \rightarrow iR', 'iR \rightarrow H+C_3H_6', 'iR+H \rightarrow C_2H_5+CH_3', 'C_2H_5+CH_3 \rightarrow iR+H', 'iR+OH \rightarrow C_3H_6+H_2O', 'C_3H_6+H_2O \rightarrow iR+OH', 'iR+O \rightarrow acetaldehyde+CH_3', 'acetaldehyde+CH_3 \rightarrow iR+O', 'nR \rightarrow CH_3+C_2H_4', 'CH_3+C_2H_4 \rightarrow nR', 'nR \rightarrow H+C_3H_6', 'H+C_3H_6 \rightarrow nR', 'propanal+nR \rightarrow propionyl+C_3H_8', 'propionyl+C_3H_8 \rightarrow propanal+nR', 'propanal+iR \rightarrow propionyl+C_3H_8', 'propionyl+C_3H_8 \rightarrow propanal+iR', 'propanal+allyl \rightarrow propionyl+C_3H_6', 'propionyl+C_3H_6 \rightarrow propanal+allyl', 'allyl+H(+M) \rightarrow C_3H_6(+M)', 'C_3H_6(+M) \rightarrow allyl+H(+M)', 'C_2H_3+CH_3(+M) \rightarrow C_3H_6(+M)', 'C_3H_6(+M) \rightarrow C_2H_3+CH_3(+M)', 'C_3H_6 \rightarrow propen1yl+H', 'propen1yl+H \rightarrow C_3H_6', 'C_3H_6 \rightarrow propen2yl+H', 'propen2yl+H \rightarrow C_3H_6', 'C_3H_6+O \rightarrow C_2H_5+HCO', 'C_2H_5+HCO \rightarrow C_3H_6+O', 'C_3H_6+O \rightarrow ketene+CH_3+H', 'ketene+CH_3+H \rightarrow C_3H_6+O', 'C_3H_6+O \rightarrow CH_3CHCO+H+H', 'CH_3CHCO+H+H \rightarrow C_3H_6+O', 'C_3H_6+O \rightarrow allyl+OH', 'allyl+OH \rightarrow C_3H_6+O', 'C_3H_6+O \rightarrow propen1yl+OH', 'propen1yl+OH \rightarrow C_3H_6+O', 'C_3H_6+O \rightarrow propen2yl+OH', 'propen2yl+OH \rightarrow C_3H_6+O', 'C_3H_6+HO_2 \rightarrow allyl+H_2O_2', 'allyl+H_2O_2 \rightarrow C_3H_6+HO_2', 'C_3H_6+HO_2 \rightarrow propen1yl+H_2O_2', 'propen1yl+H_2O_2 \rightarrow C_3H_6+HO_2', 'C_3H_6+HO_2 \rightarrow propen2yl+H_2O_2', 'propen2yl+H_2O_2 \rightarrow C_3H_6+HO_2', 'C_3H_6+H \rightarrow C_2H_4+CH_3', 'C_2H_4+CH_3 \rightarrow C_3H_6+H', 'C_3H_6+H \rightarrow allyl+H_2', 'allyl+H_2 \rightarrow C_3H_6+H', 'C_3H_6+H \rightarrow propen2yl+H_2', 'propen2yl+H_2 \rightarrow C_3H_6+H', 'C_3H_6+H \rightarrow propen1yl+H_2', 'propen1yl+H_2 \rightarrow C_3H_6+H', 'C_3H_6+O_2 \rightarrow allyl+HO_2', 'allyl+HO_2 \rightarrow C_3H_6+O_2', 'C_3H_6+O_2 \rightarrow propen1yl+HO_2', 'propen1yl+HO_2 \rightarrow C_3H_6+O_2', 'C_3H_6+O_2 \rightarrow propen2yl+HO_2', 'propen2yl+HO_2 \rightarrow C_3H_6+O_2', 'C_3H_6+CH_3 \rightarrow allyl+CH_4', 'allyl+CH_4 \rightarrow C_3H_6+CH_3', 'C_3H_6+CH_3 \rightarrow propen1yl+CH_4', 'propen1yl+CH_4 \rightarrow C_3H_6+CH_3', 'C_3H_6+CH_3 \rightarrow propen2yl+CH_4', 'propen2yl+CH_4 \rightarrow C_3H_6+CH_3', 'C_3H_6+C_2H_5 \rightarrow allyl+C_2H_6', 'allyl+C_2H_6 \rightarrow C_3H_6+C_2H_5', 'C_3H_6+acetylperoxy \rightarrow allyl+CH_3CO_3H', 'allyl+CH_3CO_3H \rightarrow C_3H_6+acetylperoxy', 'C_3H_6+CH_3OO \rightarrow allyl+CH_3OOH', 'allyl+CH_3OOH \rightarrow C_3H_6+CH_3OO', 'C_3H_6+HO_2 \rightarrow propen1ol+OH', 'propen1ol+OH \rightarrow C_3H_6+HO_2', 'C_3H_6+CH_3CH_2OO \rightarrow allyl+CH_3CH_2OOH', 'allyl+CH_3CH_2OOH \rightarrow C_3H_6+CH_3CH_2OO', 'C_3H_6+nROO \rightarrow allyl+nROOH', 'allyl+nROOH \rightarrow C_3H_6+nROO', 'C_3H_6+iROO \rightarrow allyl+iROOH', 'allyl+iROOH \rightarrow C_3H_6+iROO', 'allyl+O \rightarrow acrolein+H', 'acrolein+H \rightarrow allyl+O', 'allyl+OH \rightarrow acrolein+H+H', 'acrolein+H+H \rightarrow allyl+OH', 'allyl+O_2 \rightarrow acetyl+CH_2O', 'acetyl+CH_2O \rightarrow allyl+O_2', 'allyl+O_2 \rightarrow acrolein+OH', 'acrolein+OH \rightarrow allyl+O_2', 'allyl+HCO \rightarrow C_3H_6+CO', 'C_3H_6+CO \rightarrow allyl+HCO', 'allyl \rightarrow propen2yl', 'propen2yl \rightarrow allyl', 'allyl \rightarrow propen1yl', 'propen1yl \rightarrow allyl', 'C_2H_2+CH_3 \rightarrow allyl', 'allyl \rightarrow C_2H_2+CH_3', 'allyl+CH_3OO \rightarrow allyloxy+CH_3O', 'allyloxy+CH_3O \rightarrow allyl+CH_3OO', 'allyl+C_2H_5 \rightarrow C_2H_4+C_3H_6', 'C_2H_4+C_3H_6 \rightarrow allyl+C_2H_5', 'C_2H_2+CH_3 \rightarrow propen1yl', 'propen1yl \rightarrow C_2H_2+CH_3', 'propen1yl+O \rightarrow C_2H_4+HCO', 'C_2H_4+HCO \rightarrow propen1yl+O', 'propen1yl+OH \rightarrow C_2H_4+HCO+H', 'C_2H_4+HCO+H \rightarrow propen1yl+OH', 'propen1yl+O_2 \rightarrow acetaldehyde+HCO', 'acetaldehyde+HCO \rightarrow propen1yl+O_2', 'propen1yl+HO_2 \rightarrow C_2H_4+HCO+OH', 'C_2H_4+HCO+OH \rightarrow propen1yl+HO_2', 'propen1yl+HCO \rightarrow C_3H_6+CO', 'C_3H_6+CO \rightarrow propen1yl+HCO', 'C_2H_2+CH_3 \rightarrow propen2yl', 'propen2yl \rightarrow C_2H_2+CH_3', 'propen2yl \rightarrow propen1yl', 'propen1yl \rightarrow propen2yl', 'propen2yl+O \rightarrow CH_3+ketene', 'CH_3+ketene \rightarrow propen2yl+O', 'propen2yl+OH \rightarrow CH_3+ketene+H', 'CH_3+ketene+H \rightarrow propen2yl+OH', 'propen2yl+O_2 \rightarrow acetyl+CH_2O', 'acetyl+CH_2O \rightarrow propen2yl+O_2', 'propen2yl+HO_2 \rightarrow CH_3+ketene+OH', 'CH_3+ketene+OH \rightarrow propen2yl+HO_2', 'propen2yl+HCO \rightarrow C_3H_6+CO', 'C_3H_6+CO \rightarrow propen2yl+HCO', 'nR+HO_2 \rightarrow nRO+OH', 'nRO+OH \rightarrow nR+HO_2', 'iR+HO_2 \rightarrow iRO+OH', 'iRO+OH \rightarrow iR+HO_2', 'CH_3OO+nR \rightarrow CH_3O+nRO', 'CH_3O+nRO \rightarrow CH_3OO+nR', 'CH_3OO+iR \rightarrow CH_3O+iRO', 'CH_3O+iRO \rightarrow CH_3OO+iR', 'nROO+CH_2O \rightarrow nROOH+HCO', 'nROOH+HCO \rightarrow nROO+CH_2O', 'nROO+acetaldehyde \rightarrow nROOH+acetyl', 'nROOH+acetyl \rightarrow nROO+acetaldehyde', 'iROO+CH_2O \rightarrow iROOH+HCO', 'iROOH+HCO \rightarrow iROO+CH_2O', 'iROO+acetaldehyde \rightarrow iROOH+acetyl', 'iROOH+acetyl \rightarrow iROO+acetaldehyde', 'nROO+HO_2 \rightarrow nROOH+O_2', 'nROOH+O_2 \rightarrow nROO+HO_2', 'iROO+HO_2 \rightarrow iROOH+O_2', 'iROOH+O_2 \rightarrow iROO+HO_2', 'C_2H_4+nROO \rightarrow C_2H_3+nROOH', 'C_2H_3+nROOH \rightarrow C_2H_4+nROO', 'C_2H_4+iROO \rightarrow C_2H_3+iROOH', 'C_2H_3+iROOH \rightarrow C_2H_4+iROO', 'CH_3OH+nROO \rightarrow CH_2OH+nROOH', 'CH_2OH+nROOH \rightarrow CH_3OH+nROO', 'CH_3OH+iROO \rightarrow CH_2OH+iROOH', 'CH_2OH+iROOH \rightarrow CH_3OH+iROO', 'acrolein+nROO \rightarrow CH_2CHCO+nROOH', 'CH_2CHCO+nROOH \rightarrow acrolein+nROO', 'acrolein+iROO \rightarrow CH_2CHCO+iROOH', 'CH_2CHCO+iROOH \rightarrow acrolein+iROO', 'CH_4+nROO \rightarrow CH_3+nROOH', 'CH_3+nROOH \rightarrow CH_4+nROO', 'CH_4+iROO \rightarrow CH_3+iROOH', 'CH_3+iROOH \rightarrow CH_4+iROO', 'nROO+CH_3OO \rightarrow nRO+CH_3O+O_2', 'nRO+CH_3O+O_2 \rightarrow nROO+CH_3OO', 'iROO+CH_3OO \rightarrow iRO+CH_3O+O_2', 'iRO+CH_3O+O_2 \rightarrow iROO+CH_3OO', 'H_2+nROO \rightarrow H+nROOH', 'H+nROOH \rightarrow H_2+nROO', 'H_2+iROO \rightarrow H+iROOH', 'H+iROOH \rightarrow H_2+iROO', 'iROO+C_2H_6 \rightarrow iROOH+C_2H_5', 'iROOH+C_2H_5 \rightarrow iROO+C_2H_6', 'nROO+C_2H_6 \rightarrow nROOH+C_2H_5', 'nROOH+C_2H_5 \rightarrow nROO+C_2H_6', 'iROO+propanal \rightarrow iROOH+propionyl', 'iROOH+propionyl \rightarrow iROO+propanal', 'nROO+propanal \rightarrow nROOH+propionyl', 'nROOH+propionyl \rightarrow nROO+propanal', 'iROO+acetylperoxy \rightarrow iRO+acetyloxy+O_2', 'iRO+acetyloxy+O_2 \rightarrow iROO+acetylperoxy', 'nROO+acetylperoxy \rightarrow nRO+acetyloxy+O_2', 'nRO+acetyloxy+O_2 \rightarrow nROO+acetylperoxy', 'iROO+CH_3CH_2OO \rightarrow iRO+ethoxy+O_2', 'iRO+ethoxy+O_2 \rightarrow iROO+CH_3CH_2OO', 'nROO+CH_3CH_2OO \rightarrow nRO+ethoxy+O_2', 'nRO+ethoxy+O_2 \rightarrow nROO+CH_3CH_2OO', 'iROO+iROO \rightarrow O_2+iRO+iRO', 'O_2+iRO+iRO \rightarrow iROO+iROO', 'nROO+nROO \rightarrow O_2+nRO+nRO', 'O_2+nRO+nRO \rightarrow nROO+nROO', 'iROO+nROO \rightarrow iRO+nRO+O_2', 'iRO+nRO+O_2 \rightarrow iROO+nROO', 'iROO+CH_3 \rightarrow iRO+CH_3O', 'iRO+CH_3O \rightarrow iROO+CH_3', 'iROO+C_2H_5 \rightarrow iRO+ethoxy', 'iRO+ethoxy \rightarrow iROO+C_2H_5', 'iROO+iR \rightarrow iRO+iRO', 'iRO+iRO \rightarrow iROO+iR', 'iROO+nR \rightarrow iRO+nRO', 'iRO+nRO \rightarrow iROO+nR', 'iROO+allyl \rightarrow iRO+allyloxy', 'iRO+allyloxy \rightarrow iROO+allyl', 'nROO+CH_3 \rightarrow nRO+CH_3O', 'nRO+CH_3O \rightarrow nROO+CH_3', 'nROO+C_2H_5 \rightarrow nRO+ethoxy', 'nRO+ethoxy \rightarrow nROO+C_2H_5', 'nROO+iR \rightarrow nRO+iRO', 'nRO+iRO \rightarrow nROO+iR', 'nROO+nR \rightarrow nRO+nRO', 'nRO+nRO \rightarrow nROO+nR', 'nROO+allyl \rightarrow nRO+allyloxy', 'nRO+allyloxy \rightarrow nROO+allyl', 'C_2H_5+CH_2O \rightarrow nRO', 'nRO \rightarrow C_2H_5+CH_2O', 'propanal+H \rightarrow nRO', 'nRO \rightarrow propanal+H', 'CH_3+acetaldehyde \rightarrow iRO', 'iRO \rightarrow CH_3+acetaldehyde', 'allyloxy+O_2 \rightarrow acrolein+HO_2', 'acrolein+HO_2 \rightarrow allyloxy+O_2', 'propen1ol \rightarrow C_2H_4+CH_2O', 'C_2H_4+CH_2O \rightarrow propen1ol', 'propen1ol+OH \rightarrow CH_2O+C_2H_3+H_2O', 'CH_2O+C_2H_3+H_2O \rightarrow propen1ol+OH', 'propen1ol+H \rightarrow CH_2O+C_2H_3+H_2', 'CH_2O+C_2H_3+H_2 \rightarrow propen1ol+H', 'propen1ol+O \rightarrow CH_2O+C_2H_3+OH', 'CH_2O+C_2H_3+OH \rightarrow propen1ol+O', 'propen1ol+HO_2 \rightarrow CH_2O+C_2H_3+H_2O_2', 'CH_2O+C_2H_3+H_2O_2 \rightarrow propen1ol+HO_2', 'propen1ol+CH_3OO \rightarrow CH_2O+C_2H_3+CH_3OOH', 'CH_2O+C_2H_3+CH_3OOH \rightarrow propen1ol+CH_3OO', 'propen1ol+CH_3 \rightarrow CH_2O+C_2H_3+CH_4', 'CH_2O+C_2H_3+CH_4 \rightarrow propen1ol+CH_3', 'C_3H_6+OH \rightarrow allyl+H_2O', 'allyl+H_2O \rightarrow C_3H_6+OH', 'C_3H_6+OH \rightarrow propen1yl+H_2O', 'propen1yl+H_2O \rightarrow C_3H_6+OH', 'C_3H_6+OH \rightarrow propen2yl+H_2O', 'propen2yl+H_2O \rightarrow C_3H_6+OH', 'C_3H_6+OH \rightarrow allyl-alcohol+H', 'allyl-alcohol+H \rightarrow C_3H_6+OH', 'C_3H_6+OH \rightarrow ethenol+CH_3', 'ethenol+CH_3 \rightarrow C_3H_6+OH', 'C_3H_6+OH \rightarrow propen1ol+H', 'propen1ol+H \rightarrow C_3H_6+OH', 'C_3H_6+OH \rightarrow propen2ol+H', 'propen2ol+H \rightarrow C_3H_6+OH', 'C_3H_6+OH \rightarrow acetaldehyde+CH_3', 'acetaldehyde+CH_3 \rightarrow C_3H_6+OH', 'allyl+HO_2 \rightarrow prod_2', 'prod_2 \rightarrow allyl+HO_2', 'allyl+HO_2 \rightarrow acrolein+H_2O', 'acrolein+H_2O \rightarrow allyl+HO_2', 'allyl+HO_2 \rightarrow allyloxy+OH', 'allyloxy+OH \rightarrow allyl+HO_2', 'prod_2 \rightarrow acrolein+H_2O', 'acrolein+H_2O \rightarrow prod_2', 'prod_2 \rightarrow allyloxy+OH', 'allyloxy+OH \rightarrow prod_2', 'allyloxy \rightarrow C_2H_3+CH_2O', 'C_2H_3+CH_2O \rightarrow allyloxy', 'allyloxy \rightarrow vinoxylmethyl', 'vinoxylmethyl \rightarrow allyloxy', 'allyloxy \rightarrow formylethyl', 'formylethyl \rightarrow allyloxy', 'allyloxy \rightarrow acrolein+H', 'acrolein+H \rightarrow allyloxy', 'allyloxy \rightarrow C_2H_4+HCO', 'C_2H_4+HCO \rightarrow allyloxy', 'vinoxylmethyl \rightarrow C_2H_3+CH_2O', 'C_2H_3+CH_2O \rightarrow vinoxylmethyl', 'vinoxylmethyl \rightarrow formylethyl', 'formylethyl \rightarrow vinoxylmethyl', 'vinoxylmethyl \rightarrow acrolein+H', 'acrolein+H \rightarrow vinoxylmethyl', 'vinoxylmethyl \rightarrow C_2H_4+HCO', 'C_2H_4+HCO \rightarrow vinoxylmethyl', 'formylethyl \rightarrow C_2H_3+CH_2O', 'C_2H_3+CH_2O \rightarrow formylethyl', 'formylethyl \rightarrow acrolein+H', 'acrolein+H \rightarrow formylethyl', 'formylethyl \rightarrow C_2H_4+HCO', 'C_2H_4+HCO \rightarrow formylethyl', 'C_2H_3+CH_2O \rightarrow acrolein+H', 'acrolein+H \rightarrow C_2H_3+CH_2O', 'C_2H_3+CH_2O \rightarrow C_2H_4+HCO', 'C_2H_4+HCO \rightarrow C_2H_3+CH_2O', 'O_2+nR \rightarrow nROO', 'nROO \rightarrow O_2+nR', 'O_2+nR \rightarrow QOOH_2', 'QOOH_2 \rightarrow O_2+nR', 'O_2+nR \rightarrow QOOH_1', 'QOOH_1 \rightarrow O_2+nR', 'O_2+nR \rightarrow HO_2+C_3H_6', 'HO_2+C_3H_6 \rightarrow O_2+nR', 'O_2+nR \rightarrow OH+propoxide', 'OH+propoxide \rightarrow O_2+nR', 'nROO \rightarrow QOOH_2', 'QOOH_2 \rightarrow nROO', 'nROO \rightarrow QOOH_1', 'QOOH_1 \rightarrow nROO', 'nROO \rightarrow HO_2+C_3H_6', 'HO_2+C_3H_6 \rightarrow nROO', 'nROO \rightarrow OH+propoxide', 'OH+propoxide \rightarrow nROO', 'QOOH_2 \rightarrow QOOH_1', 'QOOH_1 \rightarrow QOOH_2', 'QOOH_2 \rightarrow HO_2+C_3H_6', 'HO_2+C_3H_6 \rightarrow QOOH_2', 'QOOH_2 \rightarrow OH+propoxide', 'OH+propoxide \rightarrow QOOH_2', 'QOOH_1 \rightarrow HO_2+C_3H_6', 'HO_2+C_3H_6 \rightarrow QOOH_1', 'QOOH_1 \rightarrow OH+propoxide', 'OH+propoxide \rightarrow QOOH_1', 'O_2+iR \rightarrow iROO', 'iROO \rightarrow O_2+iR', 'O_2+iR \rightarrow QOOH_3', 'QOOH_3 \rightarrow O_2+iR', 'O_2+iR \rightarrow HO_2+C_3H_6', 'HO_2+C_3H_6 \rightarrow O_2+iR', 'O_2+iR \rightarrow OH+propoxide', 'OH+propoxide \rightarrow O_2+iR', 'iROO \rightarrow QOOH_3', 'QOOH_3 \rightarrow iROO', 'iROO \rightarrow HO_2+C_3H_6', 'HO_2+C_3H_6 \rightarrow iROO', 'iROO \rightarrow OH+propoxide', 'OH+propoxide \rightarrow iROO', 'QOOH_3 \rightarrow HO_2+C_3H_6', 'HO_2+C_3H_6 \rightarrow QOOH_3', 'QOOH_3 \rightarrow OH+propoxide', 'OH+propoxide \rightarrow QOOH_3', 'HO_2+C_3H_6 \rightarrow OH+propoxide', 'OH+propoxide \rightarrow HO_2+C_3H_6', 'O_2+QOOH_1 \rightarrow O_2QOOH_1', 'O_2QOOH_1 \rightarrow O_2+QOOH_1', 'O_2+QOOH_1 \rightarrow OH+OH+OQ^\primeO_1', 'OH+OH+OQ^\primeO_1 \rightarrow O_2+QOOH_1', 'O_2+QOOH_1 \rightarrow HO_2+prod_2', 'HO_2+prod_2 \rightarrow O_2+QOOH_1', 'O_2+QOOH_1 \rightarrow OH+prod_3', 'OH+prod_3 \rightarrow O_2+QOOH_1', 'O_2+QOOH_2 \rightarrow well_2', 'well_2 \rightarrow O_2+QOOH_2', 'O_2+QOOH_2 \rightarrow well_3', 'well_3 \rightarrow O_2+QOOH_2', 'O_2+QOOH_2 \rightarrow well_5', 'well_5 \rightarrow O_2+QOOH_2', 'O_2+QOOH_2 \rightarrow OH+prod_3', 'OH+prod_3 \rightarrow O_2+QOOH_2', 'O_2+QOOH_2 \rightarrow OH+OH+frag_4', 'OH+OH+frag_4 \rightarrow O_2+QOOH_2', 'O_2+QOOH_2 \rightarrow OH+OH+frag_5', 'OH+OH+frag_5 \rightarrow O_2+QOOH_2', 'O_2+QOOH_2 \rightarrow HO_2+prod_2', 'HO_2+prod_2 \rightarrow O_2+QOOH_2', 'O_2+QOOH_2 \rightarrow HO_2+prod_6', 'HO_2+prod_6 \rightarrow O_2+QOOH_2', 'O_2+QOOH_2 \rightarrow HO_2+prod_7', 'HO_2+prod_7 \rightarrow O_2+QOOH_2', 'O_2+QOOH_2 \rightarrow O_2+QOOH_3', 'O_2+QOOH_3 \rightarrow O_2+QOOH_2', 'O_2+QOOH_3 \rightarrow well_2', 'well_2 \rightarrow O_2+QOOH_3', 'O_2+QOOH_3 \rightarrow well_3', 'well_3 \rightarrow O_2+QOOH_3', 'O_2+QOOH_3 \rightarrow well_5', 'well_5 \rightarrow O_2+QOOH_3', 'O_2+QOOH_3 \rightarrow OH+prod_3', 'OH+prod_3 \rightarrow O_2+QOOH_3', 'O_2+QOOH_3 \rightarrow OH+OH+frag_4', 'OH+OH+frag_4 \rightarrow O_2+QOOH_3', 'O_2+QOOH_3 \rightarrow OH+OH+frag_5', 'OH+OH+frag_5 \rightarrow O_2+QOOH_3', 'O_2+QOOH_3 \rightarrow HO_2+prod_2', 'HO_2+prod_2 \rightarrow O_2+QOOH_3', 'O_2+QOOH_3 \rightarrow HO_2+prod_6', 'HO_2+prod_6 \rightarrow O_2+QOOH_3', 'O_2+QOOH_3 \rightarrow HO_2+prod_7', 'HO_2+prod_7 \rightarrow O_2+QOOH_3', 'O_2QOOH_1 \rightarrow OH+OQ^\primeOOH_1', 'OH+OQ^\primeOOH_1 \rightarrow O_2QOOH_1', 'O_2QOOH_1 \rightarrow HO_2+prod_2', 'HO_2+prod_2 \rightarrow O_2QOOH_1', 'O_2QOOH_1 \rightarrow OH+prod_3', 'OH+prod_3 \rightarrow O_2QOOH_1', 'well_2 \rightarrow well_3', 'well_3 \rightarrow well_2', 'well_2 \rightarrow well_5', 'well_5 \rightarrow well_2', 'well_2 \rightarrow OH+prod_3', 'OH+prod_3 \rightarrow well_2', 'well_2 \rightarrow OH+prod_4', 'OH+prod_4 \rightarrow well_2', 'well_2 \rightarrow OH+prod_5', 'OH+prod_5 \rightarrow well_2', 'well_2 \rightarrow HO_2+prod_2', 'HO_2+prod_2 \rightarrow well_2', 'well_2 \rightarrow HO_2+prod_6', 'HO_2+prod_6 \rightarrow well_2', 'well_2 \rightarrow HO_2+prod_7', 'HO_2+prod_7 \rightarrow well_2', 'well_3 \rightarrow well_5', 'well_5 \rightarrow well_3', 'well_3 \rightarrow OH+prod_3', 'OH+prod_3 \rightarrow well_3', 'well_3 \rightarrow OH+prod_4', 'OH+prod_4 \rightarrow well_3', 'well_3 \rightarrow OH+prod_5', 'OH+prod_5 \rightarrow well_3', 'well_3 \rightarrow HO_2+prod_2', 'HO_2+prod_2 \rightarrow well_3', 'well_3 \rightarrow HO_2+prod_6', 'HO_2+prod_6 \rightarrow well_3', 'well_3 \rightarrow HO_2+prod_7', 'HO_2+prod_7 \rightarrow well_3', 'well_5 \rightarrow OH+prod_3', 'OH+prod_3 \rightarrow well_5', 'well_5 \rightarrow OH+prod_4', 'OH+prod_4 \rightarrow well_5', 'well_5 \rightarrow OH+prod_5', 'OH+prod_5 \rightarrow well_5', 'well_5 \rightarrow HO_2+prod_2', 'HO_2+prod_2 \rightarrow well_5', 'well_5 \rightarrow HO_2+prod_6', 'HO_2+prod_6 \rightarrow well_5', 'well_5 \rightarrow HO_2+prod_7', 'HO_2+prod_7 \rightarrow well_5', 'prod_6 \rightarrow propen1oxy+OH', 'propen1oxy+OH \rightarrow prod_6', 'prod_7 \rightarrow propen2oxy+OH', 'propen2oxy+OH \rightarrow prod_7', 'OQ^\primeOOH_1 \rightarrow OQ^\primeO_1+OH', 'OQ^\primeO_1+OH \rightarrow OQ^\primeOOH_1', 'prod_3 \rightarrow frag_3+OH', 'frag_3+OH \rightarrow prod_3', 'prod_4 \rightarrow frag_4+OH', 'frag_4+OH \rightarrow prod_4', 'prod_5 \rightarrow frag_5+OH', 'frag_5+OH \rightarrow prod_5', 'OQ^\primeO_1 \rightarrow vinoxy+CH_2O', 'vinoxy+CH_2O \rightarrow OQ^\primeO_1', 'frag_4 \rightarrow acetyl+CH_2O', 'acetyl+CH_2O \rightarrow frag_4', 'frag_5 \rightarrow CH_3+glyoxal', 'CH_3+glyoxal \rightarrow frag_5', 'frag_5 \rightarrow HCO+acetaldehyde', 'HCO+acetaldehyde \rightarrow frag_5'};

%% for plot
% markers = {'+' , 'o' , '*' , 'x' , 'square' , 'diamond' , 'v' , '^' , '>' , '<' , 'pentagram' , 'hexagram' , '.'};
markers = {'None'};

%% Current file directory
file_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', '..', '..', '..', '..', 'SOHR_DATA');
pic_dir = fullfile(fileparts(mfilename('fullpath')));

%% import time
fn_time = fullfile(file_dir, 'output', 'time_dlsode_M.csv');
delimiter = '';
formatSpec = '%f%[^\n\r]';
%% Open the text file.
fileID = fopen(fn_time,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);
time_vec = dataArray{:, 1};
%% Clear temporary variables
clearvars fn_time delimiter formatSpec fileID dataArray ans;

%% import rate_const
fn_k = fullfile(file_dir, 'output', 'drc_dlsode_M.csv');
delimiter = ',';
formatStr = '';
for i=1:N_S
    formatStr = strcat(formatStr, '%f');
end
formatStr = strcat(formatStr, '%[^\n\r]');
formatSpec = char(formatStr);
%% Open the text file.
fileID = fopen(fn_k,'r');
%% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
%% Close the text file.
fclose(fileID);
%% Create output variable
k_mat = [dataArray{1:end-1}];
%% Clear temporary variables
clearvars fn_k delimiter formatSpec fileID dataArray ans;

k_mat = 1./k_mat;

% sort by the reaction rates around 0.5 tau, idx == 3550 for example
sort_axis = round(0.42 * length(time_vec));

%% plot
fig = figure();
% https://www.mathworks.com/help/matlab/graphics_transition/why-are-plot-lines-different-colors.html
% https://www.mathworks.com/help/matlab/creating_plots/customize-graph-with-two-y-axes.html
% for hidden color order concerning axis colors
co = [    
%     0    0.4470    0.7410 % 1th plot
    1   0   0 % bl
    ]; 
set(fig,'defaultAxesColorOrder',co)

xpos = [0.15 0.95];
ypos = [0.07, 0.375, 0.375, 0.685, 0.685, 0.99];

%##########################################################################
% Panel 1
%##########################################################################
iax = 1; % Or whichever
x0=xpos(1); y0=ypos((length(ypos)/2 - iax)*2 + 1); spanx=xpos(2) - xpos(1); spany=ypos((length(ypos)/2 - iax)*2 + 1 + 1) - ypos((length(ypos)/2 - iax)*2 + 1);
%% [left bottom width height]
pos = [x0 y0 spanx spany];
subplot('Position',pos);

% target array is the species index we want to study, sort this array, not
% the whole array
% target_array = linspace(10, N_S, N_S+1-10);
target_array = [18, 60, 15, 45, 95, 87, 81]; ylim_range = [10^-6, 10^20];
% target_array = [95, 79, 91, 61, 88, 47, 102]; ylim_range = [10^-12, 10^-1];
% target_array = [12, 14, 13, 11, 4, 9]; ylim_range = [10^-8, 10^13];
[~,I] = sort(k_mat(sort_axis, target_array),'descend');
semilogy([time_vec(sort_axis), time_vec(sort_axis)], ylim_range, ...
    '--k', 'HandleVisibility','off');
hold on;

% number of lines we will plot
topN_array = linspace(1, length(target_array), length(target_array));
N_plot = length(topN_array);
colors = lines(N_plot);

% graph handler
color_idx = 1;
marker_idx = 1;
delta_n = 800;
legend_name = cell(N_plot,1);
for idx=1:N_plot
    spe_idx = target_array(I(topN_array(idx)));
    legend_name{idx, 1} = spe_name_latex{1, spe_idx};
    
    semilogy(time_vec, k_mat(:, spe_idx), ...
        'LineWidth', 2.0, ...
        'color', colors(mod(color_idx-1, length(colors))+ 1, :), ...
        'HandleVisibility','off'); hold on;
    scatter(time_vec(1:delta_n:end), k_mat(1:delta_n:end, spe_idx), ...
    'LineWidth', 2.0, ...
    'MarkerEdgeColor', colors(mod(color_idx-1, length(colors))+ 1, :), ...
    'marker', markers{1, mod(marker_idx-1, length(markers))+ 1}, ...
    'HandleVisibility','off'); hold on;
    plot(nan, nan, 'LineWidth', 2.0, 'LineStyle', '-', ...
    'color', colors(mod(color_idx-1, length(colors))+ 1, :), ...
    'marker', markers{1, mod(marker_idx-1, length(markers))+ 1});
    
    hold on;
    color_idx = color_idx + 1;
    marker_idx = marker_idx + 1;
end

%% settings
set(gca,'GridLineStyle','--');
xlabel('Time (seconds)', 'FontSize', 20);
set(gca,'Xticklabel',[]);
ylabel('Lifetime (seconds)', 'FontSize', 20);
ylim(ylim_range);

%% figure settings
grid on;
xlim([0, tau*end_t]);
leg_h = legend(legend_name);
set(leg_h, 'FontSize', 16, 'Box', 'off');
% set(leg_h, 'Location', 'West')

%##########################################################################
% Panel 2
%##########################################################################
iax = 2; % Or whichever
x0=xpos(1); y0=ypos((length(ypos)/2 - iax)*2 + 1); spanx=xpos(2) - xpos(1); spany=ypos((length(ypos)/2 - iax)*2 + 1 + 1) - ypos((length(ypos)/2 - iax)*2 + 1);
%% [left bottom width height]
pos = [x0 y0 spanx spany];
subplot('Position',pos);

% target array is the species index we want to study, sort this array, not
% the whole array
% target_array = linspace(10, N_S, N_S+1-10);
% target_array = [18, 60, 15, 45, 95, 87, 81]; ylim_range = [10^-6, 10^20];
target_array = [95, 79, 91, 61, 88, 47, 102]; ylim_range = [10^-12, 10^-1];
% target_array = [12, 14, 13, 11, 4, 9]; ylim_range = [10^-8, 10^13];
[~,I] = sort(k_mat(sort_axis, target_array),'descend');
semilogy([time_vec(sort_axis), time_vec(sort_axis)], ylim_range, ...
    '--k', 'HandleVisibility','off');
hold on;

% number of lines we will plot
topN_array = linspace(1, length(target_array), length(target_array));
N_plot = length(topN_array);
colors = lines(N_plot);

% graph handler
color_idx = 1;
marker_idx = 1;
delta_n = 800;
legend_name = cell(N_plot,1);
for idx=1:N_plot
    spe_idx = target_array(I(topN_array(idx)));
    legend_name{idx, 1} = spe_name_latex{1, spe_idx};
    
    semilogy(time_vec, k_mat(:, spe_idx), ...
        'LineWidth', 2.0, ...
        'color', colors(mod(color_idx-1, length(colors))+ 1, :), ...
        'HandleVisibility','off'); hold on;
    scatter(time_vec(1:delta_n:end), k_mat(1:delta_n:end, spe_idx), ...
    'LineWidth', 2.0, ...
    'MarkerEdgeColor', colors(mod(color_idx-1, length(colors))+ 1, :), ...
    'marker', markers{1, mod(marker_idx-1, length(markers))+ 1}, ...
    'HandleVisibility','off'); hold on;
    plot(nan, nan, 'LineWidth', 2.0, 'LineStyle', '-', ...
    'color', colors(mod(color_idx-1, length(colors))+ 1, :), ...
    'marker', markers{1, mod(marker_idx-1, length(markers))+ 1});
    
    hold on;
    color_idx = color_idx + 1;
    marker_idx = marker_idx + 1;
end

%% settings
ylim(ylim_range);
set(gca,'GridLineStyle','--');

xlabel('Time (seconds)', 'FontSize', 20);
ylabel('Lifetime (seconds)', 'FontSize', 20);
set(gca,'Xticklabel',[]);

%% figure settings
grid on;
xlim([0, tau*end_t]);
leg_h = legend(legend_name);
set(leg_h, 'FontSize', 16, 'Box', 'off');
% set(leg_h, 'Location', 'West')

%##########################################################################
% Panel 3
%##########################################################################
iax = 3; % Or whichever
x0=xpos(1); y0=ypos((length(ypos)/2 - iax)*2 + 1); spanx=xpos(2) - xpos(1); spany=ypos((length(ypos)/2 - iax)*2 + 1 + 1) - ypos((length(ypos)/2 - iax)*2 + 1);
%% [left bottom width height]
pos = [x0 y0 spanx spany];
subplot('Position',pos);

% target array is the species index we want to study, sort this array, not
% the whole array
% target_array = linspace(10, N_S, N_S+1-10);
% target_array = [18, 60, 15, 45, 95, 87, 81]; ylim_range = [10^-6, 10^20];
% target_array = [95, 79, 91, 61, 88, 47, 102]; ylim_range = [10^-12, 10^-1];
target_array = [12, 14, 13, 11, 4, 9]; ylim_range = [10^-8, 10^13];
[~,I] = sort(k_mat(sort_axis, target_array),'descend');
semilogy([time_vec(sort_axis), time_vec(sort_axis)], ylim_range, ...
    '--k', 'HandleVisibility','off');
hold on;

% number of lines we will plot
topN_array = linspace(1, length(target_array), length(target_array));
N_plot = length(topN_array);
colors = lines(N_plot);

% graph handler
color_idx = 1;
marker_idx = 1;
delta_n = 800;
legend_name = cell(N_plot,1);
for idx=1:N_plot
    spe_idx = target_array(I(topN_array(idx)));
    legend_name{idx, 1} = spe_name_latex{1, spe_idx};
    
    semilogy(time_vec, k_mat(:, spe_idx), ...
        'LineWidth', 2.0, ...
        'color', colors(mod(color_idx-1, length(colors))+ 1, :), ...
        'HandleVisibility','off'); hold on;
    scatter(time_vec(1:delta_n:end), k_mat(1:delta_n:end, spe_idx), ...
    'LineWidth', 2.0, ...
    'MarkerEdgeColor', colors(mod(color_idx-1, length(colors))+ 1, :), ...
    'marker', markers{1, mod(marker_idx-1, length(markers))+ 1}, ...
    'HandleVisibility','off'); hold on;
    plot(nan, nan, 'LineWidth', 2.0, 'LineStyle', '-', ...
    'color', colors(mod(color_idx-1, length(colors))+ 1, :), ...
    'marker', markers{1, mod(marker_idx-1, length(markers))+ 1});
    
    hold on;
    color_idx = color_idx + 1;
    marker_idx = marker_idx + 1;
end

%% settings
ylim(ylim_range);
set(gca,'GridLineStyle','--');

xlabel('Time (seconds)', 'FontSize', 20);
ylabel('Lifetime (seconds)', 'FontSize', 20);

%% figure settings
grid on;
xlim([0, tau*end_t]);
leg_h = legend(legend_name);
set(leg_h, 'FontSize', 16, 'Box', 'off');
% set(leg_h, 'Location', 'West')

%% figure size
x0=10;
y0=10;
width=400;
height=850;
set(gcf,'units','points','position',[x0,y0,width,height]);

%% save to file
figname = strcat(fig_prefix, '_3_in_1_v1.png');
print(fig, fullfile(pic_dir, figname), '-r200', '-dpng');