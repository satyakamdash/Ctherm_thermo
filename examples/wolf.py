import numpy as np
from component_contribution.thermodynamic_constants import default_RT
from component_contribution.component_contribution import ComponentContribution
from component_contribution.kegg_model import KeggModel

REACTION_FNAME = 'wolf_reactions.txt'

cc = ComponentContribution.init()
reaction_strings = open(REACTION_FNAME, 'r').readlines()
model = KeggModel.from_formulas(reaction_strings, raise_exception=True)

model.add_thermo(cc)
dG0_prime, dG0_std, sqrt_Sigma = model.get_transformed_dG0(7.5, 0.2, 298.15)

mM_conc = 1e-3 * np.matrix(np.ones((len(model.cids), 1)))
if 'C00001' in model.cids:
    mM_conc[model.cids.index('C00001'), 0] = 1.0
dGm_prime = dG0_prime + default_RT * model.S.T * np.log(mM_conc)

for i, r in enumerate(reaction_strings):
    print '-'*50
    print r.strip()
    print "dG0  = %8.1f +- %5.1f" % (model.dG0[i, 0], dG0_std[i, 0] * 1.96)
    print "dG'0 = %8.1f +- %5.1f" % (dG0_prime[i, 0], dG0_std[i, 0] * 1.96)
    print "dG'm = %8.1f +- %5.1f" % (dGm_prime[i, 0], dG0_std[i, 0] * 1.96)