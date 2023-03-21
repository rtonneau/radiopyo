
import pathlib

import numpy as np

import radiopyo_vec as rp
from radiopyo_vec.physics.beam import ConstantBeam
from radiopyo_vec.sim_env import SimEnv

file = pathlib.Path(rf"{rp.__path__[0]}\data\reactions_water.ron")

if __name__ == "__main__":
    beam = ConstantBeam("e", 5)
    env = SimEnv.from_ron_file(file)

    # Try to alter a species to check it changes everywhere
    # env.species.dyn_sp["H2O2"].label = "OOOOOOOOOOH"

    print("\nAcid Base reactions parsed: ")
    for ab in env.reactions.acid_base:
        print("\t-->", ab)

    print("\nK Reactions parsed: ")
    for kr in env.reactions.k_reaction:
        print("\t-->", kr)

    print("\nRadiolytic reactions parsed: ")
    for rad in env.reactions.radiolytic:
        print("\t-->", rad)

    print("\nEnzymatic reactions parsed: ")
    for enz in env.reactions.enzymatic:
        print("\t-->", enz)

    print(f"\nSpecies list:")
    for sp in env.species:
        print(sp)

    x = env.make_kreaction_matrix()
    np.savetxt(".//tests//output//Ri_matrix.csv", x[1], delimiter=",")
    np.savetxt(".//tests//output//Riy_matrix.csv", x[2],
               delimiter=",",
               fmt="% d")

    Gi = env.make_radiolytic_vector()
    np.savetxt(".//tests//output//Gi_matrix.csv", Gi, delimiter=",")

    Ez_subs, Ez_enz, Ez_kmki = env.make_enzymatic_matrix()
    np.savetxt(".//tests//output//Enz_subs_matrix.csv",
               Ez_subs,
               delimiter=",",
               fmt="%d")

    np.savetxt(".//tests//output//Enz_enz_matrix.csv",
               Ez_enz,
               delimiter=",",
               fmt="%d")

    np.savetxt(".//tests//output//Enz_kmki_matrix.csv",
               Ez_kmki,
               delimiter=",")
