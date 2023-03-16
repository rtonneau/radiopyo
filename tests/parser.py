
import pathlib

import radiopyo as rp
from radiopyo.chemistry.species import TrackedSpecies
from radiopyo.sim_env import SimEnv

file = pathlib.Path(rf"{rp.__path__[0]}\data\reactions.ron")

sim = SimEnv.from_ron_file(file)

print("Here are all the reactions: ")
for idx, elt in enumerate(sim.reactions):
    print(f"{idx}) {elt}")

print("\n\n")

for elt in sim.species:
    if not isinstance(elt, TrackedSpecies):
        continue
    print(f"--> {elt.label}")
    for kr in elt.kreactions:
        print(f"\t {kr}")
    print("\n\n")


print(f"Sim species map: \n{sim.cc_values}")

print(f"Number of tracked Species: {sim.number_tracked_species()}")
print(f"{list(sim.iter_TrackedSpecies())}")

print(f"Initial cc: {sim.initial_values()}")
