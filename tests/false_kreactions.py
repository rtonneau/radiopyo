
from radiopyo.chemistry.k_reactions import ReactionSpecies, KReaction
from radiopyo.chemistry import RawSpecies

new = KReaction.new_empty()

new.add_as_reactant("e_aq")
new.add_as_reactant("H2O2")
new.add_as_reactant("H2O2")

new.add_as_product("O_minus")
new.add_as_product("H_plus")
new.add_as_product("OH_r")

sp = new.reactants[0]
print(f"{new}")

print(f"Took one species: {sp}")
print(f"Is instance RawSpecies: ", isinstance(sp, RawSpecies))
print(f"Equality with 'e_aq' : ", sp == "e_aq")

print(f"Get idx of 'e_aq': ", new._get_reactant_idx("e_aq"))

