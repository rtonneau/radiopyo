
import pathlib

import radiopyo_vec as rp
from radiopyo_vec.physics.beam import ConstantBeam
from radiopyo_vec.sim_env import SimEnv

file = pathlib.Path(rf"{rp.__path__[0]}\data\reactions_bio.ron")

if __name__ == "__main__":
    beam = ConstantBeam("e", 5)
    env = SimEnv.from_ron_file(file)

    env.make_kreaction_matrix()