
import radiopyo as rp
from radiopyo.utils.ron_parser import RonFileParser
import pathlib


file = pathlib.Path(rf"{rp.__path__[0]}\data\reactions_bio.ron")

parser = RonFileParser(file)
parser.parse_file()

x, y = parser.build_species_list()

for elt in y:
    print(elt)

for elt in x:
    print(elt)
