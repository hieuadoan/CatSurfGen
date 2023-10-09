#Doping with Ni adatom
from catkit.build import molecule
from catkit.gen.surface import SlabGenerator
from catkit.gen.adsorption import Builder
from ase.io import read
import os

def generate_adsorption_structures(slab,file_base,site,baseDir):
    adsorbate_dict = { "H2O": molecule('H2O')[0]}
    ads_slabs = []
    count =0
    for key, value in adsorbate_dict.items():
        builder = Builder(slab)
        ads_slab = builder.add_adsorbate(adsorbate_dict[key], bonds=[0], index=-1)
        print("Printing number of Ads sites for slabs: {}". format(len(ads_slab)))

        for i in range(len(ads_slab)):
            print("{}_{}_{}_config_{}".format(file_base, site, key, i))
            ads_slab[i].write(baseDir+'/POSCAR_{}_{}_{}_config_{}'.format(file_base, site, key, i))
            count+=1
        ads_slabs.append(ads_slab)
    print("Total number of structures {}".format(count))
    return(ads_slabs, adsorbate_dict)


def generate_surfaces(filename,baseDir):
    slabs = []
    f = read(filename)
    file_base = filename.split('.')[0]
    indices = {"100":(1, 0, 0)}
               #, "010":(0, 1, 0), "001" :(0, 0, 1),
               #"110":(1, 1, 0),"101":(1, 0, 1), "011": (0, 1, 1),
               #"111": (1, 1, 1)}
    terminations = {"100":[0]}
                    #, "010":[0], "001" :[0],
                    #"110":[0], "101":[0], "011":[0], "111":[0]}
    for key, val in indices.items():
        gen = SlabGenerator(f, miller_index=indices[key], standardize_bulk=True,
                            primitive=True, layers=6, fixed=2, vacuum=10, tol=1e-20)
        slab= gen.get_slab(iterm=terminations[key][0], size=(3,3))
        print("{} {}".format(file_base, key))
        slab.write(baseDir+'/POSCAR_{}_{}_clean'.format(file_base,key))
        ads_slabs,_ = generate_adsorption_structures(slab, file_base, key, baseDir)
        print(len(ads_slabs))
        slabs.append(slab)
    return(slabs, indices)


def main():

    dirName = 'GenDir'
    try:
        os.mkdir(dirName)
        print("Directory " , dirName ,  " Created ")
    except FileExistsError:
        print("Directory " , dirName ,  " already exists")

    filename = "Pd.cif"
    slabs, indices = generate_surfaces(filename, dirName)

if  __name__ == "__main__":

    main()
