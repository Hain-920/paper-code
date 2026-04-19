import arcpy
import os
from glob import glob

def rescale_species_data(species_name):
    base_dir = r'D:\Rwork\connectivity\aqua_data'

    species_dir = os.path.join(base_dir, species_name)
    from_dir = os.path.join(species_dir, "{}_rasterize".format(species_name))
    to_dir = os.path.join(species_dir, "{}_rescale".format(species_name))

    if not os.path.exists(to_dir):
        os.makedirs(to_dir)

    arcpy.env.workspace = from_dir

    files = [f for f in os.listdir(from_dir) if f.endswith('.tif')]

    for f in files:
        try:
            output_path = os.path.join(to_dir, f)
            arcpy.management.Resample(f, output_path, "1 1", "NEAREST")
        except Exception as e:
            print e

rescale_species_data('Gadus')