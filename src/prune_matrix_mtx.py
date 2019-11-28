# /usr/bin/env python3

import os
import shutil
import csv
import gzip
import random
from src.utils import remove_ext
from zipfile import ZipFile
from more_itertools import first


# I (Noah Dove) have chosen to to maintain this file and it may be completely
# broken. However it is worth noting that it was already broken because for some
# reason it overwrites the actual gene expression data with the number of matrix
# entries (see line #69). So the mock matrix used in testing is completely
# corrupted, and what am I doing, obviously I need to fix this instead of just
# document it. Work for monday.

# TODO fix this shit

def prune_matrix_mtx(mtx_dir_name: str, percent_prune: float = None) -> None:
    if percent_prune is None:
        percent_prune = 0.95  # include row if larger than this value

    mtx_file = 'matrix.mtx.gz'

    # Decompress archive.
    with ZipFile(os.path.join(os.getcwd(), 'test', mtx_dir_name)) as zipObj:
        zipObj.extractall(os.path.join(os.getcwd(), 'test'))
    arc_dir = remove_ext(mtx_dir_name, '.zip')
    os.chdir(os.path.join(os.getcwd(), 'test', arc_dir))

    # Decompress matrix.mtx only.
    assert mtx_file in os.listdir('.')
    with gzip.open(mtx_file, 'rb') as f_in:
        outfile = first(os.path.splitext(f_in.name))
        with open(outfile, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(mtx_file)  # remove physical file
    mtx_file = mtx_file.rstrip('.gz')  # change value of string

    # Prune rows in MTX file according to percent_prune parameter.
    row_counter = 0
    with open(mtx_file) as infile, open('tmp.csv', 'w') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        for row in reader:
            if row_counter <= 1:
                row_counter = row_counter + 1
                writer.writerow(row)
                continue
            elif row_counter > 1:
                rand_num = random.random()
                if rand_num > percent_prune:
                    writer.writerow(row)

        os.remove(mtx_file)
        os.rename('tmp.csv', mtx_file)

    # Second (and third) pass through matrix.mtx only to correct number of entries in header.
    row_counter = 0
    with open(mtx_file) as infile, open('tmp.csv', 'w') as outfile:
        num_lines = sum(1 for line in infile) - 2  # subtract two header lines
        infile.seek(0)
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        for row in reader:
            if row_counter == 0:
                row_counter += 1
                writer.writerow(row)
            elif row_counter == 1:
                row = first(row).split(' ')
                row[2] = str(num_lines)
                writer.writerow([' '.join(row)])
            else:
                writer.writerow(row)

        os.remove(mtx_file)  # remove file so we can replace it
        os.rename('tmp.csv', mtx_file)

    # Compress file again.
    with open(mtx_file, 'rb') as f_in:
        with gzip.open(f'{mtx_file}.gz', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(mtx_file)  # remove uncompressed file

    # Create ZIP file of new mtx files, and clean up.
    os.chdir('..')
    os.remove(os.path.join(os.getcwd(), mtx_dir_name))
    with ZipFile(f'{arc_dir}.zip', 'w') as zipObj:
        for root, dirs, files in os.walk(arc_dir):
            for file in files:
                zipObj.write(os.path.join(arc_dir, file))
    shutil.rmtree(arc_dir)
