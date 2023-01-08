# -*- coding: utf-8 -*-
from os.path import abspath, join, dirname
from os import remove


class ColorizeSVG:

    def __init__(self, folder='workdir', input_file='phyml_Blosum62_test1.svg'):
        '''Create a colorized SVG object'''
        dir = abspath(join(dirname(__file__), '..', folder))
        self.input_path = join(dir, input_file)
        self.output_path = join(dir, input_file[:-4] + '_colorized.svg')
        try:
            remove(self.output_path)
        except:
            pass

    def colorize_seaview_phylo_tree(self):
        '''Colorize the id's of the phylogenetic tree SVG file from the seaview'''
        with open(self.input_path) as f1, open(self.output_path, 'a') as f2:
            # Verify each line to colorize it
            for line in f1:
                if 'Sev' in line:
                    color = 'red'
                elif 'Mod' in line:
                    color = 'blue'
                elif 'Mil' in line:
                    color = 'green'
                else:
                    color = 'black'

                line = line.replace('rgb(0,0,0)', color)
                f2.write(line)
