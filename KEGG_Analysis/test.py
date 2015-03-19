import re,numpy
from Smash.Utilities.Seq.Fasta import FastaReader,FastaWriter
from Smash.Utilities.Tree.Tree import Tree
from Smash.Utilities.Components import *
from Smash.Databases.GenomeDB.DB import GenomeDB
from Smash.Databases.Xref.DB import XrefDB

class KeggXrefDir(object):
    """
    directory containing <sp_id>_xrefall.list
    """
    def __init__(self,path):
        self.source_dir = path

class KeggKOFile(object):
    """
    ftp://ftp.genome.jp/pub/kegg/genes/ko
    """
    def __init__(self,path):
        self.source_file = path
        self.ko_to_description = {}
    def parse(self):
        ko_re = re.compile(r"""
            ENTRY\s+(?P<ko_id>K\d+)\s+KO\n
            NAME\s+(?P<ko_name>.+?)\n
            (DEFINITION\s+(?P<ko_definition>.+?)\n){0,1}""",re.X|re.DOTALL)
        txt = file(self.source_file).read()
        for ko_txt in [mod for mod in txt.split('\n///\n') if mod!='']:
            m = ko_re.match(ko_txt)
            if not m:
                sys.stderr.write('Error parsing ko contents:\n%s\n'%ko_txt)
                sys.exit(100)
            else:
                #print ko_txt
                ko_id = m.group('ko_id')
                ko_name = m.group('ko_name').replace('\n',' ')
                if m.group('ko_definition')!=None:
                    ko_definition = m.group('ko_definition').replace('\n',' ')
                else:
                    ko_definition = '' #None
                self.ko_to_description[ko_id]='%s - %s'%(ko_name,ko_definition)

class KeggModuleFile(object):
    """
    ftp://ftp.genome.jp/pub/kegg/pathway/module
    """
    def __init__(self,path):
        self.source_file = path
        self.module_to_description = {}
        self.module_to_kos = defaultdict(set)
        self.map_to_modules = defaultdict(set)
        self.module_to_maps = defaultdict(set)
        self.ko_to_modules = defaultdict(set)
        self.type_to_modules = defaultdict(set)
    def parse(self):
        mod_re = re.compile(r"""
            ENTRY\s+(?P<module_id>M\d+)\s+(?P<module_type>\w+)\s+Module\n
            NAME\s+(?P<module_name>.+?)\n
            (DEFINITION\s+(?P<definition>.+?)\n){0,1}
            (PATHWAY\s+(?P<pathways>.+?)){0,1}
            (ORTHOLOGY\s+(?P<orthologs>.+?))*
            (REACTION|COMPOUND|COMMENT|\Z)""",re.X|re.DOTALL)
        txt = file(self.source_file).read()
        for module_txt in [mod for mod in txt.split('\n///\n') if mod!='']:
            m = mod_re.match(module_txt)
            if not m:
                sys.stderr.write('Error parsing module contents:\n%s\n'%module_txt)
                sys.exit(100)
            else:
                #print module_txt
                module_id = m.group('module_id')
                module_type = m.group('module_type')
                module_name = m.group('module_name').replace('\n',' ')
                pathways = m.group('pathways')
                pathway_ids = []
                if pathways!=None:
                    for p in [path_atom for path_atom in pathways.split('\n') if path_atom!='']:
                        p_id = p.strip().split(None,1)[0]
                        if not p_id.startswith('ko'):
                            sys.stderr.write('Unrecognised pathway id: %s\n'%p)
                            sys.exit(100)
                        else:
                            pathway_ids.append(p_id.replace('ko',''))
                orthologs = set([])
                self.module_to_description[module_id]=module_name
                self.type_to_modules[module_type].add(module_id)
                for pathway_id in pathway_ids:
                    self.module_to_maps[module_id].add(pathway_id)
                    self.map_to_modules[pathway_id].add(module_id)
                #print 'definition',m.group('definition')
                #print 'pathway',pathway_id
                #print m.group('orthologs')
                if m.group('orthologs')!=None:
                    for l in [atom for atom in m.group('orthologs').split('\n') if atom!='']:
                        ko_ids = []
                        [ko_ids.extend(ko_atom.split('+')) for ko_atom in l.strip().split(None,1)[0].split(',')]
                        for ko_id in [ko_atom.replace('(','').replace(')','') for ko_atom in ko_ids]:
                            if ko_id.startswith('K'):
                                orthologs.add(ko_id)
                    self.module_to_kos[module_id].update(orthologs)
                    for ko in orthologs:
                        self.ko_to_modules[ko].add(module_id)
                #print module_id,module_type,module_name,pathway_id,','.join(orthologs)

class KeggTaxonomyFile(object):
    def __init__(self,path):
        self.source_file = path
        self.tree = None
        self.kegg_sp_id_to_short_name = {}
        self.short_to_long_name = {}
        self.long_to_short_name = {}
        self.full_to_short_name = {}
        self.short_to_full_name = {}
        self.depth_to_node_ids =defaultdict(list)
        self._short_to_ncbi_taxid = None
    def _get_short_to_ncbi_taxid(self):
        if self._short_to_ncbi_taxid == None:
            self._short_to_ncbi_taxid = {}
            fname = os.path.join(
                os.path.split(self.source_file)[0],
                'genome')
            if not os.path.exists(fname):
                sys.stderr.write('Can\'t get ncbi taxid mapping, you need to put ftp://ftp.genome.jp/pub/kegg/genes/genome in the same directory\n')
            else:
                txt = file(fname).read()
                recs = txt.split('//\n')
                for rec in recs[:-1]:
                    short_name,taxid = None,None
                    for l in rec.split('\n'):
                        if l.startswith('ENTRY'):
                            short_name = l.strip().split(None)[1]
                        elif l.startswith('TAXONOMY'):
                            taxid = l.strip().split(None)[1].split(':')[1]
                            break
                    self._short_to_ncbi_taxid[short_name]=taxid
        return self._short_to_ncbi_taxid
    short_to_ncbi_taxid = property(_get_short_to_ncbi_taxid)
    def parse(self):
        current_depth = 0
        depth_to_current_parent = {0:'root'}
        tree = Tree()
        tree.add_node('root')
        for l in [line for line in file(self.source_file) if line != '\n']:
            if l.startswith('#'): #internal node
                depth,node_id =l.split(None,1)
                depth = len(depth)
                node_id = node_id.strip()
                if current_depth==depth: #same depth
                    depth_to_current_parent[depth]=node_id
                elif current_depth<depth: #descending
                    depth_to_current_parent[depth]=node_id
                    current_depth = depth
                else: #ascending
                    depth_to_current_parent[depth]=node_id
                    current_depth = depth
                #print 'Internal Node'
                parent = depth_to_current_parent[current_depth-1]
            else:
                kegg_sp_id,short_name,long_name,full_name = l.strip().split('\t')
                node_id = short_name
                parent = depth_to_current_parent[current_depth]
                self.kegg_sp_id_to_short_name[kegg_sp_id]=short_name
                self.short_to_long_name[short_name]=long_name
                self.long_to_short_name[long_name]=short_name
                self.short_to_full_name[short_name]=full_name
                self.full_to_short_name[full_name]=short_name
                
            #print 'Lineage','; '.join([depth_to_current_parent[x] for x in range(current_depth+1)]),node_id
            parent_idx = tree.node_id_to_idx[parent][0]
            node_idx = tree.add_node(node_id)
            tree.add_branch(parent_idx,node_idx)
        tree.finalise()
        self.tree = tree

class baseTaskI(MultiStepInterface):
    def __init__(self,*args,**kwds):
        param = {
            'data_dir' : {
                'desc' : 'The directory where the KEGG flatfiles live [default: %default]',
                'metavar' : 'DIR',
                'group' :'input',
            },
            'outfile_stem' : {
                'desc' : 'A Unique id to identify the analysis',
                'default':'kegg_analysis'
            },
            'outdir' : {
                'desc' : 'The directory where results will be saved',
                'metavar' : 'DIR',
                'default':'.'
            },
            'formats' : {
                'desc' : 'The directory where results will be saved',
                'action':'append'
            },
            'overwrite':{
                'desc':'Overwrite an existing analysis',
                'default':False,
                'action':'store_true'
            },
            'species_to_kos' : {
                'desc': 'A file containing the mapping of species to the KO\'s, if none is present the data will be parsed from the <sp_id>_xrefall.list files',
                'metavar':'FILE',
                'group' :'input'
            },
            'limit_to_species': {
                'desc':'Limit analysis to a subset of species using the short species names. You can also use the full names of species groups to limit the analysis. A file may also be provided with one species/group id per line',
                'action':'append',
                'group':'subsets',
            },
            'remove_species':{
                'desc':'Remove species from analysis',
                'action':'append',
                'group':'subset'
            },
            'limit_to_maps': {
                'desc':'Limit analysis to a subset of maps using the short maps names. You can also use the full names of maps groups to limit the analysis. A file may also be provided with one maps/group id per line',
                'action':'append',
                'group':'subsets',
            },
            'limit_to_modules':{
                'desc':'Limit the analysis to a subset of KEGG modules',
                'action':'append',
                'group':'subsets',
            },
            'ko_subset':{
                'desc':'Use with run_step get_sequences to specify which orthologous groups you want to download',
                'action':'append'
            },
            'species_category_file':{
                'desc':'A file containing categorical data for the species',
            },
            'show_species_categories':{
                'desc':'The categorical data to use on the pca plots, can either be a taxonomic depth on the KEGG tree given in the format "tax_<depth>"',
                'action':'append',
            }
        }
        if kwds.has_key('param'):
            param.update(kwds['param'])
        kwds['param']=param
        if 'group_descriptions' not in kwds:
            kwds['group_descriptions']={}
        kwds['group_descriptions']['subsets']='Specify subsets of the KEGG data for analysis'
        kwds['group_descriptions']['input']='Input Data'
        super(baseTaskI,self).__init__(*args,**kwds)
    def check_values(self):
        super(baseTaskI,self).check_values()
        if self['limit_to_modules'] and self.interface['limit_to_maps']:
            self.logger.error('Can\'t specify both limit_to_modules and limit_to_maps\n')
            sys.exit(100)
        if 'get_sequences' in self['run_steps']:
            assert(self['ko_subset']!=None)
        if self['show_species_categories']==None:
            self['show_species_categories']=['tax_%d'%self['species_category_depth']]
            
class baseTask(MultiStepComponent):
    def __init__(self,*args,**kwds):
        if 'steps' not in kwds:
            kwds['steps'] = []
        kwds['steps'].append('print_maps')
        kwds['steps'].append('print_species')
        super(baseTask,self).__init__(*args,**kwds)
        self._map_to_kos = None
        self._kos_to_map = None
        self._map_tree = None
        self._species_to_kos = None
        self.map_to_depth = None
        self.leaf_maps = None #maps one level up from KOs
        self.map_to_description = None
        self.ko_to_description = None
        self._taxonomy = None
        self._target_species = None
        self._target_maps = None
        self._target_kos = None
        self._ko_data = None
        self._species_order = None #species ordered according to taxonomy
        self._map_order = None #maps ordered according to the map heirarchy
        self._ko_order = None #KO's ordered according to the map heirarchy
        self._modules = None
        self._target_modules = None
        self._species_categories = None
        species_order = []
    def _get_species_order(self):
        if self._species_order == None:
            self._species_order = []
            for idx in self.taxonomy.tree.leaf_iter(self.taxonomy.tree.root):
                node_id = self.taxonomy.tree.node_idx_to_id[idx]
                if node_id in self.target_species:
                    self._species_order.append(node_id)
        return self._species_order
    species_order = property(_get_species_order)
    def _get_map_order(self):
        if self._map_order == None:
            self._map_order = []
            counter = 0
            for idx in self.map_tree.internal_node_iter(self.map_tree.root):
                counter+=1
                node_id = self.map_tree.node_idx_to_id[idx]
                if node_id in self.target_maps:
                    self._map_order.append(node_id)
            if len(self.target_maps-set(self._map_order))>0:
                sys.stderr.write('Warning: different lenghts (due to paths with no KO\'s)\n')
        return self._map_order
    map_order = property(_get_map_order)
    def _get_ko_order(self):
        if self._ko_order == None:
            self._ko_order = []
            for map_id in self.map_order:
                #print self.map_to_kos[map_id]
                self._ko_order.extend(self.map_to_kos[map_id])
        return self._ko_order
    ko_order = property(_get_ko_order)
    def print_maps(self):
        self.map_tree.print_subtree(labels=self.map_to_description,print_leaves=False)
    def print_species(self):
        self.taxonomy.tree.print_subtree(labels=self.taxonomy.short_to_long_name)
    def _get_target_kos(self):
        """
        just remove the redundnacy from get_ko_order
        """
        if self._target_kos == None:
            self._target_kos = []
            for ko in self.ko_order:
                if ko not in self._target_kos:
                    self._target_kos.append(ko)
        return self._target_kos
    target_kos = property(_get_target_kos)
    def _get_target_maps(self):
        if self._target_maps == None:
            self.map_tree
            if self.interface['limit_to_maps']:
                self._target_maps = set([])
                for map_id in self.interface['limit_to_maps']:
                    if os.path.exists(map_id):
                        raise IOError,'not implemented'
                    else:
                        node_ids = self.map_tree.node_id_to_idx.get(map_id,None)
                        if node_ids!=None:
                            for node_id in node_ids:
                                map_id = self.map_tree.node_idx_to_id[node_id]
                                #print map_id,map_id in self.leaf_maps
                                if map_id in self.leaf_maps:
                                    self._target_maps.add(map_id)
                                else:
                                    for cn in self.map_tree.internal_node_iter(node_id):
                                        #print 'child node'
                                        cn_map_id=self.map_tree.node_idx_to_id[cn]
                                        if cn_map_id in self.leaf_maps:
                                            self._target_maps.add(cn_map_id)
                        else:
                            self.logger.error('Unrecognised Map ID: %s\n'%map_id)
                            sys.exit(100)
            else:
                self._target_maps = set([m for m in self.leaf_maps])
        return self._target_maps
    target_maps = property(_get_target_maps)
    def _get_target_species(self):
        if self._target_species == None:
            if self.interface['limit_to_species']!=None:
                self._target_species = set([])
                for t_sp in self.interface['limit_to_species']:
                    if os.path.exists(t_sp):
                        raise IOError,'not implemented'
                    else:
                        node_ids = self.taxonomy.tree.node_id_to_idx.get(t_sp,None)
                        if node_ids!=None:
                            for node_id in node_ids:
                                if node_id in self.taxonomy.tree.leaves:
                                    self._target_species.add(t_sp)
                                else:
                                    for leaf_idx in self.taxonomy.tree.leaf_iter(node_id):
                                        self._target_species.add(self.taxonomy.tree.node_idx_to_id[leaf_idx])
                        else:
                            self.logger.error('Unrecognised Species ID: %s\n'%t_sp)
                            sys.exit(1003)

            else:
                self._target_species = set(self.taxonomy.short_to_long_name.keys())
            if self.interface['remove_species']!=None:
                remove_list = set([])
                for t_sp in self.interface['remove_species']:
                    node_ids = self.taxonomy.tree.node_id_to_idx.get(t_sp,None)
                    if node_ids!=None:
                        for node_id in node_ids:
                            if node_id in self.taxonomy.tree.leaves:
                                remove_list.add(t_sp)
                            else:
                                for leaf_idx in self.taxonomy.tree.leaf_iter(node_id):
                                    remove_list.add(self.taxonomy.tree.node_idx_to_id[leaf_idx])
                    else:
                        self.logger.error('Unrecognised Species ID: %s\n'%t_sp)
                        sys.exit(1003)
                self._target_species = self._target_species-remove_list
        return self._target_species
    target_species = property(_get_target_species)
    def _get_modules(self):
        if self._modules == None:
            module_file = os.path.join(
                self.interface['data_dir'],
                'module')
            if not os.path.exists(module_file):
                sys.stderr.write('Cant\'t find KEGG module file %s\n'%module_file)
                sys.exit(1002)
            self._modules = KeggModuleFile(module_file)
            self._modules.parse()
        return self._modules
    modules = property(_get_modules)
    def _get_target_modules(self):
        if self._target_modules==None:
            tm = []
            if self.interface['limit_to_modules']==None:
                for kegg_map in self.target_maps:
                    for target_map in self.modules.map_to_modules[kegg_map]:
                        if target_map not in tm:
                            tm.append(target_map)
                self._target_modules = tm
            else:
                snag_list = []
                tm=[]
                for m in self.interface['limit_to_modules']:
                    if m not in self.modules.module_to_kos:
                        snag_list.append(m)
                    else:
                        if m not in tm: tm.append(m)
                if snag_list:
                    self.logger.error('Could\'t not find module with ids: %s\n'%','.join(snag_list))
                    sys.exit(1003)
                self._target_modules = tm
        return self._target_modules
    target_modules = property(_get_target_modules)
    def _get_taxonomy(self):
        if self._taxonomy==None:
            tax_file = os.path.join(
                self.interface['data_dir'],
                'taxonomy')
            if not os.path.exists(tax_file):
                self.logger.error('Can\'t find KEGG taxonomy file, should be at: %s\n'%tax_file)
                sys.exit(1002)
            self._taxonomy = KeggTaxonomyFile(tax_file)
            self._taxonomy.parse()
        return self._taxonomy
    taxonomy = property(_get_taxonomy)
    def _get_species_to_kos(self):
        """
        Get the species files from the xrefs files, the data_dir should contain
        an xrefs subdir with files <sp_id>_xrefall.list
        """
        if self._species_to_kos == None:
            if self.interface['species_to_kos'] != None and os.path.exists(self.interface['species_to_kos']):
                data = {}
                for l in file(self.interface['species_to_kos']):
                    sp_id,ko_data = l.split('\t')
                    try:
                        ko_data = [ko.split(':') for ko in ko_data.rstrip().split(None)]
                        data[sp_id] = dict([(ko,int(count)) for ko,count in ko_data])
                    except ValueError:
                        self.logger.warning('No KO\'s for %s (maybe KEGG hasn\'t curated it yet\n'%sp_id)
                self._species_to_kos = data
            else:
                xref_dir = os.path.join(self.interface['data_dir'],'xrefs')
                if not os.path.exists(xref_dir):
                    sys.stderr.write('Error: directory doesn\'t exist: %s\n'%xref_dir)
                    sys.exit(100)
                data = {}
                for sp_file in glob.glob('%s/*_xrefall.list'%xref_dir):
                    for l in file(sp_file):
                        try:
                            gene_id,xrefs = l.strip().split('\t',1)
                        except:
                            self.logger.error('No xrefs: %s'%l)
                            gene_id = l.strip()
                            xrefs = ''
                        sp_id,gene_id = gene_id.split(':')
                        if sp_id not in data:
                            data[sp_id] = defaultdict(int)
                        for xref in xrefs.split('\t'):
                            if xref:
                                source,xref_id = xref.split(':')
                                if source=='ko':
                                    kos = xref_id.strip().split()
                                    for k in kos:
                                        data[sp_id][k]+=1

                if self.interface['species_to_kos'] != None:
                    mapping_fh = file(self.interface['species_to_kos'],'w')
                    for sp_id in data:
                        mapping_fh.write('%s\t%s\n'%(sp_id,' '.join(['%s:%d'%(ko,count) for ko,count in data[sp_id].items()])))
                    mapping_fh.close()
                self._species_to_kos = data
        return self._species_to_kos
    species_to_kos = property(_get_species_to_kos)
    def _get_ko_data(self):
        if self._ko_data == None:
            ko_file = os.path.join(
                self.interface['data_dir'],
                'ko')
            if not os.path.exists(ko_file):
                sys.stderr.write('Cant find KO file: %s\n'%ko_file)
                sys.exit(1003)
            self._ko_data = KeggKOFile(ko_file)
            self._ko_data.parse()
        return self._ko_data
    ko_data = property(_get_ko_data)

    def _parse_kos(self):
        """
        Parse the brite file ko00001.keg to extract the KOs, the map heirarchy   
        """
        ko_file = os.path.join(self.interface['data_dir'],'ko00001.keg')
        if not os.path.exists(ko_file):
            self.logger.error('File %s doesn\'t exist, exiting\n'%ko_file)
            sys.exit(100)

        depth_to_current_parent = {0:'root'}
        #give a numeric depth to things
        letter_to_depth = dict([(l,x+1) for x,l in enumerate(['A','B','C','D','E','F','G'])])
        tree = Tree()
        tree.add_node('root')
        
        #regex to extract some data, just use spltis for the rest
        path_re = re.compile('\w\s+(\d+)\s(.+)')
        ko_re = re.compile('\w\s+<a.+?>(K\d+)</a>\s+(.+)')

        in_section =False
        self._map_to_kos = defaultdict(list) 
        self.map_to_description = {}
        self.ko_to_description = {}
        self.map_to_depth = {}
        self.leaf_maps = set([])
        for l in file(ko_file):
            if in_section and l!='#\n':
                if l=='!\n':
                    break
                if len(l)>2:
                    depth = letter_to_depth[l[0]]
                    node_id = None
                    if depth<=2:
                        try:
                            map_id,map_title = l.strip().split('<B>')[1].split('</B>')[0].split(None,1)
                        except:
                            print l,
                            raise
                        #print ' '*depth,map_id,map_title
                        node_id = map_id
                        self.map_to_description[map_id]=map_title
                    elif depth==3:
                        m = path_re.match(l)
                        if m:
                            kegg_id = m.group(1)
                            description = m.group(2)
                            try:
                                description,path_id = description.split('[PATH:ko')
                                path_id = path_id.split(']')[0]
                            except ValueError:
                                path_id = None
                            #print ' '*depth,kegg_id,description,path_id
                            node_id = kegg_id
                            self.leaf_maps.add(node_id)
                            self.map_to_description[kegg_id]=description
                        else:
                            print l,
                            raise ValueError
                    elif depth==4:
                        m = ko_re.match(l)
                        if m:
                            ko_id,ko_description = m.group(1),m.group(2)
                            node_id = ko_id
                        else:
                            print l,
                            raise ValueError
                        #print ' '*depth,ko_id,ko_description
                        self.ko_to_description[ko_id] = ko_description
                        self._map_to_kos[depth_to_current_parent[depth-1]].append(ko_id)
                    else:
                        sys.stderr.write('MAx depth is 4, what\'s up here:\n %s\n'%l)
                        sys.exit(100)
                    self.map_to_depth[node_id]=depth
                    depth_to_current_parent[depth]=node_id
                    parent_idx = tree.node_id_to_idx[depth_to_current_parent[depth-1]][0]
                    node_idx = tree.add_node(node_id)
                    tree.add_branch(parent_idx,node_idx)
            if l=='!\n':
                in_section=True


        unclassified_kos = set(self.ko_data.ko_to_description.keys())-set(self.ko_to_description.keys())
        if len(unclassified_kos)>0:
            # print 'Unclassified KOs',len(unclassified_kos)
            node_id = 'Unclassified'
            node_idx = tree.add_node('Unclassified')
            self.map_to_description[node_id]='Unclassified'
            tree.add_branch(tree.node_id_to_idx['root'][0],node_idx)
            self._map_to_kos[node_id] = list(unclassified_kos)
            self.ko_to_description.update(dict([(k,self.ko_data.ko_to_description[k]) for k in unclassified_kos]))
            self.map_to_depth[node_id]=1
            self.leaf_maps.add(node_id)
            for ko in unclassified_kos:
                ko_idx = tree.add_node(ko)
                tree.add_branch(node_idx,ko_idx)
        tree.finalise()
        self._map_tree = tree
        del_list = set([tree.node_idx_to_id[i] for i in tree.leaves]).intersection(self.leaf_maps)
        self.leaf_maps = self.leaf_maps - del_list
        #tree.print_subtree(node_id=node_idx,labels=self.map_to_description)
        
        #for leaf in tree.leaf_iter(0):
        #    print leaf,tree.node_idx_to_id[leaf],self.ko_to_description.get(tree.node_idx_to_id[leaf],None)
        #for internal_node in tree.internal_node_iter(0):
            #    print internal_node,tree.node_idx_to_id[internal_node],self.map_to_description.get(tree.node_idx_to_id[internal_node],None)
    def _get_map_to_kos(self):
        if self._map_to_kos==None:
            self._parse_kos()
        return self._map_to_kos
    map_to_kos = property(_get_map_to_kos)
    def _get_kos_to_map(self):
        if self._kos_to_map==None:
            self._parse_kos()
        return self._kos_to_map
    kos_to_map = property(_get_kos_to_map)
    def _get_map_tree(self):
        if self._map_tree==None:
            self._parse_kos()
        return self._map_tree
    map_tree = property(_get_map_tree)



class CompareKAASAnnotationI(baseTaskI):
    """
    Take the output of KAAS and compare to the 
    existing KEGG annotations
    """
    def __init__(self,*args,**kwds):
        param = {
            'kaas_annot_files':{
                'desc':'The KAAS annotation files that you want to compare. Should be given in the format <query_id>:<path_to_ko_file>',
                'action':'append',
                'required':True
            },
            'kaas_protein_files':{
                'desc':'The protein files annotated by KAAS. Should be given in the format <query_id>:<path_to_file>',
                'action':'append',
            },
            'species_category_depth':{
                'desc':'The depth in the KEGG taxonomic tree from which to take the species categories [default: %default]',
                'type':int,
                'default':3
            },
            'brite_category_depth':{
                'desc':'The depth in the KEGG Brite heirarchy of maps to use as categories [default: %default]',
                'type':int,
                'default':2,
            },
            'highlight_species':{
                'desc':'Highlight specific species on a graph, by default the query species will be highlighted',
                'action':'append'
                },
            'label_proportion_of_species':{
                'desc':'The proportion of species to label on the biplot [default: %default]',
                'default':0.1,
                'type':'float'
            },
            'label_proportion_of_maps':{
                'desc':'The proportion of maps to label on the biplot [default: %default]',
                'default':0.1,
                'type':'float'
            },
            'query':{
                'desc':'The query in a pairwise comparison'
            },
            'targets':{
                'desc':'The query in a pairwise comparison',
                'action':'append',
            },
            'label_proportion_of_modules':{
                'desc':'The proportion of modules to label on the biplot [default: %default]',
                'default':0.01,
                'type':'float'
            },
            'label_proportion_of_kos':{
                'desc':'The proportion of kos to label on the biplot [default: %default]',
                'default':0.01,
                'type':'float'
            },
            'plot_components':{
                'desc':'Which of the components to plot against each other, given in the format x:y [default: %default]',
                'default':'1:2'
            
            }
            
        }
        if kwds.has_key('param'):
            param.update(kwds['param'])
        kwds['param']=param
        super(CompareKAASAnnotationI,self).__init__(*args,**kwds)
    def check_values(self):
        super(CompareKAASAnnotationI,self).check_values()
        tmp = {}
        for af in self['kaas_annot_files']:
            qid,q = af.split(':')
            tmp[qid]=q
        self['kaas_annot_files'] = tmp
        tmp = {}
        if self['kaas_protein_files']:
            for af in self['kaas_protein_files']:
                qid,q = af.split(':')
                tmp[qid]=q
            self['kaas_protein_files'] = tmp

        if self['highlight_species']==None:
            self['highlight_species'] = []
        self['highlight_species'].extend(self['kaas_annot_files'].keys())
        self['plot_components'] = [int(i) for i in self['plot_components'].split(':')]
        assert(len(self['plot_components'])==2)


class CompareKAASAnnotation(baseTask):
    def __init__(self,*args,**kwds):
        if 'interface_class' not in kwds:
            kwds['interface_class']=CompareKAASAnnotationI
        if 'steps' not in kwds:
            kwds['steps']=[]
        kwds['exe_path'] = os.path.join(
                os.path.dirname(__file__),
                'bin',
                'analyse_kaas_annotation.py')

        kwds['steps'].append('summarise_kaas_annot')
        super(CompareKAASAnnotation,self).__init__(*args,**kwds)
        self._query_ko_mapping = None
        self.query_genes_unmapped = None
        self._query_ids = None
        self._map_total_matrix = None
        self._module_total_matrix = None
        self._module_completeness_matrix = None
        self._ko_abundance_matrix = None
        self._taxonomic_categories = None
        self._map_categories = None
    def _get_query_ids(self):
        if self._query_ids == None:
            self._get_query_ko_mapping()
        return self._query_ids
    query_ids = property(_get_query_ids)
    def _get_target_species(self):
        if self._target_species==None:
            super(CompareKAASAnnotation,self)._get_target_species()
            self._target_species.update(set(self.query_ids))
        return self._target_species
    target_species = property(_get_target_species)
    def _get_species_order(self):
        if self._species_order==None:
            super(CompareKAASAnnotation,self)._get_species_order()
            self._species_order.extend(self.query_ids)
        return self._species_order
    species_order = property(_get_species_order)
    def _get_query_ko_mapping(self):
        if not self._query_ko_mapping:
            self._query_ids = []
            self._query_ko_mapping = {}
            self.query_genes_unmapped =defaultdict(list)
            for qid,qf in self.interface['kaas_annot_files'].items():
                self._query_ids.append(qid)
                self._query_ko_mapping[qid] = defaultdict(list)
                for l in file(qf):
                    try:
                        gene_id,ko_id = l.strip().split()
                        self._query_ko_mapping[qid][ko_id].append(gene_id)
                        
                    except ValueError:
                        self.query_genes_unmapped[qid].append(l.strip())
                assert(qid not in self.species_to_kos)
                self.species_to_kos[qid] = dict([(k,len(v)) for k,v in self._query_ko_mapping[qid].items()])
        return self._query_ko_mapping
    query_ko_mapping = property(_get_query_ko_mapping)
    def _get_map_total_matrix(self):
        if self._map_total_matrix==None:
            matrix = numpy.zeros((len(self.target_species),
                                  len(self.target_maps)),
                                 dtype='d')
            for x,sp_id in enumerate(self.species_order):
                for y,map_id in enumerate(self.map_order):
                    #print x,sp_id,y,map_id
                    if sp_id in self.species_to_kos:
                        total = sum([self.species_to_kos[sp_id].get(k,0) for k in self.map_to_kos[map_id]])
                        matrix[x,y]=total
                    else:
                        sys.stderr.write('No KOs for species %s\n'%sp_id)
            self._map_total_matrix = matrix
        return self._map_total_matrix
    map_total_matrix = property(_get_map_total_matrix)
    def _get_module_total_matrix(self):
        if self._module_total_matrix==None:
            matrix = numpy.zeros((len(self.target_species),
                                  len(self.target_modules)),
                                 dtype='d')
            for x,sp_id in enumerate(self.species_order):
                for y,module_id in enumerate(self.target_modules):
                    if sp_id in self.species_to_kos:
                        total = sum([self.species_to_kos[sp_id].get(k,0) for k in self.modules.module_to_kos[module_id]])
                        matrix[x,y]=total
                    else:
                        sys.stderr.write('No KOs for species %s\n'%sp_id)
            self._module_total_matrix = matrix
        return self._module_total_matrix
    module_total_matrix = property(_get_module_total_matrix)   
    def _get_module_completeness_matrix(self):
        if self._module_completeness_matrix==None:
            matrix = numpy.zeros((len(self.target_species),
                                  len(self.target_modules)),
                                 dtype='f')
            for x,sp_id in enumerate(self.species_order):
                for y,module_id in enumerate(self.target_modules):
                    if sp_id in self.species_to_kos:
                        mod_kos = set(self.modules.module_to_kos[module_id])
                        if len(mod_kos)==0:
                            self.logger.warning('Module %s has no kos\n'%module_id)
                        else:
                            sp_kos = set([])
                            [sp_kos.add(k) for k in mod_kos if k in self.species_to_kos[sp_id]]
                            completeness = (1.0*len(sp_kos))/len(mod_kos)
                            matrix[x,y]=completeness
                    else:
                        sys.stderr.write('No KOs for species %s\n'%sp_id)
            self._module_completeness_matrix = matrix
        return self._module_completeness_matrix
    module_completeness_matrix = property(_get_module_completeness_matrix)   

    def _get_ko_abundance_matrix(self):
        if self._ko_abundance_matrix==None:
            self._ko_abundance_matrix = numpy.zeros(
                [len(self.target_species),len(self.target_kos)],dtype='i8')
            ko_to_idx = dict([(k,x) for x,k in enumerate(self.target_kos)])
            for x,sp in enumerate(self.species_order):
                if sp in self.species_to_kos:
                    for ko in self.species_to_kos[sp]:
                        y = ko_to_idx.get(ko,None)
                        if y != None:
                            self._ko_abundance_matrix[x,y]+=self.species_to_kos[sp][ko]
                else:
                    sys.stderr.write('No KOs for species %s\n'%sp)

        return self._ko_abundance_matrix
    ko_abundance_matrix = property(_get_ko_abundance_matrix)
    def _get_map_categories(self):
        if self._map_categories==None:
            depths = list(set(self.map_tree.node_data['depth']))
            depths.sort()
            self._map_categories = {}
            for depth in depths[1:-1]:
                cats = [None]*len(self.map_order)
                node_ids = numpy.where(self.map_tree.node_data['depth']==depth)[0] #all nodes at this depth
                for node_id in node_ids:
                    if node_id not in self.map_tree.leaves: #if it's not a leaf (leaf is a KO, no categories possible)
                        node_label = self.map_tree.node_idx_to_id[node_id]
                        is_terminal = True
                        for leaf_node in self.map_tree.internal_node_iter(node_id):
                            is_terminal = False
                            try:
                                #print self.map_tree.node_idx_to_id[leaf_node],self.map_order
                                sp_idx = self.map_order.index(
                                    self.map_tree.node_idx_to_id[leaf_node])
                                cats[sp_idx] = node_label
                            except ValueError:
                                pass
                        if is_terminal:
                            try:
                                sp_idx = self.map_order.index(node_label)
                                cats[sp_idx] = node_label
                            except:
                                pass
                if len([c for c in cats if c!=None]):
                    self._map_categories[depth]= cats
        return self._map_categories
    map_categories = property(_get_map_categories)
    def _get_species_categories(self):
        if self._species_categories==None:
            self._species_categories={}
            for depth,cats in self.taxonomic_categories.items():
                self._species_categories['tax_%d'%depth]=cats
            if self.interface['species_category_file']:
                fh = file(self.interface['species_category_file'])
                header = fh.readline()
                category_ids = header.strip().split('\t')[1:]
                for cat in category_ids:
                    self._species_categories[cat]=[None]*len(self.species_order)
                for l in fh:
                    fields = l.rstrip('\n').split('\t')
                    sp_idx = None
                    try:
                        sp_idx = self.species_order.index(fields[0])
                    except:
                        self.logger.warning('Species %s not member of target species list\n'%fields[0])
                    if sp_idx != None:
                        for x,c in enumerate(category_ids):
                            self._species_categories[c][sp_idx]=fields[x+1]
        return self._species_categories
    species_categories = property(_get_species_categories)
    def _get_taxonomic_categories(self):
        if self._taxonomic_categories==None:
            depths = list(set(self.taxonomy.tree.node_data['depth']))
            depths.sort()
            self._taxonomic_categories = {}
            for depth in depths[1:-1]:
                cats = [None]*len(self.species_order)
                node_ids = numpy.where(self.taxonomy.tree.node_data['depth']==depth)[0]
                for node_id in node_ids:
                    if node_id not in self.taxonomy.tree.leaves:
                        node_label = self.taxonomy.tree.node_idx_to_id[node_id]
                        for leaf_node in self.taxonomy.tree.leaf_iter(node_id):
                            try:
                                sp_idx = self.species_order.index(
                                    self.taxonomy.tree.node_idx_to_id[leaf_node])
                                cats[sp_idx] = node_label
                            except ValueError:
                                pass
                if len([c for c in cats if c!=None]):
                    self._taxonomic_categories[depth]= cats
        return self._taxonomic_categories
    taxonomic_categories = property(_get_taxonomic_categories)
    def summarise_kaas_annot(self,fh=sys.stdout):
        for qid,qf in self.interface['kaas_annot_files'].items():
            total_proteins_mapped = sum([len(self.query_ko_mapping[qid][k]) for k in self.query_ko_mapping[qid]])
            total_proteins_unmapped = len(self.query_genes_unmapped[qid])
            total_proteins = total_proteins_mapped+total_proteins_unmapped
            distinct_kos = len(self.query_ko_mapping[qid])

            fh.write('Query ID: %s\n KAAS file: %s\n'%(qid,qf))
            fh.write(' total proteins: %d\n'%(total_proteins))
            fh.write(' with KO number: %d (%.2f)\n'%(total_proteins_mapped,float(total_proteins_mapped)/total_proteins))
            fh.write('   distinct KOs: %d\n'%distinct_kos)
    def pca_on_ko_abundance(self):
        from Smash.Utilities.Graphics.statsGraphs import PCATask
        m = self.ko_abundance_matrix
        for cat in self.interface['show_species_categories']:
            pca_data = PCATask.graph_data(
                m,
                self.species_order,
                self.target_kos,
                observation_categories = [(cat,self.species_categories[cat])],
            )
            pca_g = PCATask()
            axes_dict = pca_g.run(
                outfile_stem = 'pca_ko_abundance_%s_%s'%(cat,self.interface['outfile_stem']),
                formats = self.interface['formats'],
                plot_components = self.interface['plot_components'],
                figsize = [20,20],
                bbox = [0.2,0.1,0.7,0.8],
                plot_variables = True,
                label_variables = True,
                label_observations = True,
                highlight_observations = self.interface['highlight_species'],
                outdir = self.interface['outdir'],
                label_observations_p = self.interface['label_proportion_of_species'],
                label_variables_p = self.interface['label_proportion_of_kos'],
                data = pca_data)
            annot_text = []
            for v in  pca_g.var_sorted[:int(self.interface['label_proportion_of_kos']*len(pca_g.var_sorted))]:
                annot_text.append('%s - %s'%(self.target_kos[v],self.ko_to_description[self.target_kos[v]]))
            annot_text = '\n'.join(annot_text)
            t = axes_dict['main'].text(0.001,0.01,
                                       annot_text,
                                       size = 6,
                                       bbox = {'fc':'none','alpha':0.5,'lw':2.0,'ec':'white','boxstyle':'round'},
                                       transform=pca_g.figure.transFigure,
                                       color='w')
            pca_g.save()

    def create_tables(self):
        #print self.species_order 

        all_sp_ids =  [self.interface['query']] + list(self.interface['targets'])
        [all_sp_ids.append(sp_id) for sp_id in self.species_order if sp_id not in all_sp_ids]
        sp_indices = [self.species_order.index(sp_id) for sp_id in all_sp_ids]
        print sp_indices
        sep = '\t'
        print all_sp_ids
        sorter = []
        for x,ko_id in enumerate(self.target_kos):
            abundances = self.ko_abundance_matrix[:,x].take(sp_indices)
            if abundances.any():
                #print x,ko_id,self.ko_to_description[ko_id]
                #print len(abundances)
                #print abundances
                diff = abs(abundances[0]-abundances[1])
                sorter.append((diff,ko_id,abundances))
        sorter.sort(reverse =True)
        outfh = file(os.path.join(self.interface['outdir'],
                                  '%s_ko_abudance.csv.txt'%self.interface['outfile_stem']),'w')
        outfh.write('%s\n'%sep.join(['KO','KO description','KEGG URL']+
                                     [self.taxonomy.short_to_long_name.get(sp_id,sp_id) for sp_id in all_sp_ids]))
        for diff,ko_id,abundances in sorter:
            outfh.write('%s\n'%sep.join(
                [ko_id,self.ko_to_description[ko_id],'http://www.genome.jp/dbget-bin/www_bget?ko+%s'%ko_id]+['%d'%a for a in abundances]))

    def pairwise_comparison(self):
        qid = self.interface['query']
        tids = self.interface['targets']
        all_sp_ids =  list(set([self.interface['query']]+self.interface['targets']))
        self.interface['limit_to_species'] = [sp for sp in all_sp_ids if sp not in self.query_ko_mapping]
        qidx = self.species_order.index(qid)
        #map totals
        q_map_totals = self.map_total_matrix[qidx,:]
        q_ko_totals = self.ko_abundance_matrix[qidx,:]
        for tid in tids:
            tidx = self.species_order.index(tid)
            sorter = []
            for y,ko in enumerate(self.target_kos):
                q_abund =  self.ko_abundance_matrix[qidx,y]
                t_abund =  self.ko_abundance_matrix[tidx,y]
                if q_abund>0 or t_abund>0:
                    #print ko,q_abund,t_abund,self.ko_to_description[ko]
                    sorter.append((abs(q_abund-t_abund),[ko,q_abund,t_abund,self.ko_to_description[ko]]))
            sorter.sort(reverse=True)
            print 'KOs'
            for sort_id,s in sorter[:100]:
                print ' '.join(map(str,s))
            sorter = []
            for y,map_id in enumerate(self.map_order):
                q_abund =  self.map_total_matrix[qidx,y]
                t_abund =  self.map_total_matrix[tidx,y]
                if q_abund>0 or t_abund>0:
                    #print ko,q_abund,t_abund,self.ko_to_description[ko]
                    sorter.append((abs(q_abund-t_abund),[map_id,q_abund,t_abund,self.map_to_description[map_id]]))
            sorter.sort(reverse=True)
            print 'Maps'
            for sort_id,s in sorter[:100]:
                print ' '.join(map(str,s))
            
            sorter = []
            for y,mod_id in enumerate(self.target_modules):
                q_abund =  self.module_total_matrix[qidx,y]
                t_abund =  self.module_total_matrix[tidx,y]
                if q_abund>0 or t_abund>0:
                    #print ko,q_abund,t_abund,self.ko_to_description[ko]
                    sorter.append((abs(q_abund-t_abund),[mod_id,q_abund,t_abund,self.modules.module_to_description[mod_id]]))
            sorter.sort(reverse=True)
            print 'Modules'
            for sort_id,s in sorter[:100]:
                mod_id,q_abund,t_abund,mod_desc = s
                print ' '.join(map(str,s))
                for ko in self.modules.module_to_kos[mod_id]:
                    try:
                        ko_idx = self.target_kos.index(ko)
                    except:
                        ko_idx = None
                        sys.stderr.write('Dodgy KO: %s\n'%ko)
                    if ko_idx:
                        q_abund =  self.ko_abundance_matrix[qidx,ko_idx]
                        t_abund =  self.ko_abundance_matrix[tidx,ko_idx]
                        print ' ',ko,q_abund,t_abund,self.ko_to_description.get(ko,None)
 
        print q_ko_totals
        print q_map_totals.sum()
        print q_ko_totals.sum()
        #print self.map_total_matrix 
    def pca_on_module_totals(self):
        from Smash.Utilities.Graphics.statsGraphs import PCATask
        m =  self.module_total_matrix
        if m.shape[1]==0:
            sys.stderr.write('No modules for this selection\n')
            return
        for cat in self.interface['show_species_categories']:
          
            pca_data = PCATask.graph_data(
                m,
                self.species_order,
                self.target_modules,
                observation_categories = [(cat,self.species_categories[cat])],
            )
            pca_g = PCATask()
            axes_dict = pca_g.run(
                outfile_stem = 'pca_module_abundance_%s_%s'%(cat,self.interface['outfile_stem']),
                formats = self.interface['formats'],
                figsize = [15,15],
                plot_components = self.interface['plot_components'],
                bbox = [0.2,0.1,0.7,0.8],
                plot_variables = True,
                label_variables = True,
                label_observations = True,
                highlight_observations = self.interface['highlight_species'],
                outdir = self.interface['outdir'],
                label_observations_p = self.interface['label_proportion_of_species'],
                label_variables_p = self.interface['label_proportion_of_modules'],
                data = pca_data)
            annot_text = []
            cutoff = int(self.interface['label_proportion_of_modules']*len(pca_g.var_sorted))
            for v in  pca_g.var_sorted[:cutoff]:
                annot_text.append('%s - %s'%(self.target_modules[v],self.modules.module_to_description[self.target_modules[v]]))
            annot_text = '\n'.join(annot_text)
            t = axes_dict['main'].text(0.001,0.01,
                                       annot_text,
                                       size = 6,
                                       bbox = {'fc':'none','alpha':0.5,'lw':2.0,'ec':'white','boxstyle':'round'},
                                       transform=pca_g.figure.transFigure,
                                       color='w')
            pca_g.save()
        
    def get_sequences(self):
        """
        get sequences for a particular KO
        """
        #print self.interface['ko_subset']
        from Smash.Utilities.DataSources.KEGG.Soap import KEGGQuery
        kq = KEGGQuery()
        query_fasta = {}
        for query_id,qf in self.interface['kaas_protein_files'].items():
            query_fasta[query_id] = FastaReader(qf)
            query_fasta[query_id].as_dict()
        for ko in self.interface['ko_subset']:
            out_fasta_file = os.path.join(self.interface['outdir'],
                                          '%s_%s.proteins.fa'%(self.interface['outfile_stem'],
                                                               ko))
            out_fasta_fh = file(out_fasta_file,'w')
            writer = FastaWriter(out_fasta_fh)
            target_gene_ids = []
            for org in self.interface['highlight_species']:
                gene_ids =kq.get_genes_by_ko('ko:%s'%ko,org)
                target_gene_ids.extend([g[0] for g in gene_ids])
            seqs_str = kq.get_gene_sequences(gene_list=target_gene_ids)
            out_fasta_fh.write(seqs_str)
            out_fasta_fh.flush()
            for qid in query_fasta:
                for gene_id in self.query_ko_mapping[qid][ko]:
                    writer('%s:%s'%(qid,gene_id),query_fasta[qid][gene_id])
            out_fasta_fh.close()
            print out_fasta_file

    def pca_on_map_totals(self):
        from Smash.Utilities.Graphics.statsGraphs import PCATask
        m = self.map_total_matrix
        for cat in self.interface['show_species_categories']:
            pca_data = PCATask.graph_data(
                m,
                self.species_order,
                self.map_order,
                observation_categories = [(cat,self.species_categories[cat])],
                variable_categories = [(self.interface['brite_category_depth'],
                                        ['%s (%s)'%(self.map_to_description[c],c) for c in self.map_categories[self.interface['brite_category_depth']]])],
            )
            
            pca_g = PCATask()
            axes_dict = pca_g.run(
                formats = self.interface['formats'],
                figsize = [15,15],
                plot_components = self.interface['plot_components'],
                bbox = [0.2,0.1,0.7,0.9],
                plot_variables = True,
                label_variables = True,
                label_observations = True,
                highlight_observations = self.interface['highlight_species'],
                outdir = self.interface['outdir'],
                outfile_stem = 'pca_map_totals_%s_%s'%(cat,self.interface['outfile_stem']),
                label_observations_p = self.interface['label_proportion_of_species'],
                label_variables_p = self.interface['label_proportion_of_maps'],
                data = pca_data)
            annot_text = []
            for v in  pca_g.var_sorted[:int(self.interface['label_proportion_of_maps']*len(pca_g.var_sorted))]:
                annot_text.append('%s - %s'%(self.map_order[v],self.map_to_description[self.map_order[v]]))
            annot_text = '\n'.join(annot_text)
            t = axes_dict['main'].text(0.001,0.01,
                                       annot_text,
                                       size = 6,
                                       bbox = {'fc':'none','alpha':0.5,'lw':2.0,'ec':'white','boxstyle':'round'},
                                       transform=pca_g.figure.transFigure,
                                       color='w')
            pca_g.save()
    def draw_maps(self):
        pass
    def plot_module_completeness(self):
        from Smash.Utilities.Graphics.graphTasks import MatrixGraph,CategoryGraph
        matrix = self.module_completeness_matrix
        matrix = numpy.ma.masked_equal(matrix,0.0,copy=True)
        if matrix.shape[1]==0:
            sys.stderr.write('No modules for this selection\n')
            return


        #get pathway categories
        #get species categories
        label_species_depth = self.interface['species_category_depth']
        node_ids = numpy.where(self.taxonomy.tree.node_data['depth']==label_species_depth)[0]
        sp_category_data = {'category_to_items':{},
                            'category_labels':{}}
        for species_cat in node_ids:
            leaf_nodes = set([self.taxonomy.tree.node_idx_to_id[s] for s in self.taxonomy.tree.leaf_iter(species_cat)]).intersection(self.target_species)
            if len(leaf_nodes)>0:
                sp_category_data['category_to_items'][species_cat] = [self.taxonomy.short_to_long_name.get(sp,sp) for sp in leaf_nodes]
                sp_category_data['category_labels'][species_cat] = self.taxonomy.tree.node_idx_to_id[species_cat]
        row_labels = [self.taxonomy.short_to_long_name.get(sp,sp) for sp in self.species_order]
        #cg = CategoryGraph()
        cgd = CategoryGraph.graph_data(**sp_category_data)
        num_cats = len(sp_category_data['category_to_items'])
        row_cat_width = num_cats*0.25
        row_label_width = 0.25*10
        col_labels = ['%s - %s'%(mod,self.modules.module_to_description[mod]) for mod in self.target_modules]
        figsize =  [0.25*matrix.shape[1],
                    0.25*matrix.shape[0]]
        figsize[0]+=(row_cat_width+row_label_width)
        #figsize[1]+=(col_label_height+col_cat_height)
        ag = MatrixGraph()
        ad = MatrixGraph.graph_data(matrix,
                                    #numpy.ma.masked_equal(matrix,0,copy=True),
                                    row_labels = row_labels,
                                    col_labels=col_labels)
        ag.run(
            data = ad,
            show_row_totals = False,
            figsize=figsize,
            #            row_category_data = cgd,
            #            row_category_a = row_cat_width+row_label_width,
            #            row_category_padding = float(row_label_width)/(row_cat_width+row_label_width),
            #            col_category_data = mcgd,
            #            col_category_a = col_cat_height+col_label_height,
            #            col_category_padding = float(col_label_height)/(col_cat_height+col_label_height),
            #            col_label_pos = 'top',
            #            logscale=True,
            #bgcolor='black',
            #fgcolor='white',
            logscale = False,
            #colormap = 'autumn',
            cluster_cols = 'single:euclidean',
            cluster_rows = 'single:euclidean',
            show_col_dendrogram = False,
            show_row_dendrogram = False,
            outfile_stem = 'module_completeness_%s'%self.interface['outfile_stem'],
            outdir = self.interface['outdir'],
            formats = self.interface['formats'],
            show_row_totals = True,
            show_col_totals = True)
        ag.save()

    def plot_module_abundance(self):
        from Smash.Utilities.Graphics.graphTasks import MatrixGraph,CategoryGraph
        matrix = self.module_total_matrix
        if matrix.shape[1]==0:
            sys.stderr.write('No modules for this selection\n')
            return
        #get pathway categories
        #get species categories
                
        label_species_depth = self.interface['species_category_depth']
        node_ids = numpy.where(self.taxonomy.tree.node_data['depth']==label_species_depth)[0]
        sp_category_data = {'category_to_items':{},
                            'category_labels':{}}
        for species_cat in node_ids:
            leaf_nodes = set([self.taxonomy.tree.node_idx_to_id[s] for s in self.taxonomy.tree.leaf_iter(species_cat)]).intersection(self.target_species)
            if len(leaf_nodes)>0:
                sp_category_data['category_to_items'][species_cat] = [self.taxonomy.short_to_long_name.get(sp,sp) for sp in leaf_nodes]
                sp_category_data['category_labels'][species_cat] = self.taxonomy.tree.node_idx_to_id[species_cat]
        row_labels = [self.taxonomy.short_to_long_name.get(sp,sp) for sp in self.species_order]
        #cg = CategoryGraph()
        cgd = CategoryGraph.graph_data(**sp_category_data)
        num_cats = len(sp_category_data['category_to_items'])
        row_cat_width = num_cats*0.25
        row_label_width = 0.25*10
        col_labels = ['%s - %s'%(mod,self.modules.module_to_description[mod]) for mod in self.target_modules]
        figsize =  [0.25*matrix.shape[1],
                    0.25*matrix.shape[0]]
        figsize[0]+=(row_cat_width+row_label_width)
        #figsize[1]+=(col_label_height+col_cat_height)
        ag = MatrixGraph()
        ad = MatrixGraph.graph_data(matrix,
                                    #numpy.ma.masked_equal(matrix,0,copy=True),
                                    row_labels = row_labels,
                                    col_labels=col_labels)
        ag.run(
            data = ad,
            show_row_totals = False,
            figsize=figsize,
            #            row_category_data = cgd,
            #            row_category_a = row_cat_width+row_label_width,
            #            row_category_padding = float(row_label_width)/(row_cat_width+row_label_width),
            #            col_category_data = mcgd,
            #            col_category_a = col_cat_height+col_label_height,
            #            col_category_padding = float(col_label_height)/(col_cat_height+col_label_height),
            #            col_label_pos = 'top',
            #            logscale=True,
            #bgcolor='black',
            #fgcolor='white',
            logscale = True,
            colormap = 'autumn',
            cluster_cols = 'single:euclidean',
            cluster_rows = 'single:euclidean',
            show_col_dendrogram = False,
            show_row_dendrogram = False,
            outfile_stem = 'module_abundance_%s'%self.interface['outfile_stem'],
            outdir = self.interface['outdir'],
            formats = self.interface['formats'],
            show_row_totals = True,
            show_col_totals = True)
        ag.save()
    def plot_ko_abundance(self):
        """
        plot the number of kos per organism in the selected pathways
        """
        from Smash.Utilities.Graphics.graphTasks import MatrixGraph,CategoryGraph
        matrix = self.ko_abundance_matrix
        #get pathway categories
        #get species categories
        label_species_depth = self.interface['species_category_depth']
        node_ids = numpy.where(self.taxonomy.tree.node_data['depth']==label_species_depth)[0]
        sp_category_data = {'category_to_items':{},
                            'category_labels':{}}
        for species_cat in node_ids:
            leaf_nodes = set([self.taxonomy.tree.node_idx_to_id[s] for s in self.taxonomy.tree.leaf_iter(species_cat)]).intersection(self.target_species)
            if len(leaf_nodes)>0:
                sp_category_data['category_to_items'][species_cat] = [self.taxonomy.short_to_long_name.get(sp,sp) for sp in leaf_nodes]
                sp_category_data['category_labels'][species_cat] = self.taxonomy.tree.node_idx_to_id[species_cat]
        row_labels = [self.taxonomy.short_to_long_name.get(sp,sp) for sp in self.species_order]
        #cg = CategoryGraph()
        cgd = CategoryGraph.graph_data(**sp_category_data)
        num_cats = len(sp_category_data['category_to_items'])
        row_cat_width = num_cats*0.25
        row_label_width = 0.25*10
        col_labels = ['%s - %s'%(ko,self.ko_to_description[ko]) for ko in self.target_kos]
        pass_indices = []
        for y in range(len(self.target_kos)):
            if matrix[:,y].any():
                pass_indices.append(y)

        col_labels = list(numpy.array(col_labels).take(pass_indices))
        new_matrix = numpy.zeros((matrix.shape[0],len(pass_indices)),dtype='i4')
        for x in range(len(row_labels)):
            new_matrix[x:,] = matrix[x:,].take(pass_indices)
        matrix = new_matrix

        figsize =  [15+0.25*matrix.shape[1],
                    5+0.25*matrix.shape[0]]
        figsize[0]+=(row_cat_width+row_label_width)
        #figsize[1]+=(col_label_height+col_cat_height)
        ag = MatrixGraph()
        ad = MatrixGraph.graph_data(matrix,
                                    row_labels = row_labels,
                                    col_labels=col_labels)
        ag.run(
            data = ad,
            show_row_totals = False,
            figsize=figsize,
                        row_category_data = cgd,
                        row_category_a = row_cat_width+row_label_width,
                        row_category_padding = float(row_label_width)/(row_cat_width+row_label_width),
            #            col_category_data = mcgd,
            #            col_category_a = col_cat_height+col_label_height,
            #            col_category_padding = float(col_label_height)/(col_cat_height+col_label_height),
            #            col_label_pos = 'top',
            #            logscale=True,
            #bgcolor='black',
            #fgcolor='white',
            bbox = [0.1,0.1,0.8,0.6],
            logscale = True,
            colormap = 'autumn',
            cluster_cols = 'single:euclidean',
            cluster_rows = 'single:euclidean',
            show_col_dendrogram = False,
            show_row_dendrogram = False,
            outfile_stem = 'ko_abundance_%s'%self.interface['outfile_stem'],
            outdir = self.interface['outdir'],
            formats = self.interface['formats'],
            show_row_totals = True,
            show_col_totals = True)
        ag.save()

    def graph_map_totals(self):
        """
        plot the total number of KO's for each of the target pathways
        """
        from Smash.Utilities.Graphics.graphTasks import MatrixGraph,CategoryGraph
        matrix = self.map_total_matrix

        #get pathway categories
        node_ids = numpy.where(self.map_tree.node_data['depth']==self.interface['brite_category_depth'])[0]
        map_category_data  = {'category_to_items':{},'category_labels':{}}
        for map_cat in node_ids:
            leaf_nodes = set([self.map_tree.node_idx_to_id[s] for s in self.map_tree.internal_node_iter(map_cat)]).intersection(self.target_maps)
            if len(leaf_nodes)>0:
                map_category_data['category_to_items'][map_cat] = [self.map_to_description[m] for m in leaf_nodes]
                map_category_data['category_labels'][map_cat]=self.map_to_description[self.map_tree.node_idx_to_id[map_cat]]
        #get species categories
        label_species_depth = self.interface['species_category_depth']
        node_ids = numpy.where(self.taxonomy.tree.node_data['depth']==label_species_depth)[0]
        sp_category_data = {'category_to_items':{},
                            'category_labels':{}}
        for species_cat in node_ids:
            leaf_nodes = set([self.taxonomy.tree.node_idx_to_id[s] for s in self.taxonomy.tree.leaf_iter(species_cat)]).intersection(self.target_species)
            if len(leaf_nodes)>0:
                sp_category_data['category_to_items'][species_cat] = [self.taxonomy.short_to_long_name.get(sp,sp) for sp in leaf_nodes]
                sp_category_data['category_labels'][species_cat] = self.taxonomy.tree.node_idx_to_id[species_cat]
        row_labels = [self.taxonomy.short_to_long_name.get(sp,sp) for sp in self.species_order]
        #cg = CategoryGraph()
        cgd = CategoryGraph.graph_data(**sp_category_data)
        mcgd = CategoryGraph.graph_data(**map_category_data)
        num_cats = len(sp_category_data['category_to_items'])
        row_cat_width = num_cats*0.25
        col_cat_height = len(map_category_data['category_to_items'])*0.25
        col_label_height = 0.25*20
        row_label_width = 0.25*10
        col_labels = [self.map_to_description[m] for m in self.map_order]
        figsize =  [0.25*matrix.shape[1],
                    0.25*matrix.shape[0]]
        figsize[0]+=(row_cat_width+row_label_width)
        figsize[1]+=(col_label_height+col_cat_height)
        ag = MatrixGraph()
        ad = MatrixGraph.graph_data(matrix,
                                    row_labels = row_labels,
                                    col_labels=col_labels)
        ag.run(
            data = ad,
            show_row_totals = False,
            figsize=figsize,
            #cluster_cols = 'single:euclidean',
            #cluster_rows = 'single:euclidean',
            #            row_dendrogram_padding = 0.05,
            #            row_category_data = cgd,
            #            row_category_a = row_cat_width+row_label_width,
            #            row_category_padding = float(row_label_width)/(row_cat_width+row_label_width),
            #            col_category_data = mcgd,
            #            col_category_a = col_cat_height+col_label_height,
            #            col_category_padding = float(col_label_height)/(col_cat_height+col_label_height),
            #            col_label_pos = 'top',
            #            logscale=True,
            #bgcolor='black',
            #fgcolor='white',
            outfile_stem = 'map_total_matrix_%s'%self.interface['outfile_stem'],
            outdir = self.interface['outdir'],
            formats = self.interface['formats'],
            show_row_totals = True,
            show_col_totals = True)
        ag.save()
        mask_rows = []
         


class MapGenomeDBToKEGGI(baseTaskI):
    def __init__(self,*args,**kwds):
        param = {
            'genomedb_url' : {
                'desc' : 'The url of the genome db [default: %default]',
                'required' : True,
            },
            'xrefdb_url' : {
                'desc':'The url of an xrefdb to map the ids to KEGG',
            }
        }
        if 'param' not in kwds:
            kwds['param'] = param
        else:
            kwds['param'].update(param)
        self._gdb = None
        super(MapGenomeDBToKEGGI,self).__init__(*args,**kwds)


class MapGenomeDBToKEGG(baseTask):
    def __init__(self,*args,**kwds):
        if 'interface_class' not in kwds:
            kwds['interface_class']=MapGenomeDBToKEGGI
        if 'steps' not in kwds:
            kwds['steps']=[]
        kwds['exe_path'] = os.path.join(
                os.path.dirname(__file__),
                'bin',
                'map_genomedb_to_KEGG.py')
        kwds['steps'].append('do_mapping')
        super(MapGenomeDBToKEGG,self).__init__(*args,**kwds)
        self._gdb = None
        self._xrefdb = None
    def _get_xrefdb(self):
        if self._xrefdb == None:
            self._xrefdb =XrefDB(db_url=self.interface['xrefdb_url'])
        return self._xrefdb
    xrefdb = property(_get_xrefdb)
    def _get_gdb(self):
        if not self._gdb:
            self._gdb = GenomeDB(db_url=self.interface['genomedb_url'])
        return self._gdb
    gdb = property(_get_gdb)
    def map_by_ids(self):
        q = select([self.gdb.tables['genes'].c.name]).where(self.gdb.tables['genes'].c.type=='coding')
        gene_ids = set([r[0] for r in q.execute()])
        print self.xrefdb
        protein_names = set([r[0] for r in select([self.gdb.tables['proteins'].c.name]).execute()])
        for p in protein_names:
            xrefs = self.xrefdb.get_xrefs(p)
            print p
            print xrefs
        res = {}
        return res
    def map_by_exact_sequence(self,target_taxa):
        q_fasta = self.gdb.protein_fasta[1]
        q_seq_to_ids = defaultdict(set)
        for seqid,seq in q_fasta.min_iter():
            q_seq_to_ids[seq].add(seqid)
    def summarise_mapping(self,mapping,fh=sys.stdout):
        fh.write('Mapping:\n')
        for sp in mapping:
            fh.write(' Target species: %s\n'%sp)
            fh.write('   total queries: %d\n'%mapping[sp]['summary'])
    def do_mapping(self):
        print self.gdb
        print dir(self.gdb.db_md)
        print self.gdb.db_md.organism_name
        print self.gdb.db_md.organism_taxid
        #self.create_species_files()
        id_mapping = self.map_by_ids()
        self.summarise_mapping(id_mapping)
        sp_order = [(id_mapping[sp]['overlap'],sp) for sp in id_mapping]
        sp_order.sort(reverse=True)
        #self.map_by_exact_sequence(id_mapping.keys())
