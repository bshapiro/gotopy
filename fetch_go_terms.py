from collections import Counter
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.base import get_godag
from goatools.godag.go_tasks import get_go2parents
from goatools.gosubdag.gosubdag import GoSubDag
from pickle import dump, load
from collections import defaultdict
import mygene
import pandas as pd
import numpy as np
from glob import glob

METH_ANNOTATION_FILE = '../data/hm27_hg19_annotations.txt'  # change for your purposes


def get_terms(entrez_dict, name, namespace, goto_type='0'):
    objanno = Gene2GoReader("data/gene2go", taxids=[9606])  # human
    gene2goids_human = objanno.get_id2gos(namespace=namespace, taxids=[9606])  # namespace: BP, MF, CC
    go2geneids_human = objanno.get_id2gos(namespace=namespace, go2geneids=True, taxids=[9606])
    go2geneids_counter = Counter()
    for go in go2geneids_human:
        go2geneids_counter[go] = len(go2geneids_human[go])
    godag = get_godag('data/go-basic.obo', optional_attrs='relationship')
    go_dict = {}
    for entity in entrez_dict.keys():
        entrez_ids = entrez_dict[entity]
        if entrez_ids == []:
            go_dict[entity] = ['NO ENTREZ']
            continue
        for entrez_id in entrez_ids:
            goids = gene2goids_human.get(int(entrez_id))
            if goids is None:
                if go_dict.get(entity) is None:
                    go_dict[entity] = ['NO GOIDS']
                continue
            gosubdag = GoSubDag(goids, godag, relationships={'part_of', }, prt=False)
            # go_info = [gosubdag.go2nt.get(goid) for goid in goids]
            go2parents = get_go2parents(gosubdag.go2obj, gosubdag.relationships)
            go2parents_nr = list(set(go2parents.keys()))
            go2parents_nr_filtered = go2parents_nr
            if goto_type == '0200':  # filter to GO terms with more than 0 gene products and less than 200
                go2parents_nr_filtered = [go for go in go2parents_nr if go2geneids_counter[go] > 0 and go2geneids_counter[go] < 200]
            elif goto_type == '0':  # filter to GO terms with more than 0 gene products
                go2parents_nr_filtered = [go for go in go2parents_nr if go2geneids_counter[go] > 0]
            elif goto_type == 'all':  # do not filter GO terms
                go2parents_nr_filtered = go2parents_nr
            if go_dict.get(entity) == ['NO GOIDS'] or go_dict.get(entity) == ['NO ENTREZ'] or go_dict.get(entity) is None:
                go_dict[entity] = [go2parents_nr_filtered]
            else:
                go_dict[entity].append(go2parents_nr_filtered)

    dump(go_dict, open('data/go_' + name + '_' + namespace + '_' + goto_type + '_terms.dump', 'wb'))
    return go_dict


def convert_exp_to_entrez(exp_file_location):
    filenames = glob(exp_file_location)
    all_genes = set()
    for filename in filenames:
        genes = load(open(filename, 'rb')).index
        all_genes.update(genes)
    gene_entrez_dict = convert_genes_to_entrez(all_genes)
    return gene_entrez_dict


def convert_meth_to_entrez(meth_file_location, exp_file_location=None, R_file_locations=None):
    """
    Optionally, provide the location of expression data and a matrix R that maps
    expression to methylation by index. This ensures that genes known a priori
    to be related to the methylation probe are included in the set of associated
    entrez IDs.
    """
    meth_filenames = sorted(glob(meth_file_location))
    exp_filenames = sorted(glob(exp_file_location))
    if R_file_locations is not None:
        R_filenames = sorted(glob(R_file_location))

    meth_gene_map = {}
    meth_genes = []

    for i in range(len(meth_filenames)):
        meth_filename = meth_filenames[i]
        exp_filename = exp_filenames[i]
        file_genes = load(open(exp_filename, 'rb')).index
        if R_file_locations is not None:
            R_filename = R_filenames[i]
            R = load(open(R_filename, 'rb'))
        probes = load(open(meth_filename, 'rb')).index
        hg19_meth_annotations = pd.read_csv(open(METH_ANNOTATION_FILE), keep_default_na=False, sep='\t')
        hg19_meth_map = dict(hg19_meth_annotations[hg19_meth_annotations['probeID'].isin(probes)].loc[:, ['probeID', 'gene']].values)

        for j, probe in enumerate(probes):
            if meth_gene_map.get(probe) is not None:
                continue
            comb_related_genes = hg19_meth_map[probe]
            related_genes = comb_related_genes.split(';')
            if exp_file_location is not None:
                related_genes.extend([file_genes[k] for k in np.nonzero(R.T[j, :])[0]])
                related_genes = list(set(related_genes))
            meth_gene_map[probe] = related_genes
            meth_genes.extend(related_genes)

    gene_entrez_dict = convert_genes_to_entrez(set(meth_genes))
    probe_entrez_dict = defaultdict(list)
    for probe in meth_gene_map.keys():
        for gene in meth_gene_map[probe]:
            probe_entrez_dict[probe].extend(gene_entrez_dict[gene])
        probe_entrez_dict[probe] = list(set(probe_entrez_dict[probe]))

    return probe_entrez_dict


def convert_genes_to_entrez(genes):
    mg = mygene.MyGeneInfo()
    responses = mg.querymany(genes, scopes='symbol', species='9606', entrezonly=True)
    gene_name_dict = defaultdict(list)
    for response in responses:
        gene_name = response['query']
        entrez_id = response.get('entrezgene')
        if entrez_id is not None:
            gene_name_dict[gene_name].append(entrez_id)

    for gene_name in genes:
        if gene_name_dict.get(gene_name) is None:
            print(gene_name)
            response = mg.query('alias:' + gene_name, species='9606', entrezonly=True)
            print(response)
            if response['total'] != 0:
                gene_name_dict[gene_name] = [response['hits'][i]['entrezgene'] for i in range(len(response['hits']))]
        if gene_name_dict.get(gene_name) is None:
            gene_name_dict[gene_name] = []
    return gene_name_dict


if __name__ == "__main__":
    # example use below
    view1 = 'exp_intrinsic'
    view2 = 'meth_intrinsic'

    try:
        view1_entrez_dict = load(open("data/" + view1 + "_symbol_to_entrez.dump", "rb"))
        view2_entrez_dict = load(open("data/" + view2 + "_symbol_to_entrez.dump", "rb"))
    except:
        exp_file_location = '../data/TCGA_exp_*_intrinsic.dump'
        meth_file_location = '../data/TCGA_meth_*_intrinsic.dump'
        R_file_location = '../data/TCGA_R_*_intrinsic.dump'

        view1_entrez_dict = convert_exp_to_entrez(exp_file_location)
        dump(view1_entrez_dict, open("data/exp_intrinsic_symbol_to_entrez.dump", "wb"))
        view2_entrez_dict = convert_meth_to_entrez(meth_file_location, exp_file_location, R_file_location)
        dump(view2_entrez_dict, open("data/meth_intrinsic_symbol_to_entrez.dump", "wb"))

    goto_type = '0'

    get_terms(view1_entrez_dict, view1, 'BP', goto_type)
    get_terms(view2_entrez_dict, view2, 'BP', goto_type)

    get_terms(view1_entrez_dict, view1, 'CC', goto_type)
    get_terms(view1_entrez_dict, view1, 'MF', goto_type)
    get_terms(view1_entrez_dict, view1, 'all', goto_type)

    get_terms(view2_entrez_dict, view2, 'CC', goto_type)
    get_terms(view2_entrez_dict, view2, 'MF', goto_type)
    get_terms(view2_entrez_dict, view2, 'all', goto_type)
