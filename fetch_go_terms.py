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


def get_terms(entrez_dict, name, namespace, goto_type):
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
            gosubdag = GoSubDag(goids, godag, relationships={'part_of',}, prt=False)
            # go_info = [gosubdag.go2nt.get(goid) for goid in goids]
            go2parents = get_go2parents(gosubdag.go2obj, gosubdag.relationships)
            go2parents_nr = list(set(go2parents.keys()))
            go2parents_nr_filtered = go2parents_nr
            if goto_type == '30200':
                go2parents_nr_filtered = [go for go in go2parents_nr if go2geneids_counter[go] > 30 and go2geneids_counter[go] < 200]
            elif goto_type == '0':
                go2parents_nr_filtered = [go for go in go2parents_nr if go2geneids_counter[go] > 0]
            if go_dict.get(entity) == ['NO GOIDS'] or go_dict.get(entity) == ['NO ENTREZ'] or go_dict.get(entity) is None:
                go_dict[entity] = [go2parents_nr_filtered]
            else:
                go_dict[entity].append(go2parents_nr_filtered)

    dump(go_dict, open('data/go_' + name + '_' + namespace + '_' + goto_type + '_terms.dump', 'wb'))
    return go_dict


def convert_exp_to_entrez(exp_file_location):
    genes = load(open(exp_file_location, 'rb')).index
    gene_entrez_dict = convert_genes_to_entrez(genes)
    return gene_entrez_dict


def convert_meth_to_entrez(meth_file_location, exp_file_location=None, R_file_location=None):
    if exp_file_location is not None:
        file_genes = load(open(exp_file_location, 'rb')).index
        R = load(open(R_file_location, 'rb'))
    probes = load(open(meth_file_location, 'rb')).index
    hg19_meth_annotations = pd.read_csv(open('data/TCGA/hm27_hg19_annotations.txt'), keep_default_na=False, sep='\t')
    hg19_meth_map = dict(hg19_meth_annotations[hg19_meth_annotations['probeID'].isin(probes)].loc[:, ['probeID', 'gene']].values)
    meth_gene_map = {}
    meth_genes = []
    i = 0
    for probe in probes:
        comb_related_genes = hg19_meth_map[probe]
        related_genes = comb_related_genes.split(';')
        if exp_file_location is not None:
            related_genes.extend([file_genes[j] for j in np.nonzero(R.T[i, :])[0]])
            related_genes = list(set(related_genes))
        meth_gene_map[probe] = related_genes
        meth_genes.extend(related_genes)
        i += 1
    gene_entrez_dict = convert_genes_to_entrez(set(meth_genes))
    probe_entrez_dict = defaultdict(list)
    for probe in probes:
        for gene in meth_gene_map[probe]:
            probe_entrez_dict[probe].extend(gene_entrez_dict[gene])
        probe_entrez_dict[probe] = list(set(probe_entrez_dict[probe]))
    return probe_entrez_dict


def convert_protexp_to_entrez():
    genes = pd.read_csv('/Users/benj/Documents/Research/Projects/bmvc/eval/data/pqtl/rna_norm_15_PCsremoved.txt', sep=None, index_col=0).index
    proteins = pd.read_csv('/Users/benj/Documents/Research/Projects/bmvc/eval/data/pqtl/pro_norm_9_PCsremoved.txt', sep=None, index_col=0).index
    gene_entrez_dict = convert_ensembl_to_entrez(genes)
    protein_entrez_dict = convert_ensembl_to_entrez(proteins)
    dump(gene_entrez_dict, open("data/pqtl_gene_to_entrez.dump", "wb"))
    dump(protein_entrez_dict, open("data/pqtl_protein_to_entrez.dump", "wb"))
    return gene_entrez_dict, protein_entrez_dict


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


def convert_ensembl_to_entrez(genes):
    mg = mygene.MyGeneInfo()
    responses = mg.querymany(genes, scopes='ensembl.gene', species='9606', entrezonly=True)
    gene_name_dict = defaultdict(list)
    for response in responses:
        gene_name = response['query']
        entrez_id = response.get('entrezgene')
        if entrez_id is not None:
            gene_name_dict[gene_name].append(entrez_id)

    for gene_name in genes:
        if gene_name_dict.get(gene_name) is None:
            gene_name_dict[gene_name] = []
    return gene_name_dict


if __name__ == "__main__":

    exp_file_location = '../data/TCGA_exp_intrinsic.dump'
    meth_file_location = '../data/TCGA_meth_intrinsic.dump'
    R_file_location = '../data/TCGA_R_intrinsic.dump'
    view1 = 'exp_intrinsic'
    view2 = 'meth_intrinsic'

    # view1 = 'exp_diffmeth'
    # view2 = 'meth_diffmeth'
    # view1 = 'exp_intrinsic'
    # view2 = 'meth_intrinsic'
    # view1 = 'pqtl_gene'
    # view2 = 'pqtl_protein'

    exp_entrez_dict = convert_exp_to_entrez(exp_file_location)
    dump(exp_entrez_dict, open("data/exp_symbol_to_entrez.dump", "wb"))
    meth_entrez_dict = convert_meth_to_entrez(exp_file_location, meth_file_location, R_file_location)
    dump(meth_entrez_dict, open("data/meth_symbol_to_entrez.dump", "wb"))

    # pqtl_gene_dict, pqt_protein_dict = convert_protexp_to_entrez()

    view1_entrez_dict = load(open("data/" + view1 + "_symbol_to_entrez.dump", "rb"))
    view2_entrez_dict = load(open("data/" + view2 + "_symbol_to_entrez.dump", "rb"))

    # pqtl_gene_to_entrez_dict = load(open("data/" + 'pqtl_gene' + "_to_entrez.dump", "rb"))
    # pqtl_protein_to_entrez_dict = load(open("data/" + 'pqtl_protein' + "_to_entrez.dump", "rb"))

    # get_terms(pqtl_gene_to_entrez_dict, 'pqtl_gene', 'BP', '0')
    # get_terms(pqtl_protein_to_entrez_dict, 'pqtl_protein', 'BP', '0')

    get_terms(exp_entrez_dict, view1, 'BP', '0')
    get_terms(meth_entrez_dict, view2, 'BP', '0')

    get_terms(exp_entrez_dict, view1, 'CC', '0')
    get_terms(exp_entrez_dict, view1, 'MF', '0')
    get_terms(exp_entrez_dict, view1, 'all', '0')

    get_terms(meth_entrez_dict, view2, 'CC', '0')
    get_terms(meth_entrez_dict, view2, 'MF', '0')
    get_terms(meth_entrez_dict, view2, 'all', '0')
