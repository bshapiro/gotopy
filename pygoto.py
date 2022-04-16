from collections import defaultdict
from itertools import combinations
from pickle import load
from copy import copy
import numpy as np


def compute_goto(entities, no_rel_entities, z, k, namespace, view, goto_type='0'):
    entity_go_dict = load(open('data/go_' + view + '_' + namespace + '_' + goto_type + '_terms.dump', 'rb'))
    cluster_goto_scores = []
    goto_dict = defaultdict(lambda: 0)
    cluster_sizes = []
    num_clusters = copy(k)
    for cluster in range(k):
        entity_indices = np.where(z == cluster)
        cluster_entities = entities[entity_indices]
        cluster_goto, cluster_size = compute_goto_cluster(cluster_entities, namespace, entity_go_dict)
        if cluster_goto is not None:
            goto_dict[cluster] = cluster_goto
            cluster_goto_scores.append(cluster_goto)
            cluster_sizes.append(cluster_size)
        else:
            goto_dict[cluster] = 0
            cluster_goto_scores.append(0)
            cluster_sizes.append(0)
            num_clusters = num_clusters - 1  # remove cluster without pairs
    # print(cluster_goto_scores)
    cluster_weights = [cluster_sizes[i] / float(np.sum(cluster_sizes)) for i in range(k)]
    goto_weight_dict = dict(zip([i for i in range(k)], cluster_weights))

    weighted_goto = np.dot(cluster_goto_scores, cluster_weights)
    return weighted_goto, goto_dict, goto_weight_dict


def compute_goto_cluster(entities, namespace, entity_go_dict):
    sum_intersect = 0
    entities_final = []
    for entity in entities:
        terms_list = entity_go_dict.get(entity)
        if not isinstance(terms_list[0], str):
            entities_final.append(entity)
    n_entities_final = len(entities_final)
    n_comparisons_final = n_entities_final * (n_entities_final - 1) / 2.0
    if n_comparisons_final == 0:
        return None, 0
    for (entity_i, entity_j) in combinations(entities_final, 2):
        pairwise_goto = compute_pairwise_goto(entity_i, entity_j, entity_go_dict)
        sum_intersect += pairwise_goto
    if n_comparisons_final != 0:
        return sum_intersect / n_comparisons_final, n_entities_final
    else:
        return None, 0


def compute_pairwise_goto(entity_i, entity_j, entity_go_dict):
    cur_intersections = []
    no_entrez_or_go_entities = set()
    terms_i_list = entity_go_dict.get(entity_i)
    terms_j_list = entity_go_dict.get(entity_j)

    if isinstance(terms_i_list[0], str) or isinstance(terms_j_list[0], str):
        if isinstance(terms_i_list[0], str):
            if terms_i_list[0] == 'NO ENTREZ' or terms_i_list[0] == "NO GOIDS":
                no_entrez_or_go_entities.add(entity_i)
        if isinstance(terms_j_list[0], str):
            if terms_j_list[0] == 'NO ENTREZ' or terms_j_list[0] == "NO GOIDS":
                no_entrez_or_go_entities.add(entity_j)
        return "NO ENTREZ OR GO"
    for terms_i in terms_i_list:
        for terms_j in terms_j_list:
            num_intersect = len(set(terms_i).intersection(set(terms_j)))
            # print(terms_i, terms_j, num_intersect)
            if num_intersect != 0:
                num_intersect = num_intersect
                # num_intersect = num_intersect / float(min(len(set(terms_i)), len(set(terms_j)))) # this is normalized GOTO
                # num_intersect = 1  # this is BHI
            cur_intersections.append(num_intersect)
    else:
        return np.mean(cur_intersections)
