#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 10:08:12 2018

@author: admin = Kris
Edited by Masaki
"""

#learner corpus analysis
import glob
import pickle
import spacy

nlp = spacy.load('en_core_web_md')


def spcy_process(flnm):
    out_dict = {
        "lem_bi": [],
        "v_nsubj": [],
        "v_dobj": [],
        "v_advmod": [],
        "n_amod": [],
        "n_nn": []
    }
    doc = nlp(open(flnm, "r", encoding='utf-8-sig').read())
    for sent in doc.sents:
        prev_lemm = "<START>"
        for token in sent:  #iterate through tokens
            if token.lemma_ == "-PRON-":
                lemma = token.text.lower()
            else:
                lemma = token.lemma_
            if token.pos_ == "PUNCT":
                continue

            if prev_lemm != "<START>":
                out_dict["lem_bi"].append(prev_lemm + "_" + lemma)
            prev_lemm = lemma

            if token.dep_ == "nsubj":  #if the item is a direct object:
                out_dict["v_nsubj"].append(token.head.lemma_ + "_" +
                                           lemma)  #main verb _ subject

            #verb_dobj
            if token.dep_ == "dobj":  #if the item is a direct object:
                out_dict["v_dobj"].append(token.head.lemma_ + "_" +
                                          lemma)  #main verb _ direct object

            if (token.dep_ == "advmod"
                    or token.dep_ == "npmod") and token.head.pos_ == "VERB":
                out_dict["v_advmod"].append(
                    token.head.lemma_ + "_" +
                    token.lemma_)  #main verb _ adverb modifier

            if token.dep_ == "amod":
                out_dict["n_amod"].append(
                    token.head.lemma_ + "_" +
                    token.lemma_)  #head noun _ #adjective modifier

            if token.dep_ == "compound" and token.pos_ == "NOUN" and token.head.pos_ == "NOUN":
                out_dict["n_nn"].append(token.head.lemma_ + "_" +
                                        token.lemma_)  #head noun + modifier
    return (out_dict)


def index_calc(item_l, item_d, index):
    denom = len(item_l)
    if denom == 0:
        outvar = 0
    else:
        numerator = 0
        for item in item_l:
            if item in item_d:
                numerator += item_d[item][index]

        outvar = numerator / denom

    return (str(outvar))


def item_output(item_l, item_d):  #added by Masaki for individual input
    item_dict = {}
    for item in item_l:
        index_list = []
        if item in item_d:
            #print(item)
            #print(item_d[item])
            index_list.extend(item_d[item])
        item_dict[item] = index_list
    #print(item_dict)
    return item_dict


def index_calc_strongest(item_l, item_d, index1, index2):
    denom = len(item_l)
    if denom == 0:
        outvar = 0
    else:
        numerator = 0
        for item in item_l:
            if item in item_d:
                numerator += max([item_d[item][index1], item_d[item][index2]])

        outvar = numerator / denom

    return (str(outvar))


def indexer(index_value, index_name, index_list, name_list):
    index_list.append(index_value)
    name_list.append(index_name)


lem_bi_soa = pickle.load(
    open("1_reference_norms/dep_lists-2019-11-13/lem_bi_soa.pickle", "rb"))
v_nsubj_soa = pickle.load(
    open("1_reference_norms/dep_lists-2019-11-13/v_nsubj_soa.pickle", "rb"))
v_dobj_soa = pickle.load(
    open("1_reference_norms/dep_lists-2019-11-13/v_dobj_soa.pickle", "rb"))
v_advmod_soa = pickle.load(
    open("1_reference_norms/dep_lists-2019-11-13/v_advmod_soa.pickle", "rb"))
n_amod_soa = pickle.load(
    open("1_reference_norms/dep_lists-2019-11-13/n_amod_soa.pickle", "rb"))
n_nn_soa = pickle.load(
    open("1_reference_norms/dep_lists-2019-11-13/n_nn_soa.pickle", "rb"))


def write_individual(dict_name, fname, dep_type, outname,
                     meta):  #Added by Masaki (May18 2019)
    for item in dict_name:
        fname.write("\n" + str(outname) + "\t" + "\t".join(meta[0:7]) + '\t' +
                    str(item) + "\t" + str(dep_type))
        for index in dict_name[item]:
            fname.write("\t" + str(index))


filenames = glob.glob(
    "1_corpora/2_learner_corpora/ICNALE_Edited Essays_2.1/EE_Unmerged_Unclassified/*.txt"
)
L2_edited = glob.glob(
    "1_corpora/2_learner_corpora/ICNALE_Edited Essays_2.1/EE_Unmerged_Unclassified/*.txt"
)
L2_essays = glob.glob(
    "1_corpora/2_learner_corpora/ICNALE_Written_Essays_2.3_3/Unmerged_classified/*/*.txt"
)
L1_essays = glob.glob(
    "1_corpora/2_learner_corpora/ICNALE_Written_Essays_2.3/L1/*/*.txt")

filenames = L2_edited + L1_essays

outf = open("9_collocation_prop_ICNALE/Col_data/ICNALE_L1_average_200922.csv",
            "w")
for idx, filename in enumerate(filenames):
    print("Processing " + str(idx + 1) + " of " + str(len(filenames)) +
          " files")
    outname = filename.split("/")[-1]
    meta_data = outname.split("_")
    values = [outname]
    if "ENS_" in filename:
        meta_data.append("Orig.txt")
    values.extend(meta_data)

    headers = [
        "filename", "#W", "Country", "Topic", "ID", "Level", "?", "Edited"
    ]
    ld = spcy_process(filename)

    indexer(index_calc(ld["lem_bi"], lem_bi_soa, 0), "lemm_bg_deltap_w1cue",
            values, headers)
    indexer(index_calc(ld["lem_bi"], lem_bi_soa, 1), "lemm_bg_deltap_w2cue",
            values, headers)
    indexer(index_calc_strongest(ld["lem_bi"], lem_bi_soa, 0, 1),
            "lemm_bg_deltap_strgst", values, headers)
    indexer(index_calc(ld["lem_bi"], lem_bi_soa, 2), "lemm_bg_MI", values,
            headers)
    indexer(index_calc(ld["lem_bi"], lem_bi_soa, 3), "lemm_bg_T", values,
            headers)
    indexer(index_calc(ld["lem_bi"], lem_bi_soa, 5), "lemm_bg_freq", values,
            headers)
    indexer(index_calc(ld["lem_bi"], lem_bi_soa, 6), "lemm_freq", values,
            headers)

    indexer(index_calc(ld["v_nsubj"], v_nsubj_soa, 0), "v_nsubj_deltap_govcue",
            values, headers)
    indexer(index_calc(ld["v_nsubj"], v_nsubj_soa, 1), "v_nsubj_deltap_depcue",
            values, headers)
    indexer(index_calc_strongest(ld["v_nsubj"], v_nsubj_soa, 0, 1),
            "v_nsubj_deltap_strgst", values, headers)
    indexer(index_calc(ld["v_nsubj"], v_nsubj_soa, 2), "v_nsubj_MI", values,
            headers)
    indexer(index_calc(ld["v_nsubj"], v_nsubj_soa, 3), "v_nsubj_T", values,
            headers)
    indexer(index_calc(ld["v_nsubj"], v_nsubj_soa, 5), "v_nsubj_freq", values,
            headers)
    indexer(index_calc(ld["v_nsubj"], v_nsubj_soa, 6), "v_finite_freq", values,
            headers)
    indexer(index_calc(ld["v_nsubj"], v_nsubj_soa, 7), "nsubj_freq", values,
            headers)

    indexer(index_calc(ld["v_dobj"], v_dobj_soa, 0), "v_dobj_deltap_govcue",
            values, headers)
    indexer(index_calc(ld["v_dobj"], v_dobj_soa, 1), "v_dobj_deltap_depcue",
            values, headers)
    indexer(index_calc_strongest(ld["v_dobj"], v_dobj_soa, 0, 1),
            "v_dobj_deltap_strgst", values, headers)
    indexer(index_calc(ld["v_dobj"], v_dobj_soa, 2), "v_dobj_MI", values,
            headers)
    indexer(index_calc(ld["v_dobj"], v_dobj_soa, 3), "v_dobj_T", values,
            headers)
    indexer(index_calc(ld["v_dobj"], v_dobj_soa, 5), "v_dobj_freq", values,
            headers)
    indexer(index_calc(ld["v_dobj"], v_dobj_soa, 6), "v_transitive_freq",
            values, headers)
    indexer(index_calc(ld["v_dobj"], v_dobj_soa, 7), "dobj_freq", values,
            headers)

    indexer(index_calc(ld["v_advmod"], v_advmod_soa, 0),
            "v_advmod_deltap_govcue", values, headers)
    indexer(index_calc(ld["v_advmod"], v_advmod_soa, 1),
            "v_advmod_deltap_depcue", values, headers)
    indexer(index_calc_strongest(ld["v_advmod"], v_advmod_soa, 0, 1),
            "v_advmod_deltap_strgst", values, headers)
    indexer(index_calc(ld["v_advmod"], v_advmod_soa, 2), "v_advmod_MI", values,
            headers)
    indexer(index_calc(ld["v_advmod"], v_advmod_soa, 3), "v_advmod_T", values,
            headers)
    indexer(index_calc(ld["v_advmod"], v_advmod_soa, 5), "v_advmod_freq",
            values, headers)
    indexer(index_calc(ld["v_advmod"], v_advmod_soa, 6), "v_w_advmod_freq",
            values, headers)
    indexer(index_calc(ld["v_advmod"], v_advmod_soa, 7), "advmod_freq", values,
            headers)

    indexer(index_calc(ld["n_amod"], n_amod_soa, 0), "n_amod_deltap_govcue",
            values, headers)
    indexer(index_calc(ld["n_amod"], n_amod_soa, 1), "n_amod_deltap_depcue",
            values, headers)
    indexer(index_calc_strongest(ld["n_amod"], n_amod_soa, 0, 1),
            "n_amod_deltap_strgst", values, headers)
    indexer(index_calc(ld["n_amod"], n_amod_soa, 2), "n_amod_MI", values,
            headers)
    indexer(index_calc(ld["n_amod"], n_amod_soa, 3), "n_amod_T", values,
            headers)
    indexer(index_calc(ld["n_amod"], n_amod_soa, 5), "n_amod_freq", values,
            headers)
    indexer(index_calc(ld["n_amod"], n_amod_soa, 6), "n_w_amod_freq", values,
            headers)
    indexer(index_calc(ld["n_amod"], n_amod_soa, 7), "amod_freq", values,
            headers)

    indexer(index_calc(ld["n_nn"], n_nn_soa, 0), "n_nnmod_deltap_govcue",
            values, headers)
    indexer(index_calc(ld["n_nn"], n_nn_soa, 1), "n_nnmod_deltap_depcue",
            values, headers)
    indexer(index_calc_strongest(ld["n_nn"], n_nn_soa, 0, 1),
            "n_nnmod_deltap_strgst", values, headers)
    indexer(index_calc(ld["n_nn"], n_nn_soa, 2), "n_nnmod_MI", values, headers)
    indexer(index_calc(ld["n_nn"], n_amod_soa, 3), "n_nn_T", values, headers)
    indexer(index_calc(ld["n_nn"], n_amod_soa, 5), "n_nn_freq", values,
            headers)
    indexer(index_calc(ld["n_nn"], n_amod_soa, 6), "n_w_nn_freq", values,
            headers)
    indexer(index_calc(ld["n_nn"], n_amod_soa, 7), "nn_modifier_freq", values,
            headers)

    if idx == 0:
        outf.write(",".join(headers))

    outf.write("\n" + ",".join(values))

    #individual outpud from here Added by Masaki
    indf = open("9_collocation_prop_ICNALE/itemanalysis_indout.txt", "a")
    if idx == 0:
        indf.write(
            "filename\tW\tCountry\tTopic\tID\tLevel\t?\tEdited\tItem\tdep_type\tdp_gov\tdp_dep\tMI\tT\tFred\tgov_Freq\tdep_Freq"
        )

    lem_bi_dict = item_output(ld["lem_bi"], lem_bi_soa)
    v_nsubj_dict = item_output(ld["v_nsubj"], v_nsubj_soa)
    v_dobj_dict = item_output(ld["v_dobj"], v_dobj_soa)
    v_advmod_dict = item_output(ld["v_advmod"], v_advmod_soa)
    n_amod_dict = item_output(ld["n_amod"], n_amod_soa)
    n_nn_dict = item_output(ld["n_nn"], n_nn_soa)

    write_individual(lem_bi_dict, indf, "lem_bi", outname, meta_data)
    write_individual(v_nsubj_dict, indf, "v_nsubj", outname, meta_data)
    write_individual(v_dobj_dict, indf, "v_dobj", outname, meta_data)
    write_individual(v_advmod_dict, indf, "v_advmod", outname, meta_data)
    write_individual(n_amod_dict, indf, "n_amod", outname, meta_data)
    write_individual(n_nn_dict, indf, "n_nn", outname, meta_data)

    indf.flush()
    indf.close()

outf.flush()
outf.close()

print("Finished!")
