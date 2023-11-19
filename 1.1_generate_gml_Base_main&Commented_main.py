"""
generate xml file for Gephi

Returns:
    None: create new gml file
"""
import os
from collections import defaultdict
import networkx as nx
import pandas as pd
import numpy as np


def generate_paper_mesh_dict(data:pd.DataFrame, is_major_topic=True) -> dict:
    """generate {paperID:[MeshTopicList]}

    Args:
        data (pd.DataFrame): pandas dataframe with 3 columns: PMID, Major_YN, DescriptorName_UI
        is_major_topic (bool, optional): filter major mesh headings or NOT. Defaults to True.

    Returns:
        dict: {paperID:[MeshTopicList]}
    """
    res = {} 
    if is_major_topic:
        data = data[data["DescriptorName_MajorTopicYN"] == "Y"][["PMID", "DescriptorName_UI"]]
    else:
        data = data[["PMID", "DescriptorName_UI"]]
    for name, grouped in data.groupby("PMID"):
        res[name] = grouped["DescriptorName_UI"].tolist()
    return res


def create_mesh_info(mesh_file:str) -> dict:
    """create mesh info dict, matching uid, name & index in adjacent matrix

    Args:
        mesh_file (str): mesh file path

    Returns:
        dict: {'D000001': {'index': 0, 'name': 'Calcimycin'},...}
    """
    res = {}
    df_mesh = pd.read_csv(mesh_file)
    for i, uid in enumerate(df_mesh.UI.tolist()):
        res[uid] = {"index":i, "name":df_mesh.loc[i,"NAME"]}
    return res


def create_adjacent_matrix(n:int) -> np.array:
    """create n*n adjacent matrix

    Args:
        n (int): nodes cnt

    Returns:
        np.array: numpy array
    """
    return np.zeros((n,n), dtype=np.int32)

def generate_gml(paper_meshtopic:dict, matrix:np.ndarray, mesh:dict, gml_path:str) -> None:
    """generate gml file using paper_meshTopic dict

    Args:
        paper_meshtopic (dict): {paper:[meshTopic]}
        matrix (np.ndarray): adjacent matrix
        gml_path (str): file path
    """
    # update adjacent matrix
    for one in paper_meshtopic.values():
        if len(one) == 1:
            continue
        for i in range(len(one)):
            i_word_index = mesh[one[i]]["index"]
            for j in range(i+1, len(one)):
                j_word_index = mesh[one[j]]["index"]
                matrix[i_word_index][j_word_index] += 1
    # create network according to adjacent matrix
    G = nx.from_numpy_array(matrix)
    # add nodes info using mesh_info_dict
    for uid, attr in mesh.items():
        node_index, node_name = attr["index"], attr["name"]
        G.add_node(node_index, ui=uid, name=node_name)
    nx.write_gml(G, gml_path)

if __name__ == "__main__":
    # check-outs
    def avg_topic(paper_mesh_topic:defaultdict) -> float:
        """calculate average number of main topic mesh headings

        Args:
            paper_mesh_topic (defaultdict): {paper_id:[topic1, ...], ...}

        Returns:
            float: mean value for topic counts
        """
        return np.mean([len(i) for i in paper_mesh_topic.values()])
    
    # main function
    os.chdir("Project/Mesh_Clustering")
    # with open("PKG23_A06_MeshHeadingList.csv", encoding="UTF-8") as f:
    #     next(f)  # Skip column names
    #     all_papers = f.readlines()
    df = pd.read_csv("PKG23_A06_MeshHeadingList.csv", usecols=["PMID", "DescriptorName_MajorTopicYN", "DescriptorName_UI"])
    mesh_info = create_mesh_info("MeSH_2023.csv")
    
    
    # ----------------------------Base Main-----------------------------------
    paper_meshTopic = generate_paper_mesh_dict(df, is_major_topic=True)
    empty_matrix = create_adjacent_matrix(30454)
    
    generate_gml(paper_meshTopic, empty_matrix, mesh_info, "Base_main.gml")
    print("Base_main.gml Finished")
    # -------------------------------------------------------------------------
    
    # # ----------------------------Base All-----------------------------------
    # paper_meshTopic = generate_paper_mesh_dict(df, is_major_topic=False)
    # empty_matrix = create_adjacent_matrix(30454)
    
    # generate_gml(paper_meshTopic, empty_matrix, mesh_info, "Base_all_0414.gml")
    # print("Base_all_0414.gml Finished")
    # # -------------------------------------------------------------------------
    
    
    # ------------------ merge commented papers ---------------------
    df_commented = pd.read_csv("2_PubMed_comted_pmid_220916.csv")
    df_commented.columns = ["PMID"]
    df = pd.merge(left=df, right=df_commented, on="PMID", how="right")
    # ---------------------------------------------------------------
    
    # ----------------------------Commented Main-----------------------------------
    paper_meshTopic = generate_paper_mesh_dict(df, is_major_topic=True)
    empty_matrix = create_adjacent_matrix(30454)
    
    generate_gml(paper_meshTopic, empty_matrix, mesh_info, "Commented_main.gml")
    print("Commented_main.gml Finished")
    # ------------------------------------------------------------------------------
    
    # # ----------------------------Commented All-----------------------------------
    # paper_meshTopic = generate_paper_mesh_dict(df, is_major_topic=False)
    # empty_matrix = create_adjacent_matrix(30454)
    
    # generate_gml(paper_meshTopic, empty_matrix, mesh_info, "Commented_all_0414.gml")
    # print("Commented_all_0414.gml Finished")
    # # ------------------------------------------------------------------------------