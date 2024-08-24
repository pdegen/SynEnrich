import pandas as pd
from collections import defaultdict
from typing import List, Dict, Set, Any, DefaultDict

def get_sig_dict(summary_df: pd.DataFrame, tools: List[str], metrics: List[str], qval: float) -> Dict[Any, Any]:

    sig_dict: Dict = {t: {m: set() for m in metrics} for t in tools}

    for tool in tools:
        for metric in metrics:
            terms = summary_df[(tool,metric,"qvalue")].dropna()
            sig = terms[terms<qval]
            print(tool, metric, "Terms tested:", len(terms), "Significant:",  len(sig))
            sig_dict[tool][metric] = set((sig.index + " | " + summary_df.loc[sig.index,("nan","nan","Description")]))

    return sig_dict

def count_combinations(d) -> int:
    if not isinstance(d, Dict):
        return 1
    return sum(count_combinations(v) for v in d.values())

def create_intersection_depth_df(nested_dict: Dict[Dict, Set]) -> pd.DataFrame:
    """
    Counts intersection depth of each element of a nested dict containing sets to be compared. Each level of the dict correpsodns to e.g. experimental factor.
    :param nested_dict: Nested dict of factors where the deepest level contains sets
    :returns: pd.DataFrame with index = union of set elements.

    Depth: The intersection depth (number of sets a given element appears in)
    RelativeDepth: Depth divided by total number of factor combinations
    Factors: String listing all the factors a given element appears in
    """

    intersection_depth: DefaultDict = defaultdict(int)
    factors_for_element: DefaultDict = defaultdict(list)
    def calculate_depth_factors(nested_dict, factor_list):
        for factor, sub_dict in nested_dict.items():
            new_factors = factor_list + [factor]
            if isinstance(sub_dict, set):
                for elem in sub_dict:
                    intersection_depth[elem] += 1
                    factors_for_element[elem].append(new_factors)
            else:
                calculate_depth_factors(sub_dict, new_factors)
    calculate_depth_factors(nested_dict, [])

    # Cleanup
    for f in factors_for_element:
        for i, ll in enumerate(factors_for_element[f]):
            factors_for_element[f][i] = ".".join(ll)
        factors_for_element[f] = " | ".join(factors_for_element[f])

    depth_df = pd.DataFrame(intersection_depth.values(), index=intersection_depth.keys(), columns=["Depth"])
    depth_df.sort_values(by="Depth", ascending=False, inplace=True)
    factors_df = pd.DataFrame(factors_for_element.values(), index=factors_for_element.keys(), columns=["Factors"])
    depth_df = pd.concat([depth_df,factors_df], axis=1)

    combos = count_combinations(nested_dict)
    depth_df["RelativeDepth"] = depth_df["Depth"] / combos

    if " | " in depth_df.index[0]:
        depth_df["Description"] = depth_df.index.str.split("|").str[1]
        depth_df.index = depth_df.index.str.split("|").str[0]
        return depth_df[["Description","Depth","RelativeDepth","Factors"]]
    
    return depth_df[["Depth","RelativeDepth","Factors"]]
    