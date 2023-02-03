import numpy as np
import himut.gtlib
from typing import Dict, List, Tuple
gt_lst = ["AA", "TA", "CA", "GA", "TT", "CT", "GT", "CC", "GC", "GG"]



def get_argmin_gt(
    gt_pl_lst: List[float],
    gt2gt_state: Dict[str, str]
):
    ilst = np.argsort(gt_pl_lst)
    gt = [gt_lst[i] for i in ilst][0]
    gt_pl_lst = [gt_pl_lst[i] for i in ilst]
    gq = (np.array(gt_pl_lst) - min(gt_pl_lst))[1]
    gq = int(gq) if gq < 99 else 99
    gt_state = gt2gt_state[gt]
    return gt, gq, gt_state


def get_germ_gt(
    ref: str,
    allele2bq_lst: Dict[int, List[int]],
) -> Tuple[List[str], List[float], str]:

    gt_pl_lst, gt2gt_state = himut.gtlib.get_germ_gt_pD(
        ref,
        allele2bq_lst,
    )
    gt, gq, gt_state = get_argmin_gt(gt_pl_lst, gt2gt_state) 
    if gt_state == "het":
        vgt = "0/1" 
    elif gt_state == "hetalt":
        vgt = "1/2" 
    elif gt_state == "homalt":
        vgt = "1/1" 
    elif gt_state == "homref":
        vgt = "0/0" 
        
    if gt[0] != ref and gt.count(ref) == 1:
        gt = gt[::-1]
    return gt, gq, vgt, gt_state, 
