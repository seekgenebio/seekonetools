import re
from functools import partial, lru_cache
from xopen import xopen


def generate_mismatch_sequence(sequence: str, bases: str='ATCG') -> dict:
    # 构建错配1字典
    # [A]AAAAAAA -> T + AAAAAAA, C + AAAAAAA, G + AAAAAAA, N + AAAAAAA
    tmp = {}
    for i in range(len(sequence)):
        for c in bases:
            new = sequence[:i] + c + sequence[i + 1:]
            if new != sequence:
                tmp[new] = sequence
    return tmp

def generate_indel_sequence(sequence: str, bases: str='ATCG') -> dict:
    # 构建缺失1字典
    # [A]AAAAAAA -> AAAAAAA + A, AAAAAAA + T, AAAAAAA + C, AAAAAAA + G
    tmp = {}
    for i in range(len(sequence)):
        for c in bases:
            new = sequence[:i] + sequence[i + 1:] + c
            if new != sequence:
                tmp[new] = sequence
    return tmp

def parse_structure(string: str='B8L15B8L15B8U8T10') -> tuple:
    # 拆解barcode结构
    # B表示cellbarcode碱基
    # L表示linker碱基
    # U表示UMI碱基
    # T表示T碱基
    # X表示任意碱基
    regex = re.compile(r'([BLUXT])(\d+)')
    groups = regex.findall(string)
    return tuple([(_[0], int(_[1])) for _ in groups])

def correct_seq(seq: str, perfect: dict, indel: dict=None, mis: dict=None,
                seqafter: str=None, linker: str=None) -> tuple:
    # 0: barcode不能有效纠错
    # 1: barcode来自白名单，不用纠错
    # 2: 缺失碱基纠错
    # 3: 碱基错配纠错
    # -1: 有多种纠错可能
    code = 0
    shift = 0
    # print(seq)
    if seq in perfect:
        code = 1
    elif mis and seq in mis:
        code = 3
        seq = mis[seq]
    elif indel and seq in indel:
        if (not linker) or linker.startswith(seq[-1] + seqafter):
            code = 2
            shift = -1
            seq = indel[seq]
    # print(code, shift, seq)
    return code, shift, seq


def prepare_funcs(structure: str, _barcode: list=None, _linker: list=None,
                    B: tuple=(0, 0), L:tuple=(0, 0)) -> tuple:
    r1_structure = parse_structure(structure)

    do_correct = {'L': True, 'B': True}
    if not _linker:
        do_correct['L'] = False
    if not _barcode:
        do_correct['B'] = False

    funcs = []

    barcode = dict(enumerate(_barcode))
    linker = dict(enumerate(_linker))

    barcode_idx = 0
    linker_idx = 0
    for _, (code, n) in enumerate(r1_structure):
        if code == 'B':
            if do_correct['B']:
                tmp = build_dict(barcode.get(barcode_idx, barcode[0]), *B)
                if r1_structure[_+1][0] == 'L' and do_correct['L']:
                    funcs.append(
                        partial(correct_seq, linker=linker.get(linker_idx, linker[0]),
                            **tmp
                        )
                    )
                else:
                    funcs.append(
                        partial(correct_seq,
                            **tmp
                        )
                    )
            else:
                funcs.append(None)
            barcode_idx += 1
        elif code == 'L':
            if do_correct['L']:
                tmp = build_dict(linker.get(linker_idx, linker[0]), *L)
                funcs.append(
                    partial(correct_seq,
                        **tmp
                    )
                )
            else:
                funcs.append(None)
            linker_idx += 1
        else:
            funcs.append(None)
    return r1_structure, funcs

def check_keys(d1: dict, d2: dict) -> None:
    _check = d1.keys() & d2.keys()
    if _check:
        for k in _check:
            pass
            # print(f"{k}\t{d1[k]}\t{d2[k]}")


@lru_cache(10)
def build_dict(seq: str, err: int=0, indel: int=0) -> dict:
    perfect_dict = {}
    indel_dict = {}
    mismatch_dict = {}
    searchDict = {}
    try:
        with xopen(seq, 'r') as fh:
            for line in fh:
                line = line.strip()
                if line.startswith('#'):
                    continue
                if not line:
                    continue
                perfect_dict[line] = line
                if indel == 1:
                    tmp = generate_indel_sequence(line)
                    check_keys(tmp, indel_dict)
                    indel_dict.update(tmp)
                if err == 1:
                    tmp = generate_mismatch_sequence(line)
                    check_keys(tmp, mismatch_dict)
                    mismatch_dict.update(tmp)
    except Exception:
        perfect_dict[seq] = seq
        if indel == 1:
            tmp = generate_indel_sequence(seq)
            check_keys(tmp, indel_dict)
            indel_dict.update(tmp)
        if err == 1:
            tmp = generate_mismatch_sequence(seq)
            check_keys(tmp, mismatch_dict)
            mismatch_dict.update(tmp)
    searchDict = {'perfect': perfect_dict, 'indel': indel_dict, 'mis': mismatch_dict}
    return searchDict
