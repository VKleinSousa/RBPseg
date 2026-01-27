#!/usr/bin/env python3
import argparse
from typing import List, Tuple, Optional

from Bio.PDB import PDBParser, PDBIO


def _get_chain(model, chain_id: str):
    for ch in model:
        if ch.id == chain_id:
            return ch
    return None


def _max_resseq(chain) -> int:
    """Return max residue sequence number (resseq) among standard residues in chain."""
    max_id = 0
    for res in chain.get_residues():
        # residue.id is (hetflag, resseq, icode)
        resseq = res.id[1]
        if isinstance(resseq, int) and resseq > max_id:
            max_id = resseq
    return max_id


def merge_chain_pair(
    structure,
    dst_id: str,
    src_id: str,
    gap: int = 1,
    delete_src: bool = True,
    model_index: int = 0,
):
    """
    Merge src chain into dst chain within a given model:
      - append residues from src to dst
      - renumber appended residues after the last dst residue (optionally with a gap)
      - optionally delete src chain after merge
    """
    model = list(structure.get_models())[model_index]

    dst = _get_chain(model, dst_id)
    src = _get_chain(model, src_id)

    if dst is None:
        raise ValueError(f"Destination chain '{dst_id}' not found in model {model_index}.")
    if src is None:
        raise ValueError(f"Source chain '{src_id}' not found in model {model_index}.")
    if dst_id == src_id:
        raise ValueError("Destination and source chain IDs are the same.")

    start = _max_resseq(dst) + gap
    new_resseq = start

    # Copy residues from src -> dst with new residue ids
    # Use res.copy() to avoid parent/child conflicts in Bio.PDB object tree.
    for res in list(src.get_residues()):
        res_copy = res.copy()
        hetflag, _, icode = res_copy.id

        # Set a clean, sequential residue number; keep hetflag, wipe insertion code by default.
        # If you want to preserve insertion codes, remove the "icode=' '" line.
        icode = " "
        res_copy.id = (hetflag, new_resseq, icode)
        new_resseq += 1

        dst.add(res_copy)

    # Renumber ALL residues in dst contiguously starting from 1 (optional but usually desired)
    # Comment this block out if you want to keep original dst numbering.
    reseqq = 1
    for res in list(dst.get_residues()):
        hetflag, _, icode = res.id
        res.id = (hetflag, reseqq, " ")
        reseqq += 1

    if delete_src:
        model.detach_child(src.id)


def parse_pairs(args) -> List[Tuple[str, str]]:
    pairs: List[Tuple[str, str]] = []

    if args.pairs:
        # format: "E,B;F,B;G,B"
        chunks = [c.strip() for c in args.pairs.split(";") if c.strip()]
        for c in chunks:
            a, b = [x.strip() for x in c.split(",")]
            pairs.append((a, b))

    if args.pair:
        # repeated --pair A,B
        for p in args.pair:
            a, b = [x.strip() for x in p.split(",")]
            pairs.append((a, b))

    if not pairs:
        raise ValueError("No pairs provided. Use --pair A,B (repeatable) or --pairs 'A,B;C,D'.")

    # basic validation
    for dst, src in pairs:
        if len(dst) != 1 or len(src) != 1:
            raise ValueError(f"Chain IDs must be 1 character. Got: '{dst}', '{src}'")
    return pairs


def main():
    ap = argparse.ArgumentParser(description="Merge chains in a PDB based on given (dst,src) chain pairs.")
    ap.add_argument("--pdb", required=True, help="Input PDB file")
    ap.add_argument("--out", required=True, help="Output PDB file")
    ap.add_argument("--pair", action="append",
                    help="One merge pair 'DST,SRC' (repeatable). Example: --pair B,E --pair B,F")
    ap.add_argument("--pairs",
                    help="Multiple merge pairs as 'DST,SRC;DST,SRC'. Example: --pairs 'B,E;B,F;B,G'")
    ap.add_argument("--gap", type=int, default=1,
                    help="Gap between last residue in dst and first appended residue before final renumber (default: 1)")
    ap.add_argument("--keep-src", action="store_true",
                    help="Do NOT delete source chains after merging")
    ap.add_argument("--model-index", type=int, default=0,
                    help="Which model to use if PDB has multiple MODELs (default: 0)")
    args = ap.parse_args()

    pairs = parse_pairs(args)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", args.pdb)

    # Apply merges in order
    for dst, src in pairs:
        merge_chain_pair(
            structure,
            dst_id=dst,
            src_id=src,
            gap=args.gap,
            delete_src=not args.keep_src,
            model_index=args.model_index,
        )

    io = PDBIO()
    io.set_structure(structure)
    io.save(args.out)
    print(f"Wrote merged PDB: {args.out}")


if __name__ == "__main__":
    main()

