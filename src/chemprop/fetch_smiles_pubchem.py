
from pathlib import Path
import argparse
import time
import logging
import requests
import pandas as pd
import urllib.parse
import concurrent.futures
import random
import json
import re
import sys

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
DEFAULT_DELAY = 0.25
MAX_WORKERS = 3
RETRY_LIMIT = 3
PUBCHEM_CACHE = Path("pubchem_name_cache.json")

BRD_BASE_RE = re.compile(r'^(BRD-[A-Za-z0-9]+)', re.IGNORECASE)

def canonical_brd_variants(brd):

    if brd is None:
        return []
    s = str(brd).strip()
    if not s:
        return []
    variants = set()

    variants.add(s)
    variants.add(s.upper())
    variants.add(s.lower())
    variants.add(s.replace('-', '_'))
    variants.add(s.replace('_', '-'))

    m = BRD_BASE_RE.match(s)
    if m:
        base = m.group(1)
        variants.add(base)
        variants.add(base.upper())
        variants.add(base.lower())
        variants.add(base.replace('-', '_'))

    parts = s.split('-')
    if len(parts) >= 2:
        base2 = '-'.join(parts[:2])
        variants.add(base2)
        variants.add(base2.upper())
        variants.add(base2.lower())
        variants.add(base2.replace('-', '_'))
    return list(variants)

def read_table_try(path):
    try:
        return pd.read_csv(path, sep='\t', low_memory=False)
    except Exception:
        try:
            return pd.read_csv(path, sep=',', low_memory=False)
        except Exception as e:
            raise RuntimeError(f"Failed to read {path}: {e}")

def find_repurposing_file():
    cwd = Path('.')
    candidates = list(cwd.glob('repurpos*')) + list(cwd.glob('*repurpos*.txt')) + list(cwd.glob('repurposing*')) + list(cwd.glob('*.txt')) + list(cwd.glob('*.tsv')) + list(cwd.glob('*.csv'))
    for p in candidates:
        try:
            hdr = p.open('r', encoding='utf-8', errors='ignore').readline().lower()
            if 'broad_id' in hdr or 'broadid' in hdr:
                return p
        except Exception:
            continue
    return None

def detect_smiles_col(df):
    for c in df.columns:
        if 'smiles' in c.lower() or 'canonical' in c.lower():
            return c
    return None

def build_repurposing_map(repurposing_path):

    if not repurposing_path or not Path(repurposing_path).exists():
        logging.warning("Repurposing file not found.")
        return {}, {}
    rep = read_table_try(repurposing_path)
    rep.columns = [c.strip() for c in rep.columns]

    id_col = None
    smiles_col = None
    name_col = None
    for c in rep.columns:
        cl = c.lower()
        if cl in ('broad_id','broadid','pert_id','pertid'):
            id_col = c
            break
    for c in rep.columns:
        if 'smiles' in c.lower() or 'canonical' in c.lower():
            smiles_col = c
            break
    for c in rep.columns:
        if 'pert_iname' in c.lower() or c.lower() == 'pert_iname' or c.lower() == 'pert_iname'.lower() or 'name' == c.lower():
            name_col = c
            break
    if id_col is None:
        logging.warning("Repurposing file missing expected id column (broad_id/pert_id).")
        return {}, {}
    if smiles_col is None:
        logging.warning("Repurposing file missing smiles-like column.")
    rep_map = {}
    name_map = {}
    for _, row in rep.iterrows():
        bid = str(row.get(id_col, '')).strip()
        smi = row.get(smiles_col) if smiles_col else None
        if pd.isna(smi):
            smi = None
        if smi:
            smi = str(smi).strip()
            if not smi:
                smi = None

        for v in canonical_brd_variants(bid):
            rep_map[v] = smi

        if name_col:
            nm = row.get(name_col)
            if not pd.isna(nm) and smi:
                name_map.setdefault(str(nm).strip().lower(), smi)
    logging.info(f"Built repurposing map for ~{len(rep_map)} BRD variants and {len(name_map)} names")
    return rep_map, name_map

def build_pertinfo_maps(pert_info_path):

    if not pert_info_path or not Path(pert_info_path).exists():
        logging.info("No pert_info provided.")
        return {}, {}
    pert = read_table_try(pert_info_path)
    pert.columns = [c.strip() for c in pert.columns]
    id_col = None
    name_col = None
    smiles_col = None
    for c in pert.columns:
        if c.lower() in ('pert_id','pertid','broad_id','broadid'):
            id_col = c
            break
    for c in pert.columns:
        if 'pert_iname' in c.lower() or 'iname' in c.lower() or c.lower() == 'name':
            name_col = c
            break
    for c in pert.columns:
        if 'smiles' in c.lower() or 'canonical' in c.lower():
            smiles_col = c
            break
    if id_col is None:
        logging.warning("pert_info missing id column; skipping pert_info mapping")
        return {}, {}
    pert_map = {}
    pert_name_map = {}
    for _, r in pert.iterrows():
        bid = str(r.get(id_col, '')).strip()
        if not bid:
            continue
        smi = r.get(smiles_col) if smiles_col else None
        if pd.isna(smi):
            smi = None
        if smi:
            smi = str(smi).strip()
            if not smi:
                smi = None
        name = r.get(name_col) if name_col else None
        if not pd.isna(name):
            name = str(name).strip()
        for v in canonical_brd_variants(bid):
            if smi:
                pert_map[v] = smi
            if name:
                pert_name_map[v] = name
    logging.info(f"Built pert_info maps: {len(pert_map)} id->smiles, {len(pert_name_map)} id->name")
    return pert_map, pert_name_map

def load_pubchem_cache():
    if PUBCHEM_CACHE.exists():
        try:
            return json.loads(PUBCHEM_CACHE.read_text(encoding='utf-8'))
        except Exception:
            return {}
    return {}

def save_pubchem_cache(cache):
    try:
        PUBCHEM_CACHE.write_text(json.dumps(cache, indent=2), encoding='utf-8')
    except Exception:
        pass

def fetch_smiles_by_name(name, session=None, timeout=12):

    if not name or str(name).strip() in ('-666','nan','None',''):
        return None
    s = session or requests.Session()
    q = urllib.parse.quote(str(name), safe='')

    url = f"{PUBCHEM_BASE}/compound/name/{q}/property/CanonicalSMILES/JSON"
    try:
        r = s.get(url, timeout=timeout)
        if r.status_code == 200:
            js = r.json()
            props = js.get('PropertyTable', {}).get('Properties', [])
            if props:
                smi = props[0].get('CanonicalSMILES')
                if smi:
                    return smi
    except Exception:
        pass

    try:
        url_cid = f"{PUBCHEM_BASE}/compound/name/{q}/cids/TXT"
        r2 = s.get(url_cid, timeout=timeout)
        if r2.status_code == 200 and r2.text.strip():
            first_cid = r2.text.strip().splitlines()[0].strip()
            url_by_cid = f"{PUBCHEM_BASE}/compound/cid/{first_cid}/property/CanonicalSMILES/JSON"
            r3 = s.get(url_by_cid, timeout=timeout)
            if r3.status_code == 200:
                js = r3.json()
                props = js.get('PropertyTable', {}).get('Properties', [])
                if props:
                    smi = props[0].get('CanonicalSMILES')
                    if smi:
                        return smi
    except Exception:
        pass
    return None

def pubchem_worker(name, delay, user_agent, retries=RETRY_LIMIT):
    s = requests.Session()
    s.headers.update({'User-Agent': user_agent})
    backoff = 0.5
    for attempt in range(retries):
        try:
            smi = fetch_smiles_by_name(name, session=s, timeout=12)
            if smi:
                return name, smi
        except Exception:
            pass
        time.sleep(delay + backoff * attempt + random.random() * 0.1)
    return name, None

def threaded_pubchem_query(names, delay=DEFAULT_DELAY, max_workers=MAX_WORKERS):
    uniq = [n for n in dict.fromkeys(names) if n and str(n).strip()]
    if not uniq:
        return {}
    logging.info(f"PubChem fallback: querying {len(uniq)} names with {max_workers} workers")
    results = {}
    user_agent = 'scz_connectivity/1.0 (contact: your-email@example.com)'
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = {ex.submit(pubchem_worker, name, delay, user_agent): name for name in uniq}
        for fut in concurrent.futures.as_completed(futures):
            name = futures[fut]
            try:
                n, smi = fut.result()
                results[n] = smi
                if smi:
                    logging.info(f"PubChem: found SMILES for {n}: {smi[:120]}")
            except Exception as e:
                logging.debug(f"PubChem worker failed for {name}: {e}")
                results[name] = None
    return results

def main(args):
    inp = Path(args.input)
    if not inp.exists():
        logging.error("Input dataset not found: %s", inp)
        sys.exit(1)

    logging.info(f"Loading dataset: {inp}")
    df = pd.read_csv(inp, sep='\t', low_memory=False)
    if 'BRD_ID' not in df.columns:
        logging.error("Input dataset missing 'BRD_ID' column.")
        sys.exit(1)

    rep_file = Path(args.repurposing) if args.repurposing else find_repurposing_file()
    rep_map, rep_name_map = ({}, {})
    if rep_file:
        logging.info(f"Using repurposing file: {rep_file}")
        rep_map, rep_name_map = build_repurposing_map(rep_file)
    else:
        logging.warning("No repurposing file detected.")

    pert_info_path = Path(args.pert_info) if args.pert_info else None
    pert_map, pert_name_map = build_pertinfo_maps(pert_info_path) if pert_info_path else ({}, {})

    smiles_col = None
    for c in df.columns:
        if 'smiles' == c.lower() or 'smiles' in c.lower():
            smiles_col = c
            break
    if not smiles_col:

        df['smiles'] = ''
        smiles_col = 'smiles'

    missing_mask = df[smiles_col].isna() | (df[smiles_col].astype(str).str.strip() == '')
    logging.info(f"Dataset rows: {len(df)}, currently {missing_mask.sum()} rows missing SMILES")

    if args.unmapped and Path(args.unmapped).exists():
        with open(args.unmapped, 'r', encoding='utf-8') as fh:
            unmapped_list = [l.strip() for l in fh if l.strip()]
        logging.info(f"Loaded {len(unmapped_list)} BRD IDs from {args.unmapped}")
    else:
        unmapped_list = df.loc[missing_mask, 'BRD_ID'].dropna().unique().tolist()
        logging.info(f"Detected {len(unmapped_list)} unique unmapped BRD IDs from input dataset")

    brd_to_smiles = {}
    name_candidates = {}

    pubchem_cache = load_pubchem_cache()

    mapped_offline = 0
    for brd in unmapped_list:
        if not brd:
            continue
        found = None
        for v in canonical_brd_variants(brd):
            if v in rep_map and rep_map[v]:
                found = rep_map[v]
                break
        if found:
            brd_to_smiles[brd] = found
            mapped_offline += 1
    logging.info(f"Mapped {mapped_offline} BRD IDs from repurposing (offline)")

    mapped_pert = 0
    for brd in unmapped_list:
        if brd in brd_to_smiles:
            continue
        for v in canonical_brd_variants(brd):
            if v in pert_map and pert_map[v]:
                brd_to_smiles[brd] = pert_map[v]
                mapped_pert += 1
                break
    logging.info(f"Mapped {mapped_pert} BRD IDs from pert_info (offline)")

    remaining = [b for b in unmapped_list if b not in brd_to_smiles]
    for brd in remaining:
        cand_names = []

        if 'pert_iname' in df.columns:
            vals = df.loc[df['BRD_ID'] == brd, 'pert_iname'].dropna().unique().tolist()
            cand_names.extend([str(x).strip() for x in vals if str(x).strip()])

        for v in canonical_brd_variants(brd):
            if v in pert_name_map:
                cand_names.append(pert_name_map[v])

        for v in canonical_brd_variants(brd):

            pass

        if not cand_names:
            cand_names = [brd]
        for n in cand_names:
            nk = n.strip().lower()
            if nk and nk not in name_candidates:
                name_candidates[nk] = n.strip()

    mapped_by_repname = 0
    for name_lower, sample in list(name_candidates.items()):
        if name_lower in rep_name_map and rep_name_map[name_lower]:

            for brd in remaining:

                suggested = False
                if 'pert_iname' in df.columns:
                    vals = df.loc[df['BRD_ID'] == brd, 'pert_iname'].dropna().unique().tolist()
                    if any(str(v).strip().lower() == name_lower for v in vals):
                        suggested = True
                for v in canonical_brd_variants(brd):
                    if pert_name_map.get(v, '').strip().lower() == name_lower:
                        suggested = True
                if suggested and brd not in brd_to_smiles:
                    brd_to_smiles[brd] = rep_name_map[name_lower]
                    mapped_by_repname += 1
    logging.info(f"Mapped {mapped_by_repname} BRDs by repurposing name map (offline)")

    df['fetched_smiles_offline'] = df['BRD_ID'].map(brd_to_smiles)
    mask_missing = (df[smiles_col].isna() | (df[smiles_col].astype(str).str.strip() == '')) & df['fetched_smiles_offline'].notna()
    df.loc[mask_missing, smiles_col] = df.loc[mask_missing, 'fetched_smiles_offline']
    logging.info(f"After offline mapping, missing SMILES: {(df[smiles_col].isna() | (df[smiles_col].astype(str).str.strip() == '')).sum()}")

    remaining_mask = df[smiles_col].isna() | (df[smiles_col].astype(str).str.strip() == '')
    remaining_brds = df.loc[remaining_mask, 'BRD_ID'].dropna().unique().tolist()
    logging.info(f"Remaining unmapped unique BRDs before PubChem: {len(remaining_brds)}")

    query_names = []
    for brd in remaining_brds:
        names = []
        if 'pert_iname' in df.columns:
            vals = df.loc[df['BRD_ID'] == brd, 'pert_iname'].dropna().unique().tolist()
            names.extend([str(x).strip() for x in vals if str(x).strip()])
        for v in canonical_brd_variants(brd):
            if v in pert_name_map:
                names.append(pert_name_map[v])
        if not names:
            names = [brd]
        for n in names:
            lk = str(n).strip().lower()
            if lk and (lk not in pubchem_cache):
                query_names.append(n.strip())

    uniq_query = list(dict.fromkeys(query_names))
    if args.max and int(args.max) > 0:
        uniq_query = uniq_query[:int(args.max)]
        logging.info(f"Limiting PubChem queries to first {len(uniq_query)} names (max set)")

    if uniq_query:
        pubchem_results = threaded_pubchem_query(uniq_query, delay=args.delay, max_workers=MAX_WORKERS)

        for n, s in pubchem_results.items():
            pubchem_cache[n.lower()] = s
        save_pubchem_cache(pubchem_cache)

        newly_mapped = 0
        for brd in remaining_brds:
            if brd in brd_to_smiles:
                continue

            names = []
            if 'pert_iname' in df.columns:
                vals = df.loc[df['BRD_ID'] == brd, 'pert_iname'].dropna().unique().tolist()
                names.extend([str(x).strip() for x in vals if str(x).strip()])
            for v in canonical_brd_variants(brd):
                if v in pert_name_map:
                    names.append(pert_name_map[v])
            if not names:
                names = [brd]
            found = None
            for n in names:
                s = pubchem_cache.get(str(n).strip().lower())
                if s:
                    found = s
                    break
            if found:
                brd_to_smiles[brd] = found
                newly_mapped += 1

        logging.info(f"Mapped {newly_mapped} BRDs from PubChem fallback")

        df['fetched_smiles_pubchem'] = df['BRD_ID'].map(brd_to_smiles)
        mask_missing2 = (df[smiles_col].isna() | (df[smiles_col].astype(str).str.strip() == '')) & df['fetched_smiles_pubchem'].notna()
        df.loc[mask_missing2, smiles_col] = df.loc[mask_missing2, 'fetched_smiles_pubchem']

    outpath = Path(args.output) if args.output else Path("C4A_training_dataset_filled.tsv")
    df.to_csv(outpath, sep='\t', index=False)
    logging.info(f"Wrote updated dataset to {outpath} ({len(df)} rows)")

    remaining_final = df.loc[df[smiles_col].isna() | (df[smiles_col].astype(str).str.strip() == ''), 'BRD_ID'].dropna().unique().tolist()
    rempath = Path("unmapped_brd_ids_remaining.txt")
    rempath.write_text("\n".join([str(x) for x in remaining_final]), encoding='utf-8')
    logging.info(f"Wrote {len(remaining_final)} remaining unmapped BRD IDs to {rempath}")

    mapped_now = df[smiles_col].notna().sum()
    logging.info(f"Summary: mapped SMILES count now {mapped_now} rows (of {len(df)})")

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', default='C4A_training_dataset.tsv', help='Input ML dataset (tab-separated) with BRD_ID column')
    ap.add_argument('--unmapped', default='unmapped_brd_ids.txt', help='File with unmapped BRD IDs (one per line). If missing, script will detect from input.')
    ap.add_argument('--pert_info', default='GSE92742_Broad_LINCS_pert_info.txt', help='Perturbagen metadata file to get pert_iname -> name mapping (optional but recommended)')
    ap.add_argument('--repurposing', default=None, help='Optional explicit repurposing file path (if not provided, script tries to detect it in cwd)')
    ap.add_argument('--output', default='C4A_training_dataset_filled.tsv', help='Output merged dataset path')
    ap.add_argument('--delay', type=float, default=DEFAULT_DELAY, help='Delay seconds between PubChem calls (default 0.25)')
    ap.add_argument('--max', type=int, default=0, help='Max number of names to query PubChem (0 = unlimited)')
    args = ap.parse_args()
    main(args)
