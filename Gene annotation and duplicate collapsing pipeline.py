#!/usr/bin/env python3
"""Gene annotation and duplicate collapsing pipeline."""

import argparse
import gzip
import json
import logging
import re
import sqlite3
import sys
from datetime import datetime, timezone
from pathlib import Path

import mygene
import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry
from tqdm import tqdm

ENSEMBL_REST = "https://rest.ensembl.org"
NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
INPUT_FILENAME = "top_genes_mapped.tsv"
ANNOTATED_SUFFIX = "_annotated.tsv"
COLLAPSED_SUFFIX = "_annotated_collapsed.tsv"
FAILED_MAPPINGS_FILE = "failed_mappings.tsv"
DEFAULT_CACHE_NAME = "gene_annotation_cache.sqlite3"
NCBI_TOOL = "gene-annotation-pipeline"
DEFAULT_NCBI_EMAIL = ""

def utc_now():
    return datetime.now(timezone.utc).isoformat()

def make_session():
    session = requests.Session()
    retries = Retry(
        total=5,
        backoff_factor=0.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
    )
    adapter = HTTPAdapter(max_retries=retries)
    session.mount("https://", adapter)
    session.headers.update({
        "Content-Type": "application/json",
        "User-Agent": "gene-annotation-pipeline/1.0",
    })
    return session

SESSION = make_session()
MG = mygene.MyGeneInfo()

class CacheDB:
    def __init__(self, path):
        self.path = str(path)
        self.conn = sqlite3.connect(self.path)
        self.make_tables()

    def make_tables(self):
        cur = self.conn.cursor()
        cur.execute("""
            CREATE TABLE IF NOT EXISTS mappings_v2 (
                symbol TEXT NOT NULL,
                species TEXT NOT NULL,
                entrez INTEGER,
                ensembl_gene TEXT,
                ensembl_transcripts TEXT,
                name TEXT,
                raw_json TEXT,
                retrieved_at TEXT,
                PRIMARY KEY (symbol, species)
            )
        """)
        cur.execute("""
            CREATE TABLE IF NOT EXISTS sequences_v2 (
                ensembl_gene TEXT PRIMARY KEY,
                dna_seq TEXT,
                cds_seqs TEXT,
                protein_seqs TEXT,
                retrieved_at TEXT
            )
        """)
        self.conn.commit()

    def get_mapping(self, symbol, species="human"):
        cur = self.conn.cursor()
        cur.execute(
            """
            SELECT entrez, ensembl_gene, ensembl_transcripts,
                   name, raw_json, retrieved_at
            FROM mappings_v2
            WHERE symbol = ? AND species = ?
            """,
            (symbol, species),
        )
        row = cur.fetchone()
        if not row:
            return None

        entrez, ensembl_gene, transcripts, name, raw_json, retrieved_at = row
        raw = json.loads(raw_json) if raw_json else None

        if isinstance(raw, dict) and raw.get("notfound"):
            return None

        return {
            "symbol": symbol,
            "entrez": entrez,
            "ensembl_gene": ensembl_gene,
            "ensembl_transcripts": (
                json.loads(transcripts) if transcripts else None
            ),
            "name": name,
            "raw": raw,
            "retrieved_at": retrieved_at,
            "species": species,
        }

    def set_mapping(
        self,
        symbol,
        species,
        entrez=None,
        ensembl_gene=None,
        ensembl_transcripts=None,
        name=None,
        raw_json=None,
    ):
        self.conn.execute(
            """
            INSERT OR REPLACE INTO mappings_v2
            (symbol, species, entrez, ensembl_gene, ensembl_transcripts,
             name, raw_json, retrieved_at)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                symbol,
                species,
                entrez,
                ensembl_gene,
                json.dumps(ensembl_transcripts)
                if ensembl_transcripts is not None else None,
                name,
                json.dumps(raw_json) if raw_json is not None else None,
                utc_now(),
            ),
        )
        self.conn.commit()

    def get_sequences(self, ensembl_gene):
        cur = self.conn.cursor()
        cur.execute(
            """
            SELECT dna_seq, cds_seqs, protein_seqs, retrieved_at
            FROM sequences_v2
            WHERE ensembl_gene = ?
            """,
            (ensembl_gene,),
        )
        row = cur.fetchone()
        if not row:
            return None

        dna, cds, protein, retrieved_at = row
        return {
            "ensembl_gene": ensembl_gene,
            "dna_seq": dna,
            "cds_seqs": json.loads(cds) if cds else None,
            "protein_seqs": json.loads(protein) if protein else None,
            "retrieved_at": retrieved_at,
        }

    def set_sequences(self, ensembl_gene, dna, cds, protein):
        self.conn.execute(
            """
            INSERT OR REPLACE INTO sequences_v2
            (ensembl_gene, dna_seq, cds_seqs, protein_seqs, retrieved_at)
            VALUES (?, ?, ?, ?, ?)
            """,
            (
                ensembl_gene,
                dna,
                json.dumps(cds) if cds is not None else None,
                json.dumps(protein) if protein is not None else None,
                utc_now(),
            ),
        )
        self.conn.commit()

    def close(self):
        self.conn.close()

def clean_symbol(value):
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    value = str(value).strip()
    if not value:
        return ""
    value = value.strip("'\"")
    value = re.sub(r"\.(\d+)$", "", value)
    return value

def split_candidates(raw):
    if raw is None or (isinstance(raw, float) and pd.isna(raw)):
        return []

    pieces = re.split(r"///|;|,|\||/|\\|\s+", str(raw))
    pieces = [clean_symbol(p) for p in pieces]
    pieces = [p for p in pieces if p]

    def score(token):
        upper = token.upper()
        if re.match(r"^[A-Z0-9_\-]+$", token, flags=re.I):
            return 4
        if re.match(r"^[A-Za-z]{1,5}\d{0,3}$", token):
            return 4
        if re.match(
            r"^(AK|AL|BC|CTD|CTB|AX|FLJ|DQ|AF|AY|NM|NR|XR)\d+",
            token,
            flags=re.I,
        ):
            return 3
        if re.match(r"^[A-Z]{2}\d{5,}", token, flags=re.I):
            return 3
        if upper.startswith(("RP", "LOC")):
            return 2
        return 0

    unique = []
    seen = set()
    for piece in pieces:
        if piece not in seen:
            seen.add(piece)
            unique.append(piece)

    return sorted(unique, key=lambda x: -score(x))

def is_accession(token):
    if not token:
        return False

    patterns = [
        r"^(AK|AL|BC|CTD|CTB|AX|FLJ|DQ|AF|AY|NM|NR|XR)\d+",
        r"^ENST\d+",
        r"^ENSG\d+",
        r"^[A-Z]{2}\d{5,}",
    ]

    return any(
        re.match(pattern, token, flags=re.I)
        for pattern in patterns
    ) or token.upper().startswith("AX")

def query_mygene_bulk(symbols, species="human"):
    symbols = [x for x in sorted(set(symbols)) if x]
    if not symbols:
        return {}

    results = {}
    chunk_size = 800

    for start in range(0, len(symbols), chunk_size):
        chunk = symbols[start:start + chunk_size]

        try:
            data = MG.querymany(
                chunk,
                scopes=[
                    "symbol",
                    "alias",
                    "name",
                    "accession",
                    "refseq",
                    "unigene",
                ],
                fields="entrezgene,ensembl,symbol,name",
                species=species,
                as_dataframe=False,
            )
        except Exception as exc:
            logging.warning(
                "MyGene batch query failed near item %s: %s",
                start,
                exc,
            )
            data = query_mygene_one_by_one(chunk, species)

        for item in data:
            if not isinstance(item, dict):
                continue

            key = item.get("query") or item.get("symbol")
            if not key:
                continue

            if item.get("notfound"):
                results[key] = None
            else:
                results[key] = item

    return results

def query_mygene_one_by_one(symbols, species):
    out = []

    for symbol in symbols:
        try:
            data = MG.query(
                symbol,
                species=species,
                fields="entrezgene,ensembl,symbol,name",
            )
            hits = data.get("hits", []) if data else []
            if hits:
                hit = hits[0]
                hit["query"] = symbol
                out.append(hit)
            else:
                out.append({"notfound": True, "query": symbol})
        except Exception:
            out.append({"notfound": True, "query": symbol})

    return out

def parse_ensembl_field(field):
    gene = None
    transcripts = None

    if isinstance(field, dict):
        gene = field.get("gene")
        transcripts = field.get("transcript")

    elif isinstance(field, list):
        for item in field:
            if isinstance(item, dict) and item.get("gene"):
                gene = item.get("gene")
                transcripts = item.get("transcript")
                break

        if gene is None and field:
            first = field[0]
            if isinstance(first, str) and first.startswith("ENS"):
                gene = first

    elif isinstance(field, str) and field.startswith("ENS"):
        gene = field

    return gene, transcripts

def mapping_from_mygene(item):
    if not item:
        return None

    entrez = item.get("entrezgene")
    name = item.get("name")
    ensembl_gene, transcripts = parse_ensembl_field(item.get("ensembl"))

    return {
        "entrez": entrez,
        "ensembl_gene": ensembl_gene,
        "ensembl_transcripts": transcripts,
        "name": name,
        "raw": item,
    }

def ncbi_get(path, params, timeout=20):
    try:
        response = SESSION.get(
            f"{NCBI_EUTILS_BASE}/{path}",
            params=params,
            timeout=timeout,
        )
        response.raise_for_status()
        return response.json()
    except Exception:
        return None

def ncbi_esearch_nuccore(accession, email=""):
    if not accession:
        return None

    params = {
        "db": "nuccore",
        "term": f"{accession}[Accession]",
        "retmode": "json",
        "tool": NCBI_TOOL,
        "email": email,
    }

    data = ncbi_get("esearch.fcgi", params)
    if not data:
        return None

    return data.get("esearchresult", {}).get("idlist", [])

def ncbi_elink_nuccore_to_gene(nuccore_id, email=""):
    if not nuccore_id:
        return None

    params = {
        "dbfrom": "nuccore",
        "db": "gene",
        "id": nuccore_id,
        "retmode": "json",
        "tool": NCBI_TOOL,
        "email": email,
    }

    data = ncbi_get("elink.fcgi", params)
    if not data:
        return None

    linked = []
    for linkset in data.get("linksets", []):
        for group in linkset.get("linksetdbs", []) or []:
            if group.get("dbto") == "gene":
                linked.extend(group.get("links", []) or [])

    return linked

def ncbi_esummary_gene(gene_id, email=""):
    if not gene_id:
        return None

    params = {
        "db": "gene",
        "id": gene_id,
        "retmode": "json",
        "tool": NCBI_TOOL,
        "email": email,
    }

    data = ncbi_get("esummary.fcgi", params)
    if not data:
        return None

    return data.get("result", {}).get(str(gene_id))

def ncbi_lookup(token, species="human", email=""):
    if not is_accession(token):
        return None

    ids = ncbi_esearch_nuccore(token, email=email)
    if not ids:
        return None

    for nuccore_id in ids[:3]:
        gene_ids = ncbi_elink_nuccore_to_gene(
            nuccore_id,
            email=email,
        )
        if not gene_ids:
            continue

        summary = ncbi_esummary_gene(
            gene_ids[0],
            email=email,
        )
        if not summary:
            continue

        entrez = summary.get("uid")
        try:
            entrez = int(entrez) if entrez else None
        except (TypeError, ValueError):
            entrez = None

        symbol = (
            summary.get("nomenclaturesymbol")
            or summary.get("name")
        )
        name = (
            summary.get("description")
            or summary.get("title")
            or summary.get("nomenclaturename")
            or symbol
        )

        return {
            "entrez": entrez,
            "ensembl_gene": None,
            "ensembl_transcripts": None,
            "name": name,
            "symbol": symbol,
            "raw": summary,
        }

    return None

def ensembl_request(path, accept, timeout=20, params=None):
    url = f"{ENSEMBL_REST}/{path.lstrip('/')}"

    try:
        response = SESSION.get(
            url,
            headers={"Accept": accept},
            params=params,
            timeout=timeout,
        )
        if response.status_code != 200:
            return None

        if accept == "application/json":
            return response.json()
        return response.text.strip()
    except Exception:
        return None

def ensembl_lookup_symbol(symbol, species="human"):
    return ensembl_request(
        f"lookup/symbol/{species}/{symbol}",
        "application/json",
        params={"expand": 0},
    )

def ensembl_xrefs_symbol(symbol, species="human"):
    return ensembl_request(
        f"xrefs/symbol/{species}/{symbol}",
        "application/json",
    )

def ensembl_fetch_sequence(ensembl_id, seq_type="genomic"):
    if seq_type not in {"genomic", "cdna", "cds", "protein"}:
        seq_type = "genomic"

    return ensembl_request(
        f"sequence/id/{ensembl_id}",
        "text/plain",
        timeout=30,
        params={"type": seq_type},
    )

def ensembl_mapping(token, species="human"):
    lookup = ensembl_lookup_symbol(token, species)
    if lookup and lookup.get("id"):
        return {
            "entrez": None,
            "ensembl_gene": lookup.get("id"),
            "ensembl_transcripts": None,
            "name": lookup.get("description")
            or lookup.get("display_name"),
            "raw": lookup,
            "source": "Ensembl-lookup",
        }

    xrefs = ensembl_xrefs_symbol(token, species)
    if not xrefs:
        return None

    chosen = None
    for item in xrefs:
        if item.get("type") == "gene":
            chosen = item
            break

    if chosen is None:
        chosen = xrefs[0]

    if not chosen or not chosen.get("id"):
        return None

    return {
        "entrez": None,
        "ensembl_gene": chosen.get("id"),
        "ensembl_transcripts": None,
        "name": chosen.get("display_name"),
        "raw": chosen,
        "source": "Ensembl-xrefs",
    }

def open_gtf(path):
    path = Path(path)
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("r")

def parse_gtf_attributes(text):
    attrs = {}
    for key, value in re.findall(
        r'([A-Za-z0-9_\-]+)\s+"([^"]+)"',
        text,
    ):
        attrs[key] = value
    return attrs

def build_gencode_index(path):
    index = {}
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"GENCODE file not found: {path}")

    try:
        with open_gtf(path) as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue

                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[2] != "gene":
                    continue

                attrs = parse_gtf_attributes(parts[8])
                gene_id = attrs.get("gene_id")
                gene_name = attrs.get("gene_name")

                data = {
                    "gene_id": gene_id,
                    "attrs": attrs,
                }

                if gene_name:
                    index.setdefault(gene_name, data)
                if gene_id:
                    index.setdefault(gene_id, data)

    except Exception as exc:
        logging.warning("Could not read GENCODE file: %s", exc)
        return {}

    return index

def gencode_mapping(token, index):
    if not index:
        return None

    hit = index.get(token)
    if not hit:
        return None

    attrs = hit.get("attrs", {})
    gene_id = hit.get("gene_id")
    name = attrs.get("gene_name") or attrs.get("gene_type")

    if not gene_id:
        return None

    return {
        "entrez": None,
        "ensembl_gene": gene_id,
        "ensembl_transcripts": None,
        "name": name,
        "raw": {"gencode": attrs},
        "source": "GENCODE",
    }

def choose_mapping(
    original,
    candidates,
    mygene_results,
    cache,
    species,
    gencode_index=None,
    force_refresh=False,
    use_ncbi=True,
    ncbi_email="",
):
    for token in candidates:
        if not token:
            continue

        if not force_refresh:
            cached = cache.get_mapping(token, species)
            if cached:
                cached["chosen_token"] = token
                cached["mapping_source"] = "cache"
                return cached

    for token in candidates:
        if not token:
            continue

        result = mygene_results.get(token)
        if not result:
            continue

        mapping = mapping_from_mygene(result)
        if mapping is None:
            continue

        cache.set_mapping(
            token,
            species,
            mapping.get("entrez"),
            mapping.get("ensembl_gene"),
            mapping.get("ensembl_transcripts"),
            mapping.get("name"),
            result,
        )

        saved = cache.get_mapping(token, species)
        if saved:
            saved["chosen_token"] = token
            saved["mapping_source"] = "MyGene"
            return saved

    if use_ncbi:
        for token in candidates:
            if not token or not is_accession(token):
                continue

            mapping = ncbi_lookup(
                token,
                species=species,
                email=ncbi_email,
            )
            if not mapping:
                continue

            cache.set_mapping(
                token,
                species,
                mapping.get("entrez"),
                mapping.get("ensembl_gene"),
                mapping.get("ensembl_transcripts"),
                mapping.get("name"),
                mapping.get("raw"),
            )

            saved = cache.get_mapping(token, species)
            if saved:
                saved["chosen_token"] = token
                saved["mapping_source"] = "NCBI"
                if mapping.get("symbol"):
                    saved["symbol"] = mapping["symbol"]
                return saved

    if gencode_index:
        for token in candidates:
            if not token:
                continue

            mapping = gencode_mapping(token, gencode_index)
            if not mapping:
                continue

            cache.set_mapping(
                token,
                species,
                mapping.get("entrez"),
                mapping.get("ensembl_gene"),
                mapping.get("ensembl_transcripts"),
                mapping.get("name"),
                mapping.get("raw"),
            )

            saved = cache.get_mapping(token, species)
            if saved:
                saved["chosen_token"] = token
                saved["mapping_source"] = "GENCODE"
                return saved

    for token in candidates or [original]:
        if not token:
            continue

        mapping = ensembl_mapping(token, species)
        if not mapping:
            continue

        cache.set_mapping(
            token,
            species,
            mapping.get("entrez"),
            mapping.get("ensembl_gene"),
            mapping.get("ensembl_transcripts"),
            mapping.get("name"),
            mapping.get("raw"),
        )

        saved = cache.get_mapping(token, species)
        if saved:
            saved["chosen_token"] = token
            saved["mapping_source"] = mapping.get("source", "Ensembl")
            return saved

    return None

def fetch_sequences(ensembl_ids, cache, force_refresh=False):
    for gene_id in tqdm(sorted(ensembl_ids), desc="Fetching sequences"):
        if not gene_id:
            continue

        old = cache.get_sequences(gene_id)
        if old and not force_refresh:
            continue

        dna = ensembl_fetch_sequence(gene_id, "genomic")
        cds = ensembl_fetch_sequence(gene_id, "cds")
        protein = ensembl_fetch_sequence(gene_id, "protein")

        cache.set_sequences(
            gene_id,
            dna,
            [cds] if cds else None,
            [protein] if protein else None,
        )

def fill_mapping_columns(df, mapping_results, cache):
    out = df.copy()

    out["chosen_candidate"] = None
    out["mapping_source"] = None
    out["entrez_id"] = None
    out["ensembl_gene"] = None
    out["description"] = None
    out["dna_seq_present"] = False
    out["protein_seq_present"] = False

    for index, row in out.iterrows():
        original = str(row.get("gene_symbol", ""))
        mapping = mapping_results.get(original)

        if not mapping:
            continue

        out.at[index, "chosen_candidate"] = (
            mapping.get("symbol") or mapping.get("chosen_token")
        )
        out.at[index, "mapping_source"] = (
            mapping.get("mapping_source") or "unknown"
        )
        out.at[index, "entrez_id"] = mapping.get("entrez")
        out.at[index, "ensembl_gene"] = mapping.get("ensembl_gene")
        out.at[index, "description"] = mapping.get("name")

        gene_id = mapping.get("ensembl_gene")
        if gene_id:
            sequence_data = cache.get_sequences(gene_id)
            if sequence_data:
                out.at[index, "dna_seq_present"] = bool(
                    sequence_data.get("dna_seq")
                )
                out.at[index, "protein_seq_present"] = bool(
                    sequence_data.get("protein_seqs")
                )

    return out

def canonical_id(row):
    gene_id = row.get("ensembl_gene")
    if gene_id and str(gene_id).strip().upper().startswith("ENS"):
        return str(gene_id).strip()

    entrez = row.get("entrez_id")
    if pd.notna(entrez) and str(entrez).strip() not in {
        "",
        "None",
        "nan",
    }:
        try:
            return f"ENTREZ:{int(float(entrez))}"
        except (TypeError, ValueError):
            return f"ENTREZ:{entrez}"

    chosen = row.get("chosen_candidate")
    if chosen and str(chosen).strip():
        return str(chosen).strip()

    return str(row.get("gene_symbol", "")).strip()

def add_numeric_versions(df, columns):
    numeric = []
    other = []

    for column in columns:
        values = pd.to_numeric(df[column], errors="coerce")
        if values.notna().sum() >= 1:
            numeric.append(column)
            df[column + "_numeric"] = values
        else:
            other.append(column)

    return numeric, other

def collapse_annotations(df, agg_method="median"):
    out = df.copy()
    out["_canonical_id"] = out.apply(canonical_id, axis=1)
    out["_orig_symbol"] = out["gene_symbol"].astype(str)

    metadata_columns = {
        "gene_symbol",
        "chosen_candidate",
        "mapping_source",
        "entrez_id",
        "ensembl_gene",
        "description",
        "dna_seq_present",
        "protein_seq_present",
        "_canonical_id",
        "_orig_symbol",
    }

    value_columns = [
        column
        for column in out.columns
        if column not in metadata_columns
        and not column.endswith("_numeric")
    ]

    numeric_columns, other_columns = add_numeric_versions(
        out,
        value_columns,
    )

    groups = out.groupby("_canonical_id", sort=False)
    rows = []

    mapping_fields = [
        "chosen_candidate",
        "mapping_source",
        "entrez_id",
        "ensembl_gene",
        "description",
        "dna_seq_present",
        "protein_seq_present",
    ]

    for canonical, group in groups:
        mapped = group[
            group["ensembl_gene"].notna()
            | group["entrez_id"].notna()
        ]

        base = (
            mapped.iloc[0]
            if len(mapped) else group.iloc[0]
        )

        row = {
            "canonical_id": canonical,
            "n_probes": int(len(group)),
            "orig_symbols": " | ".join(
                sorted(set(group["_orig_symbol"].astype(str)))
            ),
        }

        for field in mapping_fields:
            row[field] = base.get(field)

        for column in numeric_columns:
            values = group[column + "_numeric"]
            if agg_method == "mean":
                value = values.mean(skipna=True)
            else:
                value = values.median(skipna=True)
            row[column] = (
                float(value) if pd.notna(value) else None
            )

        for column in other_columns:
            values = [
                str(value).strip()
                for value in group[column].astype(str)
                if str(value).strip() not in {"", "nan", "None"}
            ]
            values = list(dict.fromkeys(values))
            row[column] = " | ".join(values) if values else None

        rows.append(row)

    return pd.DataFrame(rows)

def annotate_file(
    infile,
    cache,
    species="human",
    force_refresh=False,
    fetch_seq=True,
    gencode_index=None,
    agg_method="median",
    use_ncbi=True,
    ncbi_email="",
):
    infile = Path(infile)
    df = pd.read_csv(
        infile,
        sep="\t",
        dtype=str,
        low_memory=False,
    )

    if "gene_symbol" not in df.columns:
        raise ValueError(
            f"Input file {infile} does not contain a 'gene_symbol' column"
        )

    original_values = df["gene_symbol"].tolist()

    candidate_map = {}
    all_candidates = []

    for value in original_values:
        candidates = split_candidates(value)
        original_key = "" if pd.isna(value) else str(value)
        candidate_map[original_key] = candidates
        all_candidates.extend(candidates)

    mygene_results = query_mygene_bulk(
        all_candidates,
        species=species,
    )

    mappings = {}
    failed = []

    for original, candidates in candidate_map.items():
        mapping = choose_mapping(
            original,
            candidates,
            mygene_results,
            cache,
            species,
            gencode_index=gencode_index,
            force_refresh=force_refresh,
            use_ncbi=use_ncbi,
            ncbi_email=ncbi_email,
        )

        if mapping:
            mappings[original] = mapping
        else:
            mappings[original] = None
            failed.append(original)

    if fetch_seq:
        ids = {
            mapping.get("ensembl_gene")
            for mapping in mappings.values()
            if mapping and mapping.get("ensembl_gene")
        }
        fetch_sequences(
            ids,
            cache,
            force_refresh=force_refresh,
        )

    annotated = fill_mapping_columns(
        df,
        mappings,
        cache,
    )

    annotated_path = infile.with_name(
        infile.stem + ANNOTATED_SUFFIX
    )
    annotated.to_csv(
        annotated_path,
        sep="\t",
        index=False,
    )

    collapsed = collapse_annotations(
        annotated,
        agg_method=agg_method,
    )

    collapsed_path = infile.with_name(
        infile.stem + COLLAPSED_SUFFIX
    )
    collapsed.to_csv(
        collapsed_path,
        sep="\t",
        index=False,
    )

    return (
        annotated_path,
        collapsed_path,
        sorted(set(failed)),
    )

def discover_files(base_dir, process_ml_manual=False, all_tsv=False):
    base_dir = Path(base_dir)

    if all_tsv:
        return sorted(base_dir.rglob("*.tsv"))

    if process_ml_manual:
        found = []
        for folder in [base_dir / "ML", base_dir / "Manual"]:
            if not folder.exists():
                continue
            found.extend(folder.rglob(INPUT_FILENAME))
        return sorted(found)

    return sorted(base_dir.rglob(INPUT_FILENAME))

def choose_cache_path(base_dir, requested):
    if requested:
        return Path(requested).expanduser().resolve()

    cache_dir = Path(base_dir) / ".cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir / DEFAULT_CACHE_NAME

def build_parser():
    parser = argparse.ArgumentParser(
        description="Annotate gene tables and collapse duplicate gene rows."
    )

    parser.add_argument(
        "--base-dir",
        "-b",
        required=False,
        help="Folder to search for input TSV files.",
    )
    parser.add_argument(
        "--cache-db",
        default=None,
        help="SQLite cache location. Defaults to <base-dir>/.cache/.",
    )
    parser.add_argument(
        "--species",
        default="human",
        help="Species name used by MyGene and Ensembl. Default: human.",
    )
    parser.add_argument(
        "--process-ml-manual",
        action="store_true",
        help="Only search the ML and Manual folders for top_genes_mapped.tsv.",
    )
    parser.add_argument(
        "--all-tsv",
        action="store_true",
        help="Process every .tsv file under the base folder.",
    )
    parser.add_argument(
        "--force-refresh",
        action="store_true",
        help="Ignore cached mappings and sequence results.",
    )
    parser.add_argument(
        "--no-seq",
        action="store_true",
        help="Map genes but do not fetch DNA/CDS/protein sequences.",
    )
    parser.add_argument(
        "--no-ncbi",
        action="store_true",
        help="Skip the NCBI accession fallback.",
    )
    parser.add_argument(
        "--gencode-gtf",
        default=None,
        help="Optional GENCODE GTF or GTF.gz file.",
    )
    parser.add_argument(
        "--agg-method",
        choices=["median", "mean"],
        default="median",
        help="How numeric duplicate rows are combined. Default: median.",
    )
    parser.add_argument(
        "--ncbi-email",
        default=DEFAULT_NCBI_EMAIL,
        help="Email passed to NCBI E-utilities.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Only process this many files. Useful for testing.",
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Amount of logging to print.",
    )

    return parser

def configure_logging(level):
    logging.basicConfig(
        level=getattr(logging, level),
        format="[%(levelname)s] %(message)s",
    )

def load_gencode_if_requested(path):
    if not path:
        return None

    path = Path(path).expanduser().resolve()
    print(f"Reading GENCODE annotation: {path}")
    index = build_gencode_index(path)
    print(f"GENCODE entries loaded: {len(index)}")
    return index

def write_failed(base_dir, failed):
    if not failed:
        return None

    failed_path = Path(base_dir) / FAILED_MAPPINGS_FILE
    values = sorted(set(failed))

    pd.DataFrame({
        "failed_symbol": values,
    }).to_csv(
        failed_path,
        sep="\t",
        index=False,
    )

    return failed_path

def main():
    parser = build_parser()
    parser.add_argument(
        "--input",
        nargs="+",
        default=None,
        help="Process these files directly instead of searching the base folder.",
    )
    args = parser.parse_args()
    configure_logging(args.log_level)

    if not args.base_dir and not args.input:
        print("Use --base-dir or --input.", file=sys.stderr)
        return 1

    if args.base_dir:
        base_dir = Path(args.base_dir).expanduser().resolve()
    else:
        base_dir = Path(args.input[0]).expanduser().resolve().parent

    if not base_dir.exists():
        print(f"Base directory does not exist: {base_dir}", file=sys.stderr)
        return 1

    if not base_dir.is_dir():
        print(f"Base path is not a directory: {base_dir}", file=sys.stderr)
        return 1

    if args.all_tsv and args.process_ml_manual:
        print(
            "Use either --all-tsv or --process-ml-manual, not both.",
            file=sys.stderr,
        )
        return 1

    if args.input and (args.all_tsv or args.process_ml_manual):
        print(
            "Use --input by itself, not with --all-tsv or --process-ml-manual.",
            file=sys.stderr,
        )
        return 1

    if args.limit < 0:
        print("--limit cannot be negative.", file=sys.stderr)
        return 1

    cache_path = choose_cache_path(base_dir, args.cache_db)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    if args.input:
        files = []
        for value in args.input:
            path = Path(value).expanduser().resolve()
            if not path.exists():
                print(f"Input file does not exist: {path}", file=sys.stderr)
                return 1
            if not path.is_file():
                print(f"Input path is not a file: {path}", file=sys.stderr)
                return 1
            files.append(path)
        files = sorted(set(files))
    else:
        files = discover_files(
            base_dir,
            process_ml_manual=args.process_ml_manual,
            all_tsv=args.all_tsv,
        )

    if not files:
        print(
            "No input files found. Check the base directory or use --all-tsv."
        )
        return 0

    if args.limit:
        files = files[:args.limit]

    print(f"Found {len(files)} input file(s)")
    print(f"Base directory: {base_dir}")
    print(f"Species: {args.species}")
    print(f"Cache: {cache_path}")
    print(f"Sequence fetching: {'on' if not args.no_seq else 'off'}")
    print(f"NCBI fallback: {'on' if not args.no_ncbi else 'off'}")
    print(f"Duplicate aggregation: {args.agg_method}")

    gencode_index = None
    if args.gencode_gtf:
        gencode_path = Path(args.gencode_gtf).expanduser().resolve()
        if not gencode_path.exists():
            print(
                f"GENCODE file does not exist: {gencode_path}",
                file=sys.stderr,
            )
            return 1
        gencode_index = load_gencode_if_requested(gencode_path)

    if args.ncbi_email:
        logging.info("Using the supplied NCBI email")
    elif not args.no_ncbi:
        logging.info(
            "NCBI fallback is enabled. Passing --ncbi-email is recommended."
        )

    cache = CacheDB(cache_path)
    processed = []
    all_failed = []

    try:
        for number, infile in enumerate(files, start=1):
            print(f"\n[{number}/{len(files)}] Processing {infile}")

            try:
                check_input = pd.read_csv(
                    infile,
                    sep="\t",
                    nrows=0,
                )
                if "gene_symbol" not in check_input.columns:
                    raise ValueError(
                        "input file is missing the 'gene_symbol' column"
                    )

                annotated_path, collapsed_path, failed = annotate_file(
                    infile,
                    cache,
                    species=args.species,
                    force_refresh=args.force_refresh,
                    fetch_seq=not args.no_seq,
                    gencode_index=gencode_index,
                    agg_method=args.agg_method,
                    use_ncbi=not args.no_ncbi,
                    ncbi_email=args.ncbi_email,
                )

                print(f"  annotated: {annotated_path}")
                print(f"  collapsed: {collapsed_path}")

                if failed:
                    print(f"  unmapped: {len(failed)}")
                    all_failed.extend(failed)
                else:
                    print("  unmapped: 0")

                processed.append(infile)

            except Exception as exc:
                logging.exception("Failed to process %s: %s", infile, exc)

    finally:
        cache.close()

    failed_path = write_failed(base_dir, all_failed)

    print("\nDone.")
    print(f"Processed: {len(processed)}")
    print(f"Failed mappings: {len(set(all_failed))}")

    if failed_path:
        print(f"Failed mapping file: {failed_path}")

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
